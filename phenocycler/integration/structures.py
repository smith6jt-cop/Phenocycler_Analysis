"""
S3 — tissue, islets, ducts and vessels: the things serial sections actually share.

Two sections cut 5 um apart share no cells. They do share architecture: the same islets,
the same ducts, the same vessels, in almost the same places. So structures — not cells —
are what registration aligns and what matching pairs.

The single most important design decision here is that **both modalities run the same
algorithm**. If PhenoCycler islets came from QuPath annotations and Xenium islets from
DBSCAN, any difference between the two would be confounded by the method and the comparison
would mean nothing. So islets are called identically on both sides, and PhenoCycler's
existing QuPath ``islet_num`` is used to *validate* that call rather than replace it.

The seeding rule is inherited from the Xenium side, where it was established empirically:
seed DBSCAN on **hormone-positive cells**, not on the endocrine *label*. Label-seeding found
2-9 islets per section; hormone-seeding found 112+. The label is too permissive — it
includes CD99-bright and generic endocrine cells scattered through the exocrine tissue,
which bridge real islets into one giant cluster.

That rule transfers to PhenoCycler almost for free, and the coincidence is a genuinely
useful one:

    Xenium seed        any of IAPP/MAFA/PDX1/NKX6-1/SLC30A8/ARX/MAFB/GHRL/HHEX
                       at lognorm >= log1p(5)
    PhenoCycler seed   any of INS/GCG/SST at _norm >= 5

The PhenoCycler threshold is ``[lineage] hormone_min_norm``, which is already 5 — the same
value the hormone floor uses to reject noise-tail false endocrine calls. The two modalities
therefore land on structurally symmetric islet definitions without either being tuned to
match the other.

Everything here works on the cell-table contract, so it needs no imaging data and is
testable on a synthetic islet field.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, discover_partitions, read_cell_table
from .tissues import has_structures
from .tissues import structure_kinds as tissue_structure_kinds

#: Islet size classes, from `Xenium_Analysis/scripts/insulitis_analysis.py::SIZE_BINS`.
#: Kept identical so an islet's size class means the same thing in both modalities.
SIZE_BINS: tuple[tuple[float, float, str], ...] = (
    (0, 20, "tiny"), (21, 50, "small"), (51, 100, "med"),
    (101, 200, "large"), (201, np.inf, "xl"),
)
SIZE_CLASSES = tuple(b[2] for b in SIZE_BINS)

#: Insulitis grading thresholds, immune cells per 100 endocrine cells. Also from
#: insulitis_analysis.py — reused so a "peri_insulitis" islet is the same thing on both sides.
CT_PERI_THRESH = 3.0
CT_INSULITIS_THRESH = 6.0

#: Radius around an islet counted as its peri-islet zone.
PERI_ZONE_UM = 20.0

#: Common lineages whose presence defines each structure type.
ISLET_LINEAGE = "Endocrine"
DUCT_LINEAGE = "Exocrine"
VESSEL_LINEAGE = "Vascular"
IMMUNE_LINEAGES = ("T_NK", "B_Plasma", "Myeloid", "Granulocyte")

#: Structure kinds this module can extract.
STRUCTURE_KINDS = ("islet", "duct", "vessel")


@dataclass
class StructureParams:
    """DBSCAN + feature parameters, identical for both modalities by construction."""

    eps_um: float = 50.0
    min_samples: int = 10
    peri_zone_um: float = PERI_ZONE_UM
    #: Minimum cells for a duct/vessel cluster (looser than islets — these are landmarks,
    #: not analysis units).
    landmark_min_samples: int = 15
    landmark_eps_um: float = 40.0

    @classmethod
    def from_config(cls, cfg: PipelineConfig) -> "StructureParams":
        return cls(eps_um=float(cfg.islet_eps_um), min_samples=int(cfg.islet_min_samples))


def size_class(n: float) -> str:
    for lo, hi, name in SIZE_BINS:
        if lo <= n <= hi:
            return name
    return SIZE_CLASSES[-1]


def insulitis_grade(immune_per_100_endo: float) -> str:
    """``no_insulitis`` / ``peri_insulitis`` / ``insulitis`` from the per-100-endocrine rate."""
    if immune_per_100_endo >= CT_INSULITIS_THRESH:
        return "insulitis"
    if immune_per_100_endo >= CT_PERI_THRESH:
        return "peri_insulitis"
    return "no_insulitis"


def hormone_seed_mask(table: CellTable, *, hormone_min_norm: float = 5.0) -> np.ndarray:
    """Boolean mask of hormone-positive seed cells, per modality.

    PhenoCycler: any of INS/GCG/SST ``_norm >= hormone_min_norm``, read from the ``feat__*``
    columns. Falls back to ``endocrine_subtype`` (which ``export_pheno`` derived with the
    same threshold) and then to the Endocrine label.

    Xenium: cells whose fine type is a specific endocrine subtype (Beta/Alpha/Delta/Gamma/
    Epsilon). These come from the detection-gated fine tier, which requires a specific
    hormone-surrogate gene to be detected — so it is the same "specific hormone evidence"
    criterion the Xenium islet code applies directly to the counts matrix.
    """
    df = table.df
    if table.modality == "phenocycler":
        cols = [c for c in (f"feat__{m}" for m in ("INS", "GCG", "SST")) if c in df.columns]
        if cols:
            vals = df.loc[:, cols].to_numpy(np.float64)
            return np.nanmax(vals, axis=1) >= hormone_min_norm
        if "endocrine_subtype" in df.columns:
            return df["endocrine_subtype"].isin(["Beta", "Alpha", "Delta"]).to_numpy()
        return (df["lineage_common"] == ISLET_LINEAGE).to_numpy()

    specific = {"Beta", "Alpha", "Delta", "Gamma", "Epsilon"}
    fine = df["celltype_fine"].astype(str)
    mask = fine.isin(specific).to_numpy()
    if mask.sum() >= 1:
        return mask
    # Un-phenotyped or generic-only labelling: fall back to the broad endocrine class and say
    # so, because the resulting islets will be systematically over-merged.
    print("[structures] xenium: no specific endocrine subtypes found — falling back to the "
          "broad Endocrine label, which over-merges islets (see the module docstring)",
          flush=True)
    return (df["lineage_common"] == ISLET_LINEAGE).to_numpy()


def dbscan_labels(coords: np.ndarray, eps_um: float, min_samples: int) -> np.ndarray:
    """DBSCAN cluster labels (-1 = noise). Same parameters both modalities."""
    if len(coords) == 0:
        return np.zeros(0, dtype=int)
    from sklearn.cluster import DBSCAN

    return np.asarray(
        DBSCAN(eps=eps_um, min_samples=min_samples).fit_predict(np.asarray(coords, np.float64))
    )


def _hull_area(points: np.ndarray) -> float:
    """Convex-hull area in um^2; 0 for degenerate clusters."""
    if len(points) < 3:
        return 0.0
    try:
        from scipy.spatial import ConvexHull

        return float(ConvexHull(points).volume)   # 'volume' is area in 2-D
    except Exception:  # noqa: BLE001 - collinear points, scipy absent
        return 0.0


def _eccentricity(points: np.ndarray) -> float:
    """Elongation from the covariance eigenvalues; 0 = circular, ->1 = a line.

    Used as a matching feature: a long thin islet stays long and thin in the next section,
    so shape is extra evidence beyond position and size.
    """
    if len(points) < 3:
        return 0.0
    cov = np.cov(points.T)
    if not np.all(np.isfinite(cov)):
        return 0.0
    vals = np.linalg.eigvalsh(cov)
    vals = np.clip(vals, 0, None)
    if vals[-1] <= 0:
        return 0.0
    return float(np.sqrt(max(0.0, 1.0 - vals[0] / vals[-1])))


def _composition(sub: pd.DataFrame, modality: str) -> dict:
    """beta/alpha/delta counts inside one islet, from whichever labels the modality has."""
    if modality == "phenocycler" and "endocrine_subtype" in sub.columns:
        vc = sub["endocrine_subtype"].value_counts()
    else:
        vc = sub["celltype_fine"].value_counts()
    return {
        "n_beta": int(vc.get("Beta", 0)),
        "n_alpha": int(vc.get("Alpha", 0)),
        "n_delta": int(vc.get("Delta", 0)),
    }


def call_islets(
    table: CellTable,
    params: StructureParams,
    *,
    hormone_min_norm: float = 5.0,
    registered: bool = False,
) -> tuple[pd.DataFrame, pd.Series]:
    """Find islets and describe them.

    Returns ``(islets, per_cell_islet_id)`` — the second aligned to ``table.df`` so callers
    can write ``structure_id`` straight back into the cell table.
    """
    df = table.df
    xy = table.xy(registered=registered)
    seed = hormone_seed_mask(table, hormone_min_norm=hormone_min_norm)

    labels = np.full(len(df), -1, dtype=int)
    if seed.sum() >= params.min_samples:
        labels[seed] = dbscan_labels(xy[seed], params.eps_um, params.min_samples)

    n_islets = int(labels.max()) + 1 if (labels >= 0).any() else 0
    prefix = f"{table.modality[:4]}_{table.donor_id}"
    islet_id = np.full(len(df), "", dtype=object)

    immune_mask = df["lineage_common"].isin(IMMUNE_LINEAGES).to_numpy()
    immune_xy = xy[immune_mask]

    rows = []
    for k in range(n_islets):
        members = labels == k
        pts = xy[members]
        n_endo = int(members.sum())
        cx, cy = float(pts[:, 0].mean()), float(pts[:, 1].mean())
        sid = f"{prefix}_islet_{k:04d}"
        islet_id[members] = sid

        # Immune infiltration within the islet + its peri zone. This is what makes an islet
        # comparable across modalities beyond its geometry: an insulitic islet should be
        # insulitic in both sections.
        n_immune = 0
        if len(immune_xy):
            radius = float(np.max(np.linalg.norm(pts - [cx, cy], axis=1))) + params.peri_zone_um
            d = np.linalg.norm(immune_xy - [cx, cy], axis=1)
            n_immune = int((d <= radius).sum())
        per100 = 100.0 * n_immune / max(n_endo, 1)

        comp = _composition(df.loc[members], table.modality)
        rows.append({
            "structure_id": sid,
            "kind": "islet",
            "modality": table.modality,
            "donor_id": table.donor_id,
            "roi": table.roi,
            "x_um": cx,
            "y_um": cy,
            "n_cells": n_endo,
            "area_um2": _hull_area(pts),
            "eccentricity": _eccentricity(pts),
            "size_class": size_class(n_endo),
            "n_immune_zone": n_immune,
            "immune_per_100_endo": per100,
            "insulitis_grade": insulitis_grade(per100),
            **comp,
            "frac_beta": comp["n_beta"] / max(n_endo, 1),
            "frac_alpha": comp["n_alpha"] / max(n_endo, 1),
            "frac_delta": comp["n_delta"] / max(n_endo, 1),
        })

    islets = pd.DataFrame(rows)
    return islets, pd.Series(islet_id, index=df.index, dtype="object")


def call_landmarks(
    table: CellTable,
    params: StructureParams,
    kind: str,
    *,
    registered: bool = False,
) -> pd.DataFrame:
    """Duct or vessel clusters — the secondary fiducials.

    Needed because islets alone are not always enough: in advanced T1D donors beta-cell loss
    thins the islet field badly, and a section with a handful of small islets cannot
    constrain an affine transform. Ducts (Pan-CK / KRT19) and vessels (CD31 / PECAM1) persist
    regardless of disease state.
    """
    lineage = {"duct": DUCT_LINEAGE, "vessel": VESSEL_LINEAGE}[kind]
    df = table.df
    xy = table.xy(registered=registered)
    mask = (df["lineage_common"] == lineage).to_numpy()

    if mask.sum() < params.landmark_min_samples:
        return pd.DataFrame(columns=["structure_id", "kind", "modality", "donor_id", "roi",
                                     "x_um", "y_um", "n_cells", "area_um2", "eccentricity"])

    labels = dbscan_labels(xy[mask], params.landmark_eps_um, params.landmark_min_samples)
    n = int(labels.max()) + 1 if (labels >= 0).any() else 0
    prefix = f"{table.modality[:4]}_{table.donor_id}"

    sub_xy = xy[mask]
    rows = []
    for k in range(n):
        pts = sub_xy[labels == k]
        rows.append({
            "structure_id": f"{prefix}_{kind}_{k:04d}",
            "kind": kind,
            "modality": table.modality,
            "donor_id": table.donor_id,
            "roi": table.roi,
            "x_um": float(pts[:, 0].mean()),
            "y_um": float(pts[:, 1].mean()),
            "n_cells": int(len(pts)),
            "area_um2": _hull_area(pts),
            "eccentricity": _eccentricity(pts),
        })
    return pd.DataFrame(rows)


def validate_against_qupath(table: CellTable, islets: pd.DataFrame,
                            islet_id: pd.Series) -> Optional[dict]:
    """Compare the DBSCAN islet call against PhenoCycler's QuPath annotations.

    PhenoCycler carries hand-drawn islet annotations in ``structure_id`` (from
    ``islet_num``). We do not use them for the call — symmetry with Xenium matters more —
    but they are an independent check that the DBSCAN parameters are sane on this modality,
    and a silent disagreement would be worth knowing about.
    """
    if table.modality != "phenocycler":
        return None
    df = table.df
    qp = df["structure_id"].astype(str)
    annotated = qp.ne("") & qp.ne("nan")
    if not annotated.any():
        return None

    n_qupath = int(qp[annotated].nunique())
    n_dbscan = len(islets)
    # Agreement: of cells QuPath calls islet, how many did DBSCAN also assign to one?
    both = annotated & islet_id.ne("")
    recall = float(both.sum() / max(int(annotated.sum()), 1))
    return {
        "n_qupath_islets": n_qupath,
        "n_dbscan_islets": n_dbscan,
        "qupath_cell_recall": recall,
    }


def extract_donor(
    table: CellTable,
    cfg: PipelineConfig,
    *,
    kinds: Optional[Iterable[str]] = None,
    registered: bool = False,
) -> tuple[dict[str, pd.DataFrame], pd.Series]:
    """All structure kinds for one section. Returns ``({kind: frame}, per_cell_islet_id)``.

    Kinds default to what the section's TISSUE supports. Lymph node supports none — it has no
    islets, so hormone-seeded DBSCAN has nothing to seed on and islet matching would have
    nothing to match. The registration-free levels (donor, niche) and the registered grid
    still apply there; see ``phenocycler.integration.tissues``.
    """
    params = StructureParams.from_config(cfg)
    tissue = table.tissue or "pancreas"
    if kinds is None:
        kinds = tissue_structure_kinds(tissue)
    kinds = list(kinds)
    if not kinds:
        empty = pd.Series([""] * len(table.df), index=table.df.index, dtype="object")
        return {}, empty
    out: dict[str, pd.DataFrame] = {}

    islet_id = pd.Series([""] * len(table.df), index=table.df.index, dtype="object")
    if "islet" in kinds:
        islets, islet_id = call_islets(table, params,
                                       hormone_min_norm=cfg.hormone_min_norm,
                                       registered=registered)
        out["islet"] = islets
        check = validate_against_qupath(table, islets, islet_id)
        if check:
            print(f"[structures] {table.modality}/{table.donor_id}: DBSCAN found "
                  f"{check['n_dbscan_islets']} islets vs {check['n_qupath_islets']} QuPath "
                  f"annotations; {check['qupath_cell_recall']:.0%} of QuPath islet cells "
                  f"landed in a DBSCAN islet", flush=True)

    for kind in ("duct", "vessel"):
        if kind in kinds:
            out[kind] = call_landmarks(table, params, kind, registered=registered)

    return out, islet_id


def structures_path(cfg: PipelineConfig, modality: str, donor_id: str, roi: str,
                    kind: str) -> Path:
    base = cfg.structures_dir / modality / f"donor_id={donor_id}"
    if roi:
        base = base / f"roi={roi}"
    return base / f"{kind}s.parquet"


def run_structures(
    cfg: PipelineConfig,
    *,
    modalities: Iterable[str] = ("phenocycler", "xenium"),
    donors: Optional[list[str]] = None,
    kinds: Optional[Iterable[str]] = None,
    tissues: Optional[Iterable[str]] = None,
) -> pd.DataFrame:
    """Stage entry point: extract structures for every exported section that supports them."""
    roots = {"phenocycler": cfg.cells_pheno_dir, "xenium": cfg.cells_xen_dir}
    summary = []

    for modality in modalities:
        root = roots[modality]
        for donor, roi in discover_partitions(root):
            if donors and donor not in set(donors):
                continue
            table = read_cell_table(root, donor, roi, modality=modality)
            if tissues and table.tissue not in set(tissues):
                continue
            if not has_structures(table.tissue or "pancreas"):
                print(f"[structures] {modality}/{donor}/{roi or '-'}: "
                      f"tissue '{table.tissue}' has no structural unit — skipped "
                      f"(donor + niche levels still apply)", flush=True)
                continue
            frames, islet_id = extract_donor(table, cfg, kinds=kinds)

            # Write the islet assignment back so the cell table knows its structure.
            if "islet" in frames:
                table.df["structure_id"] = islet_id.to_numpy()
                from .contract import partition_path, write_cell_table
                write_cell_table(table.df, partition_path(root, donor, roi))

            for kind, frame in frames.items():
                path = structures_path(cfg, modality, donor, roi, kind)
                path.parent.mkdir(parents=True, exist_ok=True)
                frame.to_parquet(path, index=False)

            row = {"modality": modality, "donor_id": donor, "roi": roi,
                   "tissue": table.tissue, "n_cells": table.n_cells}
            row.update({f"n_{k}": len(v) for k, v in frames.items()})
            if "islet" in frames and len(frames["islet"]):
                isl = frames["islet"]
                row["n_endocrine_in_islets"] = int(isl["n_cells"].sum())
                row["insulitic_islets"] = int((isl["insulitis_grade"] == "insulitis").sum())
            summary.append(row)
            print(f"[structures] {modality}/{donor}/{roi or '-'}: "
                  + ", ".join(f"{len(v)} {k}s" for k, v in frames.items()), flush=True)

    out = pd.DataFrame(summary)
    if out.empty:
        print("[structures] nothing to do — run the export/import stages first", flush=True)
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Extract islets/ducts/vessels for both modalities")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--modality", action="append", dest="modalities",
                    choices=["phenocycler", "xenium"])
    ap.add_argument("--tissue", action="append", dest="tissues",
                    help="restrict to one tissue (repeatable); default is every tissue")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    summary = run_structures(cfg, modalities=args.modalities or ("phenocycler", "xenium"),
                             donors=args.donors, tissues=args.tissues)
    if len(summary):
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
