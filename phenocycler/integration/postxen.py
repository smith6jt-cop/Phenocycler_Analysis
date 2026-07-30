"""
S1c — the post-Xenium re-stain: protein and RNA on the *same cells*.

After the Xenium run, that section is stained for INS and GCG on the PhenoCycler. The output
is an ordinary qptiff, so it flows through the entire core pipeline unchanged — REDSEA,
RESTORE, the hormone floor, the lineage call. What makes it different is the slide underneath:
it images the **Xenium section**, so its cells are the same physical cells as the
transcriptome.

That is the only place in this cohort where protein and RNA are measured on one cell, and it
sits inside an otherwise serial-section design:

    PhenoCycler section  ──── serial ────▶  Xenium section
     (55-plex protein)                       (5K RNA)
                                                 │
                                            same slide
                                                 ▼
                                       post-Xenium re-stain
                                        (INS / GCG protein)

Two consequences worth being precise about.

**This pairing is a measurement, and the cohort mode does not change that.** ``sameslide.py``
guards on ``[integration] mode`` because it claims whole-cell correspondence for *every* cell,
which serial sections cannot support. Here the correspondence is a property of the section
pair, not of the cohort, so the guard would be answering the wrong question. It is still a
*matching* problem — QuPath segments the qptiff, Xenium segments DAPI plus its own boundary
stains, and the two will not agree on cell count or outline — but it is matching two views of
one thing rather than approximating across tissue.

**What it is good for is identity, not intensity.** Post-Xenium tissue has been through
protease digestion and decrosslinking, so epitopes are degraded and absolute intensities are
not comparable to a fresh PhenoCycler stain. INS and GCG are abundant secretory proteins and
survive it; CD3e reportedly does not, which is exactly what the mechanism predicts for a
lower-abundance surface marker. So this module carries the hormone-derived ``endocrine_subtype``
across and deliberately does **not** carry raw marker intensities into a cross-platform
comparison.

The payoff is that Xenium endocrine identity stops being an inference. INS, GCG and SST are
absent from the 5K panel, so Xenium normally calls Beta/Alpha from surrogate cores (PDX1,
ISL1, NEUROD1, ABCC8). With the re-stain those calls can be *checked against hormone protein
on the same cell* — and once stamped onto the Xenium table, the serial-section machinery
(S6b) gets endocrine types derived from protein on both sides, which is what makes
beta-to-beta matching across the pair meaningful.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from ..sections import SLIDE_XENIUM
from .contract import partition_path, read_cell_table, write_cell_table
from .manifest import load_manifest, paired_rows, rows_for_scope
from .sameslide import match_by_centroid, match_tables_by_iou

#: Centroid cap for "these two segmentations are the same cell". Tighter than the
#: serial-section cap: there is no section gap here, so the only disagreement is segmentation
#: itself, and two outlines of one cell agree on its centre to within a couple of microns.
DEFAULT_MAX_DIST_UM = 4.0

#: Columns carried from the re-stain onto the Xenium table. Identity only — see the module
#: docstring on why intensities from post-Xenium tissue do not travel.
CARRY_COLUMNS = ("endocrine_subtype",)


def _match(cfg, donor, roi, rest, xen, offset, *, method, max_dist_um, min_iou):
    """Match the re-stain to the transcriptome. Returns ``(pairs, method_used)``.

    IoU is the better criterion here for the same reason it is in same-slide mode, and more
    so: these really are one section, so two outlines of one cell overlap heavily while a
    neighbour does not — even when that neighbour's centroid is nearer, which in packed islet
    tissue it frequently is.

    The re-stain polygons are shifted by the same frame offset the centroids get, not warped:
    the two images are the same physical section on the same instrument, so the only
    difference is where the operator put the scan origin.
    """
    if method not in ("auto", "iou", "centroid"):
        raise ValueError(f"method must be auto|iou|centroid, got {method!r}")

    if method in ("auto", "iou"):
        try:
            from shapely.affinity import translate

            from .boundaries import load_pheno_polygons, load_xenium_polygons

            r_poly = load_pheno_polygons(cfg, donor, roi, base=cfg.cells_pheno_xif_dir)
            r_poly = r_poly.apply(lambda g: translate(g, float(offset[0]), float(offset[1])))

            man = rows_for_scope(load_manifest(cfg), roi=roi)
            row = man[man["donor_id"].astype(str) == str(donor)]
            bundle = str(row["xenium_bundle"].iloc[0]) if len(row) else ""
            if not bundle:
                raise RuntimeError(f"no xenium_bundle for {donor}/{roi} in the manifest")
            x_poly = load_xenium_polygons(Path(bundle), cell_ids=set(xen.df["cell_id"]))

            return match_tables_by_iou(r_poly, x_poly, rest, xen, min_iou=min_iou), "iou"
        except Exception as exc:  # noqa: BLE001 - fall back, but say so
            if method == "iou":
                raise
            print(f"[postxen] {donor}/{roi}: polygon matching unavailable ({exc}); "
                  f"falling back to centroids", flush=True)

    return match_by_centroid(rest.xy() + offset, xen.xy(), max_dist_um=max_dist_um), "centroid"


def pair_restain(
    cfg: PipelineConfig,
    donor: str,
    roi: str,
    *,
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
    method: str = "auto",
    min_iou: float = 0.25,
    n_perm: int = 25,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Match re-stain cells to Xenium cells on one section. Returns ``(pairs, stats)``.

    No registration is applied. Both are the same physical section on the same slide, so a
    residual translation is a *coordinate-frame* difference, not tissue movement — it is
    estimated here as a median offset and removed, rather than being handed to the full
    serial-section registration machinery, which is built for a problem this does not have.
    """
    rng = np.random.default_rng(seed)
    base = {"donor_id": donor, "roi": roi}
    try:
        rest = read_cell_table(cfg.cells_pheno_xif_dir, donor, roi, modality="phenocycler")
        xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")
    except FileNotFoundError as exc:
        return pd.DataFrame(), {**base, "error": str(exc)}
    if len(rest.df) < 10 or len(xen.df) < 10:
        return pd.DataFrame(), {**base, "error": "too few cells to pair"}

    r_xy, x_xy = rest.xy(), xen.xy()
    # Frame offset, not tissue motion: the two exports may be zeroed differently.
    offset = np.median(x_xy, axis=0) - np.median(r_xy, axis=0)
    r_xy = r_xy + offset

    m, used = _match(cfg, donor, roi, rest, xen, offset, method=method,
                     max_dist_um=max_dist_um, min_iou=min_iou)
    lo, hi = x_xy.min(axis=0), x_xy.max(axis=0)
    null = float(np.mean([
        len(match_by_centroid(r_xy, np.c_[rng.uniform(lo[0], hi[0], len(x_xy)),
                                          rng.uniform(lo[1], hi[1], len(x_xy))],
                              max_dist_um=max_dist_um)) for _ in range(n_perm)]))

    stats = {
        **base,
        "n_restain": len(rest.df), "n_xenium": len(xen.df), "n_pairs": len(m),
        "method": used,
        "frame_offset_um": float(np.hypot(*offset)),
        "pair_rate_restain": len(m) / len(rest.df),
        "pair_rate_xenium": len(m) / len(xen.df),
        "median_distance_um": (float(m["distance_um"].median())
                               if len(m) and "distance_um" in m.columns else float("nan")),
        "median_iou": (float(m["iou"].median())
                       if len(m) and "iou" in m.columns else float("nan")),
        "null_mean": null,
        "signal_over_null": len(m) / null if null > 0 else float("nan"),
    }
    if m.empty:
        stats["note"] = (f"no pairs within {max_dist_um} um — check that these really are the "
                         f"same section and that both exports are in microns")
        return pd.DataFrame(), stats

    pairs = pd.DataFrame({
        "donor_id": donor, "roi": roi,
        "restain_cell_id": rest.df["cell_id"].to_numpy()[m["moving_idx"].to_numpy(int)],
        "xen_cell_id": xen.df["cell_id"].to_numpy()[m["fixed_idx"].to_numpy(int)],
    })
    for col in ("distance_um", "iou"):
        if col in m.columns:
            pairs[col] = m[col].to_numpy()
    for col in CARRY_COLUMNS:
        if col in rest.df.columns:
            pairs[col] = rest.df[col].to_numpy()[m["moving_idx"].to_numpy(int)]
    return pairs, stats


def surrogate_agreement(pairs: pd.DataFrame, cfg: PipelineConfig,
                        donor: str, roi: str) -> dict:
    """How well Xenium's surrogate-panel endocrine call matches hormone protein.

    This is the number the re-stain exists to produce. Xenium infers Beta/Alpha from curated
    surrogate cores because INS/GCG/SST are off-panel; here the same cells have the hormones
    measured directly, so the inference can be scored instead of trusted.
    """
    out: dict = {"donor_id": donor, "roi": roi, "n_scored": 0}
    if pairs.empty or "endocrine_subtype" not in pairs.columns:
        return out
    try:
        xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")
    except FileNotFoundError:
        return out

    fine = dict(zip(xen.df["cell_id"], xen.df.get("celltype_fine", pd.Series(dtype=str))))
    protein = pairs["endocrine_subtype"].astype(str)
    rna = pairs["xen_cell_id"].map(lambda c: str(fine.get(c, "")))
    both = protein.isin(["Beta", "Alpha", "Delta"]) & rna.isin(["Beta", "Alpha", "Delta"])
    if both.sum() < 5:
        return out

    out["n_scored"] = int(both.sum())
    out["surrogate_agreement"] = float((protein[both] == rna[both]).mean())
    for t in ("Beta", "Alpha", "Delta"):
        sel = both & (protein == t)
        if sel.sum():
            out[f"recall_{t.lower()}"] = float((rna[sel] == t).mean())
            out[f"n_{t.lower()}"] = int(sel.sum())
    return out


def stamp_xenium(cfg: PipelineConfig, donor: str, roi: str, pairs: pd.DataFrame) -> int:
    """Write the hormone-derived ``endocrine_subtype`` onto the Xenium cell table.

    This is what makes the downstream serial-section matching (S6b) type both sides by the
    same kind of evidence — hormone protein — rather than protein on one side and surrogate
    RNA on the other. Cells with no re-stain partner keep whatever they had; the column is
    filled, never overwritten with blanks.
    """
    if pairs.empty or "endocrine_subtype" not in pairs.columns:
        return 0
    xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")
    df = xen.df
    mapping = dict(zip(pairs["xen_cell_id"], pairs["endocrine_subtype"].astype(str)))
    incoming = df["cell_id"].map(lambda c: mapping.get(c, ""))
    existing = df["endocrine_subtype"].astype(str) if "endocrine_subtype" in df.columns else ""
    df["endocrine_subtype"] = incoming.where(incoming != "", existing)
    write_cell_table(df, partition_path(cfg.cells_xen_dir, donor, roi))
    return int((incoming != "").sum())


def run_postxen(cfg: PipelineConfig, donors: Optional[list[str]] = None, *,
                roi: Optional[str] = None, tissue: Optional[str] = None,
                stamp: bool = True) -> pd.DataFrame:
    """Stage entry point: pair every re-stained section against its Xenium transcriptome."""
    if not cfg.cells_pheno_xif_dir.exists():
        print(f"[postxen] no post-Xenium re-stains under {cfg.cells_pheno_xif_dir} — nothing "
              f"to do. Sections are routed here by their section key carrying a "
              f"{SLIDE_XENIUM!r}-slide token; see phenocycler/sections.py.", flush=True)
        return pd.DataFrame()

    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]

    frames, rows = [], []
    for _, r in man.iterrows():
        donor, roi_name = r["donor_id"], r["roi"]
        pairs, stats = pair_restain(cfg, donor, roi_name)
        if "error" in stats:
            continue
        stats.update(surrogate_agreement(pairs, cfg, donor, roi_name))
        if stamp and len(pairs):
            stats["n_stamped"] = stamp_xenium(cfg, donor, roi_name, pairs)
        rows.append(stats)
        if len(pairs):
            frames.append(pairs)
        print(f"[postxen] {donor}/{roi_name}: {stats['n_pairs']:,} same-slide cell pairs "
              f"({stats['pair_rate_xenium']:.1%} of Xenium cells), median "
              f"{stats['median_distance_um']:.2f} um, {stats['signal_over_null']:.1f}x null"
              + (f", surrogate agreement {stats['surrogate_agreement']:.3f} "
                 f"over {stats['n_scored']} cells" if stats.get("n_scored") else ""),
              flush=True)

    if frames:
        out = cfg.paired_dir / "postxen_cells.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        allf = pd.concat(frames, ignore_index=True)
        allf.to_parquet(out, index=False)
        print(f"\n[postxen] wrote {out} ({len(allf):,} same-slide protein+RNA cell pairs — "
              f"MEASUREMENT, not inference)", flush=True)

    stats_df = pd.DataFrame(rows)
    if len(stats_df):
        cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(cfg.integration_qc_dir / "postxen_qc.csv", index=False)
    return stats_df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Pair the post-Xenium re-stain against the Xenium transcriptome")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None)
    ap.add_argument("--tissue", default=None)
    ap.add_argument("--no-stamp", action="store_true",
                    help="do not write endocrine_subtype onto the Xenium cell table")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_postxen(cfg, args.donors, roi=args.roi, tissue=args.tissue,
                        stamp=not args.no_stamp)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
