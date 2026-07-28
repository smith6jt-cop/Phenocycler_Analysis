"""
S1a — PhenoCycler outputs -> the common cell-table contract.

Reads the end of the core pipeline and reshapes it. Three sources, all joined on
``object_id`` (the QuPath ``PathObject.getID()`` UUID, which threads through every stage
with no centroid rounding, so the joins are exact):

    data/phenotype/broad/        broad_lineage, assign_margin, epi_default, score_*,
                                 cell_region, islet_num
    data/cells/                  X_centroid / Y_centroid in microns, cell area
    data/restore_gated_redsea/   {marker}_pos, {marker}_norm   (the validated 10)
    ..._extra/                   CD99 / B3TUBB / MPO           (the second RESTORE pass)

Nothing is recomputed. The 8-class call is taken as given; this stage only *widens* it.

Two derived columns do real work, and both exist because the 8-class call is deliberately
coarser than the data underneath it:

``immune_subclass``
    PhenoCycler collapses T, B and myeloid into one ``Immune`` lineage, while Xenium
    resolves all three as separate broad classes. Comparing them as-is would either throw
    away Xenium's resolution or invent PhenoCycler's. But ``CD3e_pos`` / ``CD20_pos`` /
    ``CD163_pos`` are still sitting in the gated parquet — the information was never lost,
    only not used by ``lineage.py``. Recovering it here costs one read and needs no
    re-running of RESTORE.

``endocrine_subtype``
    Same argument one level down: the Endocrine class lumps INS/GCG/SST/CD99 together, which
    masks beta-loss (INS down) behind alpha-persistence (GCG up) — the single most important
    signal in a T1D cohort. Xenium resolves Beta/Alpha/Delta as fine types, so we take the
    argmax of the three hormone ``_norm`` values above the floor.

Neither derivation changes ``broad_lineage``; both are additive columns beside it.
"""

from __future__ import annotations

import argparse
import glob
from functools import partial
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import EXTRA_MARKERS, HORMONE_MARKERS, LINEAGES, PipelineConfig, load_config
from ..parallel import map_donors
from .contract import as_feature, coerce_cell_table, partition_path, write_cell_table
from .vocab import map_pheno_series

#: Marker -> immune subclass. Mirrors `LINEAGES['Immune']`; a cell positive for more than one
#: is 'Mixed' rather than being silently resolved by column order.
IMMUNE_GATES: dict[str, str] = {"CD3e": "T", "CD20": "B", "CD163": "Myeloid"}

#: Hormone -> endocrine fine type, for the Beta/Alpha/Delta split.
HORMONE_SUBTYPE: dict[str, str] = {"INS": "Beta", "GCG": "Alpha", "SST": "Delta"}

#: Every marker with `_norm` / `_pos` columns across the two RESTORE passes.
ALL_GATED_MARKERS: list[str] = sorted({m for ms in LINEAGES.values() for m in ms})
MAIN_MARKERS: list[str] = [m for m in ALL_GATED_MARKERS if m not in EXTRA_MARKERS]


class ExportError(RuntimeError):
    """A PhenoCycler donor cannot be exported to the contract."""


def _read_first_parquet(directory: Path, donor: str, columns=None) -> Optional[pd.DataFrame]:
    files = sorted(glob.glob(str(Path(directory) / f"donor_id={donor}" / "*.parquet")))
    if not files:
        return None
    frames = []
    for f in files:
        try:
            frames.append(pd.read_parquet(f, columns=columns))
        except (ValueError, KeyError):
            # Column subset not present in this file — read it whole and let the caller
            # decide. Better than failing the donor over an optional column.
            frames.append(pd.read_parquet(f))
    out = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    return out


def derive_immune_subclass(gated: pd.DataFrame) -> pd.Series:
    """T / B / Myeloid / Mixed from the three immune gates; '' for immune-negative cells.

    'Mixed' is kept as its own value rather than resolved to a winner. A CD3e+CD163+ cell is
    usually a segmentation artifact or a doublet, and quietly calling it T would put a known
    artifact into a lineage comparison. ``vocab.PHENO_IMMUNE_TO_COMMON`` maps it to the
    least-specific common class.
    """
    present = [m for m in IMMUNE_GATES if f"{m}_pos" in gated.columns]
    if not present:
        return pd.Series([""] * len(gated), index=gated.index, dtype="object")

    flags = {m: gated[f"{m}_pos"].fillna(False).to_numpy(bool) for m in present}
    n_pos = np.sum(list(flags.values()), axis=0)
    out = np.full(len(gated), "", dtype=object)
    for m in present:
        out[flags[m] & (n_pos == 1)] = IMMUNE_GATES[m]
    out[n_pos > 1] = "Mixed"
    return pd.Series(out, index=gated.index, dtype="object")


def derive_endocrine_subtype(gated: pd.DataFrame, hormone_min_norm: float) -> pd.Series:
    """Beta / Alpha / Delta by hormone argmax; 'Endocrine' when only CD99 fires.

    Uses the same threshold as the hormone floor, so a cell counted as hormone-positive by
    the lineage call is exactly a cell that gets a subtype here. Cells with no hormone above
    the floor but bright CD99 keep the generic label — CD99 marks the endocrine types this
    panel has no specific marker for (PP / epsilon / EC).
    """
    cols = [m for m in HORMONE_MARKERS if f"{m}_norm" in gated.columns]
    out = np.full(len(gated), "", dtype=object)
    if not cols:
        return pd.Series(out, index=gated.index, dtype="object")

    norms = np.vstack([gated[f"{m}_norm"].fillna(0.0).to_numpy(np.float64) for m in cols]).T
    best = norms.argmax(axis=1)
    top = norms.max(axis=1)
    hot = top >= hormone_min_norm
    labels = np.array([HORMONE_SUBTYPE[m] for m in cols], dtype=object)
    out[hot] = labels[best[hot]]

    if "CD99_pos" in gated.columns:
        cd99 = gated["CD99_pos"].fillna(False).to_numpy(bool)
        out[(~hot) & cd99] = "Endocrine"
    return pd.Series(out, index=gated.index, dtype="object")


def export_donor(
    donor: str,
    cfg: PipelineConfig,
    *,
    roi: str = "panc",
    markers: Optional[list[str]] = None,
) -> pd.DataFrame:
    """Build the contract table for one PhenoCycler donor."""
    broad = _read_first_parquet(cfg.broad_dir, donor)
    if broad is None:
        raise ExportError(
            f"donor {donor}: no broad-lineage output under {cfg.broad_dir}/donor_id={donor}. "
            f"Run `python -m phenocycler.pipeline` (the 'lineage' stage) first."
        )
    broad["object_id"] = broad["object_id"].astype(str)

    cells = _read_first_parquet(cfg.cells_dir, donor)
    if cells is None:
        raise ExportError(
            f"donor {donor}: no cell table under {cfg.cells_dir}/donor_id={donor}; "
            f"the 'cells' stage must run before integration (it carries the coordinates)"
        )
    cells["object_id"] = cells["object_id"].astype(str)

    coord_cols = [c for c in ("object_id", "X_centroid", "Y_centroid") if c in cells.columns]
    missing_coords = {"X_centroid", "Y_centroid"} - set(coord_cols)
    if missing_coords:
        raise ExportError(
            f"donor {donor}: data/cells is missing {sorted(missing_coords)} — integration "
            f"needs micron centroids, which build_cells_parquet writes from QuPath's "
            f"'Centroid X/Y um' columns"
        )
    area_col = next((c for c in cells.columns
                     if c.lower().startswith("cell_area") or c == "Cell_Area_um_2"), None)
    if area_col:
        coord_cols.append(area_col)
    if "image" in cells.columns:
        coord_cols.append("image")

    df = broad.merge(cells.loc[:, coord_cols], on="object_id", how="left")

    # RESTORE gates: main pass + the extra pass, merged the same way lineage.py does it.
    markers = markers or ALL_GATED_MARKERS
    gated = _read_first_parquet(cfg.restore_gated_dir, donor)
    if gated is not None:
        gated["object_id"] = gated["object_id"].astype(str)
        keep = ["object_id"] + [f"{m}_{s}" for m in MAIN_MARKERS for s in ("pos", "norm")
                                if f"{m}_{s}" in gated.columns]
        df = df.merge(gated.loc[:, [c for c in keep if c in gated.columns]],
                      on="object_id", how="left")

    extra = _read_first_parquet(cfg.restore_gated_extra_dir, donor)
    if extra is not None:
        extra["object_id"] = extra["object_id"].astype(str)
        keep = ["object_id"] + [f"{m}_{s}" for m in EXTRA_MARKERS for s in ("pos", "norm")
                                if f"{m}_{s}" in extra.columns]
        df = df.merge(extra.loc[:, [c for c in keep if c in extra.columns]],
                      on="object_id", how="left")
        # Re-derive the bright-CD99 gate exactly as lineage.py does, so `immune_subclass` and
        # `endocrine_subtype` agree with the broad call rather than using the permissive gate.
        if "CD99_norm" in df.columns:
            df["CD99_pos"] = df["CD99_norm"].fillna(0.0) >= cfg.cd99_bright

    immune_subclass = derive_immune_subclass(df)
    endocrine_subtype = derive_endocrine_subtype(df, cfg.hormone_min_norm)

    # Finest available label: hormone subtype inside Endocrine, immune subclass inside
    # Immune, otherwise the broad class itself.
    fine = df["broad_lineage"].astype(str).copy()
    is_endo = df["broad_lineage"].astype(str) == "Endocrine"
    fine[is_endo & endocrine_subtype.ne("")] = endocrine_subtype[is_endo & endocrine_subtype.ne("")]
    is_imm = df["broad_lineage"].astype(str) == "Immune"
    fine[is_imm & immune_subclass.ne("")] = immune_subclass[is_imm & immune_subclass.ne("")]

    # `epi_default` marks cells assigned Epithelial because every positive gate failed, not
    # because keratin fired. Those cells are not evidence-positive, and pretending otherwise
    # inflates PhenoCycler epithelium against Xenium exocrine in every composition table.
    epi_default = (df["epi_default"].fillna(False).to_numpy(bool)
                   if "epi_default" in df.columns else np.zeros(len(df), bool))

    out = pd.DataFrame({
        "cell_id": df["object_id"].astype(str),
        "donor_id": str(donor),
        "roi": roi,
        "x_um": pd.to_numeric(df["X_centroid"], errors="coerce"),
        "y_um": pd.to_numeric(df["Y_centroid"], errors="coerce"),
        "area_um2": (pd.to_numeric(df[area_col], errors="coerce")
                     if area_col else np.nan),
        "lineage_native": df["broad_lineage"].astype(str),
        "celltype_fine": fine.astype(str),
        "confidence": pd.to_numeric(df.get("assign_margin", np.nan), errors="coerce"),
        "evidence_positive": ~epi_default,
        "region": df.get("cell_region", pd.Series([""] * len(df))).fillna("").astype(str),
        # QuPath's own islet annotation. structures.py recomputes islets by DBSCAN for
        # cross-modal symmetry, but keeping this lets that call be validated against it.
        "structure_id": df.get("islet_num", pd.Series([""] * len(df))).fillna("").astype(str),
    })
    out["lineage_common"] = map_pheno_series(out["lineage_native"], immune_subclass)
    out["immune_subclass"] = immune_subclass.to_numpy()
    out["endocrine_subtype"] = endocrine_subtype.to_numpy()

    # Features: the per-marker RESTORE `_norm` (raw / per-image threshold), which is already
    # the cross-image-comparable quantity — that is the whole point of RESTORE.
    for m in markers:
        col = f"{m}_norm"
        if col in df.columns:
            out[as_feature(m)] = pd.to_numeric(df[col], errors="coerce")

    n_missing = int(out["x_um"].isna().sum())
    if n_missing:
        print(f"[export_pheno] donor {donor}: {n_missing:,} cells have no centroid "
              f"(absent from data/cells) — dropped", flush=True)
        out = out[out["x_um"].notna()].reset_index(drop=True)

    return coerce_cell_table(out, modality="phenocycler", donor_id=str(donor), roi=roi)


def _export_and_write(donor: str, cfg: PipelineConfig, roi: str) -> dict:
    df = export_donor(donor, cfg, roi=roi)
    path = partition_path(cfg.cells_pheno_dir, donor, roi)
    write_cell_table(df, path)
    counts = df["lineage_common"].value_counts().to_dict()
    print(f"[export_pheno] {donor}: {len(df):,} cells -> {path.parent}", flush=True)
    return {"donor_id": donor, "n_cells": len(df), **counts}


def run_export_pheno(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
    n_jobs: Optional[int] = None,
) -> pd.DataFrame:
    """Stage entry point: export every donor with a broad-lineage output."""
    donors = donors or cfg.discover_donors(cfg.broad_dir)
    if not donors:
        raise ExportError(
            f"no donors found under {cfg.broad_dir}. Run the core pipeline "
            f"(`python -m phenocycler.pipeline`) through the 'lineage' stage first."
        )
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    rows = map_donors(partial(_export_and_write, cfg=cfg, roi=roi), donors,
                      n_jobs=n_jobs, ordered=True)
    out = pd.DataFrame([r for r in rows if r]).fillna(0)
    if len(out):
        print(f"\n[export_pheno] {len(out)} donor(s), "
              f"{int(out['n_cells'].sum()):,} cells total", flush=True)
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="PhenoCycler -> integration cell-table contract")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--jobs", type=int, default=None)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    summary = run_export_pheno(cfg, args.donors, roi=args.roi, n_jobs=args.jobs)
    if len(summary):
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
