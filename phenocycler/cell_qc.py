#!/usr/bin/env python3
"""
Step 1b — single-cell quality control, applied BEFORE any parameter is estimated.

Why this exists
---------------
Until METHOD v10 the pipeline had exactly one cell filter — ``cell_area >= 5 um^2`` in
:func:`phenocycler.cells_parquet.build_sql`, whose own comment conceded "real cells have area
>~12 um^2". Everything else (``sizes == 0 -> NaN``, the ``finite`` masks, ``valid_idx``) is a
NaN/validity guard, not QC.

That is load-bearing, because *every* RESTORE parameter — the raw Otsu/Triangle floors, the NNMF
arms, the Step-1 separator, the Step-2 divisor — is estimated over the unfiltered cloud, and the
divisor is a ``max``. One bad object therefore sets a donor-wide threshold. That is not
hypothetical: a 2-pixel object in donor 6450 raised the Vimentin Step-1 separator 4.2x, collapsed
its control set from 65,896 cells to 613, and quadrupled the Vimentin call.

Two tiers, because not every degenerate object is the same kind of problem
-------------------------------------------------------------------------
**Tier 1 — hard exclusion.** The object is not a cell. Removed from the canonical set, so it cannot
contribute to any statistic, any call, or any denominator. Measured cohort medians over the 20
eligible donors (31,454,041 cells):

===========================================  ==============  ==================================================
predicate                                    median removed  evidence it is not a cell
===========================================  ==============  ==================================================
``cell_area < min_cell_area`` (20 um^2)          4.31%       DAPI 0.67-0.88x normal *in every donor*; E-cadherin
                                                             4.5x lower, Vimentin 3x higher; 63-71% have N:C=1.00
``nucleus_area`` missing/<=0                    13.82%       DAPI **0.26-0.41x** normal in every donor. Median
                                                             area 36-61 um^2, so an area cut cannot catch these
``cell_area > 4 x donor-median area``            0.04%       above this the profile flips to E-cadherin 0.39x,
                                                             Pan-CK 0.29x, Vimentin 3.85x -- merged stromal and
                                                             vessel-wall blobs
``cell_solidity < 0.70``                         0.48%       area/convex-hull-area, the standard IMC
                                                             segmentation-artifact metric; concave/multi-lobed
rasterized area < 50% of reported               ~0%          polygon-geometry defect; confined to 6450 (1.52%)
                                                             and 6523 (1.20%), some rasterize to 2 px
``cell_area_px == 0``                            --          promotes the existing silent NaN rule to an
                                                             explicit, counted exclusion
===========================================  ==============  ==================================================

Combined Tier 1 removes a cohort median of 18.24% (15.39-27.45%).

**Tier 2 — excluded from parameter ESTIMATION only.** ``nc_area_ratio >= 0.999`` (no cytoplasmic
annulus; 8.10-13.66% of the retained population). These have nuclei and normal DAPI, so they are
real cells and still get a phenotype call — but their "whole-cell mean" is a nuclear mean, so they
must not set a membrane marker's floor, arm, separator or divisor.

What is deliberately NOT filtered
---------------------------------
**Circularity and aspect ratio.** Both are compartment-biased on this data: the cells they would
remove are real elongated stroma, endothelium and nerve, not artifacts. Measured marker medians of
the would-be-removed population, relative to the retained population:

    circularity < p1   E-cadherin 0.46x  Pan-CK 0.55x  Vimentin  4.57x
    aspect > 3         E-cadherin 0.09x  Pan-CK 0.13x  Vimentin 10.63x
    solidity < p1      E-cadherin 0.68x  Pan-CK 0.74x  Vimentin  1.56x   <- the one we use

Both are reported as diagnostics instead. Baharlou et al. (Front Immunol 2019) list circularity
among the standard morphology measures, but as something to *inspect*, not a threshold to apply
blindly -- and inspection here says do not apply it.

**A doublet cut.** Integrated DNA (``DAPI x nucleus_area``) *falls* to 0.73x in the largest area
decile while nucleus area rises only to 1.13x, so the large objects are large single cells, not
merged pairs. The published DNA x cell-length doublet gate would not flag them.

**A DAPI intensity floor.** Otsu on DAPI among nucleated cells puts 48-59% below threshold -- the
same two-class failure mode as E-cadherin. Nucleus presence is the direct criterion, and only
2.5-7.2% of nucleated cells fall below the non-nucleated DAPI median. An absolute floor would also
have to be shared across donors whose DAPI medians span 388-874, which the
no-cross-donor-normalisation rule forbids.

Honest scope
------------
This QC does NOT fix the E-cadherin under-call. Measured: the raw Otsu floor moves 1422.3 ->
1410-1429 under every filter variant and the maximum attainable E-cadherin+ fraction stays
23.2-23.6%. QC is necessary and correct on its own merits, and it removes the artifact class behind
the 6450 anomaly, but the compartment gap is a separate problem.

References
----------
Baker GJ et al., CyLinter, Nat Methods 2024 (nuclear-intensity/area artifact removal).
Windhager J et al., IMC data-analysis workflow ch.7 (min-area rule, per-marker SNR, image-level QC).
Baharlou H et al., Front Immunol 2019;10:2657 (histoCAT morphology QC set, DNA doublet gating).

    python -m phenocycler.cell_qc --all                 # compute + report, changes nothing
    python -m phenocycler.cell_qc --all --apply         # also filter data/cells to the retained set
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .cohort import ensure_eligible_donors
from .config import PipelineConfig, load_config

# Columns the QC needs from the canonical cells parquet.
QC_SOURCE_COLUMNS = [
    "object_id",
    "X_centroid",
    "Y_centroid",
    "cell_area",
    "nucleus_area",
    "nc_area_ratio",
    "cell_solidity",
    "cell_circularity",
    "cell_maxdiam",
    "cell_mindiam",
    "DAPI",
]

# Exclusion reasons, in the order they are applied. First matching reason wins, so a cell is
# attributed to exactly one cause and the per-reason counts sum to the total removed.
#
# `duplicate_detection` runs BEFORE the raster predicates on purpose. Duplicate polygons overlap, and
# the rasterizer resolves overlap by last-writer-wins, so the loser of each pair keeps only a sliver
# of its pixels and would otherwise be attributed to `raster_degenerate`. Removing the duplicate first
# attributes it to the real cause and leaves `raster_degenerate` to catch genuine geometry faults.
EXCLUSION_ORDER = (
    "no_nucleus",
    "area_below_min",
    "area_above_max",
    "low_solidity",
    "duplicate_detection",
    "raster_degenerate",
    "raster_empty",
)


def log(msg: str) -> None:
    print(msg, flush=True)


@dataclass(frozen=True)
class CellQCThresholds:
    """Frozen QC thresholds. Every value is donor-local or scale-free by construction."""

    min_cell_area: float = 20.0
    require_nucleus: bool = True
    max_area_median_multiple: float = 4.0
    min_cell_solidity: float = 0.70
    min_raster_area_ratio: float = 0.5
    no_cytoplasm_nc_ratio: float = 0.999
    pixel_area_um2: float = 0.5077663810243286 ** 2
    # Duplicate-detection removal. Two objects are the same physical cell when their centroids are
    # within `duplicate_radius_um` and their areas agree to within `duplicate_area_tol`.
    deduplicate: bool = True
    duplicate_radius_um: float = 1.0
    duplicate_area_tol: float = 0.15

    def __post_init__(self) -> None:
        if self.min_cell_area < 0:
            raise ValueError("min_cell_area must be >= 0")
        if self.duplicate_radius_um <= 0:
            raise ValueError("duplicate_radius_um must be > 0")
        if not 0 <= self.duplicate_area_tol < 1:
            raise ValueError("duplicate_area_tol must be in [0, 1)")
        if self.max_area_median_multiple <= 1:
            raise ValueError("max_area_median_multiple must be > 1")
        if not 0 <= self.min_cell_solidity <= 1:
            raise ValueError("min_cell_solidity must be in [0, 1]")
        if not 0 <= self.min_raster_area_ratio <= 1:
            raise ValueError("min_raster_area_ratio must be in [0, 1]")
        if not 0 < self.no_cytoplasm_nc_ratio <= 1:
            raise ValueError("no_cytoplasm_nc_ratio must be in (0, 1]")
        if self.pixel_area_um2 <= 0:
            raise ValueError("pixel_area_um2 must be > 0")

    @classmethod
    def from_config(cls, cfg: PipelineConfig) -> "CellQCThresholds":
        return cls(
            min_cell_area=cfg.cell_qc_min_cell_area,
            require_nucleus=cfg.cell_qc_require_nucleus,
            max_area_median_multiple=cfg.cell_qc_max_area_median_multiple,
            min_cell_solidity=cfg.cell_qc_min_cell_solidity,
            min_raster_area_ratio=cfg.cell_qc_min_raster_area_ratio,
            no_cytoplasm_nc_ratio=cfg.cell_qc_no_cytoplasm_nc_ratio,
            deduplicate=cfg.cell_qc_deduplicate,
            duplicate_radius_um=cfg.cell_qc_duplicate_radius_um,
            duplicate_area_tol=cfg.cell_qc_duplicate_area_tol,
        )


def find_duplicate_detections(
    x: np.ndarray,
    y: np.ndarray,
    area: np.ndarray,
    keep_score: np.ndarray,
    object_id: np.ndarray,
    eligible: np.ndarray,
    *,
    radius_um: float,
    area_tol: float,
) -> np.ndarray:
    """Boolean mask of the duplicate copies to DROP — one survivor is kept per cluster.

    Cause (confirmed from the QuPath workflow history): detection was run with ``selectAnnotations()``
    selecting the Tissue annotation AND every nested ``Islet_N`` child. InstanSeg's ``detectObjects()``
    runs inference per selected parent and does not deduplicate across parents, so every nucleus inside
    an islet was segmented twice — once under Tissue, once under its islet. The two passes tile from
    different ROI origins, so the copies are near-identical but not equal: centroids ~1 px apart, areas
    within a fraction of a percent, different vertex counts, and two distinct UUIDs. Both copies are
    then re-parented to the islet, which is why 100% of duplicates land in ``cell_region == 'core'``.

    Two objects are the same cell when their centroids are within ``radius_um`` AND their areas agree
    to within ``area_tol``. Clusters are connected components of that relation, so a 3-way overlap
    keeps exactly one survivor rather than dropping a pair and orphaning the third.

    The survivor is the copy with the highest ``keep_score`` — the rasterized pixel count, so the
    retained copy is the one whose polygon actually won the mask and therefore carries intact REDSEA
    values; the loser's REDSEA measurement was corrupted by the overlap. Ties break on the
    lexicographically smallest ``object_id``, which makes the choice deterministic and reproducible.

    Only cells flagged ``eligible`` participate. Running this after the other Tier-1 predicates means a
    pair whose partner already failed QC resolves for free, with no arbitrary choice made.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    drop = np.zeros(len(x), dtype=bool)
    idx = np.flatnonzero(eligible)
    if len(idx) < 2:
        return drop

    tree = cKDTree(np.column_stack([x[idx], y[idx]]))
    pairs = tree.query_pairs(r=radius_um, output_type="ndarray")
    if not len(pairs):
        return drop

    a0, a1 = area[idx][pairs[:, 0]], area[idx][pairs[:, 1]]
    denom = np.maximum(np.maximum(a0, a1), 1e-12)
    pairs = pairs[np.abs(a0 - a1) / denom < area_tol]
    if not len(pairs):
        return drop

    n = len(idx)
    graph = coo_matrix(
        (np.ones(len(pairs), dtype=np.int8), (pairs[:, 0], pairs[:, 1])), shape=(n, n)
    )
    n_comp, labels = connected_components(graph, directed=False)

    order = np.lexsort((object_id[idx], -keep_score[idx]))   # best keep_score first, then id asc
    seen = np.zeros(n_comp, dtype=bool)
    for local in order:
        comp = labels[local]
        if seen[comp]:
            drop[idx[local]] = True     # a survivor already claimed this cluster
        else:
            seen[comp] = True
    return drop


def _read_raster_px(cfg: PipelineConfig, donor: str) -> pd.DataFrame | None:
    """Rasterized pixel count per cell, from the REDSEA output. None if REDSEA has not run.

    ``cell_area_px`` is REDSEA's own polygon rasterization; comparing it to QuPath's reported
    ``cell_area`` catches degenerate geometry that no morphology column reveals.
    """
    files = sorted(glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet")))
    if not files:
        return None
    return pd.read_parquet(files[0], columns=["object_id", "cell_area_px"])


def compute_donor_qc(
    cfg: PipelineConfig, donor: str, thresholds: CellQCThresholds | None = None
) -> pd.DataFrame:
    """Per-cell QC verdicts for one donor.

    Returns one row per canonical cell with:
      ``object_id``, ``qc_pass`` (Tier 1), ``qc_exclude_reason`` ("" when passing),
      ``qc_no_cytoplasm`` (Tier 2 flag), ``qc_estimation_eligible`` (Tier 1 AND NOT Tier 2),
      plus the morphology used, so the verdict is auditable without re-reading the source.
    """
    thresholds = thresholds or CellQCThresholds.from_config(cfg)
    files = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    if not files:
        raise FileNotFoundError(f"donor {donor}: no canonical cells parquet under {cfg.cells_dir}")
    d = pd.read_parquet(files[0], columns=QC_SOURCE_COLUMNS)

    raster = _read_raster_px(cfg, donor)
    if raster is not None:
        d = d.merge(raster, on="object_id", how="left", validate="one_to_one")
    else:
        d["cell_area_px"] = np.nan
        log(f"  [{donor}] WARNING: no REDSEA output — raster-consistency predicates are skipped")

    area = d["cell_area"].to_numpy(dtype=float)
    nucleus = d["nucleus_area"].to_numpy(dtype=float)
    solidity = d["cell_solidity"].to_numpy(dtype=float)
    px = d["cell_area_px"].to_numpy(dtype=float)

    no_nucleus = ~np.isfinite(nucleus) | (nucleus <= 0)
    if not thresholds.require_nucleus:
        no_nucleus = np.zeros(len(d), dtype=bool)
    area_below = ~np.isfinite(area) | (area < thresholds.min_cell_area)

    # The upper bound is a multiple of the DONOR's own median, computed on the population that has
    # already survived the two size-independent predicates. An absolute um^2 cut would not be
    # comparable across donors whose median cell area spans 42.5-79.9 um^2.
    provisional = ~(no_nucleus | area_below)
    median_area = float(np.median(area[provisional])) if provisional.any() else np.nan
    max_area = thresholds.max_area_median_multiple * median_area
    area_above = np.isfinite(area) & np.isfinite(max_area) & (area > max_area)

    low_solidity = np.isfinite(solidity) & (solidity < thresholds.min_cell_solidity)

    expected_px = area / thresholds.pixel_area_um2
    with np.errstate(invalid="ignore", divide="ignore"):
        raster_ratio = np.where(expected_px > 0, px / expected_px, np.nan)
    raster_empty = np.isfinite(px) & (px <= 0)
    raster_degenerate = (
        np.isfinite(raster_ratio)
        & (raster_ratio < thresholds.min_raster_area_ratio)
        & ~raster_empty
    )

    # Duplicates are resolved among the cells that survive the other Tier-1 predicates, so a pair
    # whose partner already failed QC needs no arbitrary choice. Ordering matters: this must run
    # BEFORE the raster predicates, because a duplicate's loser is exactly what looks
    # "raster degenerate" once the winning polygon overwrites it.
    duplicate = np.zeros(len(d), dtype=bool)
    if thresholds.deduplicate:
        eligible = ~(no_nucleus | area_below | area_above | low_solidity)
        duplicate = find_duplicate_detections(
            d["X_centroid"].to_numpy(dtype=float),
            d["Y_centroid"].to_numpy(dtype=float),
            area,
            np.nan_to_num(px, nan=-1.0),           # rasterized px: the copy that won the mask
            d["object_id"].to_numpy(),
            eligible,
            radius_um=thresholds.duplicate_radius_um,
            area_tol=thresholds.duplicate_area_tol,
        )

    flags = {
        "no_nucleus": no_nucleus,
        "area_below_min": area_below,
        "area_above_max": area_above,
        "low_solidity": low_solidity,
        "duplicate_detection": duplicate,
        "raster_degenerate": raster_degenerate,
        "raster_empty": raster_empty,
    }

    # First reason wins, so per-reason counts partition the excluded set exactly.
    reason = np.full(len(d), "", dtype=object)
    claimed = np.zeros(len(d), dtype=bool)
    for name in EXCLUSION_ORDER:
        take = flags[name] & ~claimed
        reason[take] = name
        claimed |= take

    nc = d["nc_area_ratio"].to_numpy(dtype=float)
    no_cytoplasm = np.isfinite(nc) & (nc >= thresholds.no_cytoplasm_nc_ratio)

    out = pd.DataFrame(
        {
            "object_id": d["object_id"].to_numpy(),
            "qc_pass": ~claimed,
            "qc_exclude_reason": reason,
            "qc_no_cytoplasm": no_cytoplasm,
            "qc_estimation_eligible": (~claimed) & (~no_cytoplasm),
            "cell_area": area,
            "nucleus_area": nucleus,
            "nc_area_ratio": nc,
            "cell_solidity": solidity,
            "cell_area_px": px,
            "raster_area_ratio": raster_ratio,
        }
    )
    out.attrs["donor_median_area"] = median_area
    out.attrs["donor_max_area"] = max_area
    return out


def donor_summary(donor: str, qc: pd.DataFrame) -> dict:
    """One-row summary for the cohort QC table. Counts partition the excluded set exactly."""
    n = len(qc)
    row = {
        "donor": donor,
        "n": n,
        "retained": int(qc.qc_pass.sum()),
        "excluded": int((~qc.qc_pass).sum()),
        "excluded_pct": round(100.0 * float((~qc.qc_pass).mean()), 3),
        "estimation_eligible": int(qc.qc_estimation_eligible.sum()),
        "no_cytoplasm_of_retained_pct": round(
            100.0 * float(qc.loc[qc.qc_pass, "qc_no_cytoplasm"].mean()) if qc.qc_pass.any() else 0.0,
            3,
        ),
        "median_area": round(float(qc.attrs.get("donor_median_area", np.nan)), 3),
        "max_area_cut": round(float(qc.attrs.get("donor_max_area", np.nan)), 3),
    }
    for name in EXCLUSION_ORDER:
        row[f"n_{name}"] = int((qc.qc_exclude_reason == name).sum())
    return row


def apply_qc_to_canonical(cfg: PipelineConfig, donor: str, qc: pd.DataFrame) -> int:
    """Atomically filter this donor's canonical cells parquet to the Tier-1 retained set.

    Shrinking ``data/cells`` is sufficient and self-propagating: REDSEA's ``canonical_cell_mask``,
    the RESTORE compartment extraction and ``load_validation_sample`` all derive their cell set from
    it by join or set membership. REDSEA itself needs no recompute — it deliberately rasterizes every
    GeoJSON object for the contact graph (a fragment still donates real spillover to its neighbours)
    and only restricts rows at write time.

    The pre-QC parquet is preserved under ``<data_dir>/cells_preqc/donor_id=<id>/`` so the operation
    is reversible. It must NOT sit beside the filtered file: several stages glob
    ``cells/donor_id=*/*.parquet`` and expect exactly one match — ``redsea.reconcile_existing_donor``
    raises "expected one raw and one REDSEA parquet" the moment a second file appears there.
    """
    files = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    if len(files) != 1:
        raise ValueError(f"donor {donor}: expected exactly one canonical parquet, found {len(files)}")
    path = Path(files[0])
    backup_dir = cfg.data_dir / "cells_preqc" / f"donor_id={donor}"
    backup_dir.mkdir(parents=True, exist_ok=True)
    backup = backup_dir / path.name

    df = pd.read_parquet(path)
    keep_ids = set(qc.loc[qc.qc_pass, "object_id"])
    filtered = df[df.object_id.isin(keep_ids)].reset_index(drop=True)
    if len(filtered) != int(qc.qc_pass.sum()):
        raise AssertionError(
            f"donor {donor}: filtered {len(filtered)} rows but QC retained {int(qc.qc_pass.sum())}"
        )

    if not backup.exists():
        os.replace(path, backup)
    tmp = path.with_suffix(f".tmp{os.getpid()}.parquet")
    filtered.to_parquet(tmp, index=False, compression="zstd")
    os.replace(tmp, path)
    return len(filtered)


def run_cell_qc(
    cfg: PipelineConfig,
    donors: list[str],
    *,
    apply: bool = False,
    thresholds: CellQCThresholds | None = None,
) -> pd.DataFrame:
    """Compute (and optionally apply) cell QC for every donor. Returns the cohort summary."""
    ensure_eligible_donors(donors)
    thresholds = thresholds or CellQCThresholds.from_config(cfg)
    out_dir = cfg.data_dir / "cell_qc"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for donor in donors:
        qc = compute_donor_qc(cfg, donor, thresholds)
        part = out_dir / f"donor_id={donor}"
        part.mkdir(parents=True, exist_ok=True)
        tmp = part / f"data_0.tmp{os.getpid()}.parquet"
        qc.to_parquet(tmp, index=False, compression="zstd")
        os.replace(tmp, part / "data_0.parquet")

        row = donor_summary(donor, qc)
        if apply:
            row["applied_rows"] = apply_qc_to_canonical(cfg, donor, qc)
        rows.append(row)
        log(
            f"[{donor}] {row['n']:,} -> {row['retained']:,} retained "
            f"({row['excluded_pct']:.2f}% excluded) | estimation-eligible {row['estimation_eligible']:,} | "
            + ", ".join(f"{k}={row[f'n_{k}']:,}" for k in EXCLUSION_ORDER if row[f"n_{k}"])
        )

    summary = pd.DataFrame(rows)
    tmp = out_dir / f"cohort_summary.tmp{os.getpid()}.csv"
    summary.to_csv(tmp, index=False)
    os.replace(tmp, out_dir / "cohort_summary.csv")
    (out_dir / "thresholds.json").write_text(json.dumps(asdict(thresholds), indent=2) + "\n")

    total, retained = int(summary.n.sum()), int(summary.retained.sum())
    log(
        f"\nCOHORT: {total:,} -> {retained:,} retained "
        f"({100 * (total - retained) / total:.2f}% excluded); "
        f"estimation population {int(summary.estimation_eligible.sum()):,}"
    )
    if apply:
        log("canonical data/cells filtered in place; pre-QC copies kept as *.preqc.parquet")
    else:
        log("dry run — canonical data/cells unchanged (pass --apply to filter)")
    return summary


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Single-cell QC before any parameter estimation.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", default=None)
    ap.add_argument("--all", action="store_true", help="every eligible donor")
    ap.add_argument("--apply", action="store_true",
                    help="filter data/cells to the retained set (default: report only)")
    ap.add_argument("--min-cell-area", type=float, default=None)
    ap.add_argument("--min-solidity", type=float, default=None)
    ap.add_argument("--max-area-median-multiple", type=float, default=None)
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    donors = list(a.donor) if a.donor else (cfg.discover_donors() if a.all else [])
    if not donors:
        ap.error("pass --donor <id> (repeatable) or --all")

    base = CellQCThresholds.from_config(cfg)
    over = {}
    if a.min_cell_area is not None:
        over["min_cell_area"] = a.min_cell_area
    if a.min_solidity is not None:
        over["min_cell_solidity"] = a.min_solidity
    if a.max_area_median_multiple is not None:
        over["max_area_median_multiple"] = a.max_area_median_multiple
    thresholds = CellQCThresholds(**{**asdict(base), **over}) if over else base

    run_cell_qc(cfg, donors, apply=a.apply, thresholds=thresholds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
