"""
S6b — cell-to-cell matching across serial sections, where it is actually valid.

This module exists because an earlier version of this codebase asserted, flatly, that two
serial sections never measure the same cell. That is wrong. Pancreatic cells are 10-15 um
and sections are 5 um, so a cell cannot fit inside one section: most straddle the cut and
appear in both. The claim was never measured, only assumed.

What a measurement shows is more interesting than either "never" or "always". Simulating a
realistic 3D population (6000 cells, 800 x 800 um, adjacent 5 um sections, nucleus-seeded
detection) gives:

* **~70% of cells straddle the boundary, but only ~28% are detectable in both.** Segmentation
  is nucleus-seeded, and a nucleus (~0.6x cell diameter) needs enough volume in a slab to be
  called at all.

* **The recoverable population is strongly size-biased**, and it falls exactly where it hurts:

      8 um (lymphocyte)   0%      11 um (endocrine)   9%
      10 um               4%      13 um (acinar)     15%

  An 8 um cell has a ~4.8 um nucleus and cannot put enough of it on both sides of the cut.
  **Immune infiltrate is not recoverable cell-to-cell at any registration quality.** That is
  geometry, not method, and no amount of care changes it — which is why this module refuses
  to run on immune types rather than reporting a low match rate for them.

* **Unconstrained matching is not a measurement.** At a 3 um cap with *perfect* registration:
  precision 0.72, with the permuted null accounting for a third of all matches, because cell
  spacing (~10 um) is close to the cap. One pair in four is wrong at the ceiling.

* **Type-constrained matching is.** Restricting to one cell type cuts candidate density ~5x,
  so the chance rate collapses while true matches are untouched:

      all cells,      0 um residual   precision 0.72   signal/null  2.9
      endocrine only, 0 um residual   precision 0.96   signal/null 12.9
      endocrine only, 2 um residual   precision 0.95   signal/null  9.4
      endocrine only, 5 um residual   precision 0.40   signal/null  1.1

The binding constraint is the **local** registration residual, and it is sharp: everything
works at <=2 um and dies by 5 um. Global registration does not deliver that — the cohort gate
is islet RMSE <= 150 um. But it does not have to: matched islets (S6) already establish
correspondence, so this module refines rigidly *inside each islet* over a ~100-200 um
neighbourhood, where 2 um is achievable, and verifies it before matching anything.

Three guards, because the failure mode is a confident wrong pairing:

1. **Refinement is type-blind.** ICP runs on every cell in the islet, then matching runs
   within type. Refining on the same typed cells it then matches would manufacture the
   agreement it claims to find.
2. **The local residual is gated**, not merely reported. Above ``cellmatch_max_residual_um``
   the islet is skipped.
3. **Every islet gets a permuted null.** Residual is checkable but optimistic (ICP fits the
   points it measures); the null is not, and it is the gate that decides.

Scope, stated once and enforced in code: **endocrine cells inside matched islets, with a
verified local residual**. Everything else stays at the islet, niche and donor levels, which
need no cell correspondence at all.
"""

from __future__ import annotations

import argparse
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, read_cell_table
from .manifest import load_manifest, paired_rows
from .register import icp, registration_dir
from .sameslide import match_by_centroid
from .transform import Transform

#: Cell types this module will match. Endocrine only, and deliberately not configurable to
#: "all": the simulation above shows immune cells are never recoverable and unconstrained
#: matching is not a measurement. Widening this needs new evidence, not a config edit.
MATCHABLE_TYPES: tuple[str, ...] = ("Beta", "Alpha", "Delta")

#: Columns holding the fine type on each side. PhenoCycler derives `endocrine_subtype` from
#: the hormone argmax; Xenium carries `celltype_fine`. Post-Xenium IF for INS/GCG lands in
#: `endocrine_subtype` on the Xenium side too when available, which is the stronger source —
#: it is protein on both sides, so the two calls are made the same way.
TYPE_COLS: tuple[str, ...] = ("endocrine_subtype", "celltype_fine")


def fine_type(df: pd.DataFrame) -> pd.Series:
    """The finest endocrine call available, preferring a protein-derived one.

    Order matters. ``endocrine_subtype`` is the hormone argmax — protein on the PhenoCycler
    side, and also protein on the Xenium side when post-Xenium INS/GCG staining is present.
    Preferring it means both sides are typed by the same kind of evidence; falling back to
    ``celltype_fine`` means the Xenium side is typed by surrogate RNA panels instead, which
    is weaker but not wrong.
    """
    out = pd.Series([""] * len(df), index=df.index, dtype="object")
    for col in TYPE_COLS:
        if col in df.columns:
            vals = df[col].astype(str).str.strip()
            out = out.where(out != "", vals)
    return out.where(out.isin(MATCHABLE_TYPES), "")


def refine_local(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    *,
    max_shift_um: float = 25.0,
    max_dist_um: float = 15.0,
) -> tuple[Transform, float, int]:
    """Rigid refinement inside one islet. Returns ``(transform, residual_um, n_inliers)``.

    Seeded at identity because the caller has already applied the global transform, so the
    two point sets start roughly aligned and only a small correction is wanted. The
    correction is capped at ``max_shift_um``: a local fit that wants to move an islet 100 um
    is not refining the alignment, it is matching the wrong islet, and the cap turns that
    into a visible skip rather than a plausible-looking result.
    """
    if len(moving_xy) < 4 or len(fixed_xy) < 4:
        return Transform(), float("nan"), 0

    tf = icp(moving_xy, fixed_xy, max_dist_um=max_dist_um, similarity_only=True)
    moved = tf.apply(moving_xy)
    # How far the refinement actually moved the cells, not the raw translation term: an
    # affine's translation is measured from the origin, which for a section a centimetre
    # away is enormous even when the cells barely move.
    shift = float(np.median(np.linalg.norm(moved - moving_xy, axis=1)))
    if shift > max_shift_um:
        return Transform(), float("nan"), -1          # -1 flags "refinement rejected"

    from scipy.spatial import cKDTree

    d, j = cKDTree(fixed_xy).query(moved, k=1)
    _d2, i2 = cKDTree(moved).query(fixed_xy, k=1)
    mutual = i2[j] == np.arange(len(moved))
    keep = mutual & (d <= max_dist_um)
    if keep.sum() < 3:
        return tf, float("nan"), 0

    # Residual over the CLOSEST HALF of mutual pairs, not all of them. Only ~28% of cells are
    # detectable in both sections, so most mutual pairs are two different cells that happen to
    # be nearest — they measure inter-cell spacing (~5-10 um), not alignment error, and they
    # dominate a plain RMSE. Measured on synthetic pairs: with a good 2.5 um alignment, plain
    # RMSE reads 0.96 at 100% shared but 2.69 at 35% shared, so a 2 um gate would reject a
    # well-aligned islet for having a realistic overlap fraction. The trimmed estimator reads
    # 0.51 / 1.01 for the same two cases and 2.58 / 2.22 for a badly aligned one — it
    # separates alignment error from population mismatch, which is what the gate is asking.
    dd = np.sort(d[keep])
    inner = dd[: max(len(dd) // 2, 3)]
    return tf, float(np.sqrt(np.mean(inner ** 2))), int(keep.sum())


def match_typed(
    moving: pd.DataFrame,
    fixed: pd.DataFrame,
    *,
    max_dist_um: float = 3.0,
    n_perm: int = 50,
    seed: int = 0,
) -> tuple[pd.DataFrame, dict]:
    """Mutual-nearest matching *within* each fine type, plus a permuted null.

    The type constraint is what makes this a measurement rather than a coin flip: it cuts the
    density of candidate partners several-fold, so the chance-match rate falls while genuine
    matches are unaffected.
    """
    rng = np.random.default_rng(seed)
    out, n_null = [], 0.0

    m_type, f_type = fine_type(moving), fine_type(fixed)
    for t in MATCHABLE_TYPES:
        m = moving[m_type == t]
        f = fixed[f_type == t]
        if len(m) < 2 or len(f) < 2:
            continue
        M = m.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        F = f.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        pairs = match_by_centroid(M, F, max_dist_um=max_dist_um)
        if len(pairs):
            sub = pd.DataFrame({
                "cell_id_moving": m["cell_id"].to_numpy()[pairs["moving_idx"].to_numpy(int)],
                "cell_id_fixed": f["cell_id"].to_numpy()[pairs["fixed_idx"].to_numpy(int)],
                "celltype": t,
                "distance_um": pairs["distance_um"].to_numpy(float),
            })
            out.append(sub)

        # Null: same counts and same type constraint, positions shuffled in the islet bbox.
        lo, hi = F.min(axis=0), F.max(axis=0)
        if np.all(hi > lo):
            counts = [len(match_by_centroid(
                M, np.c_[rng.uniform(lo[0], hi[0], len(F)), rng.uniform(lo[1], hi[1], len(F))],
                max_dist_um=max_dist_um)) for _ in range(n_perm)]
            n_null += float(np.mean(counts))

    pairs = pd.concat(out, ignore_index=True) if out else pd.DataFrame(
        columns=["cell_id_moving", "cell_id_fixed", "celltype", "distance_um"])
    stats = {
        "n_matched": len(pairs),
        "null_mean": n_null,
        "signal_over_null": len(pairs) / n_null if n_null > 0 else float("inf")
        if len(pairs) else float("nan"),
        "median_distance_um": float(pairs["distance_um"].median()) if len(pairs) else float("nan"),
    }
    return pairs, stats


def match_islet_cells(
    moving: CellTable,
    fixed: CellTable,
    moving_id: str,
    fixed_id: str,
    tf: Transform,
    *,
    max_dist_um: float = 3.0,
    max_residual_um: float = 2.0,
    min_signal_over_null: float = 3.0,
) -> tuple[pd.DataFrame, dict]:
    """Match cells inside one matched islet pair. Returns ``(pairs, stats)``.

    ``stats['skipped']`` is set with a reason whenever the islet does not clear a guard, so a
    run reports *why* coverage is what it is rather than silently producing fewer pairs.
    """
    m = moving.df[moving.df["structure_id"].astype(str) == str(moving_id)]
    f = fixed.df[fixed.df["structure_id"].astype(str) == str(fixed_id)]
    base = {"moving_id": moving_id, "fixed_id": fixed_id, "n_moving": len(m), "n_fixed": len(f)}
    empty = pd.DataFrame(columns=["cell_id_moving", "cell_id_fixed", "celltype", "distance_um"])
    if len(m) < 5 or len(f) < 5:
        return empty, {**base, "skipped": True, "reason": "too few cells in the islet"}

    # Global transform first, then a type-BLIND local refinement over every cell in the islet.
    M = tf.apply(m.loc[:, ["x_um", "y_um"]].to_numpy(np.float64))
    F = f.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    local, residual, n_in = refine_local(M, F)
    if n_in == -1:
        return empty, {**base, "skipped": True, "reason": "local refinement exceeded the shift cap"}
    if not np.isfinite(residual):
        return empty, {**base, "skipped": True, "reason": "local refinement did not converge"}
    if residual > max_residual_um:
        return empty, {**base, "skipped": True, "residual_um": residual,
                       "reason": f"local residual {residual:.2f} um > {max_residual_um:.1f} um "
                                 f"— below this, matching is not better than chance"}

    m_ref = m.copy()
    m_ref[["x_um", "y_um"]] = local.apply(M)
    pairs, stats = match_typed(m_ref, f, max_dist_um=max_dist_um)
    stats.update(base)
    stats["residual_um"] = residual
    stats["n_refine_inliers"] = n_in

    if stats["n_matched"] and stats["signal_over_null"] < min_signal_over_null:
        return empty, {**stats, "skipped": True,
                       "reason": f"signal/null {stats['signal_over_null']:.1f} < "
                                 f"{min_signal_over_null:.1f} — not distinguishable from chance"}
    stats["skipped"] = False
    stats["reason"] = ""
    return pairs, stats


def match_donor_cells(cfg: PipelineConfig, donor: str, roi: str) -> tuple[pd.DataFrame, dict]:
    """Every matched islet for one section. Requires S6 (`paired/islets.parquet`)."""
    islets_path = cfg.paired_dir / "islets.parquet"
    if not islets_path.exists():
        return pd.DataFrame(), {"donor_id": donor, "roi": roi, "error":
                                "no matched islets — run the 'match' stage (S6) first"}
    islets = pd.read_parquet(islets_path)
    islets = islets[(islets["donor_id"].astype(str) == str(donor))
                    & (islets["roi"].astype(str) == str(roi))]
    if islets.empty:
        return pd.DataFrame(), {"donor_id": donor, "roi": roi,
                                "error": "no matched islets for this section"}

    reg = registration_dir(cfg, donor, roi)
    if not (reg / "transform.json").exists():
        return pd.DataFrame(), {"donor_id": donor, "roi": roi, "error": "no registration"}
    tf = Transform.load(reg)

    fixed_is_pheno = cfg.fixed_modality == "phenocycler"
    fixed_mod = "phenocycler" if fixed_is_pheno else "xenium"
    moving_mod = "xenium" if fixed_is_pheno else "phenocycler"
    fixed = read_cell_table(cfg.cells_pheno_dir if fixed_is_pheno else cfg.cells_xen_dir,
                            donor, roi, modality=fixed_mod)
    moving = read_cell_table(cfg.cells_xen_dir if fixed_is_pheno else cfg.cells_pheno_dir,
                             donor, roi, modality=moving_mod)

    all_pairs, rows = [], []
    for _, r in islets.iterrows():
        pairs, stats = match_islet_cells(
            moving, fixed, r["moving_id"], r["fixed_id"], tf,
            max_dist_um=cfg.cellmatch_max_dist_um,
            max_residual_um=cfg.cellmatch_max_residual_um,
            min_signal_over_null=cfg.cellmatch_min_signal_over_null)
        stats.update({"donor_id": donor, "roi": roi})
        rows.append(stats)
        if len(pairs):
            pairs.insert(0, "donor_id", donor)
            pairs.insert(1, "roi", roi)
            pairs.insert(2, "islet_pair", f"{r['moving_id']}__{r['fixed_id']}")
            all_pairs.append(pairs)

    per_islet = pd.DataFrame(rows)
    summary = {
        "donor_id": donor, "roi": roi,
        "n_islets": len(per_islet),
        "n_islets_used": int((~per_islet["skipped"].fillna(True)).sum()) if len(per_islet) else 0,
        "n_pairs": int(sum(len(p) for p in all_pairs)),
        "median_residual_um": float(pd.to_numeric(
            per_islet.get("residual_um", pd.Series(dtype=float)), errors="coerce").median())
        if len(per_islet) else float("nan"),
    }
    out = pd.concat(all_pairs, ignore_index=True) if all_pairs else pd.DataFrame()
    return out, {**summary, "per_islet": per_islet}


def run_cellmatch(cfg: PipelineConfig, donors: Optional[list[str]] = None, *,
                  roi: Optional[str] = None, tissue: Optional[str] = None) -> pd.DataFrame:
    """Stage entry point. Endocrine cells inside matched islets, guarded at every step."""
    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[cellmatch] no paired sections", flush=True)
        return pd.DataFrame()

    print(f"[cellmatch] matching {'/'.join(MATCHABLE_TYPES)} cells inside matched islets; "
          f"cap {cfg.cellmatch_max_dist_um} um, local residual must be "
          f"<= {cfg.cellmatch_max_residual_um} um.\n"
          f"[cellmatch] NOTE: immune cells are NOT matchable across 5 um sections at any "
          f"registration quality (a ~8 um cell cannot put enough nucleus on both sides of "
          f"the cut) — use the islet/niche levels for those.", flush=True)

    frames, rows, per_islet = [], [], []
    for _, r in man.iterrows():
        donor, roi_name = r["donor_id"], r["roi"]
        try:
            pairs, stats = match_donor_cells(cfg, donor, roi_name)
        except FileNotFoundError as exc:
            print(f"[cellmatch] {donor}/{roi_name}: {exc}", flush=True)
            continue
        if "error" in stats:
            print(f"[cellmatch] {donor}/{roi_name}: {stats['error']}", flush=True)
            continue
        pi = stats.pop("per_islet", pd.DataFrame())
        if len(pi):
            per_islet.append(pi)
        rows.append(stats)
        if len(pairs):
            frames.append(pairs)
        print(f"[cellmatch] {donor}/{roi_name}: {stats['n_pairs']:,} cell pairs from "
              f"{stats['n_islets_used']}/{stats['n_islets']} islets "
              f"(median local residual {stats['median_residual_um']:.2f} um)", flush=True)

    cfg.paired_dir.mkdir(parents=True, exist_ok=True)
    cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
    if frames:
        allf = pd.concat(frames, ignore_index=True)
        out = cfg.paired_dir / "islet_cells.parquet"
        allf.to_parquet(out, index=False)
        print(f"\n[cellmatch] wrote {out} ({len(allf):,} endocrine cell pairs)", flush=True)
    if per_islet:
        pd.concat(per_islet, ignore_index=True).to_csv(
            cfg.integration_qc_dir / "cellmatch_islets.csv", index=False)

    stats_df = pd.DataFrame(rows)
    if len(stats_df):
        stats_df.to_csv(cfg.integration_qc_dir / "cellmatch_qc.csv", index=False)
    return stats_df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Endocrine cell-to-cell matching inside matched islets (serial sections)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None)
    ap.add_argument("--tissue", default=None)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_cellmatch(cfg, args.donors, roi=args.roi, tissue=args.tissue)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
