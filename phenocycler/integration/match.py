"""
S6 — islet-to-islet matching: the primary pairing unit for serial sections.

Cells cannot be paired across serial sections. Islets can. An islet is roughly 50-200 um
across, so it persists through a 5 um section gap; there are tens to hundreds per section,
enough for statistics; and it is the unit T1D biology is actually about — beta-cell mass and
insulitis are per-islet properties.

Matching is a one-to-one assignment problem, solved with the Hungarian algorithm rather than
greedy nearest-neighbour. Greedy matching is order-dependent and produces conflicts (two
moving islets claiming one fixed islet); Hungarian gives the globally optimal assignment
under the cost function and is deterministic.

The cost combines three pieces of evidence, because position alone is not enough after an
imperfect registration:

* **distance** between registered centroids — the primary term;
* **area ratio** — an islet does not change size much between adjacent sections;
* **shape** (eccentricity) — a long thin islet stays long and thin.

Gates are applied before the assignment, not after: pairs beyond ``match_max_dist_um`` or
outside ``match_area_ratio`` are made infinite-cost so the optimiser cannot be forced into a
bad pair simply because everything else was taken.

Not every islet should match, and that is the point. Between two sections some islets appear,
some vanish, and small ones do both readily. A match fraction near 1.0 on real data would be
suspicious rather than reassuring — it would suggest the gates are too loose.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .manifest import PAIRED, load_manifest
from .register import registration_dir
from .structures import structures_path
from .transform import Transform

#: Cost weights. Distance dominates; the shape terms break ties and reject implausible pairs
#: that happen to be close. Tuned so a 100 um offset costs about the same as a 1.5x area
#: mismatch — i.e. geometry is trusted more than morphometry, but not exclusively.
W_DISTANCE = 1.0
W_AREA = 0.35
W_SHAPE = 0.15


def _pair_features(moving: pd.DataFrame, fixed: pd.DataFrame,
                   max_dist_um: float, area_ratio: float
                   ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Distance, area-ratio, shape-difference and feasibility matrices (n_moving x n_fixed)."""
    mv = moving.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    fx = fixed.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    dist = np.linalg.norm(mv[:, None, :] - fx[None, :, :], axis=2)

    def _col(df, name, default):
        return (df[name].to_numpy(np.float64) if name in df.columns
                else np.full(len(df), default, dtype=np.float64))

    m_area = np.maximum(_col(moving, "area_um2", 0.0), 1e-9)
    f_area = np.maximum(_col(fixed, "area_um2", 0.0), 1e-9)
    ratio = m_area[:, None] / f_area[None, :]
    # Symmetric log ratio: a 2x mismatch costs the same whichever way round it is.
    log_ratio = np.abs(np.log(ratio))

    m_ecc = _col(moving, "eccentricity", 0.0)
    f_ecc = _col(fixed, "eccentricity", 0.0)
    shape = np.abs(m_ecc[:, None] - f_ecc[None, :])

    feasible = (dist <= max_dist_um) & (ratio >= 1.0 / area_ratio) & (ratio <= area_ratio)
    return dist, log_ratio, shape, feasible


def match_structures(
    moving: pd.DataFrame,
    fixed: pd.DataFrame,
    *,
    max_dist_um: float = 200.0,
    area_ratio: float = 2.5,
    require_mutual: bool = True,
) -> pd.DataFrame:
    """One-to-one assignment between two structure tables (already in a common frame).

    ``moving`` must already carry registered coordinates in ``x_um`` / ``y_um`` — this
    function does no transforming, so it is equally usable for same-slide mode and for
    validating a registration.
    """
    empty = pd.DataFrame(columns=[
        "pair_id", "moving_id", "fixed_id", "distance_um", "area_ratio",
        "shape_diff", "cost", "mutual_nn",
    ])
    if moving is None or fixed is None or moving.empty or fixed.empty:
        return empty

    dist, log_ratio, shape, feasible = _pair_features(moving, fixed, max_dist_um, area_ratio)
    if not feasible.any():
        return empty

    cost = (W_DISTANCE * (dist / max(max_dist_um, 1e-9))
            + W_AREA * log_ratio
            + W_SHAPE * shape)
    big = float(cost[feasible].max() * 10 + 1e3) if feasible.any() else 1e6
    cost_masked = np.where(feasible, cost, big)

    from scipy.optimize import linear_sum_assignment

    rows, cols = linear_sum_assignment(cost_masked)

    # Mutual-nearest check, computed on raw distance rather than the composite cost. The
    # assignment is globally optimal, which is not the same as every individual pair being
    # each other's best choice; on serial sections where islets appear and vanish, the
    # difference is exactly the pairs most likely to be spurious.
    nn_fixed_for_moving = np.argmin(np.where(feasible, dist, np.inf), axis=1)
    nn_moving_for_fixed = np.argmin(np.where(feasible, dist, np.inf), axis=0)

    mv_ids = moving["structure_id"].to_numpy() if "structure_id" in moving.columns \
        else np.arange(len(moving)).astype(str)
    fx_ids = fixed["structure_id"].to_numpy() if "structure_id" in fixed.columns \
        else np.arange(len(fixed)).astype(str)

    out = []
    for i, j in zip(rows, cols):
        if not feasible[i, j]:
            continue
        mutual = bool(nn_fixed_for_moving[i] == j and nn_moving_for_fixed[j] == i)
        if require_mutual and not mutual:
            continue
        out.append({
            "pair_id": f"{mv_ids[i]}__{fx_ids[j]}",
            "moving_id": mv_ids[i],
            "fixed_id": fx_ids[j],
            "distance_um": float(dist[i, j]),
            "area_ratio": float(np.exp(log_ratio[i, j])),
            "shape_diff": float(shape[i, j]),
            "cost": float(cost[i, j]),
            "mutual_nn": mutual,
        })
    return pd.DataFrame(out) if out else empty


def permuted_null(
    moving: pd.DataFrame,
    fixed: pd.DataFrame,
    *,
    max_dist_um: float = 200.0,
    area_ratio: float = 2.5,
    n_perm: int = 200,
    seed: int = 0,
) -> dict:
    """How many matches would survive if the registration were meaningless?

    Shuffles the moving centroids within their own bounding box and re-matches. If the real
    match count is not clearly above this null, the pairing is not evidence of anything — a
    dense enough islet field will always produce *some* matches within 200 um by chance.

    Mirrors the rotation-null approach already used for insulitis statistics on the Xenium
    side, which is the right instinct: a spatial result needs a spatial null.
    """
    rng = np.random.default_rng(seed)
    if moving.empty or fixed.empty:
        return {"n_real": 0, "null_mean": 0.0, "null_p95": 0.0, "z": float("nan")}

    real = len(match_structures(moving, fixed, max_dist_um=max_dist_um,
                                area_ratio=area_ratio))
    xy = moving.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    lo, hi = xy.min(axis=0), xy.max(axis=0)

    counts = np.zeros(n_perm, dtype=int)
    for k in range(n_perm):
        shuffled = moving.copy()
        shuffled["x_um"] = rng.uniform(lo[0], hi[0], len(moving))
        shuffled["y_um"] = rng.uniform(lo[1], hi[1], len(moving))
        counts[k] = len(match_structures(shuffled, fixed, max_dist_um=max_dist_um,
                                         area_ratio=area_ratio))

    mean, sd = float(counts.mean()), float(counts.std())
    return {
        "n_real": int(real),
        "null_mean": mean,
        "null_p95": float(np.percentile(counts, 95)),
        "z": float((real - mean) / sd) if sd > 0 else float("inf") if real > mean else 0.0,
    }


def _annotate(pairs: pd.DataFrame, moving: pd.DataFrame, fixed: pd.DataFrame,
              moving_mod: str, fixed_mod: str) -> pd.DataFrame:
    """Attach both modalities' per-islet features to each matched pair."""
    if pairs.empty:
        return pairs
    keep = ["structure_id", "n_cells", "area_um2", "size_class", "eccentricity",
            "n_beta", "n_alpha", "n_delta", "frac_beta", "frac_alpha", "frac_delta",
            "n_immune_zone", "immune_per_100_endo", "insulitis_grade"]

    def _side(df, prefix, id_col):
        cols = [c for c in keep if c in df.columns]
        sub = df.loc[:, cols].rename(columns={c: f"{prefix}_{c}" for c in cols
                                              if c != "structure_id"})
        return sub.rename(columns={"structure_id": id_col})

    out = pairs.merge(_side(moving, moving_mod, "moving_id"), on="moving_id", how="left")
    out = out.merge(_side(fixed, fixed_mod, "fixed_id"), on="fixed_id", how="left")
    return out


def match_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    kind: str = "islet",
    with_null: bool = True,
) -> tuple[pd.DataFrame, dict]:
    """Match one donor's structures after applying the stored registration."""
    fixed_is_pheno = cfg.fixed_modality == "phenocycler"
    fixed_mod = "phenocycler" if fixed_is_pheno else "xenium"
    moving_mod = "xenium" if fixed_is_pheno else "phenocycler"

    fx_path = structures_path(cfg, fixed_mod, donor, roi, kind)
    mv_path = structures_path(cfg, moving_mod, donor, roi, kind)
    if not fx_path.exists() or not mv_path.exists():
        return pd.DataFrame(), {"error": f"structures missing for donor {donor}"}

    fixed = pd.read_parquet(fx_path)
    moving = pd.read_parquet(mv_path)

    reg_dir = registration_dir(cfg, donor, roi)
    if not (reg_dir / "transform.json").exists():
        return pd.DataFrame(), {"error": f"no registration for donor {donor}; run S4 first"}
    tf = Transform.load(reg_dir)

    moving_reg = moving.copy()
    xy = tf.apply(moving.loc[:, ["x_um", "y_um"]].to_numpy(np.float64))
    moving_reg["x_um"] = xy[:, 0]
    moving_reg["y_um"] = xy[:, 1]

    pairs = match_structures(moving_reg, fixed,
                             max_dist_um=cfg.match_max_dist_um,
                             area_ratio=cfg.match_area_ratio)
    pairs = _annotate(pairs, moving, fixed, moving_mod, fixed_mod)
    if not pairs.empty:
        pairs.insert(0, "donor_id", donor)
        pairs.insert(1, "roi", roi)
        pairs["moving_modality"] = moving_mod
        pairs["fixed_modality"] = fixed_mod

    stats = {
        "donor_id": donor,
        "roi": roi,
        "n_moving": len(moving),
        "n_fixed": len(fixed),
        "n_matched": len(pairs),
        "match_fraction": len(pairs) / max(min(len(moving), len(fixed)), 1),
        "median_distance_um": float(pairs["distance_um"].median()) if len(pairs) else float("nan"),
    }
    if with_null and len(moving) and len(fixed):
        stats.update(permuted_null(moving_reg, fixed,
                                   max_dist_um=cfg.match_max_dist_um,
                                   area_ratio=cfg.match_area_ratio, n_perm=100))
    return pairs, stats


def concordance(pairs: pd.DataFrame, moving_mod: str, fixed_mod: str) -> dict:
    """Cross-modal agreement on matched islets — the scientific acceptance test.

    If registration and matching worked, the *same islet* measured two ways should agree:
    protein INS-positive fraction should track RNA Beta fraction, and an islet graded
    insulitic by one modality should be graded insulitic by the other. This is the number
    that says whether any of the preceding machinery produced something real.
    """
    out: dict = {"n_pairs": len(pairs)}
    if len(pairs) < 3:
        return out

    from scipy import stats as sps

    for feature in ("frac_beta", "frac_alpha", "frac_delta", "immune_per_100_endo", "n_cells"):
        a, b = f"{moving_mod}_{feature}", f"{fixed_mod}_{feature}"
        if a not in pairs.columns or b not in pairs.columns:
            continue
        x = pd.to_numeric(pairs[a], errors="coerce")
        y = pd.to_numeric(pairs[b], errors="coerce")
        ok = x.notna() & y.notna()
        if ok.sum() < 3:
            continue
        rho, p = sps.spearmanr(x[ok], y[ok])
        out[f"spearman_{feature}"] = float(rho)
        out[f"p_{feature}"] = float(p)

    a, b = f"{moving_mod}_insulitis_grade", f"{fixed_mod}_insulitis_grade"
    if a in pairs.columns and b in pairs.columns:
        ok = pairs[a].notna() & pairs[b].notna()
        if ok.sum() >= 3:
            agree = float((pairs.loc[ok, a] == pairs.loc[ok, b]).mean())
            out["insulitis_agreement"] = agree
            try:
                from sklearn.metrics import cohen_kappa_score
                out["insulitis_kappa"] = float(cohen_kappa_score(
                    pairs.loc[ok, a].astype(str), pairs.loc[ok, b].astype(str)))
            except Exception:  # noqa: BLE001
                pass
    return out


def run_match(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
) -> pd.DataFrame:
    """Stage entry point: match every registered donor and write the paired islet table."""
    man = load_manifest(cfg)
    man = man[(man["pair_status"] == PAIRED) & (man["roi"] == roi)]
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[match] no paired donors", flush=True)
        return pd.DataFrame()

    fixed_mod = cfg.fixed_modality
    moving_mod = "xenium" if fixed_mod == "phenocycler" else "phenocycler"

    all_pairs, stats_rows = [], []
    for _, r in man.iterrows():
        donor = r["donor_id"]
        pairs, stats = match_donor(cfg, donor, roi)
        if "error" in stats:
            print(f"[match] {donor}: {stats['error']}", flush=True)
            continue
        stats.update(concordance(pairs, moving_mod, fixed_mod))
        all_pairs.append(pairs)
        stats_rows.append(stats)
        print(f"[match] {donor}: {stats['n_matched']}/{min(stats['n_moving'], stats['n_fixed'])} "
              f"islets matched (median {stats['median_distance_um']:.0f} um), "
              f"null mean {stats.get('null_mean', float('nan')):.1f} "
              f"z={stats.get('z', float('nan')):.1f}", flush=True)

    if all_pairs:
        pairs = pd.concat([p for p in all_pairs if len(p)], ignore_index=True) \
            if any(len(p) for p in all_pairs) else pd.DataFrame()
        if len(pairs):
            out = cfg.paired_dir / "islets.parquet"
            out.parent.mkdir(parents=True, exist_ok=True)
            pairs.to_parquet(out, index=False)
            print(f"\n[match] wrote {out} ({len(pairs)} matched islet pairs)", flush=True)

    stats_df = pd.DataFrame(stats_rows)
    if len(stats_df):
        out = cfg.integration_qc_dir / "match_qc.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(out, index=False)
    return stats_df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Match islets across registered serial sections")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_match(cfg, args.donors, roi=args.roi)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
