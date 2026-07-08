#!/usr/bin/env python3
"""
Step 4 — positivity `_norm` floor on the RESTORE calls (false-endocrine + false-immune fix).

Faithful port/extension of ``scripts/senior/apply_hormone_floor.py`` (Islet-Explorer-Senior), made
config-driven + per-donor parallel and given an in-process entry point so the orchestrator runs it as
a stage (which the upstream notebook OMITTED).

RESTORE's per-image mean+3σ threshold lands IN THE NOISE for sparse markers when there is no separated
bright population, producing a marginal-positive shoulder at ``_norm ~1`` that is background, not signal:
  - hormones {INS,GCG,SST} for β-loss (T1D) donors — thousands of acinar cells called INS⁺ at ~1.1
    (100% in tissue, no islet coherence) while real β sit at _norm 8–25.  Floor at K=5.
  - immune {CD3e,CD20,CD163} — e.g. donor 6476's 40% "CD3e⁺" (median _norm 1.24, diffuse across the
    acinar parenchyma, acinar-sized) was a threshold over-call, while 6539/6579's real T cells sit at
    _norm ≥ 2–6.  Floor at K=2 (validated to preserve 6579's real pancreatitis T-cell mode + the
    ND→T1D T-gain).  MPO (extra dir) is gated at the same K in lineage.py.

Requiring ``_norm >= K`` (threshold-relative, per-image-adaptive) rejects the noise-floor false
positives while keeping the real, separated positives.  Rewrites ONLY the ``{m}_pos`` boolean;
``_norm`` / ``_log2r`` and every other column pass through unchanged.  Idempotent (``_norm`` untouched).

    python -m phenocycler.hormone_floor                                    # full floor (hormones + immune)
    python -m phenocycler.hormone_floor --markers INS,GCG,SST --min-norm 5 # just those markers, that K
"""

from __future__ import annotations

import argparse
import functools
import glob
from pathlib import Path

import pandas as pd

from .config import PipelineConfig, load_config, HORMONE_MARKERS, IMMUNE_FLOOR_MARKERS
from .parallel import map_donors


def floor_donor_file(gated_f: str, out_f: str, floors: dict) -> tuple[int, list]:
    """Rewrite ``{m}_pos = (_norm >= K)`` for each ``{m: K}`` in `floors`; pass everything else through."""
    df = pd.read_parquet(gated_f)
    report = []
    for m, K in floors.items():
        pc, nc = f"{m}_pos", f"{m}_norm"
        if pc in df.columns and nc in df.columns:
            before = int(df[pc].sum())
            df[pc] = (df[nc].to_numpy() >= K)
            report.append(f"{m}>={K:g} {before:,}->{int(df[pc].sum()):,}")
    Path(out_f).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_f, index=False)
    return len(df), report


def _floor_one_donor(donor: str, gated_dir: str, out_dir: str, floors: dict):
    files = sorted(glob.glob(str(Path(gated_dir) / f"donor_id={donor}" / "*.parquet")))
    if not files:
        return None
    out_f = Path(out_dir) / f"donor_id={donor}" / "data_0.parquet"
    n, report = floor_donor_file(files[0], str(out_f), floors)
    print(f"[{donor}] {n:,} cells | " + " | ".join(report), flush=True)
    return {"donor": donor, "n": n, "report": report}


def default_floors(cfg: PipelineConfig) -> dict:
    """The pipeline floor: hormones at hormone_min_norm + main-dir immune markers at immune_min_norm.
    (MPO — extra dir — is gated at immune_min_norm inside lineage.py, not here.)"""
    return {**{m: cfg.hormone_min_norm for m in HORMONE_MARKERS},
            **{m: cfg.immune_min_norm for m in IMMUNE_FLOOR_MARKERS}}


def run_norm_floor(cfg: PipelineConfig, floors: dict, *, gated_dir=None, out_dir=None,
                   donors=None, n_jobs=None) -> list:
    """Apply a per-marker ``{m: K}`` positivity floor to every donor (parallel).

    Defaults floor ``restore_gated_redsea`` in place; the pipeline stage passes the pre-floor backup
    as the source so the operation is reproducible from a clean RESTORE.
    """
    gated_dir = Path(gated_dir) if gated_dir else cfg.restore_gated_dir
    out_dir = Path(out_dir) if out_dir else cfg.restore_gated_dir
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    donor_ids = donors or cfg.discover_sections(gated_dir)
    if not donor_ids:
        raise SystemExit(f"[err] no gated donors under {gated_dir}")

    fn = functools.partial(_floor_one_donor, gated_dir=str(gated_dir), out_dir=str(out_dir),
                           floors=dict(floors))
    results = [r for r in map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True) if r]
    print(f"[done] floored {floors} -> {out_dir}")
    return results


def run_hormone_floor(cfg: PipelineConfig, *, gated_dir=None, out_dir=None, markers=None,
                      min_norm=None, donors=None, n_jobs=None) -> list:
    """Backward-compatible entry point. No markers/min_norm -> the full pipeline floor (hormones@K_h +
    immune@K_i, via ``default_floors``); explicit markers+min_norm -> just those markers at that K."""
    if markers is None and min_norm is None:
        floors = default_floors(cfg)
    else:
        mk = list(markers) if markers else list(HORMONE_MARKERS)
        K = cfg.hormone_min_norm if min_norm is None else float(min_norm)
        floors = {m: K for m in mk}
    return run_norm_floor(cfg, floors, gated_dir=gated_dir, out_dir=out_dir,
                          donors=donors, n_jobs=n_jobs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--gated-dir", type=Path, default=None, help="source gated dir (default restore_gated_redsea)")
    ap.add_argument("--out-dir", type=Path, default=None, help="destination (default: in place)")
    ap.add_argument("--min-norm", type=float, default=None,
                    help="K for the markers in --markers (omit both for the full hormone+immune floor)")
    ap.add_argument("--markers", default=None,
                    help="comma-separated markers to floor at --min-norm (default: the full floor)")
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--donors", nargs="*", default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    markers = [m.strip() for m in a.markers.split(",") if m.strip()] if a.markers else None
    run_hormone_floor(cfg, gated_dir=a.gated_dir, out_dir=a.out_dir, markers=markers,
                      min_norm=a.min_norm, donors=a.donors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
