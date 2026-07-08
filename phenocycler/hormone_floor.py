#!/usr/bin/env python3
"""
Step 4 — hormone-strength floor on the RESTORE positivity calls (the false-endocrine fix).

Faithful port of ``scripts/senior/apply_hormone_floor.py`` (Islet-Explorer-Senior), made
config-driven + per-donor parallel and given an in-process ``run_hormone_floor(cfg)`` entry
point so the orchestrator can run it as a stage (which the upstream notebook OMITTED).

RESTORE's per-image threshold lands IN THE NOISE for β-loss (T1D) donors (no separated bright INS
population), so it calls thousands of acinar cells INS⁺ at ``_norm ~1.1`` (barely over threshold),
100% in tissue, no islet coherence — while real β-cells sit at ``_norm`` 8–25. Requiring
``_norm >= min_norm`` (threshold-relative, so per-image-adaptive) rejects the noise-floor false
positives while keeping the real, separated positives. Validated per-donor: 6380 (visually β-null)
INS⁺ 11,905 → ~16 at K=5, while 6539 (ND) and 6593 (β-retained T1D) keep their real β. See
``reassess_diag`` for the yardstick.

Rewrites ONLY the ``{m}_pos`` boolean of the listed markers (default the three endocrine hormones);
``_norm`` / ``_log2r`` and every other column are copied through unchanged. Idempotent: re-flooring an
already-floored dir at the same K is a no-op (``_norm`` is untouched).

    python -m phenocycler.hormone_floor --min-norm 5                 # in-place on restore_gated_redsea
    python -m phenocycler.hormone_floor --gated-dir <src> --out-dir <dst> --min-norm 5
"""

from __future__ import annotations

import argparse
import functools
import glob
from pathlib import Path

import pandas as pd

from .config import PipelineConfig, load_config, HORMONE_MARKERS
from .parallel import map_donors


def floor_donor_file(gated_f: str, out_f: str, markers, min_norm: float) -> tuple[int, list]:
    """Rewrite ``{m}_pos = (_norm >= min_norm)`` for `markers`; pass everything else through."""
    df = pd.read_parquet(gated_f)
    report = []
    for m in markers:
        pc, nc = f"{m}_pos", f"{m}_norm"
        if pc in df.columns and nc in df.columns:
            before = int(df[pc].sum())
            df[pc] = (df[nc].to_numpy() >= min_norm)
            report.append(f"{m} {before:,}->{int(df[pc].sum()):,}")
    Path(out_f).parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_f, index=False)
    return len(df), report


def _floor_one_donor(donor: str, gated_dir: str, out_dir: str, markers, min_norm: float):
    files = sorted(glob.glob(str(Path(gated_dir) / f"donor_id={donor}" / "*.parquet")))
    if not files:
        return None
    out_f = Path(out_dir) / f"donor_id={donor}" / "data_0.parquet"
    n, report = floor_donor_file(files[0], str(out_f), markers, min_norm)
    print(f"[{donor}] {n:,} cells | " + " | ".join(report), flush=True)
    return {"donor": donor, "n": n, "report": report}


def run_hormone_floor(cfg: PipelineConfig, *, gated_dir=None, out_dir=None, markers=None,
                      min_norm=None, donors=None, n_jobs=None) -> list:
    """Floor the endocrine hormone positivity calls for every donor (parallel).

    Defaults floor ``restore_gated_redsea`` in place; the pipeline stage passes the pre-floor backup
    as the source so the operation is reproducible from a clean RESTORE.
    """
    gated_dir = Path(gated_dir) if gated_dir else cfg.restore_gated_dir
    out_dir = Path(out_dir) if out_dir else cfg.restore_gated_dir
    markers = list(markers) if markers else list(HORMONE_MARKERS)
    min_norm = cfg.hormone_min_norm if min_norm is None else float(min_norm)
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    donor_ids = donors or cfg.discover_donors(gated_dir)
    if not donor_ids:
        raise SystemExit(f"[err] no gated donors under {gated_dir}")

    fn = functools.partial(_floor_one_donor, gated_dir=str(gated_dir), out_dir=str(out_dir),
                           markers=markers, min_norm=min_norm)
    results = [r for r in map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True) if r]
    print(f"[done] floored {markers} at _norm>={min_norm} -> {out_dir}")
    return results


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--gated-dir", type=Path, default=None, help="source gated dir (default restore_gated_redsea)")
    ap.add_argument("--out-dir", type=Path, default=None, help="destination (default: in place)")
    ap.add_argument("--min-norm", type=float, default=None,
                    help="require {m}_norm >= this for {m}_pos to stay True (default cfg.hormone_min_norm)")
    ap.add_argument("--markers", default=None,
                    help="comma-separated markers to floor (default the endocrine hormones INS,GCG,SST)")
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
