#!/usr/bin/env python3
"""
Orchestrator — run the full raw-data → broad-lineage pipeline end to end.

Mirrors ``notebooks/senior_phenotyping_redsea_restore.ipynb`` (idempotent
``run_step`` + a final status table), but in-process against the ``phenocycler``
package rather than shelling out to loose scripts.  Each step is skipped when its
outputs already exist unless ``force=True``.

    python -m phenocycler.pipeline                       # run every missing step
    python -m phenocycler.pipeline --only redsea restore
    python -m phenocycler.pipeline --force               # re-run everything
    python -m phenocycler.pipeline --status              # just print the status table
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

from .config import PipelineConfig, load_config


# step name -> (output pattern under data_dir, human label)
STAGES = [
    ("cells", "cells/donor_id=*", "Raw cells (DuckDB)"),
    ("redsea", "cells_redsea/donor_id=*", "REDSEA corrected"),
    ("restore", "restore_gated_redsea/donor_id=*", "RESTORE gated"),
    ("lineage", "phenotype/broad/donor_id=*", "Broad lineage"),
    ("qupath", "phenotype/qupath_class/pheno_class_*.csv", "QuPath CSVs"),
]


def _has_outputs(cfg: PipelineConfig, pattern: str) -> int:
    return len(glob.glob(str(cfg.data_dir / pattern)))


def run_step(cfg: PipelineConfig, name: str, func, pattern: str, force: bool) -> None:
    n = _has_outputs(cfg, pattern)
    if n and not force:
        print(f"[SKIP] {name}: {n} outputs already at {cfg.data_dir/pattern}")
        return
    print(f"\n=== {name} ===")
    func()


def run_pipeline(cfg: PipelineConfig, *, only=None, force=False) -> None:
    from . import cells_parquet, redsea, restore, lineage, qupath_export

    def _cells():
        cells_parquet.build_cells_parquet(cfg)

    def _redsea():
        params = redsea.RedseaParams.from_config(cfg)
        redsea.run_redsea(cfg, cfg.discover_donors(), params, n_jobs=cfg.n_jobs)

    def _restore():
        restore.run_restore(cfg)

    def _lineage():
        lineage.run_lineage(cfg)

    def _qupath():
        qupath_export.export_qupath_classes(cfg)

    funcs = {"cells": (_cells, "cells/donor_id=*"),
             "redsea": (_redsea, "cells_redsea/donor_id=*"),
             "restore": (_restore, "restore_gated_redsea/donor_id=*"),
             "lineage": (_lineage, "phenotype/broad/donor_id=*"),
             "qupath": (_qupath, "phenotype/qupath_class/pheno_class_*.csv")}
    order = ["cells", "redsea", "restore", "lineage", "qupath"]
    selected = only or order
    for name in order:
        if name not in selected:
            continue
        fn, pat = funcs[name]
        run_step(cfg, name, fn, pat, force)
    print_status(cfg)


def print_status(cfg: PipelineConfig) -> None:
    print("\n" + "=" * 52)
    print(f"{'stage':<22}{'outputs':>10}   status")
    print("-" * 52)
    for _, pattern, label in STAGES:
        n = _has_outputs(cfg, pattern)
        print(f"{label:<22}{n:>10}   {'OK' if n else 'MISSING'}")
    print("=" * 52)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--only", nargs="*", default=None,
                    choices=["cells", "redsea", "restore", "lineage", "qupath"])
    ap.add_argument("--jobs", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--status", action="store_true", help="print status and exit")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    if a.status:
        print_status(cfg)
        return 0
    run_pipeline(cfg, only=a.only, force=a.force)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
