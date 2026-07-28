#!/usr/bin/env python3
"""
Orchestrator — run the full raw-data → 8-class broad-lineage pipeline end to end.

Idempotent ``run_step`` (skips a stage when its outputs already exist unless
``force=True``) + a final status table, in-process against the ``phenocycler``
package.  Stage order (the two starred stages are the corrections the upstream
notebook / earlier port OMITTED):

    cells → redsea → restore → restore_extra* → hormone_floor* → lineage → qupath → figures

``hormone_floor`` (rewrites {INS,GCG,SST}_pos = _norm ≥ K) MUST run before
``lineage`` or the false-endocrine over-calling returns; ``restore_extra`` gates
B3TUBB/CD99/MPO for the Neural/Endocrine-CD99/Neutrophil classes.

    python -m phenocycler.pipeline                       # run every missing step
    python -m phenocycler.pipeline --only hormone_floor lineage figures
    python -m phenocycler.pipeline --force               # re-run everything
    python -m phenocycler.pipeline --status              # just print the status table
"""

from __future__ import annotations

import argparse
import glob
import shutil
from pathlib import Path

from .config import PipelineConfig, load_config


# step name -> (output pattern under data_dir, human label) for the status table
STAGES = [
    ("cells", "cells/donor_id=*", "Raw cells (DuckDB)"),
    ("redsea", "cells_redsea/donor_id=*", "REDSEA corrected"),
    ("restore", "restore_gated_redsea/donor_id=*", "RESTORE gated (10 markers)"),
    ("restore_extra", "restore_gated_redsea_extra/donor_id=*", "RESTORE gated (extra 3)"),
    ("hormone_floor", "restore_gated_redsea/donor_id=*", "Hormone floor (K=5)"),
    ("lineage", "phenotype/broad/donor_id=*", "Broad lineage (8-class)"),
    ("qupath", "phenotype/qupath_class/pheno_class_*.csv", "QuPath CSVs"),
    ("figures", "phenotype/celltype_marker_dotplot.png", "Identity QC figures"),
]

ORDER = ["cells", "redsea", "restore", "restore_extra", "hormone_floor",
         "lineage", "qupath", "figures"]


def _has_outputs(cfg: PipelineConfig, pattern: str) -> int:
    return len(glob.glob(str(cfg.data_dir / pattern)))


def _lineage_is_8class(cfg: PipelineConfig) -> bool:
    """True only if the broad-lineage partitions exist AND carry the 8-class score columns.

    A stale 6-class ``broad/`` therefore counts as 'not done', so the lineage step re-runs.
    """
    files = sorted(glob.glob(str(cfg.broad_dir / "donor_id=*" / "*.parquet")))
    if not files:
        return False
    try:
        import pyarrow.parquet as pq
        cols = set(pq.ParquetFile(files[0]).schema.names)
    except Exception:
        return False
    return {"score_Neural", "score_Neutrophil"}.issubset(cols)


def _run_step(cfg, name, func, pattern, checker, force) -> None:
    if not force:
        if checker is not None:
            if checker(cfg):
                print(f"[SKIP] {name}: already complete")
                return
        elif pattern is not None:
            n = _has_outputs(cfg, pattern)
            if n:
                print(f"[SKIP] {name}: {n} outputs already at {cfg.data_dir/pattern}")
                return
    print(f"\n=== {name} ===")
    func()


def run_pipeline(cfg: PipelineConfig, *, only=None, force=False) -> None:
    from . import (cells_parquet, redsea, restore, hormone_floor, lineage,
                   qupath_export, figures)

    def _cells():
        cells_parquet.build_cells_parquet(cfg)

    def _redsea():
        params = redsea.RedseaParams.from_config(cfg)
        redsea.run_redsea(cfg, cfg.discover_donors(), params, n_jobs=cfg.n_jobs)

    def _restore():
        restore.run_restore(cfg)

    def _restore_extra():
        restore.run_restore_extra(cfg)

    def _hormone_floor():
        # Guarantee a rollback point: if no pre-floor backup exists yet, snapshot the (un-floored)
        # gated dir to restore_gated_redsea.pre_hormonefloor BEFORE flooring it in place. Then always
        # floor FROM that backup so the result is reproducible from a clean RESTORE (idempotent — the
        # floor only rewrites {INS,GCG,SST}_pos; _norm is untouched).
        prefloor = cfg.restore_gated_prefloor_dir
        if not prefloor.exists() and cfg.restore_gated_dir.exists():
            shutil.copytree(cfg.restore_gated_dir, prefloor)
            print(f"[hormone_floor] backed up un-floored gates -> {prefloor}")
        src = prefloor if prefloor.exists() else cfg.restore_gated_dir
        hormone_floor.run_hormone_floor(cfg, gated_dir=src, out_dir=cfg.restore_gated_dir)

    def _lineage():
        lineage.run_lineage(cfg)

    def _qupath():
        qupath_export.export_qupath_classes(cfg)

    def _figures():
        figures.run_figures(cfg)

    # (func, existence-pattern, column-checker). hormone_floor is idempotent -> never skipped on
    # existence (pattern None, checker None) so it always re-applies when selected.
    steps = {
        "cells": (_cells, "cells/donor_id=*", None),
        "redsea": (_redsea, "cells_redsea/donor_id=*", None),
        "restore": (_restore, "restore_gated_redsea/donor_id=*", None),
        "restore_extra": (_restore_extra, "restore_gated_redsea_extra/donor_id=*", None),
        "hormone_floor": (_hormone_floor, None, None),
        "lineage": (_lineage, "phenotype/broad/donor_id=*", _lineage_is_8class),
        "qupath": (_qupath, "phenotype/qupath_class/pheno_class_*.csv", None),
        "figures": (_figures, "phenotype/celltype_marker_dotplot.png", None),
    }
    selected = only or ORDER
    for name in ORDER:
        if name not in selected:
            continue
        fn, pat, checker = steps[name]
        _run_step(cfg, name, fn, pat, checker, force)
    print_status(cfg)


def print_status(cfg: PipelineConfig) -> None:
    print("\n" + "=" * 60)
    print(f"{'stage':<30}{'outputs':>10}   status")
    print("-" * 60)
    for name, pattern, label in STAGES:
        n = _has_outputs(cfg, pattern)
        if name == "lineage":
            state = "OK (8-class)" if _lineage_is_8class(cfg) else ("STALE (6-class)" if n else "MISSING")
        elif name == "hormone_floor":
            state = "OK" if n else "MISSING"   # in-place; presence of restore_gated is the proxy
        else:
            state = "OK" if n else "MISSING"
        print(f"{label:<30}{n:>10}   {state}")
    print("=" * 60)

    # The optional PhenoCycler <-> Xenium integration layer picks up from `phenotype/broad/`.
    # It is a separate orchestrator with its own (heavier) environment, so it is only
    # mentioned here — never run as part of this pipeline, whose stage list and dependencies
    # are deliberately unchanged by it.
    n_integration = _has_outputs(cfg, "integration/*")
    if n_integration:
        print(f"[integration] {n_integration} artifact(s) under {cfg.data_dir/'integration'} — "
              f"`python -m phenocycler.integration.pipeline --status` for detail")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--only", nargs="*", default=None, choices=ORDER)
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
