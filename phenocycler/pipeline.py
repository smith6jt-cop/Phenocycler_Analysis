#!/usr/bin/env python3
"""
Orchestrator — run the full raw-data → 5-compartment broad-phenotyping pipeline end to end.

Idempotent ``run_step`` (skips a stage when its outputs already exist unless
``force=True``) + a final status table, in-process against the ``phenocycler``
package.  Stage order:

    cells → redsea → restore → lineage → qupath → figures

``restore`` thresholds every marker in ONE pass (curated directional pairs); ``lineage`` is the ordered
residual gating tree (5 compartments + explicit ``Other``). The old ``restore_extra`` + ``hormone_floor``
stages were removed with the curated multi-pair redesign.

    python -m phenocycler.pipeline                       # run every missing step
    python -m phenocycler.pipeline --only lineage figures
    python -m phenocycler.pipeline --force               # re-run everything
    python -m phenocycler.pipeline --status              # just print the status table
"""

from __future__ import annotations

import argparse
import glob
import shutil
import subprocess
from pathlib import Path

from .config import PipelineConfig, load_config


# step name -> (output pattern under data_dir, human label) for the status table
STAGES = [
    ("cells", "cells/donor_id=*", "Raw cells (DuckDB)"),
    ("redsea", "cells_redsea/donor_id=*", "REDSEA corrected"),
    ("restore", "restore_gated_redsea/donor_id=*", "RESTORE gated (accepted pairs)"),
    ("proliferation", "restore_gated_proliferation/donor_id=*", "Proliferation (Ki67/PCNA bimodal)"),
    ("restore_diag", "restore_redsea/mxnorm/restore_mxnorm_efficacy.csv", "RESTORE efficacy (mxnorm)"),
    ("lineage", "phenotype/broad/donor_id=*", "Compartments (5 + Other)"),
    ("qupath", "phenotype/qupath_class/pheno_class_*.csv", "QuPath CSVs"),
    ("figures", "phenotype/celltype_marker_dotplot.png", "Identity QC figures"),
]

ORDER = ["cells", "redsea", "restore", "proliferation", "restore_diag", "lineage", "qupath", "figures"]


def _has_outputs(cfg: PipelineConfig, pattern: str) -> int:
    return len(glob.glob(str(cfg.data_dir / pattern)))


def _lineage_is_current(cfg: PipelineConfig) -> bool:
    """True only if broad partitions exist AND carry the CURRENT ordered-residual schema
    (``compartment``/``cell_type``). A stale ``broad_lineage``/score_* ``broad/`` counts as 'not done',
    so the lineage step re-runs.
    """
    files = sorted(glob.glob(str(cfg.broad_dir / "donor_id=*" / "*.parquet")))
    if not files:
        return False
    try:
        import pyarrow.parquet as pq
        cols = set(pq.ParquetFile(files[0]).schema.names)
    except Exception:
        return False
    return "compartment" in cols and "cell_type" in cols


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


def _run_restore_diagnostic(cfg: PipelineConfig) -> None:
    """Shell out to the mxnorm RESTORE before/after efficacy diagnostic (an R script).

    Runs under **R 4.5.1** (``/usr/local/bin/Rscript``): the shared ``/usr/local/lib/R/site-library``
    is built for R 4.6 and its ``rlang.so`` will not load under 4.5.1, whereas the 4.5 personal lib has
    a working rlang + arrow/lme4/uwot + mxnorm (see ``scripts/R/install_mxnorm.R``). Reads the applied
    per-(donor,marker) thresholds from ``restore_redsea/positive_fractions.csv``; writes the efficacy
    CSV + figures under ``cfg.restore_mxnorm_dir``.
    """
    repo = cfg.data_dir.parent
    # The R diagnostic ships INSIDE this submodule (self-contained); cwd stays at the data-dir parent so
    # the script's relative --data fallbacks resolve there (all real inputs are passed as absolute paths).
    script = Path(__file__).resolve().parents[1] / "scripts" / "R" / "restore_mxnorm_diagnostics.R"
    rscript = "/usr/local/bin/Rscript" if Path("/usr/local/bin/Rscript").exists() else (shutil.which("Rscript") or "Rscript")
    cmd = [
        rscript, str(script),
        "--threshold-csv", str(cfg.restore_dir / "positive_fractions.csv"),
        "--cells-redsea-dir", str(cfg.cells_redsea_dir),
        "--cells-dir", str(cfg.cells_dir),
        "--donor-metadata", str(cfg.donor_metadata),
        "--outdir", str(cfg.restore_mxnorm_dir),
        "--subsample", str(cfg.restore_diag_subsample),
        "--seed", str(cfg.restore_diag_seed),
    ]
    print("[restore_diag] " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=str(repo))


def run_pipeline(cfg: PipelineConfig, *, only=None, force=False) -> None:
    from . import (cells_parquet, redsea, restore, restore_apply, lineage, qupath_export, figures)

    def _cells():
        cells_parquet.build_cells_parquet(cfg)

    def _redsea():
        params = redsea.RedseaParams.from_config(cfg)
        redsea.run_redsea(cfg, cfg.discover_donors(), params, n_jobs=cfg.n_jobs)

    def _restore():
        # Gate 4 (2026-07-24): the production apply runs the frozen ACCEPTED_PAIRS (one reference per
        # target, maintainer sign-off) on all cells per donor and writes the tri-state gated parquet that
        # `lineage` consumes. It replaces the retired paired-SSC driver (deleted with its CD20 pair web).
        restore_apply.run_apply(cfg)

    def _proliferation():
        restore.run_restore_proliferation(cfg)

    def _restore_diag():
        _run_restore_diagnostic(cfg)

    def _lineage():
        lineage.run_lineage(cfg)

    def _qupath():
        qupath_export.export_qupath_classes(cfg)

    def _figures():
        figures.run_figures(cfg)

    # (func, existence-pattern, column-checker).
    steps = {
        "cells": (_cells, "cells/donor_id=*", None),
        "redsea": (_redsea, "cells_redsea/donor_id=*", None),
        "restore": (_restore, "restore_gated_redsea/donor_id=*", None),
        "proliferation": (_proliferation, "restore_gated_proliferation/donor_id=*", None),
        "restore_diag": (_restore_diag, "restore_redsea/mxnorm/restore_mxnorm_efficacy.csv", None),
        "lineage": (_lineage, "phenotype/broad/donor_id=*", _lineage_is_current),
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
            state = "OK (5+Other)" if _lineage_is_current(cfg) else ("STALE (re-run)" if n else "MISSING")
        else:
            state = "OK" if n else "MISSING"
        print(f"{label:<30}{n:>10}   {state}")
    print("=" * 60)


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
