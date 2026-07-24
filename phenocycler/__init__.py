"""
phenocycler — Phenocycler / CODEX spatial-proteomics analysis pipeline.

This package ports the *Islet-Explorer-Senior* "Senior phenotyping pipeline"
(raw QuPath outputs → broad lineage classification) into a reusable, testable,
notebook-friendly Python package:

    build_cells_parquet   (DuckDB)          QuPath CSV        -> data/cells/
    redsea                (pixel spillover) qptiff + GeoJSON  -> data/cells_redsea/
    restore               (RESTORE norm)    cells_redsea      -> data/restore_gated_redsea/
    lineage               (5 compartments)  restore_gated     -> data/phenotype/broad/
    qupath_export         (QC round-trip)   broad             -> data/phenotype/qupath_class/
    figures               (identity QC)     restore_gated     -> data/phenotype/celltype_marker_*.png

The scientific algorithms are faithful ports of the upstream scripts; the
additions are (a) config-driven paths (no hardcoded absolute paths), (b)
per-donor multiprocessing, and (c) an optional CuPy/RAPIDS GPU backend with a
transparent CPU fallback.  See ``config.ini`` and the notebooks/ directory.
"""

from __future__ import annotations

import os as _os
# Avoid the BLAS/OpenMP threadpool deadlock observed in RESTORE's spams/sklearn fit (a `futex_wait`
# hang — 0 markers in 15 h). Default the numeric libraries to single-threaded BEFORE numpy/spams load;
# pipeline parallelism is per-donor (`--jobs`), so this does not cost throughput. Override by exporting
# the vars yourself.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
    _os.environ.setdefault(_v, "1")

from .config import (
    PipelineConfig, load_config,
    COMPARTMENT_ORDER, COMPARTMENT_GATES, COMPARTMENT_COLORS, COMPARTMENT_ABBR, OTHER_LABEL,
    UNRESOLVED_LABEL, ENDOCRINE_MARKERS, STATUS_ORDER,
)
from . import marker_taxonomy

__all__ = [
    "PipelineConfig", "load_config",
    "COMPARTMENT_ORDER", "COMPARTMENT_GATES", "COMPARTMENT_COLORS", "COMPARTMENT_ABBR", "OTHER_LABEL",
    "UNRESOLVED_LABEL", "ENDOCRINE_MARKERS", "STATUS_ORDER",
    "marker_taxonomy",
]

__version__ = "0.4.0"
