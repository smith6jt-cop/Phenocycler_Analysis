"""
phenocycler — Phenocycler / CODEX spatial-proteomics analysis pipeline.

This package ports the *Islet-Explorer-Senior* "Senior phenotyping pipeline"
(raw QuPath outputs → broad lineage classification) into a reusable, testable,
notebook-friendly Python package:

    build_cells_parquet   (DuckDB)          QuPath CSV        -> data/cells/
    redsea                (pixel spillover) qptiff + GeoJSON  -> data/cells_redsea/
    restore               (RESTORE norm)    cells_redsea      -> data/restore_gated_redsea/
    lineage               (broad lineage)   restore_gated     -> data/phenotype/broad/
    qupath_export         (QC round-trip)   broad             -> data/phenotype/qupath_class/

The scientific algorithms are faithful ports of the upstream scripts; the
additions are (a) config-driven paths (no hardcoded absolute paths), (b)
per-donor multiprocessing, and (c) an optional CuPy/RAPIDS GPU backend with a
transparent CPU fallback.  See ``config.ini`` and the notebooks/ directory.
"""

from __future__ import annotations

from .config import PipelineConfig, load_config

__all__ = ["PipelineConfig", "load_config"]

__version__ = "0.2.0"
