"""
phenocycler — Phenocycler / CODEX spatial-proteomics analysis pipeline.

This package ports the *Islet-Explorer-Senior* "Senior phenotyping pipeline"
(raw QuPath outputs → broad lineage classification) into a reusable, testable,
notebook-friendly Python package:

    build_cells_parquet   (DuckDB)          QuPath CSV        -> data/cells/
    redsea                (pixel spillover) qptiff + GeoJSON  -> data/cells_redsea/
    restore               (RESTORE norm)    cells_redsea      -> data/restore_gated_redsea/
    restore  (extra pass) (RESTORE norm)    cells_redsea      -> data/restore_gated_redsea_extra/
    hormone_floor         (endo+immune fix) restore_gated     -> data/restore_gated_redsea/ (floored)
    lineage               (7 broad classes) restore_gated     -> data/phenotype/broad/
    qupath_export         (QC round-trip)   broad             -> data/phenotype/qupath_class/
    figures               (identity QC)     restore_gated     -> data/phenotype/celltype_marker_*.png

The scientific algorithms are faithful ports of the upstream scripts; the
additions are (a) config-driven paths (no hardcoded absolute paths), (b)
per-donor multiprocessing, and (c) an optional CuPy/RAPIDS GPU backend with a
transparent CPU fallback.  See ``config.ini`` and the notebooks/ directory.
"""

from __future__ import annotations

from .config import (
    PipelineConfig, load_config,
    LINEAGES, STRUCT_LINEAGES, LINEAGE_COLORS, LINEAGE_ABBR, STATUS_ORDER,
    DEFAULT_MARKER_PAIRS, EXTRA_MARKER_PAIRS, EXTRA_MARKERS, HORMONE_MARKERS,
)
from . import marker_taxonomy

__all__ = [
    "PipelineConfig", "load_config",
    "LINEAGES", "STRUCT_LINEAGES", "LINEAGE_COLORS", "LINEAGE_ABBR", "STATUS_ORDER",
    "DEFAULT_MARKER_PAIRS", "EXTRA_MARKER_PAIRS", "EXTRA_MARKERS", "HORMONE_MARKERS",
    "marker_taxonomy", "__version__",
]


def _resolve_version() -> str:
    """Version from installed metadata, falling back to the setuptools-scm generated file.

    Deliberately not a literal. setuptools-scm derives it from the nearest git tag, so a
    tagged commit yields `0.4.0` and anything after it yields `0.4.1.dev3+g2fc0002` — which
    means a set of results can always be traced to the exact tree that produced it. A
    hardcoded string drifts from the tag the first time someone bumps one and not the other.

    Three sources in order of trustworthiness: installed distribution metadata (correct for
    both regular and editable installs), the generated `_version.py` (present in a build
    tree), then a marker that says plainly that neither was available rather than inventing
    a number.
    """
    try:
        from importlib.metadata import PackageNotFoundError, version

        return version("phenocycler-analysis")
    except (ImportError, PackageNotFoundError):
        pass
    try:
        from ._version import __version__ as scm_version

        return scm_version
    except ImportError:
        return "0.0.0+unknown"


__version__ = _resolve_version()
