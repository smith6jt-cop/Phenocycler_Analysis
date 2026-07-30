"""Donor-local Phenocycler expression calibration and hierarchical typing."""

from .config import PipelineConfig, load_config
from .marker_registry import MarkerRegistry, load_registry

__all__ = [
    "MarkerRegistry",
    "PipelineConfig",
    "load_config",
    "load_registry",
    "__version__",
]


def _resolve_version() -> str:
    """Version from installed metadata, falling back to the setuptools-scm generated file.

    Deliberately not a literal. setuptools-scm derives it from the nearest git tag, so a
    tagged commit yields `0.4.0` and anything after it yields `0.4.1.dev3+g2fc0002` — which
    means a set of results can always be traced back to the exact tree that produced it. A
    hardcoded string drifts from the tag the first time someone bumps one and not the other.

    That is the same argument this branch already makes one level down, where run_id hashes
    the source tree so a result can be tied to the code that produced it; a literal
    ``__version__`` would be the one provenance claim in the package that is maintained by
    hand.

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
