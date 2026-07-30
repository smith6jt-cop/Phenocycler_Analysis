"""Compatibility views derived from the marker registry.

No marker roster is authored here.  The versioned JSON registry is the only
source of marker role, selected compartment and spillover policy.
"""

from __future__ import annotations

import numpy as np

from .marker_registry import load_registry


def _registry():
    return load_registry()


TYPE = [marker.name for marker in _registry().markers if marker.kind == "type"]
PROCESS = [marker.name for marker in _registry().markers if marker.kind == "state"]
EXCLUDED = [
    marker.name
    for marker in _registry().markers
    if marker.kind == "excluded" or "structural_qc" in marker.uses
]
MASK = [
    marker.name
    for marker in _registry().markers
    if "microenvironment" in marker.uses
]
PROLIFERATION = [
    marker.name
    for marker in _registry().markers
    if "proliferation" in marker.uses
]
NO_SPILLOVER = [
    marker.name
    for marker in _registry().markers
    if marker.spillover_policy != "redsea_corrected"
]


def spillover_corrected_mask(cols):
    """Return the registry-declared REDSEA compensation mask."""

    policies = {
        marker.name: marker.spillover_policy for marker in _registry().markers
    }
    unknown = sorted(set(cols).difference(policies))
    if unknown:
        raise ValueError(f"qptiff channel(s) absent from marker registry: {unknown}")
    return np.asarray(
        [policies[column] == "redsea_corrected" for column in cols],
        dtype=bool,
    )


def heatmap_markers(available):
    """Identity markers present in ``available``, in registry order."""

    available = set(available)
    return [marker for marker in TYPE if marker in available]


__all__ = [
    "EXCLUDED",
    "MASK",
    "NO_SPILLOVER",
    "PROCESS",
    "PROLIFERATION",
    "TYPE",
    "heatmap_markers",
    "spillover_corrected_mask",
]
