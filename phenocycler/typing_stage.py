"""Production wrapper for simultaneous broad and parent-constrained typing."""

from __future__ import annotations

import pandas as pd

from .hierarchical_typing import (
    AssignmentThresholds,
    TypingRegistry,
    UNAVAILABLE,
    UNAVAILABLE_LABEL,
    assign_hierarchical_types,
)
from .marker_registry import MarkerRegistry, load_registry


METHOD_VERSION = "hierarchical-typing-stage-v1"


def type_calibrated_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    typing_registry: TypingRegistry | None = None,
    thresholds: AssignmentThresholds | None = None,
) -> pd.DataFrame:
    """Assign every cell while retaining uncertainty and ranked alternatives.

    Curated anchors can establish a label directly. Non-anchor inference remains
    ambiguous unless a separately versioned whole-donor stability result is
    supplied through the lower-level APIs; this wrapper never invents a
    stability score from the same cells being classified.
    """

    markers = load_registry() if marker_registry is None else marker_registry
    rules = (
        TypingRegistry.from_marker_registry(markers)
        if typing_registry is None
        else typing_registry
    )
    assignments = assign_hierarchical_types(
        evidence_wide,
        rules,
        thresholds=thresholds or AssignmentThresholds(),
    )
    keys = evidence_wide.loc[:, ["donor_id", "object_id"]].copy()
    keys["donor_id"] = keys["donor_id"].astype(str)
    keys["object_id"] = keys["object_id"].astype(str)
    if keys.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("calibrated evidence contains duplicate cells")
    if len(assignments) != len(keys):
        raise AssertionError(
            f"typing returned {len(assignments):,}/{len(keys):,} cells"
        )
    observed = set(map(tuple, assignments[["donor_id", "object_id"]].to_numpy()))
    expected = set(map(tuple, keys.to_numpy()))
    if observed != expected:
        raise AssertionError("typing changed the donor/object universe")

    if "qc_analysis_eligible" not in evidence_wide:
        raise ValueError(
            "calibrated evidence is missing qc_analysis_eligible"
        )
    eligibility = evidence_wide.loc[
        :, ["donor_id", "object_id", "qc_analysis_eligible"]
    ].copy()
    eligibility["donor_id"] = eligibility["donor_id"].astype(str)
    eligibility["object_id"] = eligibility["object_id"].astype(str)
    if eligibility["qc_analysis_eligible"].isna().any() or not eligibility[
        "qc_analysis_eligible"
    ].isin([True, False, 0, 1]).all():
        raise ValueError("qc_analysis_eligible must contain explicit booleans")
    eligibility = eligibility.set_index(["donor_id", "object_id"])[
        "qc_analysis_eligible"
    ].astype(bool)
    assignment_index = pd.MultiIndex.from_frame(
        assignments[["donor_id", "object_id"]].astype(str)
    )
    analysis_eligible = eligibility.reindex(assignment_index).to_numpy()
    if pd.isna(analysis_eligible).any():
        raise AssertionError("typing could not align geometry-QC eligibility")
    assignments["qc_analysis_eligible"] = analysis_eligible.astype(bool)
    excluded = ~assignments["qc_analysis_eligible"]
    assignments.loc[excluded, "broad_label"] = UNAVAILABLE_LABEL
    assignments.loc[excluded, "broad_type"] = UNAVAILABLE_LABEL
    assignments.loc[excluded, "broad_assignment_status"] = UNAVAILABLE
    assignments.loc[excluded, "broad_reason"] = "geometry_qc_ineligible"
    assignments.loc[excluded, "subtype_label"] = UNAVAILABLE_LABEL
    assignments.loc[excluded, "subtype_assignment_status"] = UNAVAILABLE
    assignments.loc[excluded, "subtype_reason"] = "geometry_qc_ineligible"
    assignments.loc[excluded, "specific_type"] = UNAVAILABLE_LABEL
    assignments.loc[excluded, "cell_type"] = UNAVAILABLE_LABEL
    assignments.loc[excluded, "assignment_status"] = UNAVAILABLE
    assignments.loc[excluded, "specific_assignment_status"] = UNAVAILABLE
    assignments.loc[excluded, "confidence"] = float("nan")
    assignments.loc[excluded, "broad_authoritative_markers"] = ""
    assignments.loc[excluded, "broad_anchor_classes"] = ""
    assignments["typing_rules_version"] = rules.rules_version
    assignments["typing_rules_fingerprint"] = rules.rules_fingerprint
    assignments["marker_registry_fingerprint"] = markers.fingerprint
    return assignments


__all__ = ["METHOD_VERSION", "type_calibrated_cells"]
