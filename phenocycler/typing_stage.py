"""Production wrapper for simultaneous broad and parent-constrained typing."""

from __future__ import annotations

import pandas as pd

from .hierarchical_typing import (
    AUTHORITATIVE_ANCHOR,
    INFERRED,
    UNAVAILABLE,
    UNAVAILABLE_LABEL,
    AssignmentThresholds,
    TypingRegistry,
    assign_broad_types,
    assign_hierarchical_types,
    assign_subtypes,
)
from .marker_registry import MarkerRegistry, load_registry
from .typing_model_bundle import (
    TypingModelBundle,
    TypingModelPromotionManifest,
    ensemble_stability_wide,
)

METHOD_VERSION = "hierarchical-typing-stage-v3"
RULE_SCORE_SEMANTICS = "uncalibrated_rule_weight_softmax_score"
PRODUCTION_TYPING_MODE = "probabilistic_bundle"
CANDIDATE_TYPING_MODE = "probabilistic_candidate"


def type_calibrated_cells(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None = None,
    typing_registry: TypingRegistry | None = None,
    thresholds: AssignmentThresholds | None = None,
    model_bundle: TypingModelBundle | None = None,
    model_promotion_manifest: TypingModelPromotionManifest | None = None,
) -> pd.DataFrame:
    """Assign production cells while retaining uncertainty and alternatives.

    Rules-only typing requires no release artifact. Probabilistic production
    typing requires the exact :class:`TypingModelPromotionManifest` that was
    loaded and signature-verified for ``model_bundle``. The keyless validation
    here deliberately consumes that loader-issued, nonserialized capability;
    deserializing or constructing a manifest directly cannot authorize this
    entry point.
    """

    if model_bundle is None:
        if model_promotion_manifest is not None:
            raise ValueError(
                "a typing promotion manifest cannot be used without its model bundle"
            )
        mode = "rules_only"
    else:
        if not isinstance(model_promotion_manifest, TypingModelPromotionManifest):
            raise TypeError(
                "production probabilistic typing requires a verified "
                "TypingModelPromotionManifest"
            )
        model_promotion_manifest.validate_for_bundle(model_bundle)
        mode = PRODUCTION_TYPING_MODE
    return _type_calibrated_cells_impl(
        evidence_wide,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        thresholds=thresholds,
        model_bundle=model_bundle,
        typing_mode=mode,
    )


def type_candidate_cells(
    evidence_wide: pd.DataFrame,
    *,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry | None = None,
    typing_registry: TypingRegistry | None = None,
    thresholds: AssignmentThresholds | None = None,
) -> pd.DataFrame:
    """Run an unpromoted bundle in the isolated development/evaluation lane.

    This explicit entry point is not a production authorization. Its output is
    stamped ``typing_mode='probabilistic_candidate'`` and is accepted only by
    the candidate review, threshold-selection, and held-out evaluation
    contracts.
    """

    if model_bundle is None:
        raise TypeError("candidate probabilistic typing requires a model bundle")
    return _type_calibrated_cells_impl(
        evidence_wide,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        thresholds=thresholds,
        model_bundle=model_bundle,
        typing_mode=CANDIDATE_TYPING_MODE,
    )


def _type_calibrated_cells_impl(
    evidence_wide: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry | None,
    typing_registry: TypingRegistry | None,
    thresholds: AssignmentThresholds | None,
    model_bundle: TypingModelBundle | None,
    typing_mode: str,
) -> pd.DataFrame:
    """Shared deterministic implementation behind the gated public lanes."""

    if model_bundle is None and typing_mode != "rules_only":
        raise ValueError("rules-only inference must use typing_mode='rules_only'")
    if model_bundle is not None and typing_mode not in {
        PRODUCTION_TYPING_MODE,
        CANDIDATE_TYPING_MODE,
    }:
        raise ValueError("unknown probabilistic typing mode")

    markers = load_registry() if marker_registry is None else marker_registry
    rules = (
        TypingRegistry.from_marker_registry(markers)
        if typing_registry is None
        else typing_registry
    )
    broad_model = None
    subtype_models = None
    broad_stability = None
    subtype_stability = None
    effective_thresholds = thresholds or AssignmentThresholds()
    if model_bundle is not None:
        model_bundle.validate_current(
            marker_registry=markers,
            typing_registry=rules,
        )
        # The input fingerprint captured when the bundle was loaded is part of
        # the type-stage contract.  Validate bytes here, immediately before
        # inference, rather than relying only on the pipeline's fast check.
        if model_bundle.input_artifact is not None:
            model_bundle.input_artifact.validate_current(mode="content")
        if thresholds is not None and thresholds != model_bundle.thresholds:
            raise ValueError(
                "configured typing model bundle owns its frozen thresholds"
            )
        effective_thresholds = model_bundle.thresholds
        broad_model, subtype_models = model_bundle.materialize_models(rules)
        broad_stability = ensemble_stability_wide(
            evidence_wide,
            ensemble=model_bundle.broad,
        )
        broad_assignments = assign_broad_types(
            evidence_wide,
            rules,
            model=broad_model,
            stability=broad_stability,
            thresholds=effective_thresholds,
        )
        evidence_keys = evidence_wide.loc[:, ["donor_id", "object_id"]].astype(str)
        broad_keys = broad_assignments.loc[:, ["donor_id", "object_id"]].astype(str)
        if not evidence_keys.reset_index(drop=True).equals(
            broad_keys.reset_index(drop=True)
        ):
            raise AssertionError(
                "broad typing changed cell order before parent-constrained stability"
            )
        accepted_parent = broad_assignments["broad_assignment_status"].isin(
            {AUTHORITATIVE_ANCHOR, INFERRED}
        )
        subtype_stability = {}
        for ensemble in model_bundle.subtypes:
            parent = str(ensemble.parent)
            parent_mask = accepted_parent & broad_assignments["broad_label"].eq(
                parent
            )
            if not parent_mask.any():
                continue
            # A subtype model can only affect cells already accepted into its
            # broad parent.  Restricting both the feature transform and every
            # bootstrap multiplication here avoids one full-donor stability
            # table per subtype parent on million-cell images.
            subtype_stability[parent] = ensemble_stability_wide(
                evidence_wide.loc[parent_mask.to_numpy()],
                ensemble=ensemble,
            )
        assignments = assign_subtypes(
            evidence_wide,
            broad_assignments,
            rules,
            models=subtype_models,
            stability=subtype_stability,
            thresholds=effective_thresholds,
        )
    else:
        assignments = assign_hierarchical_types(
            evidence_wide,
            rules,
            broad_model=broad_model,
            subtype_models=subtype_models,
            broad_stability=broad_stability,
            subtype_stability=subtype_stability,
            thresholds=effective_thresholds,
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
    inference_pass_columns = (
        "unique_winner_pass",
        "probability_pass",
        "margin_pass",
        "stability_pass",
        "valid_evidence_pass",
        "threshold_eligible",
    )
    for prefix in ("broad", "subtype"):
        for suffix in inference_pass_columns:
            column = f"{prefix}_{suffix}"
            if column in assignments:
                assignments.loc[excluded, column] = False
        failure_column = f"{prefix}_inference_failure_reasons"
        if failure_column in assignments:
            assignments.loc[excluded, failure_column] = "geometry_qc_ineligible"
    # Required-marker summaries are inference gates, not valid statements
    # about cells removed by the upstream geometry/QC boundary.
    for prefix in ("broad", "subtype"):
        for suffix in ("all_required_valid", "all_required_negative"):
            column = f"{prefix}_{suffix}"
            if column in assignments:
                assignments.loc[excluded, column] = False
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
    assignments["typing_mode"] = typing_mode
    assignments["typing_model_bundle_content_id"] = (
        "" if model_bundle is None else model_bundle.content_id
    )
    assignments["typing_model_bundle_version"] = (
        "" if model_bundle is None else model_bundle.bundle_version
    )
    assignments["typing_model_method_version"] = (
        "" if model_bundle is None else model_bundle.method_version
    )
    assignments["typing_feature_schema_version"] = (
        "" if model_bundle is None else model_bundle.feature_schema_version
    )
    assignments["typing_score_semantics"] = (
        RULE_SCORE_SEMANTICS
        if model_bundle is None
        else model_bundle.score_semantics
    )
    # The legacy/global field describes broad scores. Subtype bundles may be
    # partial by policy, so stamp the actual parent-level source rather than
    # mislabeling a registry fallback as a calibrated probability.
    assignments["typing_broad_score_semantics"] = assignments[
        "typing_score_semantics"
    ]
    configured_subtype_parents = {
        str(rule.parent) for rule in rules.subtype_rules
    }
    accepted_parent = assignments["broad_assignment_status"].isin(
        {AUTHORITATIVE_ANCHOR, INFERRED}
    )
    subtype_applicable = accepted_parent & assignments["broad_label"].isin(
        configured_subtype_parents
    )
    assignments["typing_subtype_score_semantics"] = ""
    assignments["typing_subtype_model_source"] = "not_applicable"
    if model_bundle is None:
        assignments.loc[
            subtype_applicable, "typing_subtype_score_semantics"
        ] = RULE_SCORE_SEMANTICS
        assignments.loc[
            subtype_applicable, "typing_subtype_model_source"
        ] = "rules_only"
    else:
        bundled_parents = set(model_bundle.subtype_by_parent)
        bundled = subtype_applicable & assignments["broad_label"].isin(
            bundled_parents
        )
        fallback = subtype_applicable & ~bundled
        assignments.loc[
            bundled, "typing_subtype_score_semantics"
        ] = model_bundle.score_semantics
        assignments.loc[
            bundled, "typing_subtype_model_source"
        ] = "probabilistic_bundle"
        assignments.loc[
            fallback, "typing_subtype_score_semantics"
        ] = RULE_SCORE_SEMANTICS
        assignments.loc[
            fallback, "typing_subtype_model_source"
        ] = "rules_fallback"
    assignments["typing_broad_bootstrap_replicates"] = (
        0 if model_bundle is None else model_bundle.broad.bootstrap_replicates
    )
    assignments["typing_broad_valid_bootstrap_replicates"] = (
        0
        if model_bundle is None
        else model_bundle.broad.valid_bootstrap_replicates
    )
    if model_bundle is None:
        assignments["typing_subtype_bootstrap_replicates"] = 0
        assignments["typing_subtype_valid_bootstrap_replicates"] = 0
    else:
        subtype_replicates = {
            str(ensemble.parent): ensemble.bootstrap_replicates
            for ensemble in model_bundle.subtypes
        }
        valid_subtype_replicates = {
            str(ensemble.parent): ensemble.valid_bootstrap_replicates
            for ensemble in model_bundle.subtypes
        }
        assignments["typing_subtype_bootstrap_replicates"] = (
            assignments["broad_label"]
            .map(subtype_replicates)
            .fillna(0)
            .astype(int)
        )
        assignments["typing_subtype_valid_bootstrap_replicates"] = (
            assignments["broad_label"]
            .map(valid_subtype_replicates)
            .fillna(0)
            .astype(int)
        )
    return assignments


__all__ = [
    "CANDIDATE_TYPING_MODE",
    "METHOD_VERSION",
    "PRODUCTION_TYPING_MODE",
    "RULE_SCORE_SEMANTICS",
    "type_calibrated_cells",
    "type_candidate_cells",
]
