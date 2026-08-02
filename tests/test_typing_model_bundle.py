from __future__ import annotations

import copy
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from phenocycler.hierarchical_typing import (
    FEATURE_SCHEMA_VERSION,
    AssignmentThresholds,
    TypingRegistry,
    classifier_from_registry,
)
from phenocycler.marker_registry import load_registry
from phenocycler.refinement_contracts import (
    DEVELOPMENT_LABEL_PURPOSE,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
)
from phenocycler.typing_model_bundle import (
    FrozenModelEnsemble,
    PromotionGatePolicy,
    PromotionReplayReceipt,
    TypingModelBundle,
    TypingModelPromotionManifest,
    ensemble_stability_wide,
    load_typing_model_bundle,
    load_typing_model_promotion_manifest,
    promotion_replay_environment_sha256,
)


def _splits() -> DonorSplitManifest:
    return DonorSplitManifest(
        split_version="split-v1",
        source_run_id="source-run",
        cohort_content_id="cohort",
        development=("D1", "D2"),
        calibration=("D3",),
        validation=("D4",),
        locked_test=("D5",),
    )


def _label_provenance(splits: DonorSplitManifest) -> DevelopmentLabelProvenance:
    return DevelopmentLabelProvenance(
        bundle_id="independent-development-v1",
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
        fingerprint="a" * 64,
    )


def _evidence_provenance(
    splits: DonorSplitManifest,
) -> DevelopmentEvidenceProvenance:
    return DevelopmentEvidenceProvenance(
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        source_run_manifest_content_id="b" * 64,
        source_stage_manifest_content_id="c" * 64,
        source_artifact_content_id="d" * 64,
        labeled_frame_fingerprint="e" * 64,
        feature_schema_version=FEATURE_SCHEMA_VERSION,
        row_count=20,
    )


def _ensemble(
    rules: TypingRegistry,
    *,
    level: str = "broad",
    parent: str | None = None,
    n_bootstrap: int = 2,
) -> FrozenModelEnsemble:
    model = classifier_from_registry(rules, level=level, parent=parent)
    coefficients = tuple(
        tuple(float(item) for item in row)
        for row in np.asarray(model.coefficients, dtype=float)
    )
    return FrozenModelEnsemble(
        level=level,
        parent=parent,
        classes=model.classes,
        markers=model.markers,
        coefficients=coefficients,
        bootstrap_coefficients=tuple(coefficients for _ in range(n_bootstrap)),
        training_donors=("D1", "D2"),
        calibration_donors=("D3",),
        bootstrap_draws=tuple(
            ("D1", "D2") if index % 2 == 0 else ("D2", "D2")
            for index in range(n_bootstrap)
        ),
        bootstrap_valid=tuple(True for _ in range(n_bootstrap)),
        bootstrap_failures=tuple("" for _ in range(n_bootstrap)),
        temperature=1.25,
        n_development_labels=20,
        n_calibration_labels=10,
    )


def _final_specific_classes(rules: TypingRegistry) -> list[str]:
    return sorted(
        {
            leaf
            for broad_class in classifier_from_registry(rules).classes
            for leaf in (
                tuple(
                    rule.name
                    for rule in rules.rules_for(level="subtype", parent=broad_class)
                )
                or (broad_class,)
            )
        }
    )


def _fit_parameters(
    rules: TypingRegistry,
    *,
    modeled_subtype_parents: tuple[str, ...] = (),
    n_bootstrap: int = 2,
) -> dict[str, object]:
    configured_parents = {
        str(rule.parent) for rule in rules.subtype_rules if rule.parent is not None
    }
    structural = {
        parent
        for parent in configured_parents
        if len(rules.rules_for(level="subtype", parent=parent)) == 1
    }
    missing = configured_parents - structural - set(modeled_subtype_parents)
    return {
        "n_bootstrap": n_bootstrap,
        "random_state": 7,
        "l2": 0.02,
        "max_abs_weight": 8.0,
        "maxiter": 500,
        "require_all_subtype_parents": False,
        "structural_rules_only_parents": sorted(structural),
        "missing_probabilistic_subtype_parents": sorted(missing),
        "final_specific_classes": _final_specific_classes(rules),
        "minimum_valid_bootstrap_fraction": 0.8,
        "class_priors": "uniform_no_intercept",
        "donor_weighting": "class_then_donor_balanced",
    }


def _candidate_bundle(
    subtype_parents: tuple[str, ...] = (),
) -> tuple[TypingModelBundle, object, TypingRegistry]:
    markers = load_registry()
    rules = TypingRegistry.from_marker_registry(markers)
    splits = _splits()
    bundle = TypingModelBundle.create(
        bundle_version="candidate-v1",
        marker_registry_fingerprint=markers.fingerprint,
        typing_rules_fingerprint=rules.rules_fingerprint,
        split_manifest=splits,
        label_provenance=_label_provenance(splits),
        evidence_provenance=_evidence_provenance(splits),
        thresholds=AssignmentThresholds(),
        broad=_ensemble(rules),
        subtypes=tuple(
            _ensemble(rules, level="subtype", parent=parent)
            for parent in subtype_parents
        ),
        fit_parameters=_fit_parameters(
            rules,
            modeled_subtype_parents=subtype_parents,
        ),
    )
    bundle.validate_current(marker_registry=markers, typing_registry=rules)
    return bundle, markers, rules


def _promotion_policy(**updates: object) -> PromotionGatePolicy:
    values: dict[str, object] = {
        "minimum_donors": 5,
        "minimum_rows": 20,
        "minimum_coverage_lower": 0.70,
        "maximum_selective_error_upper": 0.10,
        "minimum_rescue_coverage_lower": 0.40,
        "maximum_rescue_selective_error_upper": 0.10,
        "maximum_class_false_positive_rate_upper": 0.10,
        "maximum_class_false_negative_rate_upper": 0.20,
        "minimum_class_positive_rows": 3,
        "minimum_class_negative_rows": 3,
        "minimum_class_positive_donors": 3,
        "minimum_class_negative_donors": 3,
        "confidence_level": 0.95,
        "bootstrap_replicates": 1000,
    }
    values.update(updates)
    return PromotionGatePolicy(**values)


def test_candidate_bundle_roundtrip_is_content_addressed_and_registry_validated(
    tmp_path,
):
    bundle, markers, rules = _candidate_bundle()
    path = tmp_path / "typing-model.json"
    bundle.write_json(path)

    loaded = load_typing_model_bundle(
        path,
        marker_registry=markers,
        typing_registry=rules,
    )

    assert loaded.content_id == bundle.content_id
    assert loaded.created_at_utc == ""
    repeated, _, _ = _candidate_bundle()
    assert repeated.content_id == bundle.content_id
    timestamped = replace(
        bundle,
        created_at_utc="2026-08-01T00:00:00+00:00",
        content_id="",
    )
    assert timestamped.content_id == bundle.content_id
    assert loaded.input_artifact is not None
    assert loaded.to_dict() == bundle.to_dict()
    broad, subtypes = loaded.materialize_models(rules)
    assert broad.classes == bundle.broad.classes
    assert np.asarray(broad.coefficients) == pytest.approx(
        np.asarray(bundle.broad.coefficients) / bundle.broad.temperature
    )
    assert subtypes == {}
    with pytest.raises(FileExistsError, match="immutable"):
        bundle.write_json(path)


def test_model_fit_identity_ignores_release_thresholds_and_timestamp():
    candidate, _markers, _rules = _candidate_bundle()
    changed = replace(
        candidate,
        bundle_version="candidate-v2",
        thresholds=replace(
            candidate.thresholds,
            inferred_probability=0.91,
            inferred_margin=0.17,
            inferred_stability=0.88,
        ),
        created_at_utc="2026-08-01T00:00:00+00:00",
        content_id="",
    )

    assert changed.model_fit_content_id == candidate.model_fit_content_id
    assert changed.content_id != candidate.content_id
    assert changed.model_fit_identity_payload() == candidate.model_fit_identity_payload()


def test_model_fit_identity_rejects_stored_id_and_scientific_fit_tampering():
    candidate, _markers, _rules = _candidate_bundle()
    payload = candidate.to_dict()
    payload["model_fit_content_id"] = "f" * 64
    with pytest.raises(ValueError, match="model_fit_content_id is invalid"):
        TypingModelBundle.from_dict(payload)

    changed_model = replace(candidate.broad, temperature=1.5)
    with pytest.raises(ValueError, match="model_fit_content_id is invalid"):
        replace(candidate, broad=changed_model, content_id="")


def test_bundle_rejects_tampering_unknown_fields_and_registry_drift():
    bundle, markers, _rules = _candidate_bundle()
    payload = bundle.to_dict()
    tampered = copy.deepcopy(payload)
    tampered["thresholds"]["inferred_probability"] = 0.01
    with pytest.raises(ValueError, match="content_id is invalid"):
        TypingModelBundle.from_dict(tampered)

    unknown = copy.deepcopy(payload)
    unknown["prevalence_target"] = {"Immune": 0.2}
    with pytest.raises(ValueError, match="unknown fields"):
        TypingModelBundle.from_dict(unknown)

    missing_identity = copy.deepcopy(payload)
    missing_identity.pop("content_id")
    with pytest.raises(ValueError, match="missing content_id"):
        TypingModelBundle.from_dict(missing_identity)

    wrong_boolean = copy.deepcopy(payload)
    wrong_boolean["broad"]["fit_converged"] = "false"
    with pytest.raises(TypeError, match="JSON boolean"):
        TypingModelBundle.from_dict(wrong_boolean)

    drifted_rules = TypingRegistry.from_mapping(
        {
            "broad": {"Only": {"anchors": ["M"], "required": ["M"]}},
            "rules_fingerprint": "f" * 64,
        }
    )
    with pytest.raises(ValueError, match="rules fingerprint"):
        bundle.validate_current(
            marker_registry=markers,
            typing_registry=drifted_rules,
        )


def test_bundle_rejects_wrong_shape_nonfinite_and_constraint_violations():
    bundle, _markers, rules = _candidate_bundle()
    wrong_shape = bundle.to_dict()
    wrong_shape["broad"]["coefficients"] = [[0.0]]
    with pytest.raises(ValueError, match="shape"):
        TypingModelBundle.from_dict(wrong_shape)

    coefficients = np.asarray(bundle.broad.coefficients, dtype=float)
    template = classifier_from_registry(rules)
    positive_location = next(
        (
            class_index,
            template.markers.index((rule.anchor_markers + rule.support_markers)[0]),
        )
        for class_index, rule in enumerate(template.rules)
        if rule.anchor_markers or rule.support_markers
    )
    coefficients[positive_location] = -1.0
    broken = FrozenModelEnsemble(
        **{
            **bundle.broad.identity_payload(),
            "coefficients": coefficients,
        }
    )
    with pytest.raises(ValueError, match="positive sign"):
        broken.validate_registry(rules)

    nonfinite = bundle.broad.identity_payload()
    nonfinite["coefficients"][0][0] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        FrozenModelEnsemble.from_dict(nonfinite)

    outside_domain = bundle.broad.identity_payload()
    outside_domain["temperature"] = 0.01
    with pytest.raises(ValueError, match="lie in"):
        FrozenModelEnsemble.from_dict(outside_domain)

    too_large = bundle.broad.identity_payload()
    too_large["coefficients"][0][0] = 8.01
    with pytest.raises(ValueError, match="coefficient cap"):
        FrozenModelEnsemble.from_dict(too_large)

    malformed_fit = bundle.to_dict()
    malformed_fit["fit_parameters"].pop("max_abs_weight")
    with pytest.raises(ValueError, match="missing required fields"):
        TypingModelBundle.from_dict(malformed_fit)

    excessive_fit = bundle.to_dict()
    excessive_fit["fit_parameters"]["max_abs_weight"] = 9.0
    with pytest.raises(ValueError, match="max_abs_weight"):
        TypingModelBundle.from_dict(excessive_fit)


def test_registry_validation_checks_declared_subtype_inventory_and_final_ontology():
    bundle, markers, rules = _candidate_bundle(("Immune",))
    assert bundle.fit_parameters["missing_probabilistic_subtype_parents"] == [
        "Endothelial",
        "Epithelial",
    ]

    wrong_missing = dict(bundle.fit_parameters)
    wrong_missing["missing_probabilistic_subtype_parents"] = ["Endothelial"]
    inconsistent = replace(
        bundle,
        fit_parameters=wrong_missing,
        model_fit_content_id="",
        content_id="",
    )
    with pytest.raises(ValueError, match="declared missing probabilistic subtype"):
        inconsistent.validate_current(marker_registry=markers, typing_registry=rules)

    wrong_ontology = dict(bundle.fit_parameters)
    wrong_ontology["final_specific_classes"] = ["Alpha", "Beta"]
    inconsistent = replace(
        bundle,
        fit_parameters=wrong_ontology,
        model_fit_content_id="",
        content_id="",
    )
    with pytest.raises(ValueError, match="final specific classes"):
        inconsistent.validate_current(marker_registry=markers, typing_registry=rules)


def test_chunked_stability_matches_naive_and_is_row_order_invariant():
    ensemble = FrozenModelEnsemble(
        level="broad",
        parent=None,
        classes=("A", "B"),
        markers=("A_marker", "B_marker"),
        coefficients=((4.0, 0.0), (0.0, 4.0)),
        bootstrap_coefficients=(
            ((4.0, 0.0), (0.0, 4.0)),
            ((4.0, 0.0), (0.0, 0.0)),
        ),
        training_donors=("D1", "D2"),
        calibration_donors=("D3",),
        bootstrap_draws=(("D1", "D2"), ("D2", "D1")),
        bootstrap_valid=(True, True),
        bootstrap_failures=("", ""),
        temperature=1.0,
        n_development_labels=4,
        n_calibration_labels=2,
    )
    evidence = pd.DataFrame(
        {
            "donor_id": ["D", "D", "D"],
            "object_id": ["a", "b", "invalid"],
            "A_marker__state": ["positive", "negative", "indeterminate"],
            "A_marker__log2_threshold_ratio": [3.0, -2.0, 100.0],
            "A_marker__model_valid": [True, True, False],
            "B_marker__state": ["negative", "positive", "indeterminate"],
            "B_marker__log2_threshold_ratio": [-2.0, 3.0, 100.0],
            "B_marker__model_valid": [True, True, False],
        }
    )
    observed = ensemble_stability_wide(
        evidence,
        ensemble=ensemble,
        chunk_size=1,
    ).set_index("object_id")

    assert observed.loc["a", "best_candidate"] == "A"
    assert observed.loc["a", "stability"] == 1.0
    assert observed.loc["b", "best_candidate"] == "B"
    assert observed.loc["b", "stability"] == 0.5
    # Invalid high-valued evidence is masked, leaving a deterministic tie.
    assert observed.loc["invalid", "best_candidate"] == "A"
    assert observed.loc["invalid", "stability"] == 0.0

    shuffled = ensemble_stability_wide(
        evidence.sample(frac=1, random_state=4),
        ensemble=ensemble,
        chunk_size=2,
    ).set_index("object_id")
    pd.testing.assert_frame_equal(observed.sort_index(), shuffled.sort_index())


def test_invalid_bootstrap_draws_are_conservative_denominator_members():
    ensemble = FrozenModelEnsemble(
        level="broad",
        parent=None,
        classes=("A", "B"),
        markers=("A_marker", "B_marker"),
        coefficients=((4.0, 0.0), (0.0, 4.0)),
        bootstrap_coefficients=(
            ((4.0, 0.0), (0.0, 4.0)),
            ((4.0, 0.0), (0.0, 4.0)),
        ),
        training_donors=("D1", "D2"),
        calibration_donors=("D3",),
        bootstrap_draws=(("D1", "D2"), ("D1", "D1")),
        bootstrap_valid=(True, False),
        bootstrap_failures=("", "missing_classes:B"),
        n_development_labels=4,
        n_calibration_labels=2,
    )
    evidence = pd.DataFrame(
        {
            "donor_id": ["D"],
            "object_id": ["a"],
            "A_marker__state": ["positive"],
            "B_marker__state": ["negative"],
        }
    )

    observed = ensemble_stability_wide(evidence, ensemble=ensemble).iloc[0]

    assert observed["stability"] == 0.5
    assert observed["bootstrap_replicates"] == 2
    assert observed["valid_bootstrap_replicates"] == 1


def test_promotion_policy_has_content_identity_strict_io_and_scientific_floors(
    tmp_path,
):
    policy = _promotion_policy()
    path = tmp_path / "promotion-policy.json"
    policy.write_json(path)

    assert PromotionGatePolicy.read_json(path) == policy
    with pytest.raises(FileExistsError, match="immutable"):
        policy.write_json(path)
    with pytest.raises(ValueError, match="absolute production floors"):
        _promotion_policy(minimum_donors=4)
    with pytest.raises(ValueError, match="must be > 0"):
        _promotion_policy(minimum_coverage_lower=0.0)
    with pytest.raises(ValueError, match="must be < 1"):
        _promotion_policy(maximum_selective_error_upper=1.0)


def test_promotion_loader_grants_capability_only_after_key_validation(monkeypatch):
    promotion = object.__new__(TypingModelPromotionManifest)
    object.__setattr__(promotion, "_signature_verified", False)
    seen: dict[str, object] = {}

    monkeypatch.setattr(
        TypingModelPromotionManifest,
        "read_json",
        classmethod(lambda cls, path: promotion),
    )

    def validate_for_bundle(self, model_bundle, **kwargs):
        assert self._signature_verified is False
        seen["bundle"] = model_bundle
        seen.update(kwargs)

    monkeypatch.setattr(
        TypingModelPromotionManifest,
        "validate_for_bundle",
        validate_for_bundle,
    )
    bundle = SimpleNamespace(content_id="bundle")

    loaded = load_typing_model_promotion_manifest(
        "promotion.json",
        model_bundle=bundle,
        release_signing_key="external-release-secret",
    )

    assert loaded is promotion
    assert loaded._signature_verified is True
    assert seen["bundle"] is bundle
    assert seen["release_signing_key"] == "external-release-secret"
    capability_field = TypingModelPromotionManifest.__dataclass_fields__[
        "_signature_verified"
    ]
    assert capability_field.init is False
    assert capability_field.compare is False


def test_failed_promotion_validation_never_grants_capability(monkeypatch):
    promotion = object.__new__(TypingModelPromotionManifest)
    object.__setattr__(promotion, "_signature_verified", False)
    monkeypatch.setattr(
        TypingModelPromotionManifest,
        "read_json",
        classmethod(lambda cls, path: promotion),
    )

    def reject(self, model_bundle, **kwargs):
        raise ValueError("receipt signature is invalid")

    monkeypatch.setattr(TypingModelPromotionManifest, "validate_for_bundle", reject)

    with pytest.raises(ValueError, match="signature is invalid"):
        load_typing_model_promotion_manifest(
            "promotion.json",
            model_bundle=SimpleNamespace(),
            release_signing_key="wrong-release-secret",
        )
    assert promotion._signature_verified is False


def test_promotion_receipt_timestamp_must_be_utc():
    digest = "a" * 64
    fields = {
        name: digest
        for name in (
            "typing_model_bundle_content_id",
            "model_fit_content_id",
            "split_manifest_content_id",
            "threshold_selection_provenance_content_id",
            "preliminary_bundle_content_id",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
            "calibration_assignment_run_id",
            "calibration_source_run_manifest_content_id",
            "calibration_source_stage_manifest_content_id",
            "calibration_source_artifact_content_id",
            "calibration_review_context_artifact_content_id",
            "promotion_policy_content_id",
            "validation_report_content_id",
            "locked_test_report_content_id",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
            "replay_code_sha256",
            "replay_environment_sha256",
        )
    }

    with pytest.raises(ValueError, match="must use UTC"):
        PromotionReplayReceipt(
            **fields,
            issuer="release-operator",
            issued_at_utc="2026-08-01T12:00:00-04:00",
            signature_hmac_sha256=digest,
        )

    with pytest.raises(ValueError, match="canonical ISO UTC form"):
        PromotionReplayReceipt(
            **fields,
            issuer="release-operator",
            issued_at_utc="2026-08-01T16:00:00Z",
            signature_hmac_sha256=digest,
        )


def test_promotion_environment_digest_binds_optional_image(monkeypatch):
    try:
        promotion_replay_environment_sha256.cache_clear()
        monkeypatch.delenv("PHENOCYCLER_RUNTIME_IMAGE_DIGEST", raising=False)
        local_digest = promotion_replay_environment_sha256()
        promotion_replay_environment_sha256.cache_clear()
        monkeypatch.setenv("PHENOCYCLER_RUNTIME_IMAGE_DIGEST", "b" * 64)
        image_bound_digest = promotion_replay_environment_sha256()

        assert len(local_digest) == 64
        assert len(image_bound_digest) == 64
        assert image_bound_digest != local_digest
    finally:
        promotion_replay_environment_sha256.cache_clear()


def test_bundle_and_promotion_readers_reject_ambiguous_json(tmp_path):
    duplicate = tmp_path / "duplicate-bundle.json"
    duplicate.write_text('{"schema_version": 2, "schema_version": 2}\n')
    with pytest.raises(ValueError, match="duplicate key"):
        TypingModelBundle.read_json(duplicate)

    nonfinite = tmp_path / "nonfinite-promotion.json"
    nonfinite.write_text('{"schema_version": NaN}\n')
    with pytest.raises(ValueError, match="non-finite"):
        TypingModelPromotionManifest.read_json(nonfinite)
