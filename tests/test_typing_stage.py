from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
import pytest

from phenocycler.hierarchical_typing import AssignmentThresholds
from phenocycler.marker_registry import load_registry
from phenocycler.typing_model_bundle import TypingModelPromotionManifest
from phenocycler.typing_stage import (
    RULE_SCORE_SEMANTICS,
    type_calibrated_cells,
    type_candidate_cells,
)


def test_typing_stage_keeps_all_cells_and_provenance():
    registry = load_registry()
    frame = pd.DataFrame(
        {
            "donor_id": ["D", "D"],
            "object_id": ["immune", "uncertain"],
            "qc_analysis_eligible": [True, False],
        }
    )
    for marker in registry.active_markers:
        frame[f"{marker.name}__state"] = "negative"
        frame[f"{marker.name}__expression_probability"] = 0.0
        frame[f"{marker.name}__model_valid"] = True
    frame.loc[0, "CD3e__state"] = "positive"
    frame.loc[0, "CD3e__expression_probability"] = 0.999
    # This row would otherwise be a confident broad/subtype anchor.  Geometry
    # exclusion must invalidate every downstream inference gate without
    # erasing its ranked audit candidates.
    frame.loc[1, "CD3e__state"] = "positive"
    frame.loc[1, "CD3e__expression_probability"] = 0.999

    out = type_calibrated_cells(frame, marker_registry=registry).set_index("object_id")
    assert set(out.index) == {"immune", "uncertain"}
    assert out.loc["immune", "broad_type"] == "Immune"
    assert out.loc["immune", "assignment_status"] == "anchor"
    assert out.loc["uncertain", "best_broad_type"]
    assert out.loc["uncertain", "subtype_best_candidate"]
    assert out.loc["uncertain", "assignment_status"] == "unavailable"
    assert out.loc["uncertain", "broad_reason"] == "geometry_qc_ineligible"
    assert not out.loc["uncertain", "qc_analysis_eligible"]
    for prefix in ("broad", "subtype"):
        for suffix in (
            "unique_winner_pass",
            "probability_pass",
            "margin_pass",
            "stability_pass",
            "valid_evidence_pass",
            "threshold_eligible",
        ):
            assert not out.loc["uncertain", f"{prefix}_{suffix}"]
        assert (
            out.loc["uncertain", f"{prefix}_inference_failure_reasons"]
            == "geometry_qc_ineligible"
        )
    assert not out.loc["uncertain", "broad_all_required_valid"]
    assert not out.loc["uncertain", "broad_all_required_negative"]
    assert out["typing_rules_fingerprint"].str.len().eq(64).all()
    assert out["typing_mode"].eq("rules_only").all()
    assert out["typing_model_bundle_content_id"].eq("").all()
    assert out["typing_score_semantics"].eq(RULE_SCORE_SEMANTICS).all()
    assert out["typing_broad_score_semantics"].eq(RULE_SCORE_SEMANTICS).all()
    assert out["typing_broad_bootstrap_replicates"].eq(0).all()
    assert out["typing_broad_valid_bootstrap_replicates"].eq(0).all()
    assert out["typing_subtype_valid_bootstrap_replicates"].eq(0).all()


def test_typing_stage_uses_bundle_models_stability_and_frozen_provenance(
    monkeypatch,
):
    from phenocycler import typing_stage

    registry = load_registry()
    evidence = pd.DataFrame(
        {
            "donor_id": ["D", "D"],
            "object_id": ["cell", "outside-parent"],
            "qc_analysis_eligible": [True, True],
        }
    )
    calls: dict[str, object] = {}

    class Artifact:
        def validate_current(self, *, mode):
            calls["artifact_mode"] = mode

    class Bundle:
        thresholds = AssignmentThresholds(inferred_probability=0.75)
        broad = SimpleNamespace(
            parent=None,
            bootstrap_replicates=17,
            valid_bootstrap_replicates=15,
        )
        subtypes = (
            SimpleNamespace(
                parent="Immune",
                bootstrap_replicates=11,
                valid_bootstrap_replicates=9,
            ),
        )
        subtype_by_parent = {"Immune": subtypes[0]}
        input_artifact = Artifact()
        content_id = "b" * 64
        bundle_version = "typing-v1"
        method_version = "bundle-method-v1"
        feature_schema_version = "features-v1"
        score_semantics = "calibrated_probability"

        def validate_current(self, *, marker_registry, typing_registry):
            calls["validated"] = (marker_registry, typing_registry)

        def materialize_models(self, typing_registry):
            calls["materialized"] = typing_registry
            return "broad-model", {"Immune": "subtype-model"}

    bundle = Bundle()

    def fake_stability(frame, *, ensemble):
        calls.setdefault("stability_calls", []).append(
            (ensemble, frame["object_id"].tolist())
        )
        return frame[["donor_id", "object_id"]].assign(stability=1.0)

    def fake_assign_broad(frame, rules, **kwargs):
        calls["broad_assign"] = kwargs
        return pd.DataFrame(
            {
                "donor_id": ["D", "D"],
                "object_id": ["cell", "outside-parent"],
                "broad_label": ["Immune", "Epithelial"],
                "broad_assignment_status": ["inferred", "anchor"],
            }
        )

    def fake_assign_subtypes(frame, broad, rules, **kwargs):
        calls["subtype_assign"] = kwargs
        out = broad.copy()
        out["broad_type"] = out["broad_label"]
        out["broad_reason"] = [
            "posterior_margin_stability_pass",
            "single_authoritative_anchor",
        ]
        out["subtype_label"] = ["T_cell", ""]
        out["subtype_assignment_status"] = ["inferred", "not_applicable"]
        out["subtype_reason"] = [
            "posterior_margin_stability_pass",
            "terminal_parent",
        ]
        out["specific_type"] = ["T_cell", "Epithelial"]
        out["cell_type"] = out["specific_type"]
        out["assignment_status"] = ["inferred", "anchor"]
        out["specific_assignment_status"] = out["assignment_status"]
        out["confidence"] = [0.9, 0.99]
        out["broad_authoritative_markers"] = ["", "E_cadherin"]
        out["broad_anchor_classes"] = ["", "Epithelial"]
        return out

    monkeypatch.setattr(typing_stage, "ensemble_stability_wide", fake_stability)
    monkeypatch.setattr(typing_stage, "assign_broad_types", fake_assign_broad)
    monkeypatch.setattr(typing_stage, "assign_subtypes", fake_assign_subtypes)
    out = type_candidate_cells(
        evidence,
        model_bundle=bundle,
        marker_registry=registry,
    )

    assert calls["artifact_mode"] == "content"
    assert calls["stability_calls"] == [
        (bundle.broad, ["cell", "outside-parent"]),
        (bundle.subtypes[0], ["cell"]),
    ]
    broad_assigned = calls["broad_assign"]
    assert broad_assigned["model"] == "broad-model"
    assert broad_assigned["thresholds"] == bundle.thresholds
    expected_broad_stability = pd.DataFrame(
        {
            "donor_id": ["D", "D"],
            "object_id": ["cell", "outside-parent"],
            "stability": [1.0, 1.0],
        }
    )
    pd.testing.assert_frame_equal(
        broad_assigned["stability"], expected_broad_stability
    )
    subtype_assigned = calls["subtype_assign"]
    assert subtype_assigned["models"] == {"Immune": "subtype-model"}
    assert subtype_assigned["thresholds"] == bundle.thresholds
    expected_subtype_stability = pd.DataFrame(
        {"donor_id": ["D"], "object_id": ["cell"], "stability": [1.0]}
    )
    pd.testing.assert_frame_equal(
        subtype_assigned["stability"]["Immune"], expected_subtype_stability
    )
    assert out.loc[0, "typing_mode"] == "probabilistic_candidate"
    assert out.loc[0, "typing_model_bundle_content_id"] == bundle.content_id
    assert out.loc[0, "typing_model_bundle_version"] == bundle.bundle_version
    assert out.loc[0, "typing_model_method_version"] == bundle.method_version
    assert out.loc[0, "typing_feature_schema_version"] == bundle.feature_schema_version
    assert out.loc[0, "typing_score_semantics"] == bundle.score_semantics
    assert out.loc[0, "typing_broad_score_semantics"] == bundle.score_semantics
    assert out.loc[0, "typing_subtype_score_semantics"] == bundle.score_semantics
    assert out.loc[0, "typing_subtype_model_source"] == "probabilistic_bundle"
    assert out.loc[1, "typing_subtype_score_semantics"] == RULE_SCORE_SEMANTICS
    assert out.loc[1, "typing_subtype_model_source"] == "rules_fallback"
    assert out.loc[0, "typing_broad_bootstrap_replicates"] == 17
    assert out.loc[0, "typing_broad_valid_bootstrap_replicates"] == 15
    assert out.loc[0, "typing_subtype_bootstrap_replicates"] == 11
    assert out.loc[0, "typing_subtype_valid_bootstrap_replicates"] == 9
    assert out.loc[1, "typing_subtype_bootstrap_replicates"] == 0
    assert out.loc[1, "typing_subtype_valid_bootstrap_replicates"] == 0

    with pytest.raises(ValueError, match="owns its frozen thresholds"):
        type_candidate_cells(
            evidence,
            model_bundle=bundle,
            marker_registry=registry,
            thresholds=AssignmentThresholds(inferred_probability=0.70),
        )


def test_public_probabilistic_typing_rejects_candidate_bypass(monkeypatch):
    from phenocycler import typing_stage

    evidence = pd.DataFrame(
        {
            "donor_id": ["D"],
            "object_id": ["cell"],
            "qc_analysis_eligible": [True],
        }
    )
    bundle = SimpleNamespace(content_id="b" * 64)

    with pytest.raises(TypeError, match="verified TypingModelPromotionManifest"):
        type_calibrated_cells(evidence, model_bundle=bundle)
    with pytest.raises(TypeError, match="verified TypingModelPromotionManifest"):
        type_calibrated_cells(
            evidence,
            model_bundle=bundle,
            model_promotion_manifest=SimpleNamespace(
                validate_for_bundle=lambda _bundle: None
            ),
        )

    promotion = object.__new__(TypingModelPromotionManifest)
    object.__setattr__(promotion, "_signature_verified", True)
    calls: dict[str, object] = {}

    def validate_for_bundle(self, observed_bundle):
        calls["validated"] = (self, observed_bundle)

    def implementation(frame, **kwargs):
        calls["implementation"] = (frame, kwargs)
        return pd.DataFrame({"ok": [True]})

    monkeypatch.setattr(
        TypingModelPromotionManifest,
        "validate_for_bundle",
        validate_for_bundle,
    )
    monkeypatch.setattr(
        typing_stage,
        "_type_calibrated_cells_impl",
        implementation,
    )

    result = type_calibrated_cells(
        evidence,
        model_bundle=bundle,
        model_promotion_manifest=promotion,
    )

    assert result["ok"].tolist() == [True]
    assert calls["validated"] == (promotion, bundle)
    assert calls["implementation"][1]["typing_mode"] == "probabilistic_bundle"


def test_rules_only_rejects_orphan_promotion_manifest():
    promotion = object.__new__(TypingModelPromotionManifest)
    with pytest.raises(ValueError, match="cannot be used without its model bundle"):
        type_calibrated_cells(
            pd.DataFrame(),
            model_promotion_manifest=promotion,
        )
