from __future__ import annotations

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
from phenocycler.probabilistic_training import (
    _fit_ensemble,
    acceptance_operating_points,
    choose_operating_point,
    fit_probabilistic_typing_bundle,
    fit_temperature_scaling,
    load_labeled_evidence,
    main,
)
from phenocycler.refinement_contracts import (
    DEVELOPMENT_LABEL_PURPOSE,
    DEVELOPMENT_LABEL_SOURCE,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
    evidence_frame_fingerprint,
)
from phenocycler.typing_stage import type_candidate_cells

REGISTRY = TypingRegistry.from_mapping(
    {
        "marker_roles": {"A": "identity", "B": "identity"},
        "broad": {
            "TypeA": {
                "anchors": {"A": 4.0},
                "exclusions": {"B": 2.0},
                "required": ["A", "B"],
            },
            "TypeB": {
                "anchors": {"B": 4.0},
                "exclusions": {"A": 2.0},
                "required": ["A", "B"],
            },
        },
    }
)


def _fit_cli_required_arguments() -> list[str]:
    return [
        "fit",
        "--evidence",
        "evidence",
        "--labels",
        "labels.parquet",
        "--splits",
        "splits.json",
        "--source-run-manifest",
        "run.json",
        "--label-provenance",
        "labels.json",
        "--evidence-provenance",
        "evidence.json",
        "--bundle-version",
        "candidate-v1",
        "--out",
        "candidate.json",
    ]


def test_fit_cli_only_builds_an_explicit_unfrozen_threshold_candidate(capsys):
    required = _fit_cli_required_arguments()
    with pytest.raises(SystemExit):
        main([*required, "--unfrozen-candidate", "--probability-threshold", "0.8"])
    assert "requires --probability-threshold" in capsys.readouterr().err

    with pytest.raises(SystemExit):
        main(
            [
                *required,
                "--threshold-provenance",
                "thresholds.json",
            ]
        )
    error = capsys.readouterr().err
    assert "--unfrozen-candidate" in error or "unrecognized arguments" in error


def test_probabilistic_cli_dispatches_freeze_without_entering_fit(monkeypatch):
    from phenocycler import probabilistic_training

    observed = []
    monkeypatch.setattr(
        probabilistic_training,
        "_freeze_main",
        lambda arguments: observed.append(arguments) or 23,
    )

    assert probabilistic_training.main(["freeze", "--model-bundle", "x"]) == 23
    assert observed == [["--model-bundle", "x"]]


def _make_evidence_provenance(evidence, splits):
    return DevelopmentEvidenceProvenance(
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        source_run_manifest_content_id="r" * 64,
        source_stage_manifest_content_id="s" * 64,
        source_artifact_content_id="e" * 64,
        labeled_frame_fingerprint=evidence_frame_fingerprint(evidence),
        feature_schema_version=FEATURE_SCHEMA_VERSION,
        row_count=len(evidence),
    )


def _fixture():
    evidence_rows = []
    label_rows = []
    for donor in ("D1", "D2", "D3"):
        samples = (
            (
                ("TypeA", "good", 0.80, 0.20),
                ("TypeA", "challenge", 0.35, 0.65),
                ("TypeB", "good", 0.20, 0.80),
                ("TypeB", "moderate", 0.30, 0.70),
            )
            if donor == "D3"
            else (
                ("TypeA", "seed", 0.98, 0.0),
                ("TypeB", "seed", 0.0, 0.98),
            )
        )
        for label, suffix, a, b in samples:
            object_id = f"{donor}_{label}_{suffix}"
            evidence_rows.extend(
                [
                    {
                        "donor_id": donor,
                        "object_id": object_id,
                        "marker": "A",
                        "evidence_probability": a,
                    },
                    {
                        "donor_id": donor,
                        "object_id": object_id,
                        "marker": "B",
                        "evidence_probability": b,
                    },
                ]
            )
            label_rows.append(
                {
                    "donor_id": donor,
                    "object_id": object_id,
                    "level": "broad",
                    "parent": "",
                    "label": label,
                    "evidence_sufficient": True,
                    "audit_only": False,
                    "label_purpose": DEVELOPMENT_LABEL_PURPOSE,
                    "label_source": DEVELOPMENT_LABEL_SOURCE,
                    "reviewer_1": "reviewer-A",
                    "reviewer_2": "reviewer-B",
                    "adjudicated": True,
                    "root_cause": "independent_image_review",
                }
            )
    evidence = pd.DataFrame(evidence_rows)
    labels = pd.DataFrame(label_rows)
    splits = DonorSplitManifest(
        split_version="split-v1",
        source_run_id="source-run",
        cohort_content_id="cohort",
        development=("D1", "D2"),
        calibration=("D3",),
    )
    provenance = DevelopmentLabelProvenance.from_labels(
        labels,
        bundle_id="new-development-v1",
        splits=splits,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    evidence_provenance = _make_evidence_provenance(evidence, splits)
    marker_registry = SimpleNamespace(fingerprint="a" * 64)
    return (
        evidence,
        labels,
        splits,
        provenance,
        evidence_provenance,
        marker_registry,
    )


def test_donor_disjoint_fit_builds_calibrated_bootstrap_bundle():
    evidence, labels, splits, provenance, evidence_provenance, marker_registry = (
        _fixture()
    )
    bundle = fit_probabilistic_typing_bundle(
        evidence,
        labels,
        marker_registry=marker_registry,
        typing_registry=REGISTRY,
        splits=splits,
        label_provenance=provenance,
        evidence_provenance=evidence_provenance,
        bundle_version="candidate-v1",
        thresholds=AssignmentThresholds(
            inferred_probability=0.75,
            inferred_margin=0.20,
            inferred_stability=0.90,
        ),
        n_bootstrap=3,
        random_state=9,
        maxiter=200,
        created_at_utc="2026-08-01T00:00:00+00:00",
    )
    assert bundle.broad.bootstrap_replicates == 3
    assert bundle.broad.training_donors == ("D1", "D2")
    assert bundle.broad.calibration_donors == ("D3",)
    assert np.isfinite(bundle.broad.temperature)
    assert bundle.broad.temperature > 0
    assert bundle.thresholds.inferred_probability == 0.75
    assert bundle.label_provenance.fingerprint == provenance.fingerprint
    assert bundle.fit_parameters["class_priors"] == "uniform_no_intercept"

    production = pd.DataFrame(
        {
            "donor_id": ["NEW"],
            "object_id": ["rescue"],
            "qc_analysis_eligible": [True],
            # Strong enough for the fitted model, but below the separate 0.95
            # authoritative numeric-pilot anchor gate.
            "A__evidence_probability": [0.90],
            "A__model_valid": [True],
            "B__evidence_probability": [0.0],
            "B__model_valid": [True],
        }
    )
    assignment = type_candidate_cells(
        production,
        model_bundle=bundle,
        marker_registry=marker_registry,
        typing_registry=REGISTRY,
    ).iloc[0]
    assert assignment["broad_label"] == "TypeA"
    assert assignment["broad_assignment_status"] == "inferred"
    assert assignment["typing_mode"] == "probabilistic_candidate"
    assert assignment["broad_valid_evidence_pass"]


def test_locked_or_validation_labels_are_rejected_before_fit():
    evidence, labels, splits, _provenance, _evidence_provenance, marker_registry = (
        _fixture()
    )
    split_with_locked = DonorSplitManifest(
        split_version="split-v1",
        source_run_id="source-run",
        cohort_content_id="cohort",
        development=("D1", "D2"),
        calibration=("D3",),
        locked_test=("D4",),
    )
    leaked = pd.concat(
        [
            labels,
            pd.DataFrame(
                {
                    "donor_id": ["D4"],
                    "object_id": ["locked"],
                    "level": ["broad"],
                    "parent": [""],
                    "label": ["TypeA"],
                    "evidence_sufficient": [True],
                    "audit_only": [False],
                    "label_purpose": [DEVELOPMENT_LABEL_PURPOSE],
                    "label_source": [DEVELOPMENT_LABEL_SOURCE],
                    "reviewer_1": ["reviewer-A"],
                    "reviewer_2": ["reviewer-B"],
                    "adjudicated": [True],
                    "root_cause": ["independent_image_review"],
                }
            ),
        ],
        ignore_index=True,
    )
    provenance = DevelopmentLabelProvenance.from_labels(
        leaked,
        bundle_id="leaked",
        splits=split_with_locked,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    evidence_provenance = _make_evidence_provenance(evidence, split_with_locked)
    with pytest.raises(ValueError, match="locked-test"):
        fit_probabilistic_typing_bundle(
            evidence,
            leaked,
            marker_registry=marker_registry,
            typing_registry=REGISTRY,
            splits=split_with_locked,
            label_provenance=provenance,
            evidence_provenance=evidence_provenance,
            bundle_version="bad",
            thresholds=AssignmentThresholds(),
            n_bootstrap=1,
        )


def test_each_development_class_requires_two_donor_sources():
    evidence, labels, splits, _labels_provenance, _evidence_provenance, registry = (
        _fixture()
    )
    sparse = labels.loc[
        ~(
            labels["donor_id"].eq("D2")
            & labels["label"].eq("TypeB")
        )
    ].copy()
    sparse_keys = sparse[["donor_id", "object_id"]].drop_duplicates()
    sparse_evidence = evidence.merge(
        sparse_keys,
        on=["donor_id", "object_id"],
        how="inner",
        validate="m:1",
    )
    label_provenance = DevelopmentLabelProvenance.from_labels(
        sparse,
        bundle_id="under-supported",
        splits=splits,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    evidence_provenance = _make_evidence_provenance(sparse_evidence, splits)
    with pytest.raises(ValueError, match="at least two donors"):
        fit_probabilistic_typing_bundle(
            sparse_evidence,
            sparse,
            marker_registry=registry,
            typing_registry=REGISTRY,
            splits=splits,
            label_provenance=label_provenance,
            evidence_provenance=evidence_provenance,
            bundle_version="bad-support",
            thresholds=AssignmentThresholds(),
            n_bootstrap=1,
        )


def test_whole_donor_bootstrap_keeps_incomplete_draws_as_instability(monkeypatch):
    splits = DonorSplitManifest(
        split_version="bootstrap-v1",
        source_run_id="source-run",
        cohort_content_id="cohort",
        development=("D1", "D2", "D3"),
        calibration=("D4",),
    )
    labels = pd.DataFrame(
        [
            ("D1", "D1_A", "TypeA"),
            ("D1", "D1_B", "TypeB"),
            ("D2", "D2_A", "TypeA"),
            ("D2", "D2_B", "TypeB"),
            ("D3", "D3_A", "TypeA"),
            ("D4", "D4_A", "TypeA"),
            ("D4", "D4_B", "TypeB"),
        ],
        columns=["donor_id", "object_id", "label"],
    )
    fitted = classifier_from_registry(REGISTRY)
    monkeypatch.setattr(
        "phenocycler.probabilistic_training.fit_constrained_classifier",
        lambda *args, **kwargs: fitted,
    )
    monkeypatch.setattr(
        "phenocycler.probabilistic_training.fit_temperature_scaling",
        lambda *args, **kwargs: 1.0,
    )

    class PredeterminedDraws:
        def __init__(self):
            self.draws = iter(
                [
                    ["D3", "D3", "D3"],
                    ["D1", "D2", "D3"],
                    ["D2", "D1", "D3"],
                    ["D1", "D1", "D3"],
                    ["D2", "D2", "D3"],
                ]
            )

        def choice(self, *_args, **_kwargs):
            return np.asarray(next(self.draws), dtype=object)

    ensemble = _fit_ensemble(
        pd.DataFrame(),
        REGISTRY,
        labels,
        splits,
        level="broad",
        parent=None,
        n_bootstrap=5,
        rng=PredeterminedDraws(),
        l2=0.02,
        max_abs_weight=8.0,
        maxiter=100,
    )
    assert ensemble.bootstrap_replicates == 5
    assert ensemble.valid_bootstrap_replicates == 4
    assert ensemble.bootstrap_valid == (False, True, True, True, True)
    assert ensemble.bootstrap_failures[0] == "missing_classes:TypeB"


def test_operating_grid_balances_fp_fn_without_prevalence_target():
    probabilities = pd.DataFrame(
        {
            "TypeA": [0.95, 0.55, 0.10, 0.45],
            "TypeB": [0.05, 0.45, 0.90, 0.55],
        }
    )
    reference = ["TypeA", "TypeA", "TypeB", "TypeA"]
    stability = [1.0, 0.95, 1.0, 0.95]
    table = acceptance_operating_points(
        reference,
        probabilities,
        ("TypeA", "TypeB"),
        stability,
        probability_thresholds=(0.5, 0.8),
        margin_thresholds=(0.0, 0.2),
        stability_thresholds=(0.9,),
        false_positive_cost=2.0,
        false_negative_cost=1.0,
    )
    assert len(table) == 4
    permissive = table.query(
        "probability_threshold == 0.5 and margin_threshold == 0.0"
    ).iloc[0]
    strict = table.query(
        "probability_threshold == 0.8 and margin_threshold == 0.2"
    ).iloc[0]
    assert permissive["coverage"] > strict["coverage"]
    assert permissive["false_positive_weight"] > strict["false_positive_weight"]
    chosen = choose_operating_point(
        table,
        maximum_selective_error=0.01,
        minimum_coverage=0.4,
    )
    assert chosen["selective_error"] <= 0.01


def test_operating_grid_reproduces_fixed_production_dispositions_and_other():
    table = acceptance_operating_points(
        reference=["TypeA", "Other", "TypeB", "TypeB"],
        probabilities=pd.DataFrame(
            {
                "TypeA": [0.95, 0.85, 0.20, 0.10],
                "TypeB": [0.05, 0.15, 0.80, 0.90],
            }
        ),
        classes=("TypeA", "TypeB"),
        stability=[1.0, 1.0, 1.0, 1.0],
        probability_thresholds=(0.75,),
        margin_thresholds=(0.20,),
        stability_thresholds=(0.90,),
        current_prediction=["TypeA", "Other", "Ambiguous", "Ambiguous"],
        current_assignment_status=["anchor", "other", "ambiguous", "ambiguous"],
        current_reason=[
            "single_authoritative_anchor",
            "all_required_models_valid_and_negative",
            "posterior_margin_or_stability_failed",
            "posterior_margin_or_stability_failed",
        ],
        candidate_prediction=["TypeA", "TypeA", "TypeB", "TypeB"],
        threshold_eligible=[False, False, True, True],
        valid_evidence=[True, True, False, True],
    ).iloc[0]

    # The anchor is fixed, Other remains a non-call/negative, one otherwise
    # eligible cell fails the invariant evidence gate, and one is rescued.
    assert table["coverage"] == pytest.approx(2 / 4)
    assert table["false_positive_weight"] == 0.0
    assert table["false_negative_weight"] == 1.0
    assert table["rescue_eligible_weight"] == 2.0
    assert table["rescue_called_weight"] == 1.0
    assert table["rescue_coverage"] == 0.5


def test_operating_grid_rejects_a_gamed_fixed_inference_lane():
    kwargs = {
        "reference": ["TypeA", "TypeB"],
        "probabilities": pd.DataFrame(
            {"TypeA": [0.9, 0.2], "TypeB": [0.1, 0.8]}
        ),
        "classes": ("TypeA", "TypeB"),
        "stability": [1.0, 1.0],
        "current_prediction": ["TypeA", "Ambiguous"],
        "current_assignment_status": ["inferred", "ambiguous"],
        "current_reason": [
            "posterior_margin_stability_pass",
            "posterior_margin_or_stability_failed",
        ],
        "candidate_prediction": ["TypeA", "TypeB"],
        "valid_evidence": [True, True],
    }

    with pytest.raises(ValueError, match="status/reason contract"):
        acceptance_operating_points(
            **kwargs,
            threshold_eligible=[False, True],
        )


def test_operating_point_can_require_rescue_lane_performance():
    table = pd.DataFrame(
        {
            "probability_threshold": [0.9, 0.8],
            "margin_threshold": [0.2, 0.2],
            "stability_threshold": [0.9, 0.9],
            "coverage": [0.5, 0.6],
            "selective_error": [0.0, 0.05],
            "rescue_coverage": [0.0, 0.4],
            "rescue_selective_error": [np.nan, 0.1],
            "error_cost": [0.0, 0.1],
        }
    )

    selected = choose_operating_point(
        table,
        maximum_selective_error=0.1,
        minimum_coverage=0.4,
        maximum_rescue_selective_error=0.15,
        minimum_rescue_coverage=0.25,
    )

    assert selected["probability_threshold"] == 0.8
    with pytest.raises(ValueError, match="no operating point"):
        choose_operating_point(
            table,
            maximum_selective_error=0.1,
            maximum_rescue_selective_error=0.05,
            minimum_rescue_coverage=0.25,
        )


def test_probability_matrix_must_be_normalized():
    with pytest.raises(ValueError, match="sum to one"):
        acceptance_operating_points(
            ["TypeA"],
            np.array([[0.9, 0.9]]),
            ("TypeA", "TypeB"),
            [1.0],
        )


def test_temperature_scaling_rejects_unidentified_uniform_model():
    evidence, labels, *_ = _fixture()
    calibration = labels.loc[
        labels["donor_id"].eq("D3"),
        ["donor_id", "object_id", "label"],
    ]
    uniform_model = classifier_from_registry(REGISTRY)
    uniform_model.coefficients[:] = 0.0
    with pytest.raises(ValueError, match="unidentified"):
        fit_temperature_scaling(uniform_model, evidence, calibration)


def test_temperature_scaling_uses_neutral_value_for_one_class_parent():
    registry = TypingRegistry.from_mapping(
        {
            "marker_roles": {"A": "identity"},
            "broad": {"Only": {"anchors": {"A": 4.0}, "required": ["A"]}},
        }
    )
    model = classifier_from_registry(registry)
    evidence = pd.DataFrame(
        {
            "donor_id": ["D3"],
            "object_id": ["only"],
            "marker": ["A"],
            "evidence_probability": [0.9],
        }
    )
    labels = pd.DataFrame(
        {"donor_id": ["D3"], "object_id": ["only"], "label": ["Only"]}
    )
    assert fit_temperature_scaling(model, evidence, labels) == 1.0


def test_temperature_scaling_rejects_separation_without_finite_optimum():
    evidence, labels, splits, *_ = _fixture()
    development = labels.loc[
        labels["donor_id"].isin(splits.development),
        ["donor_id", "object_id", "label"],
    ]
    model = classifier_from_registry(REGISTRY)
    # Registry coefficients already rank these two perfectly. A calibration
    # set containing no countervailing errors has an optimum at T -> 0.
    separated_evidence = pd.DataFrame(
        [
            {"donor_id": "D3", "object_id": "A", "marker": "A", "evidence_probability": 0.99},
            {"donor_id": "D3", "object_id": "A", "marker": "B", "evidence_probability": 0.0},
            {"donor_id": "D3", "object_id": "B", "marker": "A", "evidence_probability": 0.0},
            {"donor_id": "D3", "object_id": "B", "marker": "B", "evidence_probability": 0.99},
        ]
    )
    separated_labels = pd.DataFrame(
        {
            "donor_id": ["D3", "D3"],
            "object_id": ["A", "B"],
            "label": ["TypeA", "TypeB"],
        }
    )
    assert len(development) > 0
    with pytest.raises(ValueError, match="boundary|unidentified"):
        fit_temperature_scaling(model, separated_evidence, separated_labels)


def test_labeled_evidence_loader_projects_partitioned_parquet(tmp_path):
    root = tmp_path / "marker_evidence"
    for donor in ("D1", "D2"):
        partition = root / f"donor_id={donor}"
        partition.mkdir(parents=True)
        pd.DataFrame(
            {
                "object_id": ["keep", "drop"],
                "A__state": ["positive", "negative"],
            }
        ).to_parquet(partition / "part.parquet", index=False)
    labels = pd.DataFrame(
        {"donor_id": ["D1", "D2"], "object_id": ["keep", "keep"]}
    )
    observed = load_labeled_evidence(root, labels)
    assert set(map(tuple, observed[["donor_id", "object_id"]].to_numpy())) == {
        ("D1", "keep"),
        ("D2", "keep"),
    }
    assert len(observed) == 2

    missing = pd.concat(
        [labels, pd.DataFrame({"donor_id": ["D2"], "object_id": ["absent"]})],
        ignore_index=True,
    )
    with pytest.raises(ValueError, match="no marker evidence"):
        load_labeled_evidence(root, missing)
