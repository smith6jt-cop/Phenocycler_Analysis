"""Tests for deterministic audit sampling and error-aware validation."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from refinement import (
    POSTERIOR_PROBABILITY_SEMANTICS,
    challenge_flags,
    coverage_risk_curve,
    deterministic_stratified_sample,
    grouped_bootstrap_coverage_risk,
    grouped_bootstrap_metrics,
    one_vs_rest_reliability_table,
    per_class_metrics,
)


def _sampling_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": ["D1"] * 5 + ["D2"] * 3,
            "object_id": [f"c{index}" for index in range(8)],
            "review_group": ["accepted"] * 5 + ["ambiguous"] * 3,
            "value": np.arange(8),
        }
    )


def test_stratified_sample_is_deterministic_independent_of_input_order():
    source = _sampling_frame()
    expected = deterministic_stratified_sample(
        source, "review_group", 2, salt="locked-audit-v1"
    )
    shuffled = source.sample(frac=1, random_state=77).reset_index(drop=True)
    observed = deterministic_stratified_sample(
        shuffled, "review_group", 2, salt="locked-audit-v1"
    )

    pd.testing.assert_frame_equal(observed, expected)
    assert expected.groupby("sampling_stratum").size().to_dict() == {
        "accepted": 2,
        "ambiguous": 2,
    }
    audit = expected.groupby("sampling_stratum").first()
    assert audit.loc["accepted", "stratum_population"] == 5
    assert audit.loc["accepted", "stratum_sampled"] == 2
    assert audit.loc["accepted", "sampling_weight"] == pytest.approx(2.5)
    assert audit.loc["ambiguous", "sampling_weight"] == pytest.approx(1.5)


def test_stratified_sample_supports_per_stratum_sizes_and_complete_strata():
    sampled = deterministic_stratified_sample(
        _sampling_frame(),
        "review_group",
        {"accepted": 1, "ambiguous": 10},
        salt="locked-audit-v1",
    )

    assert sampled.groupby("sampling_stratum").size().to_dict() == {
        "accepted": 1,
        "ambiguous": 3,
    }
    ambiguous = sampled[sampled["sampling_stratum"] == "ambiguous"]
    assert (ambiguous["stratum_population"] == 3).all()
    assert (ambiguous["stratum_sampled"] == 3).all()
    assert (ambiguous["sampling_weight"] == 1.0).all()


@pytest.mark.parametrize(
    ("mutator", "error"),
    [
        (
            lambda frame: pd.concat([frame, frame.iloc[[0]]], ignore_index=True),
            "must be unique",
        ),
        (
            lambda frame: frame.assign(
                donor_id=frame["donor_id"].mask(frame.index == 0)
            ),
            "missing values",
        ),
        (
            lambda frame: frame.assign(
                review_group=frame["review_group"].mask(frame.index == 0, "")
            ),
            "blank values",
        ),
    ],
)
def test_stratified_sample_rejects_unstable_keys(mutator, error):
    with pytest.raises(ValueError, match=error):
        deterministic_stratified_sample(mutator(_sampling_frame()), "review_group", 2)


def test_per_class_metrics_counts_abstentions_as_false_negatives():
    metrics = per_class_metrics(
        reference=["A", "A", "B", "B"],
        prediction=["A", "Ambiguous", "A", "Unavailable"],
        weights=[1.0, 2.0, 3.0, 4.0],
    ).set_index("label")

    assert set(metrics.index) == {"A", "B"}
    assert "Ambiguous" not in metrics.index
    assert metrics.loc["A", ["TP", "FP", "FN", "TN"]].tolist() == [
        1.0,
        3.0,
        2.0,
        4.0,
    ]
    assert metrics.loc["A", "precision"] == pytest.approx(1 / 4)
    assert metrics.loc["A", "recall"] == pytest.approx(1 / 3)
    assert metrics.loc["A", "FPR"] == pytest.approx(3 / 7)
    assert metrics.loc["A", "FNR"] == pytest.approx(2 / 3)
    assert metrics.loc["A", "NPV"] == pytest.approx(4 / 6)
    assert metrics.loc["A", "F1"] == pytest.approx(2 / 7)
    assert metrics.loc["A", "support"] == 3.0
    assert metrics.loc["A", "calls"] == 4.0
    assert metrics.loc["B", "TP"] == 0.0
    assert metrics.loc["B", "FN"] == 7.0
    assert np.isnan(metrics.loc["B", "precision"])
    assert metrics.loc["B", "recall"] == 0.0
    assert metrics.loc["B", "F1"] == 0.0


def test_reliability_table_reports_weighted_bins_brier_and_ece():
    probabilities = pd.DataFrame(
        {
            "A": [0.9, 0.8, 0.4, 0.1],
            "B": [0.1, 0.2, 0.6, 0.9],
        }
    )
    table = one_vs_rest_reliability_table(
        reference=["A", "B", "A", "B"],
        probabilities=probabilities,
        weights=[1.0, 2.0, 1.0, 2.0],
        score_semantics=POSTERIOR_PROBABILITY_SEMANTICS,
        n_bins=2,
    )
    a_bins = table[table["label"] == "A"].set_index("bin_index")

    assert table.groupby("label").size().to_dict() == {"A": 2, "B": 2}
    assert a_bins.loc[0, "bin_weight"] == 3.0
    assert a_bins.loc[0, "mean_probability"] == pytest.approx(0.2)
    assert a_bins.loc[0, "observed_positive_fraction"] == pytest.approx(1 / 3)
    assert a_bins.loc[1, "mean_probability"] == pytest.approx(2.5 / 3)
    assert a_bins.loc[1, "observed_positive_fraction"] == pytest.approx(1 / 3)
    np.testing.assert_allclose(a_bins["brier_score"], 1.67 / 6)
    expected_ece = 0.5 * abs(1 / 3 - 0.2) + 0.5 * abs(1 / 3 - 2.5 / 3)
    np.testing.assert_allclose(a_bins["expected_calibration_error"], expected_ece)

    with_empty_bins = one_vs_rest_reliability_table(
        reference=["A", "B", "A", "B"],
        probabilities=probabilities,
        score_semantics=POSTERIOR_PROBABILITY_SEMANTICS,
        n_bins=4,
    )
    empty = with_empty_bins[
        (with_empty_bins["label"] == "A")
        & (with_empty_bins["bin_index"] == 2)
    ].iloc[0]
    assert empty["bin_count"] == 0
    assert empty["bin_weight"] == 0.0
    assert np.isnan(empty["mean_probability"])
    assert np.isnan(empty["observed_positive_fraction"])

    # A resolved reference outside the softmax vocabulary is a negative for
    # every modeled class, not a row to discard. This makes overconfident
    # forced-class calls on true Other cells visible in Brier/reliability.
    with_other = one_vs_rest_reliability_table(
        reference=["A", "Other"],
        probabilities={"A": [0.9, 0.8], "B": [0.1, 0.2]},
        score_semantics=POSTERIOR_PROBABILITY_SEMANTICS,
        n_bins=2,
    )
    assert with_other.loc[with_other["label"].eq("A"), "brier_score"].iloc[
        0
    ] == pytest.approx((0.1**2 + 0.8**2) / 2)


@pytest.mark.parametrize(
    "semantics",
    ["empirical_tail_probability", "margin", "logit", None],
)
def test_reliability_table_rejects_non_probability_score_semantics(semantics):
    with pytest.raises(
        ValueError,
        match="donor_class_balanced_temperature_scaled_softmax_probability",
    ):
        one_vs_rest_reliability_table(
            reference=["A", "B"],
            probabilities={"A": [0.8, 0.2], "B": [0.2, 0.8]},
            score_semantics=semantics,
        )


def test_reliability_table_requires_normalized_softmax_rows():
    with pytest.raises(ValueError, match="sum to one"):
        one_vs_rest_reliability_table(
            reference=["A", "Other"],
            probabilities={"A": [0.8, 0.7], "B": [0.8, 0.2]},
            score_semantics=POSTERIOR_PROBABILITY_SEMANTICS,
        )


def test_grouped_bootstrap_metrics_is_clustered_and_deterministic():
    kwargs = dict(
        reference=["A", "B", "A", "B", "A", "B"],
        prediction=["A", "A", "A", "B", "Ambiguous", "B"],
        groups=["D1", "D1", "D2", "D2", "D3", "D3"],
        weights=[1, 1, 2, 2, 3, 3],
        replicates=80,
        confidence=0.90,
        seed=42,
    )
    first = grouped_bootstrap_metrics(**kwargs)
    second = grouped_bootstrap_metrics(**kwargs)
    order = np.array([4, 1, 5, 0, 3, 2])
    shuffled_kwargs = {
        key: [value[index] for index in order]
        if key in {"reference", "prediction", "groups", "weights"}
        else value
        for key, value in kwargs.items()
    }
    shuffled = grouped_bootstrap_metrics(**shuffled_kwargs)

    pd.testing.assert_frame_equal(first, second)
    pd.testing.assert_frame_equal(first, shuffled)
    assert set(first["metric"]) == {
        "precision", "recall", "specificity", "NPV", "FPR", "FNR", "F1"
    }
    assert (first["groups"] == 3).all()
    assert (first["replicates"] == 80).all()
    assert first["valid_replicates"].between(1, 80).all()
    assert (first["ci_low"] <= first["ci_high"]).all()

    with pytest.raises(ValueError, match="at least two groups"):
        grouped_bootstrap_metrics(["A"], ["A"], ["D1"], replicates=2)


def test_per_class_metrics_never_creates_an_abstention_pseudo_class():
    metrics = per_class_metrics(
        ["A", "B"],
        ["Ambiguous", None],
        labels=["A", "B", "Ambiguous", "Unavailable"],
    )

    assert metrics["label"].tolist() == ["A", "B"]
    assert metrics.set_index("label")["FN"].to_dict() == {"A": 1.0, "B": 1.0}
    with pytest.raises(ValueError, match="reference contains abstention"):
        per_class_metrics(["A", "Ambiguous"], ["A", "A"])


def test_coverage_risk_curve_excludes_abstentions_and_tracks_error():
    curve = coverage_risk_curve(
        reference=["A", "A", "B", "B"],
        prediction=["A", "Ambiguous", "A", "B"],
        confidence=[0.90, 0.99, 0.85, 0.70],
        thresholds=[0.0, 0.8, 0.95],
    ).set_index("threshold")

    assert curve.loc[0.0, "coverage"] == pytest.approx(3 / 4)
    assert curve.loc[0.0, "correct_weight"] == 2.0
    assert curve.loc[0.0, "selective_error"] == pytest.approx(1 / 3)
    assert curve.loc[0.8, "coverage"] == pytest.approx(1 / 2)
    assert curve.loc[0.8, "selective_error"] == pytest.approx(1 / 2)
    assert curve.loc[0.95, "called_weight"] == 0.0
    assert np.isnan(curve.loc[0.95, "selective_error"])
    assert (curve["total_weight"] == 4.0).all()


def test_coverage_risk_default_thresholds_include_endpoints():
    curve = coverage_risk_curve(
        ["A", "B"],
        ["A", "B"],
        [0.2, np.nan],
    )

    assert curve["threshold"].tolist() == [0.0, 0.2, 1.0]
    assert curve.loc[curve["threshold"] == 0.0, "coverage"].iloc[0] == 0.5
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        coverage_risk_curve(["A"], ["A"], [1.2])


def test_coverage_risk_rethresholds_only_the_inference_lane():
    curve = coverage_risk_curve(
        reference=["A", "Other", "B", "B"],
        prediction=["A", "Other", "Ambiguous", "B"],
        confidence=[0.10, 0.05, 0.70, 0.80],
        thresholds=[0.65, 0.75],
        candidate_prediction=["A", "A", "B", "B"],
        threshold_eligible=[False, False, True, True],
    ).set_index("threshold")

    # Anchors and Other retain their dispositions, but Other is a reference
    # negative rather than a named call. Lowering the gate can rescue an
    # otherwise eligible candidate; raising it can abstain an inferred candidate.
    assert curve.loc[0.65, "coverage"] == pytest.approx(3 / 4)
    assert curve.loc[0.65, "selective_error"] == 0.0
    assert curve.loc[0.65, "rescue_called_weight"] == 2.0
    assert curve.loc[0.65, "rescue_coverage"] == 1.0
    assert curve.loc[0.75, "coverage"] == pytest.approx(1 / 2)
    assert curve.loc[0.75, "selective_error"] == 0.0
    assert curve.loc[0.75, "rescue_called_weight"] == 1.0
    assert curve.loc[0.75, "rescue_coverage"] == 0.5

    with pytest.raises(ValueError, match="provided together"):
        coverage_risk_curve(
            ["A"],
            ["Ambiguous"],
            [0.8],
            candidate_prediction=["A"],
        )


def test_coverage_risk_exposes_zero_rescue_calls_despite_fixed_anchor():
    curve = coverage_risk_curve(
        reference=["A", "Other", "B"],
        prediction=["A", "Other", "Ambiguous"],
        confidence=[0.99, 0.99, 0.40],
        thresholds=[0.80],
        candidate_prediction=["A", "A", "B"],
        threshold_eligible=[False, False, True],
    ).iloc[0]

    assert curve["called_weight"] == 1.0
    assert curve["coverage"] == pytest.approx(1 / 3)
    assert curve["rescue_eligible_weight"] == 1.0
    assert curve["rescue_called_weight"] == 0.0
    assert curve["rescue_coverage"] == 0.0
    assert np.isnan(curve["rescue_selective_error"])


def test_grouped_bootstrap_coverage_risk_is_clustered_and_row_order_invariant():
    kwargs = dict(
        reference=["A", "B", "A", "B", "A", "B", "A", "B"],
        prediction=["A", "A", "A", "B", "Ambiguous", "B", "B", "B"],
        score=[0.95, 0.90, 0.80, 0.75, 0.99, 0.60, 0.55, np.nan],
        groups=["D1", "D1", "D2", "D2", "D3", "D3", "D4", "D4"],
        thresholds=[0.0, 0.7, 0.9, 1.0],
        weights=[1, 1, 2, 2, 3, 3, 4, 4],
        replicates=60,
        confidence=0.90,
        seed=17,
    )
    first = grouped_bootstrap_coverage_risk(**kwargs)
    second = grouped_bootstrap_coverage_risk(**kwargs)
    order = np.array([6, 1, 4, 7, 2, 0, 5, 3])
    shuffled_kwargs = {
        key: [value[index] for index in order]
        if key in {"reference", "prediction", "score", "groups", "weights"}
        else value
        for key, value in kwargs.items()
    }
    shuffled = grouped_bootstrap_coverage_risk(**shuffled_kwargs)

    pd.testing.assert_frame_equal(first, second)
    pd.testing.assert_frame_equal(first, shuffled)
    point = coverage_risk_curve(
        kwargs["reference"],
        kwargs["prediction"],
        kwargs["score"],
        thresholds=kwargs["thresholds"],
        weights=kwargs["weights"],
    )
    pd.testing.assert_frame_equal(first[point.columns], point)
    assert (first["coverage_valid_replicates"] == 60).all()
    assert (first["groups"] == 4).all()
    assert (first["replicates"] == 60).all()
    assert (first["coverage_ci_low"] <= first["coverage_ci_high"]).all()
    no_calls = first[first["threshold"] == 1.0].iloc[0]
    assert no_calls["selective_error_valid_replicates"] == 0
    assert np.isnan(no_calls["selective_error_ci_low"])
    assert np.isnan(no_calls["selective_error_ci_high"])

    with pytest.raises(ValueError, match="at least two groups"):
        grouped_bootstrap_coverage_risk(
            ["A", "B"],
            ["A", "B"],
            [0.9, 0.8],
            ["D1", "D1"],
            replicates=2,
        )

    lane_kwargs = dict(
        reference=["A", "Other", "B", "B"],
        prediction=["A", "Other", "Ambiguous", "B"],
        score=[0.1, 0.1, 0.7, 0.8],
        groups=["D1", "D1", "D2", "D2"],
        thresholds=[0.65, 0.75],
        candidate_prediction=["A", "A", "B", "B"],
        threshold_eligible=[False, False, True, True],
        replicates=10,
        seed=5,
    )
    lane = grouped_bootstrap_coverage_risk(**lane_kwargs)
    lane_point = coverage_risk_curve(
        lane_kwargs["reference"],
        lane_kwargs["prediction"],
        lane_kwargs["score"],
        thresholds=lane_kwargs["thresholds"],
        candidate_prediction=lane_kwargs["candidate_prediction"],
        threshold_eligible=lane_kwargs["threshold_eligible"],
    )
    pd.testing.assert_frame_equal(lane[lane_point.columns], lane_point)


def test_challenge_flags_cover_fp_fn_modes_and_keep_random_baseline():
    assignments = pd.DataFrame(
        {
            "object_id": [
                "low-margin",
                "contradiction",
                "high-abstention",
                "other-candidate",
                "clean",
                "qc-excluded",
            ],
            "broad_assignment_status": [
                "inferred",
                "anchor",
                "ambiguous",
                "other",
                "anchor",
                "unavailable",
            ],
            "broad_margin": [0.10, 0.50, 0.15, 0.0, 0.60, np.nan],
            "broad_best_probability": [0.90, 0.95, 0.91, 0.20, 0.99, 0.30],
            "broad_best_candidate": ["A", "A", "B", "B", "A", ""],
            "broad_reason": [
                "posterior_margin_stability_pass",
                "anchor_exclusion_conflict",
                "posterior_margin_or_stability_failed",
                "all_required_models_valid_and_negative",
                "single_authoritative_anchor",
                "geometry_qc_ineligible",
            ],
            "broad_anchor_classes": ["", "A|B", "", "", "A", ""],
            "qc_analysis_eligible": [True, True, True, True, True, False],
        }
    ).set_index("object_id")

    flagged = challenge_flags(assignments)

    assert flagged["review_strata"].str.startswith("random_baseline").all()
    assert flagged.loc["low-margin", "candidate_fp"]
    assert "accepted_low_margin" in flagged.loc["low-margin", "review_strata"]
    assert flagged.loc["contradiction", "candidate_fp"]
    assert flagged.loc["contradiction", "candidate_fn"]
    assert "status_reason_contradiction" in flagged.loc[
        "contradiction", "review_strata"
    ]
    assert "anchor_conflict" in flagged.loc["contradiction", "review_strata"]
    assert flagged.loc["high-abstention", "candidate_fn"]
    assert "abstained_high_probability" in flagged.loc[
        "high-abstention", "review_strata"
    ]
    assert flagged.loc["other-candidate", "candidate_fn"]
    assert "other_with_candidate" in flagged.loc[
        "other-candidate", "review_strata"
    ]
    assert not flagged.loc["clean", "candidate_fp"]
    assert not flagged.loc["clean", "candidate_fn"]
    assert flagged.loc["clean", "review_strata"] == "random_baseline"
    assert "geometry_qc_ineligible" in flagged.loc[
        "qc-excluded", "review_strata"
    ]
    assert flagged.loc["qc-excluded", "candidate_fn"]


def test_challenge_flags_detect_parent_constrained_unclassified_cells():
    assignments = pd.DataFrame(
        {
            "assignment_status": ["other"],
            "confidence": [0.85],
            "best_specific_type": ["Beta"],
            "specific_type": ["Epithelial-unclassified"],
        }
    )
    flagged = challenge_flags(assignments)

    assert flagged.loc[0, "candidate_fn"]
    assert flagged.loc[0, "review_strata"] == (
        "random_baseline;other_with_candidate;subtype_unresolved"
    )


def test_combined_specific_status_does_not_misread_terminal_parent_reason():
    flagged = challenge_flags(
        pd.DataFrame(
            {
                "specific_assignment_status": ["anchor"],
                "subtype_reason": ["terminal_parent"],
                "specific_type": ["Neural"],
            }
        )
    )

    assert not flagged.loc[0, "candidate_fp"]
    assert not flagged.loc[0, "candidate_fn"]
    assert flagged.loc[0, "review_strata"] == "random_baseline"


def test_challenge_flags_do_not_mutate_assignments():
    assignments = pd.DataFrame(
        {"assignment_status": ["anchor"], "margin": [0.5]}
    )
    original = assignments.copy(deep=True)
    challenge_flags(assignments)
    pd.testing.assert_frame_equal(assignments, original)
