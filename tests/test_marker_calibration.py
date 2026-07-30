from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.marker_calibration import (
    MAX_CONTROL_CAP,
    CalibrationConfig,
    CalibrationStatus,
    EvidenceState,
    apply_marker_calibration,
    calibrate_donor,
    calibrate_marker_frame,
    empirical_tail_probability,
    evidence_from_wide,
    evidence_to_wide,
    fit_marker_calibration,
    models_to_long,
    spatially_balanced_control_indices,
)
from phenocycler.marker_registry import load_registry


def _config(**overrides):
    values = {
        "control_cap": 220,
        "min_controls": 100,
        "spatial_bins": 4,
        "bootstrap_replicates": 40,
        "random_seed": 17,
    }
    values.update(overrides)
    return CalibrationConfig(**values)


def _cells(*, contamination_n=0, noncontrol_positive_n=50, seed=4):
    rng = np.random.default_rng(seed)
    n_controls = 200
    control_log = rng.normal(np.log1p(5.0), 0.08, n_controls)
    if contamination_n:
        control_log[:contamination_n] = np.log1p(80.0)
    controls = np.expm1(control_log)
    positives = np.expm1(
        rng.normal(np.log1p(60.0), 0.08, noncontrol_positive_n)
    )
    values = np.r_[controls, positives]
    n = len(values)
    return {
        "values": values,
        "reference": np.arange(n) < n_controls,
        "eligible": np.ones(n, dtype=bool),
        "x": (np.arange(n) % 25).astype(float),
        "y": (np.arange(n) // 25).astype(float),
        "ids": np.array([f"cell-{index:04d}" for index in range(n)]),
    }


def _fit(data, **kwargs):
    return fit_marker_calibration(
        data["values"],
        data["reference"],
        data["eligible"],
        data["x"],
        data["y"],
        data["ids"],
        donor_id="donor-A",
        marker="B3TUBB",
        registry=load_registry(),
        config=_config(),
        **kwargs,
    )


def test_config_locks_default_control_minimum_and_hard_cap():
    config = CalibrationConfig()
    assert config.min_controls == 10_000
    assert config.control_cap == MAX_CONTROL_CAP == 20_000
    with pytest.raises(ValueError, match="control_cap"):
        CalibrationConfig(control_cap=MAX_CONTROL_CAP + 1)
    with pytest.raises(ValueError, match="min_controls"):
        CalibrationConfig(control_cap=100, min_controls=101)


def test_spatial_selection_is_deterministic_row_order_invariant_and_balanced():
    # Four tiles, with tile 0 much denser than the others.
    x = np.r_[np.linspace(0, 0.4, 40), np.linspace(9.6, 10, 8),
              np.linspace(0, 0.4, 8), np.linspace(9.6, 10, 8)]
    y = np.r_[np.linspace(0, 0.4, 40), np.linspace(0, 0.4, 8),
              np.linspace(9.6, 10, 8), np.linspace(9.6, 10, 8)]
    ids = np.array([f"id-{i:03d}" for i in range(len(x))])
    candidates = np.ones(len(x), bool)
    selected = spatially_balanced_control_indices(
        x, y, candidates, ids, cap=12, spatial_bins=2, salt="fixed"
    )
    selected_ids = ids[selected]
    tiles = (x[selected] >= 5).astype(int) + 2 * (y[selected] >= 5)
    assert np.bincount(tiles, minlength=4).tolist() == [3, 3, 3, 3]

    permutation = np.random.default_rng(9).permutation(len(x))
    permuted = spatially_balanced_control_indices(
        x[permutation],
        y[permutation],
        candidates[permutation],
        ids[permutation],
        cap=12,
        spatial_bins=2,
        salt="fixed",
    )
    assert ids[permutation][permuted].tolist() == selected_ids.tolist()


def test_fit_uses_log_median_qn_empirical_tail_and_bootstrap_ci():
    data = _cells(contamination_n=4)
    model = _fit(data)
    assert model.status is CalibrationStatus.VALID
    assert model.n_controls_selected == 200
    assert model.n_controls_contaminated == 4
    assert model.contamination_fraction == pytest.approx(0.02)
    assert model.n_controls_clean == 196
    assert model.background_log1p_median == pytest.approx(
        np.median(np.log1p(data["values"][4:200]))
    )
    assert model.background_log1p_qn > 0
    assert model.threshold_ci_low_raw <= model.threshold_raw
    assert model.threshold_raw <= model.threshold_ci_high_raw
    assert model.marker_alpha == pytest.approx(0.01)

    evidence = apply_marker_calibration(
        data["values"], model, object_ids=data["ids"]
    )
    assert set(evidence.loc[200:, "state"]) == {EvidenceState.POSITIVE.value}
    assert (
        evidence.loc[200:, "empirical_tail_probability"] <= model.marker_alpha
    ).all()
    assert evidence.loc[200:, "log2_threshold_ratio"].gt(0).all()


def test_more_than_five_percent_contamination_invalidates_model_and_calls():
    data = _cells(contamination_n=11)
    model = _fit(data)
    assert model.status is CalibrationStatus.EXCESS_CONTROL_CONTAMINATION
    assert model.contamination_fraction == pytest.approx(11 / 200)

    values = data["values"].copy()
    values[-1] = np.nan
    evidence = apply_marker_calibration(values, model, object_ids=data["ids"])
    assert set(evidence.iloc[:-1]["state"]) == {
        EvidenceState.INDETERMINATE.value
    }
    assert evidence.iloc[-1]["state"] == EvidenceState.UNAVAILABLE.value


def test_estimation_eligibility_is_explicit_but_application_uses_all_cells():
    data = _cells()
    # A bright target value is reference-positive but explicitly ineligible.
    data["values"][0] = 500.0
    data["eligible"][0] = False
    model = _fit(data)
    assert model.status is CalibrationStatus.VALID
    assert data["ids"][0] not in model.control_object_ids
    assert model.n_estimation_eligible == len(data["values"]) - 1

    evidence = apply_marker_calibration(
        data["values"], model, object_ids=data["ids"]
    ).set_index("object_id")
    assert evidence.loc[data["ids"][0], "state"] == EvidenceState.POSITIVE.value


def test_noncontrol_target_frequency_cannot_change_calibration():
    sparse = _cells(noncontrol_positive_n=5)
    abundant = _cells(noncontrol_positive_n=300)
    first = _fit(sparse)
    second = _fit(abundant)
    assert first.status is second.status is CalibrationStatus.VALID
    assert first.control_object_ids == second.control_object_ids
    assert first.background_log1p_median == second.background_log1p_median
    assert first.background_log1p_qn == second.background_log1p_qn
    assert first.threshold_log1p == second.threshold_log1p
    assert first.threshold_ci_low_log1p == second.threshold_ci_low_log1p
    assert first.threshold_ci_high_log1p == second.threshold_ci_high_log1p


def test_log_transform_is_equivariant_to_multiplicative_measurement_scale():
    data = _cells()
    original = _fit(data)
    scaled = {key: value.copy() for key, value in data.items()}
    scaled["values"] = 9.0 * (data["values"] + 1.0) - 1.0
    transformed = _fit(scaled)
    expected_shift = np.log(9.0)
    assert transformed.threshold_log1p == pytest.approx(
        original.threshold_log1p + expected_shift
    )
    assert transformed.background_log1p_median == pytest.approx(
        original.background_log1p_median + expected_shift
    )

    first = apply_marker_calibration(data["values"], original)
    second = apply_marker_calibration(scaled["values"], transformed)
    assert first["state"].tolist() == second["state"].tolist()
    assert np.allclose(
        first["log2_threshold_ratio"],
        second["log2_threshold_ratio"],
        equal_nan=True,
    )


def test_insufficient_controls_are_status_not_biological_negative():
    data = _cells()
    data["reference"][:] = False
    data["reference"][:20] = True
    model = _fit(data)
    assert model.status is CalibrationStatus.INSUFFICIENT_CONTROLS
    evidence = apply_marker_calibration(data["values"], model)
    assert set(evidence["state"]) == {EvidenceState.INDETERMINATE.value}


def test_frame_api_and_long_wide_helpers_preserve_evidence():
    data = _cells()
    registry = load_registry()
    frame = pd.DataFrame(
        {
            "object_id": data["ids"],
            "X_centroid": data["x"],
            "Y_centroid": data["y"],
            "qc_estimation_eligible": data["eligible"],
            "reference_positive": data["reference"],
            "Cell__B3TUBB": data["values"],
            # A second marker with the same synthetic distribution exercises
            # donor-long and wide products without changing scientific policy.
            "Cell__Vimentin": data["values"],
        }
    )
    model, evidence = calibrate_marker_frame(
        frame,
        donor_id="donor-A",
        marker="B3TUBB",
        reference_positive="reference_positive",
        registry=registry,
        config=_config(),
    )
    model_table = models_to_long([model])
    assert model_table[["donor_id", "marker"]].to_dict("records") == [
        {"donor_id": "donor-A", "marker": "B3TUBB"}
    ]

    wide = evidence_to_wide(evidence)
    assert "B3TUBB__state" in wide
    round_trip = evidence_from_wide(wide)
    merged = evidence[
        [
            "donor_id",
            "object_id",
            "marker",
            "state",
            "log2_threshold_ratio",
            "empirical_tail_probability",
            "expression_probability",
            "empirical_tail_evidence",
        ]
    ].sort_values(["donor_id", "object_id", "marker"], ignore_index=True)
    pd.testing.assert_frame_equal(
        round_trip[merged.columns], merged, check_dtype=False
    )

    result = calibrate_donor(
        frame,
        donor_id="donor-A",
        markers=["B3TUBB", "Vimentin"],
        reference_positive={
            "EpCAM": "reference_positive",
            "E_cadherin": "reference_positive",
        },
        registry=registry,
        config=_config(),
    )
    assert set(result.model_long["marker"]) == {"B3TUBB", "Vimentin"}
    assert len(result.evidence_long) == 2 * len(frame)
    assert {"B3TUBB__state", "Vimentin__state"} <= set(
        result.evidence_wide().columns
    )


def test_empirical_tail_probability_has_finite_sample_floor():
    controls = np.array([1.0, 2.0, 3.0, 4.0])
    probability = empirical_tail_probability(
        np.array([0.0, 2.5, 10.0]), controls
    )
    assert probability.tolist() == pytest.approx([1.0, 0.6, 0.2])
