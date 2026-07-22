from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler import restore_validation as rv
from phenocycler import restore_validation_figures as rvf
from phenocycler import restore_pair_review as rpr


def test_candidate_pairs_match_revised_first_pass():
    assert rv.BASELINE_PAIRS == (
        ("E_cadherin", "CD31"),
        ("CD31", "E_cadherin"),
        ("CD3e", "E_cadherin"),
        ("CD20", "E_cadherin"),
        ("CD68", "E_cadherin"),
        ("CD11b", "E_cadherin"),
        ("B3TUBB", "EpCAM"),
        ("Vimentin", "E_cadherin"),
    )
    assert ("CD3e", "CD163") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert ("CD20", "CD68") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert ("CD68", "CD3e") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert ("CD11b", "CD20") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert ("CD20", "Pan_Cytokeratin") in rv.CANDIDATE_PAIRS
    assert rv.QUPATH_MARKERS["Pan_Cytokeratin"] == "Pan-Cytokeratin"
    assert rv.QUPATH_MARKERS["Ker8_18"] == "Ker8-18"
    assert set(rv.MARKER_INPUT_FLOOR_METHODS) == set(rv.QUPATH_MARKERS)


def test_directional_groups_selects_reference_arm_not_cluster_number():
    target = np.array([10, 11, 9, 1, 2, 1, 1, 1], dtype=float)
    reference = np.array([1, 2, 1, 10, 9, 11, 1, 2], dtype=float)
    labels = np.array([7, 7, 7, 3, 3, 3, 5, 5])
    result = rv.directional_groups(target, reference, labels)
    assert result.reference_group == 3
    assert result.target_group == 7
    assert result.directional_mean
    assert result.directional_median


def test_directional_groups_rejects_same_direction_clusters():
    target = np.array([1, 2, 8, 9], dtype=float)
    reference = np.array([1, 2, 8, 9], dtype=float)
    labels = np.array([0, 0, 1, 1])
    result = rv.directional_groups(target, reference, labels)
    assert result.reference_group == result.target_group
    assert not result.directional_mean
    assert not result.directional_median


def test_manuscript_maximum_and_target_component_rescue_call():
    stats = rv.negative_control_statistics(np.array([1.0, 2.0, 4.0]))
    assert stats["maximum"] == 4.0
    normalized, positive = rv.normalize_and_call(
        np.array([0.0, 4.0, 4.0001, 8.0]),
        stats["maximum"],
        target_component_supported=np.array([True, False, False, False]),
    )
    assert np.allclose(normalized, [0.0, 1.0, 1.000025, 2.0])
    assert positive.tolist() == [True, False, True, True]


def test_target_component_rescue_requires_shape_matched_explicit_mask():
    with pytest.raises(ValueError, match="must match"):
        rv.normalize_and_call(
            np.array([1.0, 2.0]),
            2.0,
            target_component_supported=np.array([True]),
        )


def test_negative_control_maximum_requires_positive_nonnegative_values():
    with pytest.raises(ValueError, match="nonnegative"):
        rv.negative_control_statistics(np.array([-1.0, 2.0]))
    with pytest.raises(ValueError, match="positive"):
        rv.negative_control_statistics(np.zeros(3))


def test_nnmf_recovers_compartment_rich_mutually_exclusive_arms():
    rng = np.random.default_rng(1)
    target = np.c_[
        rng.normal(10, 0.2, 80),
        rng.normal(12, 0.2, 80),
        rng.normal(11, 0.2, 80),
    ]
    reference = np.c_[
        rng.normal(1, 0.05, 80),
        rng.normal(1.2, 0.05, 80),
        rng.normal(1.1, 0.05, 80),
    ]
    target_arm = np.column_stack(
        [target[:, 0], reference[:, 0], target[:, 1], reference[:, 1],
         target[:, 2], reference[:, 2]]
    )
    reference_arm = target_arm[:, [1, 0, 3, 2, 5, 4]]
    features = np.vstack([target_arm, reference_arm])
    labels, components, metadata = rv.fit_nnmf(
        features, n_components=2, seed=0
    )
    result = rv.directional_groups(
        features[:, 4], features[:, 5], labels
    )
    assert components.shape == (2, 6)
    assert metadata["converged"]
    assert result.directional_mean
    assert result.directional_median


def test_nnmf_predict_accepts_explicit_balanced_fit_indices():
    features = np.array(
        [
            [10.0, 1.0],
            [9.0, 1.0],
            [8.0, 1.0],
            [1.0, 8.0],
            [1.0, 9.0],
            [1.0, 10.0],
        ]
    )
    labels, _, metadata = rv.fit_nnmf_predict(
        features,
        n_components=2,
        seed=0,
        fit_cap=None,
        fit_indices=np.array([0, 1, 4, 5]),
    )
    assert metadata["fit_n"] == 4
    assert len(labels) == len(features)
    with pytest.raises(ValueError, match="cannot both"):
        rv.fit_nnmf_predict(
            features,
            n_components=2,
            seed=0,
            fit_cap=4,
            fit_indices=np.array([0, 1, 4, 5]),
        )


def test_extract_compartment_sample_matches_canonical_ids(tmp_path):
    donors = ["1", "2"]
    canonical = tmp_path / "cells"
    for donor, ids in {"1": ["a", "b"], "2": ["c"]}.items():
        part = canonical / f"donor_id={donor}"
        part.mkdir(parents=True)
        pd.DataFrame({"object_id": ids}).to_parquet(part / "data_0.parquet")

    columns = {
        "Image": ["1_scan", "1_scan", "1_scan", "2_scan"],
        "Object ID": ["a", "b", "extra", "c"],
        "Cell: Area µm^2": [10, 11, 1, 12],
        "Centroid X µm": [1, 2, 3, 4],
        "Centroid Y µm": [5, 6, 7, 8],
    }
    for compartment in rv.COMPARTMENTS:
        for marker in rv.QUPATH_MARKERS:
            columns[rv.qupath_measurement(compartment, marker)] = [1, 2, 3, 4]
    csv_path = tmp_path / "measurements.csv"
    pd.DataFrame(columns).to_csv(csv_path, index=False)

    output = tmp_path / "sample"
    audit = rv.extract_compartment_sample(
        [csv_path],
        canonical,
        output,
        donors=donors,
        sample_modulus=1,
        threads=1,
    )
    assert audit["observed_n"].sum() == 3
    result = pd.read_parquet(output)
    assert set(result["object_id"]) == {"a", "b", "c"}
    assert rv.compartment_column("Nucleus", "E_cadherin") in result
    assert "extra" not in set(result["object_id"])


def test_load_validation_sample_accepts_present_redsea_row_with_all_nan_values(tmp_path):
    donor = "1"
    sample_root = tmp_path / "sample"
    raw_root = tmp_path / "raw"
    redsea_root = tmp_path / "redsea"
    for root in (sample_root, raw_root, redsea_root):
        (root / f"donor_id={donor}").mkdir(parents=True)

    sample = {
        "donor_id": [donor],
        "image": ["1_scan"],
        "object_id": ["cell"],
        "cell_area": [10.0],
        "x_um": [1.0],
        "y_um": [2.0],
        "source_csv": ["source.csv"],
    }
    raw = {
        "object_id": ["cell"],
        "X_centroid": [1.0],
        "Y_centroid": [2.0],
    }
    redsea = {"object_id": ["cell"]}
    for compartment in rv.COMPARTMENTS:
        for marker in rv.QUPATH_MARKERS:
            sample[rv.compartment_column(compartment, marker)] = [3.0]
    for marker in rv.QUPATH_MARKERS:
        raw[marker] = [3.0]
        redsea[marker] = [np.nan]

    pd.DataFrame(sample).to_parquet(
        sample_root / f"donor_id={donor}" / "data_0.parquet"
    )
    pd.DataFrame(raw).to_parquet(
        raw_root / f"donor_id={donor}" / "data_0.parquet"
    )
    pd.DataFrame(redsea).to_parquet(
        redsea_root / f"donor_id={donor}" / "data_0.parquet"
    )

    result, audit = rv.load_validation_sample(
        sample_root, raw_root, redsea_root, donor
    )
    assert len(result) == 1
    assert audit["raw_redsea_key_match"] is True


def test_spatial_agreement_uses_nearest_nonself_neighbor():
    coordinates = np.array([[0.0, 0.0], [1.0, 0.0], [100.0, 0.0], [101.0, 0.0]])
    control = np.array([False, False, True, True])
    agreement, excess = rv._spatial_agreement(coordinates, control, cap=4)
    assert agreement == 1.0
    assert excess == 1.0


def _assessment(**overrides):
    values = {
        "reference_n": 100,
        "target_n": 50,
        "distinct": True,
        "directional_mean": True,
        "directional_median": True,
        "converged": True,
        "stable": True,
        "target_arm_fold": 3.0,
        "reference_arm_fold": 3.0,
        "pair_review": rv.PairReview.ACCEPTED,
    }
    values.update(overrides)
    return rv.assess_pair_state(**values)


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({"technical_error": "missing values"}, rv.PairState.TECHNICAL_FAILURE),
        ({"reference_n": 1}, rv.PairState.NO_REFERENCE_POPULATION),
        ({"reference_arm_fold": 1.5}, rv.PairState.NO_REFERENCE_POPULATION),
        ({"directional_median": False}, rv.PairState.INVALID_PAIR),
        ({"stable": False}, rv.PairState.MODEL_UNSTABLE),
        ({"target_arm_fold": 1.99}, rv.PairState.LOW_SIGNAL_UNRESOLVED),
        (
            {
                "target_arm_fold": 1.99,
                "expression_review": rv.ExpressionReview.CONFIRMED_ABSENT,
            },
            rv.PairState.NO_TARGET_EXPRESSION_CONFIRMED,
        ),
        ({"pair_review": rv.PairReview.UNREVIEWED}, rv.PairState.PAIR_REVIEW_PENDING),
        ({}, rv.PairState.VALID),
    ],
)
def test_pair_state_machine_is_explicit(overrides, expected):
    assert _assessment(**overrides).state is expected


def test_confirmed_presence_with_weak_separation_invalidates_pair():
    result = _assessment(
        target_arm_fold=1.5,
        expression_review=rv.ExpressionReview.CONFIRMED_PRESENT,
    )
    assert result.state is rv.PairState.INVALID_PAIR


def test_calls_are_blocked_except_valid_or_confirmed_absence():
    values = np.array([0.0, 2.0, 2.0001])
    normalized, positive = rv.normalize_for_assessment(
        values,
        _assessment(),
        divisor=2.0,
        target_component_supported=np.array([True, True, False]),
    )
    assert np.allclose(normalized, [0.0, 1.0, 1.00005])
    assert positive.tolist() == [True, True, True]

    absent = _assessment(
        target_arm_fold=1.0,
        expression_review=rv.ExpressionReview.CONFIRMED_ABSENT,
    )
    normalized, positive = rv.normalize_for_assessment(
        values,
        absent,
        divisor=None,
        target_component_supported=np.ones(3, dtype=bool),
    )
    assert np.isnan(normalized).all()
    assert not positive.any()
    with pytest.raises(ValueError, match="must not fabricate"):
        rv.normalize_for_assessment(
            values,
            absent,
            divisor=2.0,
            target_component_supported=np.ones(3, dtype=bool),
        )
    with pytest.raises(ValueError, match="cannot produce"):
        rv.normalize_for_assessment(
            values,
            _assessment(target_arm_fold=1.0),
            divisor=None,
            target_component_supported=np.ones(3, dtype=bool),
        )


def test_ordered_controls_apply_vertical_separator_before_target_divisor():
    result = rv.ordered_control_groups(
        target_raw=np.array([100, 90, 100, 1, 2, 1], dtype=float),
        reference_raw=np.array([1, 2, 100, 100, 80, 100], dtype=float),
        target_redsea=np.array([90, 80, 95, 1, 2, 50], dtype=float),
        reference_redsea=np.array([1, 2, 90, 8, 9, 1.5], dtype=float),
        labels=np.array([1, 1, 1, 0, 0, 0]),
        qc_retained=np.ones(6, dtype=bool),
        target_group=1,
        reference_group=0,
        target_floor=10,
        reference_floor=10,
    )
    assert result.target_population.tolist() == [
        True,
        True,
        False,
        False,
        False,
        False,
    ]
    assert result.target_supported.tolist() == [
        True,
        True,
        False,
        False,
        False,
        False,
    ]
    assert result.reference_separator == 2
    assert result.reference_candidates.tolist() == [
        False,
        False,
        False,
        True,
        True,
        True,
    ]
    assert result.reference_control.tolist() == [
        False,
        False,
        False,
        True,
        True,
        False,
    ]
    target_divisor = np.max(
        np.array([90, 80, 95, 1, 2, 50], dtype=float)[
            result.reference_control
        ]
    )
    assert target_divisor == 2


def test_reference_negative_target_component_requires_input_qc_retention():
    result = rv.ordered_control_groups(
        target_raw=np.array([100, 2, 100, 1], dtype=float),
        reference_raw=np.array([1, 1, 100, 100], dtype=float),
        target_redsea=np.array([80, 0.5, 40, 1], dtype=float),
        reference_redsea=np.array([2, 1.5, 40, 8], dtype=float),
        labels=np.array([1, 1, 1, 0]),
        qc_retained=np.array([True, False, True, True]),
        target_group=1,
        reference_group=0,
        target_floor=10,
        reference_floor=10,
    )
    assert result.reference_separator == 2
    assert result.target_population.tolist() == [True, False, False, False]
    assert result.target_supported.tolist() == [True, False, False, False]
    assert not result.target_supported[~np.array([True, False, True, True])].any()
    assert result.reference_control.tolist() == [False, False, False, True]


def test_density_is_finite_and_respects_limits():
    centers, density = rvf._density(np.linspace(0, 2, 100), (0, 2))
    assert len(centers) == len(density) == 100
    assert np.isfinite(density).all()
    assert centers.min() >= 0
    assert centers.max() <= 2


def test_figure_coordinates_follow_manuscript_reference_x_target_y():
    target = np.array([9.0, 54.0, 99.0])
    reference = np.array([1.0, 2.0, 3.0])
    x, y = rvf._manuscript_coordinates(
        target,
        reference,
        target_bounds=(9.0, 99.0),
        reference_bounds=(1.0, 3.0),
    )
    assert np.allclose(x, [0.0, 0.5, 1.0])
    assert np.allclose(y, [0.0, 0.5, 1.0])


def test_raw_otsu_floor_separates_low_and_high_intensity_modes():
    values = np.r_[np.linspace(1, 3, 100), np.linspace(90, 110, 100)]
    floor = rv.raw_otsu_floor(values)
    assert 2.5 < floor < 90


def test_pair_input_qc_drops_only_double_low_cells():
    target = np.array([1.0, 2.0, 100.0, 1.0])
    reference = np.array([1.0, 2.0, 1.0, 100.0])
    retained, target_floor, reference_floor = rv.pair_input_qc(
        target,
        reference,
        target_marker="E_cadherin",
        reference_marker="CD31",
    )
    assert target_floor.value < 100
    assert reference_floor.value < 100
    assert retained.tolist() == [False, False, True, True]


def test_sparse_marker_floor_requires_both_otsu_and_triangle():
    values = np.r_[np.linspace(1, 3, 1_000), np.linspace(20, 100, 20)]
    result = rv.estimate_marker_input_floor(values, "CD3e")
    assert result.method == rv.MAX_LINEAR_OTSU_TRIANGLE_FLOOR
    assert result.linear_triangle is not None
    assert result.value == max(result.linear_otsu, result.linear_triangle)
    assert result.value >= result.linear_otsu


def test_structural_marker_floor_remains_linear_otsu():
    values = np.r_[np.linspace(1, 3, 100), np.linspace(90, 110, 100)]
    result = rv.estimate_marker_input_floor(values, "CD31")
    assert result.method == rv.LINEAR_OTSU_FLOOR
    assert result.value == result.linear_otsu
    assert result.linear_triangle is None


def test_marker_input_floor_fails_closed_without_explicit_policy():
    with pytest.raises(ValueError, match="no input-floor method"):
        rv.estimate_marker_input_floor(np.array([1.0, 2.0, 3.0]), "Unknown")


def test_balanced_arm_fit_sampling_equalizes_exclusive_arms():
    target = np.array([10.0] * 6 + [1.0] * 5 + [10.0])
    reference = np.array([1.0] * 6 + [10.0] * 5 + [10.0])
    indices, counts = rv.balanced_arm_fit_indices(
        target,
        reference,
        target_floor=5.0,
        reference_floor=5.0,
        fit_cap=8,
        seed=0,
    )
    assert counts == {
        "candidate_target_n": 6,
        "candidate_reference_n": 5,
        "candidate_double_high_n": 1,
        "fit_per_arm": 4,
    }
    assert len(indices) == 8
    assert 11 not in indices
    assert ((target[indices] > 5) & (reference[indices] <= 5)).sum() == 4
    assert ((target[indices] <= 5) & (reference[indices] > 5)).sum() == 4


def test_robust_minmax_scaling_uses_balanced_fit_sample_bounds():
    features = np.array([[0.0, 10.0], [5.0, 15.0], [10.0, 20.0], [1_000.0, 2_000.0]])
    scaled, lower, upper = rv.robust_minmax_scale(
        features,
        np.array([0, 1, 2]),
        lower_quantile=0,
        upper_quantile=1,
    )
    assert np.array_equal(lower, [0.0, 10.0])
    assert np.array_equal(upper, [10.0, 20.0])
    assert np.allclose(scaled[:3], [[0, 0], [0.5, 0.5], [1, 1]])
    assert np.array_equal(scaled[3], [1.0, 1.0])


def test_display_sampling_uses_equal_component_counts():
    labels = np.array([0] * 100 + [1] * 7)
    selected = rvf._display_indices(labels, cap_per_group=20, seed=0)
    assert np.bincount(labels[selected]).tolist() == [7, 7]


def test_resampled_background_statistics_reports_sampling_uncertainty():
    stats = rv.resampled_negative_control_statistics(
        np.arange(1.0, 101.0),
        sample_n=10,
        repeats=100,
        seed=0,
    )
    assert stats["n_available"] == 100
    assert stats["sample_n"] == 10
    assert stats["full_maximum"] == 100
    assert stats["maximum_q05"] <= stats["maximum"] <= stats["maximum_q95"]
    assert stats["maximum"] < stats["full_maximum"]


def test_sampled_maximum_grid_preserves_full_maximum_as_diagnostic_reference():
    grid = rv.sampled_maximum_grid(
        np.arange(1.0, 101.0),
        sample_sizes=(5, 10, 200),
        repeats=20,
        seed=0,
    )
    assert grid["sample_n"].tolist() == [5, 10]
    assert (grid["full_maximum"] == 100.0).all()
    assert (grid["maximum_to_full"] <= 1.0).all()


def test_method_specification_is_json_native_and_declares_data_roles():
    config = rv.PairValidationConfig()
    specification = config.specification()
    assert specification["method_version"] == rv.METHOD_VERSION
    assert specification["parameters"]["seeds"] == [0, 1, 2]
    assert "REDSEA" in specification["data_roles"]["divisor_and_calls"]
    assert specification["marker_input_floor_methods"]["CD3e"] == (
        rv.MAX_LINEAR_OTSU_TRIANGLE_FLOOR
    )
    assert specification["marker_input_floor_methods"]["CD20"] == (
        rv.MAX_LINEAR_OTSU_TRIANGLE_FLOOR
    )
    assert specification["marker_input_floor_methods"]["CD31"] == rv.LINEAR_OTSU_FLOOR
    assert "full maximum" in specification["divisor_policy"]
    assert "vertical reference separator" in specification["ordered_control_policy"]
    assert "input-QC-retained" in specification["ordered_control_policy"]
    assert "double-low" in specification["calling_policy"]
    assert "even below the horizontal divisor" in specification["calling_policy"]
    assert "ACCEPTED" in specification["pair_acceptance_policy"]


@pytest.mark.parametrize("donor", ["6457", "6579"])
def test_restore_pair_evaluation_rejects_excluded_donor_before_data_access(donor):
    with pytest.raises(ValueError, match=rf"excluded donor.*{donor}"):
        rv.evaluate_locked_pair(
            pd.DataFrame(),
            donor,
            "E_cadherin",
            "CD31",
        )


def test_locked_pair_assigns_all_cells_but_uses_double_low_qc_only_for_fitting():
    rng = np.random.default_rng(4)
    target_n = 40
    reference_n = 40
    double_low_n = 20
    target = np.r_[
        rng.normal(100, 2, target_n),
        rng.normal(1, 0.05, reference_n + double_low_n),
    ]
    reference = np.r_[
        rng.normal(1, 0.05, target_n),
        rng.normal(100, 2, reference_n),
        rng.normal(1, 0.05, double_low_n),
    ]
    frame = pd.DataFrame(
        {
            "raw__E_cadherin": target,
            "raw__CD31": reference,
            "redsea__E_cadherin": np.clip(target - 0.2, 0, None),
            "redsea__CD31": np.clip(reference - 0.2, 0, None),
            rv.compartment_column("Membrane", "E_cadherin"): target * 0.9,
            rv.compartment_column("Membrane", "CD31"): reference * 0.9,
            rv.compartment_column("Cell", "E_cadherin"): target,
            rv.compartment_column("Cell", "CD31"): reference,
        }
    )
    config = rv.PairValidationConfig(
        seeds=(0, 1, 2),
        fit_cap=60,
        threshold_sample_sizes=(5, 10),
        threshold_resamples=10,
        min_control_jaccard=0.8,
    )
    result = rv.evaluate_locked_pair(
        frame,
        "synthetic",
        "E_cadherin",
        "CD31",
        config=config,
        pair_review=rv.PairReview.ACCEPTED,
    )
    assert result["state"] == rv.PairState.VALID.value
    assert result["n_assigned"] == len(frame)
    assert (result["full_labels"] >= 0).all()
    assert result["n_qc_retained"] < len(frame)
    assert result["target_input_floor_method"] == rv.LINEAR_OTSU_FLOOR
    assert result["target_input_floor_linear_triangle"] is None
    assert result["reference_input_floor_method"] == rv.LINEAR_OTSU_FLOOR
    assert result["fit_per_arm"] == 30
    assert result["canonical_divisor"] == result["candidate_full_maximum"]
    assert result["target_supported_n"] >= result["target_anchor_n"]
    assert not result["target_supported"][~result["qc_retained"]].any()
    assert result["lineage_positive_n"] >= result["target_supported_n"]
    assert result["target_threshold_grid"]["sample_n"].tolist() == [5, 10]

    with_missing = frame.copy()
    with_missing.loc[len(with_missing) - 1, "redsea__CD31"] = np.nan
    missing_result = rv.evaluate_locked_pair(
        with_missing,
        "synthetic",
        "E_cadherin",
        "CD31",
        config=config,
        pair_review=rv.PairReview.ACCEPTED,
    )
    assert missing_result["state"] == rv.PairState.VALID.value
    assert missing_result["n_unavailable"] == 1
    assert missing_result["full_labels"][-1] == -1


def test_green_interval_is_hidden_only_when_it_has_no_visible_width():
    assert rvf._visible_scaled_interval(2.0, 4.0, (0.0, 10.0)) == (0.2, 0.4)
    assert rvf._visible_scaled_interval(4.0, 4.0, (0.0, 10.0)) is None
    assert rvf._visible_scaled_interval(12.0, 14.0, (0.0, 10.0)) is None


def test_figure_legend_explains_axis_specific_green_band_without_qc_box():
    labels = [handle.get_label() for handle in rvf._legend_handles()]
    assert not any("double-low" in label for label in labels)
    assert any("Step 1 vertical" in label for label in labels)
    assert any("Step 2 horizontal" in label for label in labels)
    assert any("axis-specific" in label for label in labels)
    assert any("absent when zero-width" in label for label in labels)


def test_bounded_footer_stays_within_subplot_lateral_span():
    fig, axes = rvf.plt.subplots(1, 2)
    rvf._add_bounded_footer(
        fig,
        list(axes),
        "Long explanatory text " * 30,
        bottom=0.01,
        height=0.05,
        fontsize=8,
    )
    footer = fig.axes[-1]
    assert footer.get_position().x0 >= min(axis.get_position().x0 for axis in axes)
    assert footer.get_position().x1 <= max(axis.get_position().x1 for axis in axes)
    assert "\n" in footer.texts[0].get_text()
    rvf.plt.close(fig)


def test_figure_thresholds_follow_their_manuscript_marker_axes():
    panel = {
        "distinct": True,
        "target_plot_bounds": (0.0, 100.0),
        "reference_plot_bounds": (0.0, 4.0),
        "target_background": {
            "maximum": 99.0,
            "maximum_q05": 90.0,
            "maximum_q95": 105.0,
            "mean_plus_3sd": 9.0,
        },
        "reference_background": {
            "maximum": 3.0,
            "maximum_q05": 2.0,
            "maximum_q95": 4.0,
            "mean_plus_3sd": 1.0,
        },
        "target_full_stats": {
            "maximum": 100.0,
            "mean_plus_3sd": 10.0,
        },
        "reference_full_stats": {
            "maximum": 4.0,
            "mean_plus_3sd": 2.0,
        },
    }
    fig, ax = rvf.plt.subplots()
    rvf._draw_thresholds(ax, panel)
    assert np.allclose(ax.lines[0].get_xdata(), 3.0 / 4.0)
    assert np.allclose(ax.lines[1].get_ydata(), 99.0 / 100.0)
    assert np.allclose(ax.lines[2].get_xdata(), 1.0)
    assert np.allclose(ax.lines[3].get_ydata(), 1.0)
    assert np.allclose(ax.lines[4].get_xdata(), 0.5)
    assert np.allclose(ax.lines[5].get_ydata(), 0.1)
    rvf.plt.close(fig)


def test_review_queue_selects_extremes_without_frequency_acceptance():
    rows = []
    for donor, target_fold, reference_fold, jaccard, double_high, ratio in (
        ("1", 1.2, 5.0, 0.99, 0.01, 0.8),
        ("2", 3.0, 1.5, 0.98, 0.02, 0.7),
        ("6476", 2.0, 4.0, 0.70, 0.50, 0.9),
        ("6539", 1.8, 3.0, 0.95, 0.10, 0.6),
    ):
        rows.append(
            {
                "donor": donor,
                "target": "CD20",
                "reference": "CD68",
                "target_fold": target_fold,
                "reference_fold": reference_fold,
                "min_control_jaccard": jaccard,
                "double_high_fraction": double_high,
                "target_sample_maximum": ratio * 100,
                "candidate_full_maximum": 100,
                "target_n": 10,
                "n_assigned": 100,
                "state": rv.PairState.PAIR_REVIEW_PENDING.value,
            }
        )
    queue = rpr.select_review_queue(
        pd.DataFrame(rows), pairs=(("CD20", "CD68"),)
    )
    by_donor = queue.set_index("donor")["selection_reasons"]
    assert "lowest seed stability" in by_donor["6476"]
    assert "pancreatitis" not in "; ".join(by_donor)
    assert "maintainer-confirmed target absence" in by_donor["6539"]
    assert "lowest target-arm fold" in by_donor["1"]
    assert "lowest reference-arm fold" in by_donor["2"]
    shuffled = pd.DataFrame(rows).sample(frac=1, random_state=19)
    pd.testing.assert_frame_equal(
        queue,
        rpr.select_review_queue(
            shuffled, pairs=(("CD20", "CD68"),)
        ),
    )


def test_review_categories_keep_unavailable_and_double_low_distinct():
    evaluation = {
        "valid_idx": np.array([0, 1, 2, 3, 4]),
        "qc_retained": np.array([True, True, True, False, True]),
        "reference_control": np.array([True, False, False, False, False]),
        "target_population": np.array([False, True, False, False, False]),
        "target_supported": np.array([False, True, False, True, True]),
        "target_raw": np.array([1.0, 10.0, 10.0, 1.0, 10.0]),
        "reference_raw": np.array([10.0, 1.0, 10.0, 1.0, 1.0]),
        "target_input_floor": 5.0,
        "reference_input_floor": 5.0,
    }
    categories = rpr.review_categories(evaluation, n_rows=6)
    assert categories.tolist() == [
        "Reference control arm",
        "High-confidence target anchor (sets Step 1)",
        "Double-high",
        "Double-low (fit excluded)",
        "Lower-confidence target component (reference-negative; retained)",
        "Unavailable",
    ]


def test_representative_zoom_is_deterministic_and_target_focused():
    donor_df = pd.DataFrame(
        {
            "x_canonical": np.r_[np.linspace(0, 100, 20), 800, 810, 820],
            "y_canonical": np.r_[np.linspace(0, 100, 20), 800, 810, 820],
        }
    )
    categories = np.array(
        ["Other retained (not a divisor control)"] * 20
        + ["Lower-confidence target component (reference-negative; retained)"] * 3,
        dtype=object,
    )
    first = rpr.representative_zoom_bounds(donor_df, categories)
    second = rpr.representative_zoom_bounds(donor_df, categories)
    assert first == second
    assert first[0][0] <= 810 <= first[0][1]
    assert first[1][0] <= 810 <= first[1][1]
    width = first[0][1] - first[0][0]
    height = first[1][1] - first[1][0]
    assert np.isclose(width / height, 2.2)
    assert np.isclose(
        width * height,
        (rpr.REPRESENTATIVE_ZOOM_SPAN_FRACTION * 820) ** 2,
    )
    with pytest.raises(ValueError, match="display_aspect"):
        rpr.representative_zoom_bounds(
            donor_df,
            categories,
            display_aspect=0,
        )


def test_marker_channel_display_is_deterministic_and_bounded():
    pixels = np.array(
        [[0.0, 1.0, 10.0, 100.0], [2.0, 20.0, 200.0, np.nan]]
    )
    first, first_window = rpr._display_marker_channel(pixels)
    second, second_window = rpr._display_marker_channel(pixels)
    assert np.array_equal(first, second)
    assert first_window == second_window
    assert np.nanmin(first) >= 0
    assert np.nanmax(first) <= 1
    assert first[0, 0] == 0
    assert rpr._normalized_marker_name("E_cadherin") == "ecadherin"
    assert rpr._normalized_marker_name("E-cadherin") == "ecadherin"
    assert rpr._normalized_marker_name("Ker8_18") == "ker818"


def test_full_cell_projection_labels_every_crop_cell_and_matches_sample():
    feature_columns = [
        "Membrane__Target",
        "Membrane__Reference",
        "Cell__Target",
        "Cell__Reference",
    ]
    components = np.array(
        [[1.0, 0.0, 1.0, 0.0], [0.0, 1.0, 0.0, 1.0]]
    ) / np.sqrt(2.0)
    full_df = pd.DataFrame(
        {
            "object_id": [f"cell-{index}" for index in range(5)],
            "image": ["synthetic.qptiff - resolution #1"] * 5,
            "x_canonical": np.arange(5, dtype=float),
            "y_canonical": np.arange(5, dtype=float),
            "Membrane__Target": [10.0, 10.0, 0.0, 10.0, np.nan],
            "Membrane__Reference": [0.0, 0.0, 10.0, 0.0, 0.0],
            "Cell__Target": [10.0, 10.0, 0.0, 10.0, 1.0],
            "Cell__Reference": [0.0, 0.0, 10.0, 0.0, 1.0],
            "raw__Target": [10.0, 2.0, 1.0, 10.0, 1.0],
            "raw__Reference": [1.0, 1.0, 10.0, 10.0, 1.0],
            "redsea__Target": [10.0, 2.0, 1.0, 10.0, 1.0],
            "redsea__Reference": [1.0, 1.0, 5.0, 5.0, 1.0],
        }
    )
    evaluation = {
        "target": "Target",
        "reference": "Reference",
        "feature_columns": feature_columns,
        "scale_lower": np.zeros(4),
        "scale_upper": np.full(4, 10.0),
        "components": components,
        "target_group": 0,
        "reference_group": 1,
        "target_input_floor": 5.0,
        "reference_input_floor": 5.0,
        "reference_separator": 2.0,
        "candidate_full_maximum": 5.0,
    }
    projected = rpr.project_review_crop(
        full_df, evaluation, (0.0, 4.0), (0.0, 4.0)
    )
    assert projected["review_category"].tolist() == [
        "High-confidence target anchor (sets Step 1)",
        "Double-low (fit excluded)",
        "Reference control arm",
        "Double-high",
        "Unavailable",
    ]
    assert projected["input_qc_retained"].tolist() == [
        True,
        False,
        True,
        True,
        False,
    ]
    assert projected["target_anchor"].tolist() == [True, False, False, False, False]
    assert projected["target_component_supported"].tolist() == [
        True,
        False,
        False,
        False,
        False,
    ]
    sample_categories = projected["review_category"].to_numpy(object)
    assert (
        rpr._validate_projection_overlap(
            full_df,
            sample_categories,
            projected,
            (0.0, 4.0),
            (0.0, 4.0),
        )
        == 5
    )


def test_review_row_combines_marker_images_overlay_and_merged_view():
    donor_df = pd.DataFrame(
        {
            "x_canonical": np.linspace(0, 100, 12),
            "y_canonical": np.linspace(0, 100, 12),
        }
    )
    categories = np.array(
        ["Other retained (not a divisor control)"] * 8
        + ["High-confidence target anchor (sets Step 1)"] * 2
        + ["Lower-confidence target component (reference-negative; retained)"] * 2,
        dtype=object,
    )
    evaluation = {
        "donor": "synthetic",
        "target": "Target",
        "reference": "Reference",
        "full_labels": np.zeros(12, dtype=int),
        "target_redsea": np.linspace(0, 10, 12),
        "reference_redsea": np.linspace(0, 10, 12),
        "display_target_bounds": (0.0, 10.0),
        "display_reference_bounds": (0.0, 10.0),
        "target_full_stats": {"maximum": 2.0},
        "reference_full_stats": {"maximum": 2.0},
        "state": rv.PairState.PAIR_REVIEW_PENDING.value,
        "target_fold": 2.0,
        "reference_fold": 3.0,
        "min_control_jaccard": 0.95,
        "target_anchor_n": 2,
        "target_supported_n": 4,
        "target_supported_additional_n": 2,
        "target_supported_below_divisor_n": 2,
        "target_supported_additional_below_divisor_n": 1,
    }
    pixels = np.arange(40 * 88, dtype=np.uint16).reshape(40, 88)
    marker_zoom = rpr.MarkerZoom(
        target_pixels=pixels,
        reference_pixels=np.flipud(pixels),
        x_bounds_um=(0.0, 100.0),
        y_bounds_um=(0.0, 100.0),
        target_channel="Target",
        reference_channel="Reference",
        source_image="synthetic.qptiff - resolution #1",
        source_path=Path("synthetic.qptiff"),
        pyramid_level=2,
        downsample_x=4.0,
        downsample_y=4.0,
        pixel_size_x_um=0.5,
        pixel_size_y_um=0.5,
    )
    projected_crop = donor_df.copy()
    projected_crop["review_category"] = categories
    fig, axes = rpr.plt.subplots(1, 2, figsize=(12, 5))
    rpr._draw_review_row(
        axes,
        donor_df,
        evaluation,
        categories,
        "synthetic review case",
        marker_zoom,
        projected_crop,
    )
    assert len(axes[1].child_axes) == 4
    target_axis, reference_axis, overlay_axis, merged_axis = axes[1].child_axes
    assert len(target_axis.images) == 1
    assert len(reference_axis.images) == 1
    assert len(overlay_axis.images) == 0
    assert len(merged_axis.images) == 1
    assert len(overlay_axis.child_axes) == 1
    assert overlay_axis.child_axes[0].get_title() == ""
    assert "Representative QPTIFF close-up" in axes[1].get_title()
    assert "All 12 canonical crop cells labeled" in axes[1].get_title()
    assert "Raw target marker pixels" in target_axis.get_title()
    assert "Raw reference marker pixels" in reference_axis.get_title()
    assert "Merged markers" in merged_axis.get_title()
    assert len(overlay_axis.collections[0].get_offsets()) == 12
    assert len(merged_axis.collections[0].get_offsets()) == 12
    assert not axes[0].texts
    assert "High-confidence anchor n=2" in axes[0].get_title()
    assert "input-QC-retained lower-confidence target n=2" in axes[0].get_title()
    assert "Lower-confidence at/below Step 2 n=1" in axes[0].get_title()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    assert not axes[0].title.get_window_extent(renderer).overlaps(
        axes[0].get_window_extent(renderer)
    )
    rpr.plt.close(fig)


def test_target_confidence_colors_are_distinct_and_high_contrast():
    high = np.asarray(
        rpr.to_rgba(
            rpr.CATEGORY_COLOR[
                "High-confidence target anchor (sets Step 1)"
            ]
        )[:3]
    )
    lower = np.asarray(
        rpr.to_rgba(
            rpr.CATEGORY_COLOR[
                "Lower-confidence target component (reference-negative; retained)"
            ]
        )[:3]
    )
    assert not np.array_equal(high, lower)
    assert np.linalg.norm(high - lower) > 0.6


def test_review_figure_uses_large_fonts_and_nonoverlapping_legend_band():
    fig, axes, legend_axis = rpr._create_review_figure(6)
    axes[0, 0].set_title(
        "LOW_SIGNAL_UNRESOLVED\n"
        "Target/reference fold 1.56/3.04 | control Jaccard 0.98\n"
        "High-confidence anchor n=28,685 | input-QC-retained lower-confidence target n=4,531\n"
        "Lower-confidence at/below Step 2 n=3,884",
        fontsize=rpr.PDF_ROW_TITLE_FONT_SIZE,
        pad=10,
        linespacing=1.3,
    )
    heading = rpr._add_review_suptitle(fig, "CD3e", "CD163")
    assert "CD3e = max linear otsu triangle" in heading.get_text()
    assert "CD163 = linear otsu" in heading.get_text()
    legend = rpr._add_review_legend(legend_axis)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_box = legend.get_window_extent(renderer)
    heading_box = heading.get_window_extent(renderer)
    first_row_title_box = axes[0, 0].title.get_window_extent(renderer)
    data_boxes = [
        axis.get_window_extent(renderer)
        for axis in axes.flat
    ]
    assert all(not legend_box.overlaps(box) for box in data_boxes)
    assert heading_box.y0 > first_row_title_box.y1
    assert legend_axis.get_position().y1 < min(
        axis.get_position().y0 for axis in axes[-1]
    )
    assert rpr.PDF_TITLE_HEIGHT_IN >= 2.0
    assert rpr.PDF_SUPTITLE_FONT_SIZE >= 20
    assert rpr.PDF_AXIS_LABEL_FONT_SIZE >= 14
    assert rpr.PDF_ROW_TITLE_FONT_SIZE >= 12
    assert rpr.PDF_LEGEND_FONT_SIZE >= 11
    rpr.plt.close(fig)


def test_review_shortlist_has_one_rationale_per_unique_pair():
    assert len(rpr.SHORTLIST_PAIRS) == 13
    assert len(set(rpr.SHORTLIST_PAIRS)) == len(rpr.SHORTLIST_PAIRS)
    assert set(rpr.SHORTLIST_PAIRS) == set(rpr.PAIR_RATIONALE)


def test_qupath_importer_is_rendered_with_bundle_csv_path(tmp_path):
    review_dir = tmp_path / "review bundle" / "qupath_review"
    output = tmp_path / rpr.QUPATH_IMPORTER_NAME
    rpr._write_qupath_importer(output, review_dir)
    text = output.read_text()
    assert rpr.QUPATH_REVIEW_DIR_TOKEN not in text
    assert str(review_dir.resolve()) in text
    assert "Dialogs.showChoiceDialog" in text
    assert "Choose a target <- reference overlay" in text
    assert "[Clear temporary RESTORE review classes]" in text
    assert "CSV/image mismatch; nothing was imported." in text
    assert "def setTemporaryClasses = true" in text
    assert "classifyNotSampled" not in text
    assert '"RESTORE review - Not sampled"' not in text
    assert '"RESTORE review - High-confidence target"' in text
    assert '"RESTORE review - Lower-confidence target"' in text
    assert '"5" : "RESTORE review - Lower-confidence target"' in text
    assert "[106, 61, 154]" in text
    assert "outside the representative crop retained" in text
    assert 'def originalClassMetadataKey = "RESTORE review original PathClass"' in text
    assert "PathClass.fromString(originalClass)" in text
    assert "Restored original PathClasses" in text
    assert "have non-RESTORE PathClasses" not in text
    assert "allowClassOverwrite" not in text
    assert "<=>" not in text


def test_review_workflow_is_copied_into_bundle(tmp_path):
    output = tmp_path / rpr.REVIEW_WORKFLOW_NAME
    rpr._write_review_workflow(output)
    text = output.read_text()
    assert "Use only these four surfaces" in text
    assert "review_decisions.csv" in text
    assert "No filename editing is needed" in text
    assert "Never reject or adjust a pair because a donor differs" in text
    assert "input-QC-retained lower-confidence" in text
    assert "double-low and cannot be overwritten" in text
    assert "below Step 2" in text
    assert "Every canonical cell" in text
    assert "Donor 6476 is not a pancreatitis donor" in text
    assert "Donor 6457 is excluded" in text
    assert "Donor 6579 is excluded" in text


def test_review_checklist_is_minimal_and_follows_pair_order():
    queue = pd.DataFrame(
        [
            {
                "donor": "6476",
                "target": "CD31",
                "reference": "E_cadherin",
                "selection_reasons": "stress case",
                "state": "PAIR_REVIEW_PENDING",
                "state_reason": "pending",
            },
            {
                "donor": "6450",
                "target": "E_cadherin",
                "reference": "CD31",
                "selection_reasons": "weak case",
                "state": "PAIR_REVIEW_PENDING",
                "state_reason": "pending",
            },
        ]
    )
    exports = [
        {
            "donor": "6476",
            "target": "CD31",
            "reference": "E_cadherin",
            "image": "6476 image",
            "file": "cd31_from_e_cadherin__6476.csv",
        },
        {
            "donor": "6450",
            "target": "E_cadherin",
            "reference": "CD31",
            "image": "6450 image",
            "file": "e_cadherin_from_cd31__6450.csv",
        },
    ]
    checklist = rpr.build_review_checklist(queue, exports)
    assert checklist["pair"].tolist() == [
        "E_cadherin <- CD31",
        "CD31 <- E_cadherin",
    ]
    assert checklist["pair_order"].tolist() == [1, 2]
    assert checklist["review_order"].tolist() == [1, 2]
    assert checklist["pdf_row"].tolist() == [1, 1]
    assert checklist["donor_disposition"].tolist() == ["", ""]
    assert checklist["pdf_file"].tolist() == [
        "spatial_review/E_cadherin_from_CD31.pdf",
        "spatial_review/CD31_from_E_cadherin.pdf",
    ]


def test_review_bundle_removes_temporary_output_after_failure(
    tmp_path, monkeypatch
):
    results = tmp_path / "results"
    results.mkdir()
    (results / "pair_evaluations.csv").write_text(
        f"method_version,donor,target,reference\n{rv.METHOD_VERSION},1,T,R\n"
    )
    (results / "method_spec.json").write_text(
        f'{{"method_version": "{rv.METHOD_VERSION}"}}\n'
    )
    (results / "expression_reviews.csv").write_text(
        "donor,target,review,evidence\n"
    )
    (results / "pair_reviews.csv").write_text(
        "target,reference,review,evidence\n"
    )
    monkeypatch.setattr(
        rpr, "_config_from_specification", lambda specification: object()
    )
    monkeypatch.setattr(
        rpr, "load_expression_reviews", lambda path: ({}, pd.DataFrame())
    )
    monkeypatch.setattr(
        rpr, "load_pair_reviews", lambda path: ({}, pd.DataFrame())
    )
    monkeypatch.setattr(
        rpr, "summarize_pairs", lambda evaluations: pd.DataFrame({"x": [1]})
    )
    monkeypatch.setattr(
        rpr, "select_review_queue", lambda evaluations: pd.DataFrame()
    )

    def fail_write(*args, **kwargs):
        raise RuntimeError("simulated write failure")

    monkeypatch.setattr(pd.DataFrame, "to_csv", fail_write)
    output = tmp_path / "review"
    images = tmp_path / "images"
    images.mkdir()
    full_cells = tmp_path / "full_cells"
    full_cells.mkdir()
    with pytest.raises(RuntimeError, match="simulated write failure"):
        rpr.build_review_bundle(
            tmp_path / "sample",
            full_cells,
            tmp_path / "raw",
            tmp_path / "redsea",
            images,
            results,
            output,
        )
    assert not output.exists()
    assert not output.with_name("review.tmp").exists()
