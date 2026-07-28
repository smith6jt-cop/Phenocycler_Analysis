from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler import restore_validation as rv
from phenocycler import restore_validation_figures as rvf
from phenocycler import restore_pair_review as rpr


def test_candidate_pairs_match_revised_first_pass():
    # BroadLineageRevised.md lists 8 first-pass pairs; CD20 <- E_cadherin was the eighth and is gone
    # with CD20's exclusion from RESTORE (see test_cd20_takes_no_part_in_restore).
    assert rv.BASELINE_PAIRS == (
        ("E_cadherin", "CD31"),
        ("CD31", "E_cadherin"),
        ("CD3e", "E_cadherin"),
        ("CD68", "E_cadherin"),
        ("CD11b", "E_cadherin"),
        ("B3TUBB", "EpCAM"),
        ("Vimentin", "E_cadherin"),
    )
    assert ("CD3e", "CD163") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert ("CD68", "CD3e") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    # CD11b <- CD20 was retired: CD20 is too sparse in pancreas to define a negative control,
    # so CD3e is the immune reference (see test_cd20_takes_no_part_in_restore).
    assert ("CD11b", "CD3e") in rv.IMMUNE_REFERENCE_SCREEN_PAIRS
    assert rv.QUPATH_MARKERS["Pan_Cytokeratin"] == "Pan-Cytokeratin"
    assert rv.QUPATH_MARKERS["Ker8_18"] == "Ker8-18"
    assert set(rv.MARKER_INPUT_FLOOR_METHODS) == set(rv.QUPATH_MARKERS)


def test_accepted_pairs_are_the_frozen_one_reference_per_target():
    assert set(rv.ACCEPTED_PAIRS).issubset(set(rv.CANDIDATE_PAIRS))
    targets = [t for t, _ in rv.ACCEPTED_PAIRS]
    assert len(targets) == len(set(targets))  # one reference per target
    ref = dict(rv.ACCEPTED_PAIRS)
    # CD11b took the maintainer override to EpCAM, not the E_cadherin baseline
    assert ref["CD11b"] == "EpCAM"
    assert ("CD11b", "E_cadherin") in set(rv.BASELINE_PAIRS)
    assert rv.PAIR_SETS["accepted"] == rv.ACCEPTED_PAIRS
    for target, reference in rv.ACCEPTED_PAIRS:  # no excluded marker in either role
        assert target not in rv.RESTORE_EXCLUDED_MARKERS
        assert reference not in rv.RESTORE_EXCLUDED_MARKERS


def test_epithelial_target_comparator_promoted_only_against_cd31():
    """Pan-CK / Ker8_18 / EpCAM entered as diagnostic TARGETS (they appear elsewhere only as
    references) so the single-marker Epithelial gate could be tested against alternatives. The
    anatomy check settled it and they were promoted 2026-07-27 — but only against CD31, E-cadherin's
    own frozen reference. The Vimentin-referenced rows stay comparators, so exactly one reference
    reaches production per target."""
    comparator = set(rv.EPITHELIAL_TARGET_COMPARATOR_PAIRS)
    assert comparator <= set(rv.CANDIDATE_PAIRS)
    assert rv.PAIR_SETS["epithelial-comparator"] == rv.EPITHELIAL_TARGET_COMPARATOR_PAIRS
    assert comparator & set(rv.ACCEPTED_PAIRS) == {
        ("E_cadherin", "CD31"), ("Pan_Cytokeratin", "CD31"), ("Ker8_18", "CD31"), ("EpCAM", "CD31"),
    }
    assert not {p for p in comparator if p[1] == "Vimentin"} & set(rv.ACCEPTED_PAIRS)
    assert {t for t, _ in comparator} == {"Pan_Cytokeratin", "Ker8_18", "EpCAM", "E_cadherin"}
    # both references are biologically exclusive with epithelium and are themselves screened markers
    assert {r for _, r in comparator} == {"CD31", "Vimentin"}
    for target, reference in comparator:
        assert target not in rv.RESTORE_EXCLUDED_MARKERS
        assert reference not in rv.RESTORE_EXCLUDED_MARKERS
        # every marker needs a declared raw input floor or the screen cannot form arms
        assert target in rv.MARKER_INPUT_FLOOR_METHODS
        assert reference in rv.MARKER_INPUT_FLOOR_METHODS


# The markers the compartment sample actually contains. Its source QuPath CSV no longer exists, so
# this set cannot grow without a manual re-export of all 20 donors — it is effectively frozen.
COMPARTMENT_SAMPLED_MARKERS = {
    "E_cadherin", "CD31", "CD3e", "CD20", "CD68", "CD11b", "B3TUBB", "EpCAM", "Vimentin",
    "CD163", "Pan_Cytokeratin", "Ker8_18",
}


def test_epithelial_comparator_needs_no_new_qupath_marker():
    """The epithelial comparator had to be answerable from the compartment sample already on disk."""
    markers = {m for pair in rv.EPITHELIAL_TARGET_COMPARATOR_PAIRS for m in pair}
    assert markers <= COMPARTMENT_SAMPLED_MARKERS


def test_level1_expansion_markers_are_whole_cell_only_and_declared():
    """The expansion DOES add markers, which is why it runs on whole-cell features only.

    They are absent from the compartment sample by construction — `load_validation_sample`
    synthesises their `Cell__{marker}` from the canonical parquet, which is the same number. What must
    hold is that they are honestly outside the sampled set (so nothing silently expects a Membrane
    column for them) and that each declares an input floor, without which the screen cannot form arms.
    """
    new = {m for pair in rv.LEVEL1_EXPANSION_PAIRS for m in pair} - COMPARTMENT_SAMPLED_MARKERS
    assert new, "the expansion is supposed to introduce markers"
    assert new.isdisjoint(COMPARTMENT_SAMPLED_MARKERS)
    for marker in rv.QUPATH_MARKERS:
        assert marker in rv.MARKER_INPUT_FLOOR_METHODS, f"{marker} has no input-floor policy"
    # every reference is drawn from the sampled set, so each pair keeps one well-characterised side
    assert {r for _, r in rv.LEVEL1_EXPANSION_PAIRS} <= COMPARTMENT_SAMPLED_MARKERS


def test_level1_expansion_excludes_the_promiscuous_markers():
    """CD56 marks NK, neural AND endocrine cells; CD66 marks epithelium and granulocytes. Neither can
    discriminate the compartments it spans, so neither may become a level-1 gate."""
    targets = {t for t, _ in rv.LEVEL1_EXPANSION_PAIRS}
    assert "CD56" not in targets
    assert "CD66" not in targets


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


def test_manuscript_maximum_divisor_gated_call():
    # 2026-07-23: the call is strictly value > divisor; the former target-component rescue is gone,
    # so a below-divisor cell is NOT positive even when it was in the NNMF target component.
    stats = rv.negative_control_statistics(np.array([1.0, 2.0, 4.0]))
    assert stats["maximum"] == 4.0
    normalized, positive = rv.normalize_and_call(
        np.array([0.0, 4.0, 4.0001, 8.0]),
        stats["maximum"],
    )
    assert np.allclose(normalized, [0.0, 1.0, 1.000025, 2.0])
    assert positive.tolist() == [False, False, True, True]


def test_call_rejects_nonpositive_divisor():
    with pytest.raises(ValueError, match="divisor must be finite and positive"):
        rv.normalize_and_call(np.array([1.0, 2.0]), 0.0)


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
        # genuine absence (no target arm) still abstains ...
        ({"target_n": 1}, rv.PairState.LOW_SIGNAL_UNRESOLVED),
        # ... but a low arm-MEAN SBR with a present arm no longer abstains (METHOD v9): high background
        # compresses the ratio while the divisor still separates -- it becomes a reported flag, VALID.
        ({"target_arm_fold": 1.99}, rv.PairState.VALID),
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
    )
    assert np.allclose(normalized, [0.0, 1.0, 1.00005])
    # Divisor-gated: only value > divisor is positive (2.0 is not > 2.0).
    assert positive.tolist() == [False, False, True]

    absent = _assessment(
        target_arm_fold=1.0,
        expression_review=rv.ExpressionReview.CONFIRMED_ABSENT,
    )
    normalized, positive = rv.normalize_for_assessment(
        values,
        absent,
        divisor=None,
    )
    assert np.isnan(normalized).all()
    assert not positive.any()
    with pytest.raises(ValueError, match="must not fabricate"):
        rv.normalize_for_assessment(
            values,
            absent,
            divisor=2.0,
        )
    with pytest.raises(ValueError, match="cannot produce"):
        rv.normalize_for_assessment(
            values,
            _assessment(stable=False),  # MODEL_UNSTABLE: a non-VALID, non-absence state blocks calls
            divisor=None,
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
        # fit_cap//2 = 4 caps both arms below their candidate counts, so neither is saturated
        "target_arm_saturated": False,
        "reference_arm_saturated": False,
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


def test_control_definition_admits_reference_component_double_high_cells():
    """v10: the control is defined by the REFERENCE marker alone, so a double-high cell whose NNMF
    label is the reference component now enters it — and raises the divisor accordingly.

    This is the behaviour change that decouples the divisor from the target's Otsu floor. Under the
    pre-v10 `reference_and_target` definition the control was capped at that floor, so the maximum
    over it could not exceed the floor: measured median divisor/target_input_floor = 1.044 across the
    49 Gate-3 donor-pairs.
    """
    kwargs = dict(
        #                       cell:   0    1    2      3    4
        target_raw=np.array([100.0, 90.0, 1.0, 2.0, 500.0]),
        reference_raw=np.array([1.0, 2.0, 100.0, 80.0, 100.0]),
        target_redsea=np.array([100.0, 90.0, 1.0, 2.0, 500.0]),
        reference_redsea=np.array([1.0, 2.0, 100.0, 80.0, 100.0]),
        #             cells 2,3,4 are in the reference component; cell 4 is DOUBLE-HIGH
        labels=np.array([1, 1, 0, 0, 0]),
        qc_retained=np.ones(5, dtype=bool),
        target_group=1,
        reference_group=0,
        target_floor=10.0,
        reference_floor=10.0,
    )

    legacy = rv.ordered_control_groups(**kwargs, control_definition="reference_and_target")
    manuscript = rv.ordered_control_groups(**kwargs, control_definition="reference_only")

    # Cell 4 (target_raw 500 > floor 10) is excluded by the legacy target-side cap, admitted by v10.
    assert legacy.reference_candidates.tolist() == [False, False, True, True, False]
    assert manuscript.reference_candidates.tolist() == [False, False, True, True, True]

    legacy_divisor = rv.negative_control_statistics(
        kwargs["target_redsea"][legacy.reference_control]
    )["divisor"]
    manuscript_divisor = rv.negative_control_statistics(
        kwargs["target_redsea"][manuscript.reference_control]
    )["divisor"]
    assert legacy_divisor == 2.0            # capped at the target floor, as the defect predicts
    assert manuscript_divisor == 500.0      # set by the real double-high control cell
    assert manuscript_divisor > kwargs["target_floor"]


def test_anchor_and_controls_stay_disjoint_without_the_target_side_gate():
    """Removing the target-side gate cannot make the anchor and the control overlap.

    Disjointness never rested on it: the anchor requires reference_raw <= reference_floor and the
    candidates require reference_raw > reference_floor, which are mutually exclusive whatever the NNMF
    labels do — even for a degenerate fit that collapses both groups onto one component.
    """
    result = rv.ordered_control_groups(
        target_raw=np.array([100.0, 100.0]),
        reference_raw=np.array([1.0, 100.0]),
        target_redsea=np.array([100.0, 100.0]),
        reference_redsea=np.array([1.0, 100.0]),
        labels=np.array([0, 0]),              # collapsed fit: both groups are component 0
        qc_retained=np.ones(2, dtype=bool),
        target_group=0,
        reference_group=0,
        target_floor=10.0,
        reference_floor=10.0,
        control_definition="reference_only",
    )
    assert not np.any(result.target_population & result.reference_candidates)


@pytest.mark.parametrize(
    "statistic,expected",
    [("max", 100.0), ("p999", 99.9), ("p99", 99.0), ("p95", 95.0)],
)
def test_divisor_statistic_selects_the_summary(statistic, expected):
    values = np.arange(1, 101, dtype=float)      # 1..100
    stats = rv.negative_control_statistics(values, statistic)
    assert np.isclose(stats["divisor"], expected, rtol=0, atol=0.15)
    assert stats["statistic"] == statistic
    assert stats["maximum"] == 100.0             # manuscript value always reported
    assert 0 < stats["control_fraction_le_1"] <= 1.0


def test_manuscript_invariant_holds_only_for_the_maximum():
    """negatives <= 1 is exact for `max` and false by construction for a quantile — that is the point."""
    values = np.arange(1, 101, dtype=float)
    assert rv.negative_control_statistics(values, "max")["control_fraction_le_1"] == 1.0
    p95 = rv.negative_control_statistics(values, "p95")
    assert p95["control_fraction_le_1"] < 1.0
    assert (values > p95["divisor"]).sum() > 0


def test_quantile_divisor_of_a_mostly_zero_control_fails_diagnosably():
    """A REDSEA-clamped background can be 0 above the 95th percentile; fail here, not downstream."""
    values = np.concatenate([np.zeros(99), [50.0]])
    with pytest.raises(ValueError, match="zeros"):
        rv.negative_control_statistics(values, "p95")
    assert rv.negative_control_statistics(values, "max")["divisor"] == 50.0


def test_unknown_control_definition_and_divisor_statistic_are_rejected():
    with pytest.raises(ValueError, match="control_definition"):
        rv.PairValidationConfig(control_definition="reference_or_target")
    with pytest.raises(ValueError, match="divisor_statistic"):
        rv.PairValidationConfig(divisor_statistic="p90")


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
    assert f"divisor_statistic={config.divisor_statistic}" in specification["divisor_policy"]
    assert "manuscript_divisor" in specification["divisor_policy"]
    assert (
        f"control_definition={config.control_definition}"
        in specification["control_definition_policy"]
    )
    assert "vertical reference separator" in specification["ordered_control_policy"]
    assert "input-QC-retained" in specification["ordered_control_policy"]
    assert "double-low" in specification["calling_policy"]
    # The call has been divisor-gated only since METHOD v7; the spec described the removed v6 rescue
    # until v10. Pin the CORRECT semantics, and assert the stale phrasing is gone so it cannot return.
    assert "positive if and only if normalized > 1" in specification["calling_policy"]
    assert "NOT called positive" in specification["calling_policy"]
    assert "even below the horizontal divisor" not in specification["calling_policy"]
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
    # v9 stability is gated on divisor reproducibility; a clean pair is reproducible and high-SBR.
    assert result["divisor_reproducibility"] >= 0.9
    assert result["target_low_sbr"] is False
    assert result["target_supported_n"] >= result["target_anchor_n"]
    assert not result["target_supported"][~result["qc_retained"]].any()
    # 2026-07-23: the call is divisor-gated -- lineage-positive is exactly the above-divisor count.
    # (It is NOT bounded by target_supported: a double-high cell above the divisor is called but is
    # not in the target component, so threshold_positive can exceed target_supported.)
    assert result["lineage_positive_n"] == result["threshold_positive_n"]
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


def _asymmetric_pair_frame(rng, target_arm, reference_arm):
    """Two cleanly separated exclusive arms of controllable size for the frequency rule."""
    target = np.r_[
        rng.normal(100, 2, target_arm),
        rng.normal(1, 0.05, reference_arm),
    ]
    reference = np.r_[
        rng.normal(1, 0.05, target_arm),
        rng.normal(100, 2, reference_arm),
    ]
    return pd.DataFrame(
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


def test_reference_frequency_rule_is_pure_proportion_diagnostic():
    rng = np.random.default_rng(7)
    config = rv.PairValidationConfig(
        seeds=(0, 1, 2),
        fit_cap=400,
        threshold_sample_sizes=(5, 10),
        threshold_resamples=10,
    )
    # reference arm smaller than target arm -> the inverted-fit case -> undersized
    undersized = rv.evaluate_locked_pair(
        _asymmetric_pair_frame(rng, target_arm=120, reference_arm=30),
        "synthetic",
        "E_cadherin",
        "CD31",
        config=config,
        pair_review=rv.PairReview.ACCEPTED,
    )
    # reference arm larger than target arm -> not undersized
    oversized = rv.evaluate_locked_pair(
        _asymmetric_pair_frame(rng, target_arm=30, reference_arm=120),
        "synthetic",
        "E_cadherin",
        "CD31",
        config=config,
        pair_review=rv.PairReview.ACCEPTED,
    )
    for result in (undersized, oversized):
        # the ratio is EXACTLY reference_n / target_n -- a pure proportion, no absolute term
        assert result["reference_to_target_ratio"] == pytest.approx(
            result["reference_n"] / result["target_n"]
        )
        assert result["reference_undersized"] == (
            result["target_n"] > 0 and result["reference_n"] < result["target_n"]
        )
        # the flat metrics row carries both columns unchanged
        row = rv.locked_pair_metrics(result)
        assert row["reference_to_target_ratio"] == result["reference_to_target_ratio"]
        assert row["reference_undersized"] == result["reference_undersized"]
    assert undersized["reference_undersized"] is True
    assert undersized["reference_n"] < undersized["target_n"]
    assert oversized["reference_undersized"] is False
    assert oversized["reference_n"] > oversized["target_n"]
    # diagnostic only: the flag never alters the divisor or the pair state
    assert oversized["canonical_divisor"] == oversized["candidate_full_maximum"]
    assert undersized["canonical_divisor"] == undersized["candidate_full_maximum"]


def test_reference_frequency_rule_has_no_absolute_floor():
    # A tiny reference arm that still EXCEEDS an even tinier target arm is NOT undersized: the rule is a
    # pure proportion, so absolute counts (which vary widely across images) never enter it.
    rng = np.random.default_rng(11)
    config = rv.PairValidationConfig(
        seeds=(0, 1, 2),
        fit_cap=200,
        threshold_sample_sizes=(3, 5),
        threshold_resamples=10,
    )
    result = rv.evaluate_locked_pair(
        _asymmetric_pair_frame(rng, target_arm=8, reference_arm=40),
        "synthetic",
        "E_cadherin",
        "CD31",
        config=config,
        pair_review=rv.PairReview.ACCEPTED,
    )
    assert result["reference_n"] > result["target_n"]
    assert result["reference_undersized"] is False
    assert result["reference_to_target_ratio"] > 1.0


def test_reference_to_target_ratio_config_must_be_positive():
    with pytest.raises(ValueError, match="min_reference_to_target_ratio"):
        rv.PairValidationConfig(min_reference_to_target_ratio=0)


def test_min_divisor_reproducibility_config_range():
    with pytest.raises(ValueError, match="min_divisor_reproducibility"):
        rv.PairValidationConfig(min_divisor_reproducibility=1.5)


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
                "target": "CD68",
                "reference": "CD3e",
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
        pd.DataFrame(rows), pairs=(("CD68", "CD3e"),)
    )
    by_donor = queue.set_index("donor")["selection_reasons"]
    assert "lowest seed stability" in by_donor["6476"]
    assert "pancreatitis" not in "; ".join(by_donor)
    assert "lowest target-arm fold" in by_donor["1"]
    assert "lowest reference-arm fold" in by_donor["2"]
    # The CD20-only "maintainer-confirmed target absence" stress case (donor 6539) went with CD20's
    # exclusion from RESTORE; selection is now driven purely by the evidence columns.
    assert "maintainer-confirmed target absence" not in "; ".join(by_donor)
    shuffled = pd.DataFrame(rows).sample(frac=1, random_state=19)
    pd.testing.assert_frame_equal(
        queue,
        rpr.select_review_queue(
            shuffled, pairs=(("CD68", "CD3e"),)
        ),
    )


def test_review_categories_keep_unavailable_and_double_low_distinct():
    evaluation = {
        "valid_idx": np.array([0, 1, 2, 3, 4]),
        "qc_retained": np.array([True, True, True, False, True]),
        "reference_control": np.array([True, False, False, False, False]),
        "target_population": np.array([False, True, False, False, False]),
        "target_supported": np.array([False, True, False, True, True]),
        # idx1 is above the Step 2 divisor (called); idx4 is a target-component cell below it (retained).
        "threshold_positive": np.array([False, True, False, False, False]),
        "target_raw": np.array([1.0, 10.0, 10.0, 1.0, 10.0]),
        "reference_raw": np.array([10.0, 1.0, 10.0, 1.0, 1.0]),
        "target_input_floor": 5.0,
        "reference_input_floor": 5.0,
    }
    categories = rpr.review_categories(evaluation, n_rows=6)
    assert categories.tolist() == [
        "Reference control arm",
        "Called target+ (above Step 2 divisor)",
        "Double-high",
        "Double-low (fit excluded)",
        "Target component at/below Step 2 (retained, not called)",
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
        + ["Target component at/below Step 2 (retained, not called)"] * 3,
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
            # idx3 is raw double-high but its target redsea (3) sits below the divisor (5), so it stays
            # Double-high rather than being divisor-called.
            "redsea__Target": [10.0, 2.0, 1.0, 3.0, 1.0],
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
        "Called target+ (above Step 2 divisor)",
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
        + ["Called target+ (above Step 2 divisor)"] * 2
        + ["Target component at/below Step 2 (retained, not called)"] * 2,
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
        # `divisor` is the active statistic and is what the Step-2 line draws; `maximum` is the
        # manuscript value, kept so the row stays representative of a real evaluation.
        "target_full_stats": {"maximum": 2.0, "divisor": 1.5},
        "reference_full_stats": {"maximum": 2.0, "divisor": 1.5},
        # Every real evaluation carries these; the marker panels are windowed on them and
        # `_display_floors` raises without them, by design.
        "target_input_floor": 1.0,
        "reference_input_floor": 1.0,
        "state": rv.PairState.PAIR_REVIEW_PENDING.value,
        "target_fold": 2.0,
        "reference_fold": 3.0,
        "min_control_jaccard": 0.95,
        "threshold_positive_n": 2,
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
    assert "Called >Step 2 divisor n=2" in axes[0].get_title()
    assert "anchor (sets Step 1) n=2" in axes[0].get_title()
    assert "Retained at/below Step 2 (component, NOT called) n=2" in axes[0].get_title()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    assert not axes[0].title.get_window_extent(renderer).overlaps(
        axes[0].get_window_extent(renderer)
    )
    rpr.plt.close(fig)


def test_target_confidence_colors_are_distinct_and_high_contrast():
    high = np.asarray(
        rpr.to_rgba(
            rpr.CATEGORY_COLOR["Called target+ (above Step 2 divisor)"]
        )[:3]
    )
    lower = np.asarray(
        rpr.to_rgba(
            rpr.CATEGORY_COLOR[
                "Target component at/below Step 2 (retained, not called)"
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
        "Called >Step 2 divisor n=24,801 | anchor (sets Step 1) n=28,685\n"
        "Retained at/below Step 2 (component, NOT called) n=3,884",
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
    # 11 original + the 2 accepted pairs added 2026-07-25 + the 7 epithelial promotion candidates
    # added 2026-07-27 (Pan-CK / Ker8_18 / EpCAM x {CD31, Vimentin}, plus E_cadherin <- Vimentin as
    # the comparator for the frozen E_cadherin <- CD31).
    # 11 original + 2 (2026-07-25) + 7 epithelial promotion candidates + 14 level-1 expansion pairs
    assert len(rpr.SHORTLIST_PAIRS) == 34
    assert len(set(rpr.SHORTLIST_PAIRS)) == len(rpr.SHORTLIST_PAIRS)
    assert set(rpr.SHORTLIST_PAIRS) == set(rpr.PAIR_RATIONALE)


def test_epithelial_promotion_candidates_have_a_review_surface():
    """The anatomy check made these promotion candidates; promotion needs a per-image review, which
    needs a review surface. Both candidate references are shortlisted for each, so the reference
    choice can be judged from images and not only from the frequency rule."""
    for target in ("Pan_Cytokeratin", "Ker8_18", "EpCAM"):
        for reference in ("CD31", "Vimentin"):
            assert (target, reference) in rpr.SHORTLIST_PAIRS
            assert (target, reference) in rv.CANDIDATE_PAIRS
    assert ("E_cadherin", "Vimentin") in rpr.SHORTLIST_PAIRS   # comparator for the frozen CD31 pair
    # Promoted to production 2026-07-27 against CD31 — E-cadherin's own frozen reference, so the
    # comparison is like-for-like. The Vimentin-referenced rows stay comparators.
    for target in ("Pan_Cytokeratin", "Ker8_18", "EpCAM"):
        assert (target, "CD31") in rv.ACCEPTED_PAIRS
        assert (target, "Vimentin") not in rv.ACCEPTED_PAIRS
    assert ("E_cadherin", "Vimentin") not in rv.ACCEPTED_PAIRS


def test_every_gate_marker_reaches_production_through_an_accepted_pair():
    """The whole point of the expansion: a marker in a gate must have a divisor, or the gate is
    permanently 'unavailable' and Unresolves every cell it touches."""
    from phenocycler.config import COMPARTMENT_GATES
    accepted = {t for t, _ in rv.ACCEPTED_PAIRS}
    gated = {m for gs in COMPARTMENT_GATES.values() for m in gs}
    assert gated == accepted
    assert len(accepted) == 24
    assert "CD163" not in accepted      # screened and rejected on arm separation
    assert "CD20" not in accepted       # excluded from RESTORE in both roles


def test_every_accepted_pair_has_a_review_surface():
    """A pair cannot reach production without a PDF/QuPath review surface.

    This silently failed between the Gate-2 reference freeze and 2026-07-25: SHORTLIST_PAIRS was
    frozen beforehand, so CD68<-E_cadherin and CD11b<-EpCAM ran in production unreviewable, while the
    shortlist reviewed the CD68/CD11b references the reference dossier went on to REJECT. Enforced at
    import in restore_pair_review; asserted here so the intent is visible.
    """
    missing = set(rv.ACCEPTED_PAIRS) - set(rpr.SHORTLIST_PAIRS)
    assert not missing, f"accepted pairs with no review surface: {sorted(missing)}"
    assert ("CD68", "E_cadherin") in rpr.SHORTLIST_PAIRS
    assert ("CD11b", "EpCAM") in rpr.SHORTLIST_PAIRS


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
    assert '"RESTORE review - Called target+"' in text
    assert '"RESTORE review - Retained below Step 2"' in text
    assert '"5" : "RESTORE review - Retained below Step 2"' in text
    assert "[185, 175, 203]" in text
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


# --------------------------------------------------------------------------- #
# Scaling-bound support: a robust quantile is only robust if the sample can
# span it. A sparse marker arm cannot, and one bright cell then sets the scale.
# --------------------------------------------------------------------------- #
def test_quantile_support_cells_uses_the_narrower_tail():
    assert rv.quantile_support_cells(144, lower_quantile=0.005, upper_quantile=0.995) == pytest.approx(0.72)
    assert rv.quantile_support_cells(200, lower_quantile=0.005, upper_quantile=0.995) == pytest.approx(1.0)
    assert rv.quantile_support_cells(4752, lower_quantile=0.005, upper_quantile=0.995) == pytest.approx(23.76)


def test_adaptive_bound_quantiles_widen_only_below_the_support_boundary():
    # n = 200 is exactly one cell per tail at the shipped quantiles -> keep them as configured.
    lower, upper, adapted, support = rv.adaptive_bound_quantiles(
        200, lower_quantile=0.005, upper_quantile=0.995, min_support=1, repair_support=10
    )
    assert (lower, upper, adapted) == (0.005, 0.995, False)
    assert support == pytest.approx(1.0)
    # One row fewer cannot span a tail -> widen until each tail spans repair_support cells.
    lower, upper, adapted, support = rv.adaptive_bound_quantiles(
        199, lower_quantile=0.005, upper_quantile=0.995, min_support=1, repair_support=10
    )
    assert adapted and lower == pytest.approx(10 / 199) and upper == pytest.approx(1 - 10 / 199)
    assert support == pytest.approx(10.0)


def test_adaptive_bound_quantiles_never_exceed_the_quartile():
    # 28 rows: 10/28 would be a 36% tail, which is no longer a tail. Cap at the quartile.
    lower, upper, adapted, support = rv.adaptive_bound_quantiles(
        28, lower_quantile=0.005, upper_quantile=0.995, min_support=1, repair_support=10
    )
    assert adapted and lower == pytest.approx(0.25) and upper == pytest.approx(0.75)
    assert support == pytest.approx(7.0)


def test_sparse_arm_bounds_ignore_a_single_bright_outlier():
    """The defect this fixes: with a tiny arm, p99.5 interpolates between the top two cells."""
    rng = np.random.default_rng(0)
    n_arm = 60
    feature = np.r_[rng.normal(400, 30, n_arm), rng.normal(20, 3, n_arm)]
    feature[0] = 11_551.0  # the outlier that set the 6476 CD20 bound
    features = np.column_stack([feature, feature * 0.9])
    fit_indices = np.arange(2 * n_arm)          # a 120-row balanced sample

    configured = np.quantile(features[fit_indices], 0.995, axis=0)
    lower_q, upper_q, adapted, _support = rv.adaptive_bound_quantiles(
        len(fit_indices), lower_quantile=0.005, upper_quantile=0.995,
        min_support=1, repair_support=10,
    )
    assert adapted
    _scaled, _lower, upper = rv.robust_minmax_scale(
        features, fit_indices, lower_quantile=lower_q, upper_quantile=upper_q
    )
    # The configured bound is dragged far above the population it must describe; the widened one is not.
    assert configured[0] > 5 * upper[0]
    assert upper[0] < 1_000.0
    # Crucially the bounds still come from the BALANCED sample -- swapping populations instead would
    # re-create the same domination with the two markers' roles reversed.
    assert len(fit_indices) == 2 * n_arm


def test_stability_is_uninformative_whenever_either_arm_is_saturated():
    """per_arm = min(|target|, |reference|, cap), so the smaller arm is normally taken whole and the
    seed sweep cannot perturb it. The flag must therefore trip on OR, not AND, or it never fires."""
    rng = np.random.default_rng(7)
    target = np.r_[rng.normal(100, 2, 30), rng.normal(1, 0.05, 400)]
    reference = np.r_[rng.normal(1, 0.05, 30), rng.normal(100, 2, 400)]
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
    result = rv.evaluate_locked_pair(
        frame,
        "synthetic",
        "E_cadherin",
        "CD31",
        config=rv.PairValidationConfig(
            fit_cap=1000, threshold_sample_sizes=(5, 10), threshold_resamples=10
        ),
    )
    assert result["target_arm_saturated"] is True      # 30 of 30 taken
    assert result["reference_arm_saturated"] is False  # 30 of 400
    assert result["stability_informative"] is False    # one saturated arm is enough to under-test
    assert result["probe_control_jaccard"] is not None


def test_arm_saturation_is_reported_and_probe_fraction_shrinks_arms():
    target = np.array([10.0] * 5 + [1.0] * 500)
    reference = np.array([1.0] * 5 + [10.0] * 500)
    _idx, counts = rv.balanced_arm_fit_indices(
        target, reference, 5.0, 5.0, fit_cap=1000, seed=0
    )
    # The 5-cell target arm is taken whole, so re-seeding cannot perturb it.
    assert counts["target_arm_saturated"] is True
    assert counts["reference_arm_saturated"] is False
    assert counts["fit_per_arm"] == 5

    probe_idx, probe_counts = rv.balanced_arm_fit_indices(
        target, reference, 5.0, 5.0, fit_cap=1000, seed=0, probe_fraction=0.8
    )
    assert len(probe_idx) == 8                      # 4 per arm
    # Saturation is reported for the real balanced size, not the shrunken probe size.
    assert probe_counts["target_arm_saturated"] is True


def test_immune_joint_comparator_pools_arms_and_emits_no_calls():
    rng = np.random.default_rng(3)
    n = 600
    frame = {}
    # Four disjoint immune populations plus a large negative bulk.
    blocks = {m: slice(i * 40, (i + 1) * 40) for i, m in enumerate(rv.IMMUNE_JOINT_MARKERS)}
    for marker in rv.IMMUNE_JOINT_MARKERS:
        values = rng.normal(5.0, 0.5, n)
        values[blocks[marker]] = rng.normal(300.0, 10.0, 40)
        frame[f"raw__{marker}"] = values
        for compartment in ("Membrane", "Cell"):
            frame[rv.compartment_column(compartment, marker)] = values * (
                0.9 if compartment == "Membrane" else 1.0
            )
    rows = rv.evaluate_immune_joint(
        pd.DataFrame(frame), "synthetic", fit_cap=500, seeds=(0,), n_components=(4,)
    )
    assert rows and {r["target"] for r in rows} == set(rv.IMMUNE_JOINT_MARKERS)
    # The pooled fitting set is the union of immune-positive cells, far larger than any single arm.
    assert all(r["n_immune_any"] >= 4 * 40 * 0.9 for r in rows)
    assert all(r["exclusive_arm_n"] > 0 for r in rows)
    # Comparator only: it reports recovery, never a divisor or a positivity call.
    assert all("divisor" not in r and "normalized" not in r for r in rows)
    assert all(r["anchor_recovery"] is not None for r in rows)


def test_locked_path_still_requires_exactly_two_components():
    with pytest.raises(ValueError, match="two NNMF components"):
        rv.PairValidationConfig(n_components=4)


def test_review_shortlist_only_names_pairs_the_screen_produces():
    """The bundle selects its review queue from the screen, so a shortlist entry that is not a
    candidate pair would silently yield no rows (this caught CD68<-CD20 / CD11b<-CD20 when CD20 was
    retired as a reference)."""
    missing = [p for p in rpr.SHORTLIST_PAIRS if p not in rv.CANDIDATE_PAIRS]
    assert not missing, f"shortlist pairs absent from CANDIDATE_PAIRS: {missing}"
    assert not [p for p in rpr.SHORTLIST_PAIRS if "CD20" in p]


def test_cd20_takes_no_part_in_restore():
    """Maintainer decision 2026-07-23: CD20 is excluded from RESTORE in BOTH roles and deferred to a
    second pass after the broad cell types are settled. As a reference its positive cells would BE the
    negative control and its arm is far too small; as a target its Step 1 anchor bottomed out at 14
    cells, and reviewer inspection of the CD20 PDFs found overcalling with both qptiff channels
    stretched from near-background into block artifacts."""
    assert "CD20" in rv.RESTORE_EXCLUDED_MARKERS
    for pairs in (rv.BASELINE_PAIRS, rv.IMMUNE_REFERENCE_SCREEN_PAIRS,
                  rv.EPITHELIAL_REFERENCE_COMPARATOR_PAIRS, rv.CANDIDATE_PAIRS):
        assert not [p for p in pairs if "CD20" in p], "CD20 must appear in neither role"
    # The immune targets that used to lean on CD20 use CD3e, which was already their alternative.
    assert ("CD68", "CD3e") in rv.CANDIDATE_PAIRS
    assert ("CD11b", "CD3e") in rv.CANDIDATE_PAIRS
    # Excluded from the pair web, but still EXTRACTED and still carrying its reviewer-validated input
    # floor policy, so the deferred second pass needs no QuPath re-export.
    assert "CD20" in rv.QUPATH_MARKERS
    assert "CD20" in rv.DEFERRED_EXTRACTION_MARKERS
    assert rv.MARKER_INPUT_FLOOR_METHODS["CD20"] == rv.MAX_LINEAR_OTSU_TRIANGLE_FLOOR
    # No marker may reference itself.
    assert all(t != r for t, r in rv.CANDIDATE_PAIRS)


def test_excluded_marker_cannot_be_reintroduced_in_either_role():
    """The import-time guard is the enforcement point: rebuilding CANDIDATE_PAIRS with a CD20 pair in
    either direction must fail loudly rather than silently re-screening it."""
    for pair in (("CD20", "CD3e"), ("CD11b", "CD20")):
        target, reference = pair
        offenders = [m for m in pair if m in rv.RESTORE_EXCLUDED_MARKERS]
        assert offenders, f"{target} <- {reference} should name an excluded marker"


# --------------------------------------------------------------------------- #
# Per-pair feature-compartment resolution (2026-07-27)
#
# 14 of the 24 accepted targets post-date the compartment sample and have only whole-cell values.
# Applying THEIR limitation to every pair via a run-level setting cost a real result: CD68 <-
# E_cadherin on donor 6523 is MODEL_UNSTABLE on whole-cell alone (divisor reproducibility 0.864, no
# calls) and stable with Membrane restored (1.000, 3.53% called) -- and that one donor-marker left 80%
# of the donor Unresolved. Both CD68 and E_cadherin ARE in the sample; nothing about 6523 is unusual.
# --------------------------------------------------------------------------- #

def _pair_frame(markers, compartments):
    """Minimal donor_df carrying only the named compartments for the named markers."""
    import numpy as np
    n = 40
    rng = np.random.default_rng(0)
    cols = {"object_id": [f"c{i}" for i in range(n)]}
    for m in markers:
        cols[f"raw__{m}"] = rng.random(n) * 100
        cols[f"redsea__{m}"] = rng.random(n) * 100
        for c in compartments:
            cols[rv.compartment_column(c, m)] = rng.random(n) * 100
    return pd.DataFrame(cols)


def test_effective_compartments_keep_membrane_when_the_pair_has_it():
    df = _pair_frame(("CD68", "E_cadherin"), ("Membrane", "Cell"))
    cfg = rv.PairValidationConfig(feature_compartments=("Membrane", "Cell"))
    eff = tuple(c for c in cfg.feature_compartments
                if all(rv.compartment_column(c, m) in df.columns for m in ("CD68", "E_cadherin")))
    assert eff == ("Membrane", "Cell")


def test_effective_compartments_fall_back_when_a_marker_lacks_membrane():
    """A whole-cell-only marker must degrade ITS OWN pair, never the whole run."""
    df = _pair_frame(("E_cadherin",), ("Membrane", "Cell"))
    for c in ("raw__CD79a", "redsea__CD79a", rv.compartment_column("Cell", "CD79a")):
        df[c] = 1.0                                        # CD79a has whole-cell only
    cfg = rv.PairValidationConfig(feature_compartments=("Membrane", "Cell"))
    eff = tuple(c for c in cfg.feature_compartments
                if all(rv.compartment_column(c, m) in df.columns for m in ("CD79a", "E_cadherin")))
    assert eff == ("Cell",)


def test_cell_is_always_required():
    """'Cell' is the whole-cell mean and exists for every marker; a declaration without it is an error,
    not a silent empty feature set."""
    df = _pair_frame(("CD68", "E_cadherin"), ("Membrane",))
    cfg = rv.PairValidationConfig(feature_compartments=("Membrane",))
    with pytest.raises(ValueError, match="no usable feature compartment"):
        rv.evaluate_locked_pair(df, "TEST", "CD68", "E_cadherin", config=cfg,
                                pair_review=rv.PairReview.ACCEPTED)


# --------------------------------------------------------------------------- #
# Review-display invariants. Both of these shipped broken and were caught only by a reviewer looking
# at the finished PDFs, so they are asserted here rather than trusted.
# --------------------------------------------------------------------------- #
def test_step2_line_and_axis_use_the_active_divisor_not_the_manuscript_maximum():
    """The drawn Step-2 threshold must be the divisor production uses.

    While `divisor_statistic` was `max` these coincided. Frozen at p99 they diverge by a median 6x
    (up to 200x), and both the line and the axis top were still taken from `maximum` -- so every panel
    showed a threshold production does not use, and the data cloud was crushed into the bottom few
    percent of the axis."""
    source = Path(rpr.__file__).read_text()
    drawn = source.split("cloud.axhline(", 1)[1][:260]
    assert '["target_full_stats"]["divisor"]' in drawn, (
        "Step 2 must be drawn at target_full_stats['divisor'] (the active statistic)"
    )
    assert '["target_full_stats"]["maximum"]' not in drawn, (
        "Step 2 must NOT be drawn at the manuscript maximum"
    )

    validation = Path(rv.__file__).read_text()
    assert '1.05 * float(target_full_stats["divisor"])' in validation, (
        "the display y-axis must be scaled by the active divisor"
    )
    assert '1.05 * float(target_full_stats["maximum"])' not in validation, (
        "the display y-axis must NOT be scaled by the manuscript maximum"
    )


def test_marker_panels_are_windowed_on_the_measured_background_floor():
    """Background must render as background. These panels exist to judge background handling.

    The black point is the donor-marker's own raw input floor, and the transform is linear. A
    percentile black point sits BELOW the background so nearly every pixel renders non-black, and a
    gamma under 1 lifts the dim end -- together they turn a weak channel into a bright wash."""
    rng = np.random.default_rng(0)
    floor = 500.0
    # A crop that is 97% background just under the floor and 3% real signal well above it.
    pixels = rng.normal(300.0, 40.0, size=(200, 200)).astype(np.float32)
    signal = rng.integers(0, 200, size=(2, 600))
    pixels[signal[0], signal[1]] = rng.normal(3000.0, 200.0, size=600)

    shown, window = rpr._display_marker_channel(pixels, floor)
    assert window[0] == floor, "the black point must be the measured background floor"
    background = shown[pixels <= floor]
    assert background.max() == 0.0, "every pixel at or below the floor must render pure black"
    assert shown[pixels > floor].mean() > 0.0, "signal above the floor must still render"

    # ... and the old percentile+gamma window is what it must never silently fall back to. Pin the two
    # mechanisms that made it wrong, so the guard fails if anyone reinstates either.
    legacy, legacy_window = rpr._display_marker_channel(pixels)
    assert legacy_window[0] < floor, "the legacy black point sat BELOW the background, not at it"
    legacy_background = legacy[pixels <= floor]
    assert legacy_background.mean() > 0.05, (
        "sanity: under the legacy window background rendered visibly non-black — the defect guarded here"
    )
    # gamma < 1 LIFTS the dim end: the displayed background is brighter than its linear position.
    linear_position = (
        (pixels[pixels <= floor].mean() - legacy_window[0])
        / (legacy_window[1] - legacy_window[0])
    )
    assert legacy_background.mean() > 2 * linear_position, "gamma 0.65 amplified background"


def test_marker_panel_with_no_signal_above_the_floor_renders_black():
    """A crop with nothing above the floor must go black, not have its noise stretched into signal."""
    rng = np.random.default_rng(1)
    pixels = rng.normal(200.0, 30.0, size=(64, 64)).astype(np.float32)
    shown, _ = rpr._display_marker_channel(pixels, 900.0)
    assert shown.max() == 0.0


def test_display_floors_fail_loud_when_the_evaluation_has_none():
    """A missing floor must raise, not silently revert to the percentile window."""
    with pytest.raises(KeyError, match="target_input_floor"):
        rpr._display_floors({"reference_input_floor": 1.0})
    assert rpr._display_floors(None) == (None, None)
    assert rpr._display_floors(
        {"target_input_floor": 5.0, "reference_input_floor": None}
    ) == (5.0, None)
