from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.marker_calibration import MAX_CONTROL_CAP
from phenocycler.marker_registry import load_registry
from phenocycler.reference_controls import (
    MASK_COLUMNS,
    METHOD_VERSION,
    ReferenceControlConfig,
    ReferenceControlStatus,
    build_reference_control,
    build_reference_controls,
    load_reference_masks,
    reference_masks_as_arrays,
    reference_masks_from_wide,
    reference_masks_to_wide,
    registry_references,
    spatially_balanced_rank_indices,
)


def _config(**overrides):
    values = {
        "seed_cap": 40,
        "min_controls": 20,
        "spatial_bins": 2,
        "random_seed": 5,
    }
    values.update(overrides)
    return ReferenceControlConfig(**values)


def _cells(n=100):
    index = np.arange(n)
    frame = pd.DataFrame(
        {
            "object_id": [f"cell-{value:04d}" for value in index],
            "X_centroid": (index % 10).astype(float),
            "Y_centroid": (index // 10).astype(float),
            "qc_estimation_eligible": np.ones(n, dtype=bool),
            "Cell__CD31": 1.0 + index,
            "Cell__E_cadherin": 2.0 + index[::-1],
            "Cell__Pan_Cytokeratin": 3.0 + (index * 7) % n,
            "Cell__EpCAM": 4.0 + (index * 11) % n,
            # An unrelated target column must never enter seed selection.
            "Cell__CD3e": np.linspace(0, 1_000, n),
        }
    )
    return frame


def test_config_defaults_match_calibrator_minimum_and_hard_cap():
    config = ReferenceControlConfig()
    assert config.seed_cap == MAX_CONTROL_CAP == 20_000
    assert config.min_controls == 10_000
    with pytest.raises(ValueError, match="seed_cap"):
        ReferenceControlConfig(seed_cap=MAX_CONTROL_CAP + 1)
    with pytest.raises(ValueError, match="min_controls"):
        ReferenceControlConfig(seed_cap=10, min_controls=11)


def test_registry_reference_roster_is_small_explicit_graph_projection():
    registry = load_registry()
    references = registry_references(registry)
    assert set(references) == {
        "CD31",
        "E_cadherin",
        "Pan_Cytokeratin",
        "EpCAM",
    }
    assert all(reference in registry.by_name for reference in references)


def test_single_spatial_bin_is_exact_fixed_top_rank():
    n = 50
    intensity = np.arange(1, n + 1, dtype=float)
    ids = np.array([f"id-{value:03d}" for value in range(n)])
    selected, tiles = spatially_balanced_rank_indices(
        intensity,
        np.ones(n, bool),
        np.arange(n, dtype=float),
        np.zeros(n),
        ids,
        cap=10,
        spatial_bins=1,
        salt="fixed",
    )
    assert set(ids[selected]) == set(ids[-10:])
    assert intensity[selected].tolist() == list(range(50, 40, -1))
    assert set(tiles) == {0}


def test_spatial_round_robin_takes_highest_ranks_from_every_tile():
    # Four equal tiles, eight cells per tile. Intensities are intentionally
    # much higher in tile 0; spatial balance must still take 3 from each.
    x = np.repeat([0.0, 10.0, 0.0, 10.0], 8)
    y = np.repeat([0.0, 0.0, 10.0, 10.0], 8)
    tile = (x >= 5).astype(int) + 2 * (y >= 5)
    within = np.tile(np.arange(8), 4)
    intensity = 1_000.0 * (4 - tile) + within
    ids = np.array([f"id-{value:03d}" for value in range(len(x))])
    selected, selected_tiles = spatially_balanced_rank_indices(
        intensity,
        np.ones(len(x), bool),
        x,
        y,
        ids,
        cap=12,
        spatial_bins=2,
        salt="fixed",
    )
    assert np.bincount(selected_tiles, minlength=4).tolist() == [3, 3, 3, 3]
    for tile_id in range(4):
        expected = set(np.flatnonzero((tile == tile_id) & (within >= 5)))
        assert set(selected[tile[selected] == tile_id]) == expected


def test_build_all_reference_masks_and_audits_selected_policy():
    cells = _cells()
    registry = load_registry()
    result = build_reference_controls(
        cells,
        donor_id="D1",
        registry=registry,
        config=_config(),
    )
    references = registry_references(registry)
    assert len(result.masks_long) == len(cells) * len(references)
    assert tuple(result.masks_long.columns) == MASK_COLUMNS
    assert set(result.audit["reference"]) == set(references)
    assert set(result.audit["status"]) == {ReferenceControlStatus.VALID.value}
    assert set(result.audit["n_selected"]) == {40}
    assert set(result.audit["method_version"]) == {METHOD_VERSION}
    assert set(result.audit["spillover_policy"]) == {"redsea_corrected"}
    assert set(result.audit["compartment"]) == {"Cell"}
    assert (
        result.masks_long.groupby("reference")["reference_positive"].sum()
        == 40
    ).all()
    arrays = result.as_arrays(cells["object_id"].to_numpy())
    assert set(arrays) == set(references)
    assert all(mask.dtype == bool and mask.sum() == 40 for mask in arrays.values())

    wide = reference_masks_to_wide(result.masks_long, registry=registry)
    assert len(wide) == len(cells)
    assert not wide["object_id"].duplicated().any()
    roundtrip = reference_masks_from_wide(wide, registry=registry)
    left = result.masks_long.sort_values(
        ["reference", "object_id"], ignore_index=True
    )
    right = roundtrip.sort_values(
        ["reference", "object_id"], ignore_index=True
    )
    pd.testing.assert_frame_equal(left, right, check_dtype=False)


def test_selection_is_target_independent_and_row_order_invariant():
    cells = _cells()
    first = build_reference_controls(
        cells,
        donor_id="D1",
        references=["E_cadherin"],
        config=_config(),
    )
    changed = cells.copy()
    changed["Cell__CD3e"] = np.random.default_rng(44).normal(size=len(changed))
    permutation = np.random.default_rng(8).permutation(len(changed))
    second = build_reference_controls(
        changed.iloc[permutation].reset_index(drop=True),
        donor_id="D1",
        references=["E_cadherin"],
        config=_config(),
    )
    first_selected = first.masks_long.loc[
        first.masks_long["reference_positive"], ["object_id", "seed_rank"]
    ].sort_values("seed_rank", ignore_index=True)
    second_selected = second.masks_long.loc[
        second.masks_long["reference_positive"], ["object_id", "seed_rank"]
    ].sort_values("seed_rank", ignore_index=True)
    pd.testing.assert_frame_equal(first_selected, second_selected)
    assert (
        first.audit.loc[0, "selection_fingerprint"]
        == second.audit.loc[0, "selection_fingerprint"]
    )


def test_qc_zero_missing_and_coordinates_exclude_candidates():
    cells = _cells()
    # The brightest values would otherwise win.
    cells.loc[99, "qc_estimation_eligible"] = False
    cells.loc[98, "Cell__CD31"] = 0.0
    cells.loc[97, "Cell__CD31"] = np.nan
    cells.loc[96, "X_centroid"] = np.nan
    result = build_reference_controls(
        cells,
        donor_id="D1",
        references=["CD31"],
        config=_config(),
    )
    selected = set(
        result.masks_long.loc[
            result.masks_long["reference_positive"], "object_id"
        ]
    )
    assert not {
        "cell-0099",
        "cell-0098",
        "cell-0097",
        "cell-0096",
    } & selected
    assert result.audit.loc[0, "n_control_candidates"] == 96


def test_fewer_than_minimum_controls_fails_closed():
    cells = _cells(n=30)
    cells.loc[19:, "Cell__EpCAM"] = 0.0
    mask, audit = build_reference_control(
        cells["Cell__EpCAM"],
        cells["qc_estimation_eligible"],
        cells["X_centroid"],
        cells["Y_centroid"],
        cells["object_id"],
        donor_id="D1",
        reference="EpCAM",
        config=_config(seed_cap=25, min_controls=20),
    )
    assert audit["status"] == ReferenceControlStatus.INSUFFICIENT_CONTROLS.value
    assert audit["n_control_candidates"] == 19
    assert audit["n_selected"] == 0
    assert not mask["reference_positive"].any()
    assert mask["seed_rank"].isna().all()


def test_external_masks_load_only_through_versioned_schema():
    cells = _cells()
    generated = build_reference_controls(
        cells,
        donor_id="D1",
        references=["CD31"],
        config=_config(),
    ).masks_long
    external = generated.copy()
    external["source"] = "external"
    external["method_version"] = "external-reference-controls-v1"
    external = external.sample(frac=1, random_state=3).reset_index(drop=True)
    loaded = load_reference_masks(
        external,
        donor_id="D1",
        object_ids=cells["object_id"],
        references=["CD31"],
        config=_config(),
    )
    assert loaded.audit.loc[0, "status"] == ReferenceControlStatus.VALID.value
    assert loaded.audit.loc[0, "source"] == "external"
    assert loaded.audit.loc[0, "n_selected"] == 40

    bad_version = external.copy()
    bad_version["method_version"] = "manual"
    with pytest.raises(ValueError, match="explicit method version"):
        load_reference_masks(
            bad_version,
            donor_id="D1",
            object_ids=cells["object_id"],
            references=["CD31"],
            config=_config(),
        )
    annotated = external.assign(cell_type="Endothelial")
    with pytest.raises(ValueError, match="annotation"):
        load_reference_masks(
            annotated,
            donor_id="D1",
            object_ids=cells["object_id"],
            references=["CD31"],
            config=_config(),
        )


def test_external_masks_require_exact_object_universe_and_fail_small_set_closed():
    cells = _cells()
    generated = build_reference_controls(
        cells,
        donor_id="D1",
        references=["CD31"],
        config=_config(),
    ).masks_long
    external = generated.copy()
    external["source"] = "external"
    external["method_version"] = "external-reference-controls-v2"
    with pytest.raises(ValueError, match="object universe"):
        load_reference_masks(
            external.iloc[:-1],
            donor_id="D1",
            object_ids=cells["object_id"],
            references=["CD31"],
            config=_config(),
        )

    positive_rows = external.index[external["reference_positive"]]
    keep = positive_rows[:10]
    external["reference_positive"] = False
    external["seed_rank"] = pd.NA
    external.loc[keep, "reference_positive"] = True
    external.loc[keep, "seed_rank"] = np.arange(1, len(keep) + 1)
    loaded = load_reference_masks(
        external,
        donor_id="D1",
        object_ids=cells["object_id"],
        references=["CD31"],
        config=_config(),
    )
    assert (
        loaded.audit.loc[0, "status"]
        == ReferenceControlStatus.INSUFFICIENT_CONTROLS.value
    )
    assert loaded.audit.loc[0, "n_selected"] == 0
    assert not loaded.masks_long["reference_positive"].any()


def test_array_adapter_aligns_to_requested_object_order():
    cells = _cells()
    result = build_reference_controls(
        cells,
        donor_id="D1",
        references=["EpCAM"],
        config=_config(),
    )
    reversed_ids = cells["object_id"].to_numpy()[::-1]
    arrays = reference_masks_as_arrays(
        result.masks_long,
        object_ids=reversed_ids,
        donor_id="D1",
    )
    expected = (
        result.masks_long.set_index("object_id")["reference_positive"]
        .reindex(reversed_ids)
        .to_numpy(bool)
    )
    assert np.array_equal(arrays["EpCAM"], expected)
