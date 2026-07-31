"""Focused tests for marker-independent, pre-REDSEA geometry eligibility."""

from __future__ import annotations

from dataclasses import FrozenInstanceError

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from phenocycler.artifacts import GeometryContractError
from phenocycler.geometry_qc import (
    DuplicateGeometryError,
    GeometryQCConfig,
    evaluate_geometry_metrics,
)


def _metrics(**overrides) -> pd.DataFrame:
    n = 6
    values = {
        "object_id": [f"c{i}" for i in range(n)],
        "cell_geometry_present": [True] * n,
        "nucleus_geometry_present": [True] * n,
        "centroid_x_um": [10.0 * i for i in range(n)],
        "centroid_y_um": [0.0] * n,
        "cell_area_um2": [60.0] * n,
        "cell_area_px": [60.0] * n,
        "nucleus_area_px": [30.0] * n,
        "cell_solidity": [0.92] * n,
    }
    values.update(overrides)
    return pd.DataFrame(values)


def _config(**overrides) -> GeometryQCConfig:
    return GeometryQCConfig(
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        **overrides,
    )


def test_clean_geometry_has_three_explicit_eligibility_flags():
    result = evaluate_geometry_metrics(_metrics(), donor_id="1", config=_config())
    output = result.to_frame()

    assert output.analysis_eligible.all()
    assert output.estimation_eligible.all()
    assert output.spillover_context_eligible.all()
    assert (output.analysis_reason == "").all()
    assert result.donor_median_cell_area_um2 == 60.0
    assert result.donor_max_cell_area_um2 == 240.0
    result.verify_identity()


def test_result_records_are_immutable_and_frame_is_only_a_copy():
    result = evaluate_geometry_metrics(_metrics(), donor_id="1", config=_config())
    with pytest.raises(FrozenInstanceError):
        result.records[0].analysis_eligible = False

    frame = result.to_frame()
    frame.loc[0, "analysis_eligible"] = False
    assert result.records[0].analysis_eligible is True


def test_missing_required_columns_and_geometry_are_fatal_by_default():
    with pytest.raises(GeometryContractError, match="missing columns"):
        evaluate_geometry_metrics(
            _metrics().drop(columns=["cell_area_px"]),
            donor_id="1",
            config=_config(),
        )

    missing_nucleus = _metrics(
        nucleus_geometry_present=[False, True, True, True, True, True]
    )
    with pytest.raises(GeometryContractError, match="lack nucleus geometry"):
        evaluate_geometry_metrics(missing_nucleus, donor_id="1", config=_config())

    missing_cell = _metrics(
        cell_geometry_present=[False, True, True, True, True, True]
    )
    with pytest.raises(GeometryContractError, match="lack cell geometry"):
        evaluate_geometry_metrics(missing_cell, donor_id="1", config=_config())


def test_explicit_exclude_policy_rejects_missing_nucleus_from_every_qc_role():
    metrics = _metrics(
        nucleus_geometry_present=[False, True, True, True, True, True],
        nucleus_area_px=[float("nan"), 30, 30, 30, 30, 30],
    )
    result = evaluate_geometry_metrics(
        metrics,
        donor_id="1",
        config=_config(missing_geometry_policy="exclude"),
    ).to_frame()

    first = result.loc[result.object_id == "c0"].iloc[0]
    assert not first.analysis_eligible
    assert not first.estimation_eligible
    assert not first.spillover_context_eligible
    assert first.analysis_reason == "missing_nucleus_geometry"
    assert first.spillover_context_reason == "missing_nucleus_geometry"


def test_empty_nucleus_raster_is_rejected_from_every_qc_role():
    metrics = _metrics(
        nucleus_geometry_present=[True, True, True, True, True, True],
        nucleus_area_px=[0, 30, 30, 30, 30, 30],
    )
    first = evaluate_geometry_metrics(
        metrics,
        donor_id="1",
        config=_config(missing_geometry_policy="exclude"),
    ).to_frame().loc[0]

    assert not first.analysis_eligible
    assert not first.estimation_eligible
    assert not first.spillover_context_eligible
    assert first.analysis_reason == "nucleus_raster_empty"
    assert first.spillover_context_reason == "nucleus_raster_empty"


def test_duplicate_policy_is_fatal_by_default():
    metrics = _metrics(
        centroid_x_um=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
        cell_area_um2=[60.0, 59.8, 60.0, 60.0, 60.0, 60.0],
    )
    with pytest.raises(DuplicateGeometryError, match="duplicate geometry groups"):
        evaluate_geometry_metrics(metrics, donor_id="1", config=_config())


def test_explicit_duplicate_recovery_keeps_the_raster_winner_deterministically():
    metrics = _metrics(
        centroid_x_um=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
        cell_area_um2=[60.0, 59.8, 60.0, 60.0, 60.0, 60.0],
        cell_area_px=[50.0, 59.8, 60.0, 60.0, 60.0, 60.0],
    )
    config = _config(duplicate_policy="deterministic_keep")
    first = evaluate_geometry_metrics(metrics, donor_id="1", config=config)
    second = evaluate_geometry_metrics(
        metrics.sample(frac=1.0, random_state=7),
        donor_id="1",
        config=config,
    )
    output = first.to_frame().set_index("object_id")

    assert output.loc["c1", "duplicate_survivor"]
    assert output.loc["c1", "analysis_eligible"]
    assert not output.loc["c0", "analysis_eligible"]
    assert not output.loc["c0", "spillover_context_eligible"]
    assert output.loc["c0", "analysis_reason"] == "duplicate_detection"
    assert first.content_id == second.content_id


def test_duplicate_tie_break_is_lexicographic():
    metrics = _metrics(
        object_id=["z", "a", "c2", "c3", "c4", "c5"],
        centroid_x_um=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
    )
    output = evaluate_geometry_metrics(
        metrics,
        donor_id="1",
        config=_config(duplicate_policy="deterministic_keep"),
    ).to_frame().set_index("object_id")

    assert output.loc["a", "duplicate_survivor"]
    assert not output.loc["z", "duplicate_survivor"]


def test_duplicate_recovery_prefers_the_nucleus_bearing_copy():
    metrics = _metrics(
        centroid_x_um=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
        cell_area_um2=[60.0, 59.8, 60.0, 60.0, 60.0, 60.0],
        # The nucleus-less duplicate owns more raster pixels, so viability
        # must outrank raster ownership when choosing the survivor.
        cell_area_px=[70.0, 59.8, 60.0, 60.0, 60.0, 60.0],
        nucleus_geometry_present=[False, True, True, True, True, True],
        nucleus_area_px=[0.0, 30.0, 30.0, 30.0, 30.0, 30.0],
    )
    output = evaluate_geometry_metrics(
        metrics,
        donor_id="1",
        config=_config(
            missing_geometry_policy="exclude",
            duplicate_policy="deterministic_keep",
        ),
    ).to_frame().set_index("object_id")

    assert output.loc["c1", "duplicate_survivor"]
    assert output.loc["c1", "analysis_eligible"]
    assert not output.loc["c0", "duplicate_survivor"]
    assert not output.loc["c0", "analysis_eligible"]
    assert not output.loc["c0", "spillover_context_eligible"]
    assert "missing_nucleus_geometry" in output.loc["c0", "analysis_reasons"]
    assert "duplicate_detection" in output.loc["c0", "analysis_reasons"]


def test_duplicate_columns_have_one_parquet_schema_with_or_without_groups(
    tmp_path,
):
    clean = evaluate_geometry_metrics(
        _metrics(),
        donor_id="1",
        config=_config(duplicate_policy="deterministic_keep"),
    ).to_frame()
    duplicated = evaluate_geometry_metrics(
        _metrics(
            centroid_x_um=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
            cell_area_um2=[60.0, 59.8, 60.0, 60.0, 60.0, 60.0],
        ),
        donor_id="2",
        config=_config(duplicate_policy="deterministic_keep"),
    ).to_frame()
    clean_path = tmp_path / "clean.parquet"
    duplicate_path = tmp_path / "duplicate.parquet"
    clean.to_parquet(clean_path, index=False)
    duplicated.to_parquet(duplicate_path, index=False)

    clean_schema = pq.ParquetFile(clean_path).schema_arrow
    duplicate_schema = pq.ParquetFile(duplicate_path).schema_arrow
    assert clean_schema.equals(duplicate_schema, check_metadata=True)
    duplicate_group_type = clean_schema.field("duplicate_group").type
    assert pa.types.is_string(duplicate_group_type) or pa.types.is_large_string(
        duplicate_group_type
    )
    assert str(clean_schema.field("duplicate_survivor").type) == "bool"


def test_analysis_estimation_and_spillover_context_are_not_conflated():
    metrics = _metrics(
        cell_area_um2=[10.0, 60.0, 60.0, 60.0, 60.0, 60.0],
        cell_area_px=[10.0, 60.0, 60.0, 60.0, 0.0, 60.0],
        nucleus_area_px=[5.0, 60.0, 30.0, 30.0, 30.0, 30.0],
        cell_solidity=[0.92, 0.92, 0.50, 0.92, 0.92, 0.92],
    )
    output = evaluate_geometry_metrics(
        metrics, donor_id="1", config=_config()
    ).to_frame().set_index("object_id")

    # Tiny and low-solidity polygons are not analysis cells, but remain physical
    # boundary context by the explicit default policy.
    assert not output.loc["c0", "analysis_eligible"]
    assert output.loc["c0", "spillover_context_eligible"]
    assert not output.loc["c2", "analysis_eligible"]
    assert output.loc["c2", "spillover_context_eligible"]

    # A nucleus-fills-cell object is still called, but cannot estimate a
    # compartment model.
    assert output.loc["c1", "analysis_eligible"]
    assert not output.loc["c1", "estimation_eligible"]
    assert output.loc["c1", "estimation_reason"] == "no_cytoplasm"

    # An empty raster has no physical footprint and is eligible for nothing.
    assert not output.loc["c4", "analysis_eligible"]
    assert not output.loc["c4", "spillover_context_eligible"]
    assert output.loc["c4", "spillover_context_reason"] == "raster_empty"


def test_donor_relative_area_ceiling_and_raster_scale_are_geometry_only():
    metrics = _metrics(
        cell_area_um2=[60.0, 60.0, 60.0, 60.0, 60.0, 300.0],
        cell_area_px=[60.0, 60.0, 60.0, 60.0, 10.0, 300.0],
    )
    output = evaluate_geometry_metrics(
        metrics, donor_id="1", config=_config()
    ).to_frame().set_index("object_id")

    assert output.loc["c4", "analysis_reason"] == "raster_area_mismatch"
    assert "cell_area_above_max" in output.loc["c5", "analysis_reasons"]
