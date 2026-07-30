from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
import pytest

from phenocycler.artifacts import (
    PartialArtifactError,
    StaleArtifactError,
    hash_identifiers,
)
from phenocycler.geometry_qc import GeometryQCConfig, evaluate_geometry_metrics
from phenocycler.resume import (
    _hardlink_once,
    _validate_geometry_partition,
    _validate_geometry_policy_transition,
)


def _geometry_values(policy: str) -> dict:
    return {
        "min_cell_area_um2": 20.0,
        "max_cell_area_median_multiple": 4.0,
        "min_cell_solidity": 0.7,
        "min_raster_area_ratio": 0.5,
        "max_raster_area_ratio": 2.0,
        "no_cytoplasm_nc_ratio": 0.999,
        "duplicate_radius_um": 1.0,
        "duplicate_area_tolerance": 0.15,
        "missing_geometry_policy": "exclude",
        "duplicate_policy": policy,
    }


def test_hardlink_reuse_is_idempotent_and_rejects_foreign_content(tmp_path):
    source = tmp_path / "source.parquet"
    source.write_bytes(b"same bytes")
    target = tmp_path / "target" / "data.parquet"

    _hardlink_once(source, target)
    _hardlink_once(source, target)
    assert source.stat().st_ino == target.stat().st_ino

    foreign = tmp_path / "foreign.parquet"
    foreign.write_bytes(b"different bytes")
    with pytest.raises(PartialArtifactError, match="not the source hard link"):
        _hardlink_once(foreign, target)


def test_geometry_prefix_reuse_allows_only_the_policy_transition():
    source = SimpleNamespace(
        stage="ingest",
        config=SimpleNamespace(
            canonical=json.dumps(
                {"configuration": {"geometry_qc": _geometry_values("fatal")}}
            )
        ),
    )
    context = SimpleNamespace(
        identity_config={
            "geometry_qc": _geometry_values("deterministic_keep")
        }
    )

    _validate_geometry_policy_transition(context, source)
    context.identity_config["geometry_qc"]["min_cell_solidity"] = 0.8
    with pytest.raises(
        StaleArtifactError, match="changed beyond duplicate_policy"
    ):
        _validate_geometry_policy_transition(context, source)


def test_geometry_partition_validation_checks_ingest_universe(tmp_path):
    metrics = pd.DataFrame(
        {
            "object_id": [f"c{i}" for i in range(4)],
            "cell_geometry_present": [True] * 4,
            "nucleus_geometry_present": [True] * 4,
            "centroid_x_um": [0.0, 10.0, 20.0, 30.0],
            "centroid_y_um": [0.0] * 4,
            "cell_area_um2": [60.0] * 4,
            "cell_area_px": [60.0] * 4,
            "nucleus_area_px": [30.0] * 4,
            "cell_solidity": [0.92] * 4,
        }
    )
    output = evaluate_geometry_metrics(
        metrics,
        donor_id="1",
        config=GeometryQCConfig(
            pixel_size_um_x=1.0,
            pixel_size_um_y=1.0,
            duplicate_policy="deterministic_keep",
        ),
    ).to_frame()
    path = tmp_path / "data.parquet"
    output.to_parquet(path, index=False)
    expected = SimpleNamespace(
        donor_id="1",
        row_count=len(output),
        object_id_sha256=hash_identifiers(output["object_id"]),
    )

    check, _ = _validate_geometry_partition(
        path,
        expected=expected,
        require_no_duplicates=True,
    )
    assert check.row_count == 4
    assert check.duplicate_group_count == 0

    changed = SimpleNamespace(
        donor_id="1",
        row_count=len(output),
        object_id_sha256=hash_identifiers(["other", "ids"]),
    )
    with pytest.raises(StaleArtifactError, match="object universe differs"):
        _validate_geometry_partition(
            path,
            expected=changed,
            require_no_duplicates=True,
        )
