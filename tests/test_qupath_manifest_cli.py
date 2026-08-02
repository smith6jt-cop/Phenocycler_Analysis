from __future__ import annotations

import pytest

from phenocycler.artifacts import ChannelMapping
from phenocycler import qupath_manifest
from phenocycler.qupath_manifest import example_specification


def test_manifest_template_names_every_required_contract_field():
    spec = example_specification()
    assert spec["spec_version"] == 1
    row = spec["images"][0]
    assert {
        "donor_id",
        "image_id",
        "qptiff",
        "cell_geojson",
        "nucleus_geojson",
        "pixel_size_um_x",
        "pixel_size_um_y",
        "panel_id",
        "channel_map",
    } <= set(row)
    assert "measurement_csv" in spec["defaults"]
    assert len(spec["defaults"]["measurement_csv"]) == 2
    assert row["measurement_csv_index"] == 0
    assert "segmentation_version" in spec["defaults"]


def test_measurement_csv_list_selects_the_declared_source(tmp_path):
    row = {
        "measurement_csv": ["batch1.csv", "batch2.csv"],
        "measurement_csv_index": 1,
    }
    assert qupath_manifest._measurement_path(tmp_path, row) == (
        tmp_path / "batch2.csv"
    ).resolve()


def test_measurement_csv_list_requires_a_valid_index(tmp_path):
    with pytest.raises(ValueError, match="require an integer"):
        qupath_manifest._measurement_path(
            tmp_path, {"measurement_csv": ["batch1.csv", "batch2.csv"]}
        )
    with pytest.raises(ValueError, match="outside 0..1"):
        qupath_manifest._measurement_path(
            tmp_path,
            {
                "measurement_csv": ["batch1.csv", "batch2.csv"],
                "measurement_csv_index": 2,
            },
        )
    with pytest.raises(ValueError, match="only valid"):
        qupath_manifest._measurement_path(
            tmp_path,
            {
                "measurement_csv": "batch1.csv",
                "measurement_csv_index": 0,
            },
        )


def test_explicit_channel_map_must_match_qptiff_names_and_indices(
    tmp_path, monkeypatch
):
    qptiff = tmp_path / "image.qptiff"
    qptiff.write_bytes(b"test")
    observed = (
        ChannelMapping("DAPI", "DAPI", 0),
        ChannelMapping("CD3e", "CD3e", 1),
    )
    monkeypatch.setattr(
        qupath_manifest, "_channels_from_qptiff", lambda _path: observed
    )
    row = {
        "channel_map": [
            {
                "marker": "DAPI",
                "channel_name": "DAPI",
                "channel_index": 0,
            },
            {
                "marker": "CD3e",
                "channel_name": "wrong",
                "channel_index": 1,
            },
        ]
    }
    with pytest.raises(ValueError, match="not declared"):
        qupath_manifest._channel_mapping(row, qptiff)
