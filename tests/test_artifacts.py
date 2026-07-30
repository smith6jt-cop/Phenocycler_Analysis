"""Focused tests for content-derived input and stage/run contracts."""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import pytest

import phenocycler.artifacts as artifacts
from phenocycler.artifacts import (
    ChannelMapping,
    CodeMetadata,
    ContractError,
    DuplicateObjectError,
    GeometryContractError,
    InputArtifact,
    PartialArtifactError,
    QuPathCohortManifest,
    QuPathImageManifest,
    RunManifest,
    SchemaContractError,
    StageManifest,
    StaleArtifactError,
    hash_identifiers,
    image_id_to_donor_id,
    snapshot_parquet_dataset,
)


def _geometry(x: float = 0.0):
    return {
        "type": "Polygon",
        "coordinates": [
            [
                [x, 0.0],
                [x + 1.0, 0.0],
                [x + 1.0, 1.0],
                [x, 0.0],
            ]
        ],
    }


def _write_geojson(
    path: Path,
    object_ids=("a", "b"),
    *,
    include_nucleus: bool = True,
    missing_cell_index: int | None = None,
):
    features = []
    for index, object_id in enumerate(object_ids):
        feature = {
            "type": "Feature",
            "id": object_id,
            "geometry": None if index == missing_cell_index else _geometry(float(index)),
            "properties": {"objectType": "cell"},
        }
        if include_nucleus:
            feature["nucleusGeometry"] = _geometry(float(index) + 0.2)
        features.append(feature)
    path.write_text(json.dumps({"type": "FeatureCollection", "features": features}))


def _qupath_inputs(tmp_path: Path):
    csv_path = tmp_path / "cells.csv"
    csv_path.write_text(
        "Image,Object ID,Cell: DAPI: Mean,Cell: E-cadherin: Mean\n"
        "image-1,a,10,20\nimage-1,b,11,21\n"
    )
    qptiff = tmp_path / "image.qptiff"
    qptiff.write_bytes(b"small synthetic qptiff")
    geojson = tmp_path / "cells.geojson"
    _write_geojson(geojson)
    return csv_path, qptiff, geojson


def _qupath_manifest(tmp_path: Path) -> QuPathImageManifest:
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    return QuPathImageManifest.create(
        donor_id="6374",
        image_id="6374_Scan1",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=0.5,
        pixel_size_um_y=0.5,
        panel_id="panel-v1",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("E_cadherin", "E-cadherin", 1),
        ),
        segmentation_version="instanseg-0.1.1+script-v2",
    )


def _cohort_image(
    tmp_path: Path,
    *,
    donor_id: str,
    image_id: str,
    panel_id: str = "panel-v1",
    channel_map=(
        ChannelMapping("DAPI", "DAPI", 0),
        ChannelMapping("E_cadherin", "E-cadherin", 1),
    ),
) -> QuPathImageManifest:
    image_dir = tmp_path / f"{donor_id}_{image_id.replace('/', '_')}"
    image_dir.mkdir(parents=True, exist_ok=True)
    csv_path, qptiff, geojson = _qupath_inputs(image_dir)
    return QuPathImageManifest.create(
        donor_id=donor_id,
        image_id=image_id,
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=0.5,
        pixel_size_um_y=0.5,
        panel_id=panel_id,
        channel_map=channel_map,
        segmentation_version="instanseg-0.1.1+script-v2",
    )


def test_qupath_manifest_roundtrip_and_current_validation(tmp_path):
    manifest = _qupath_manifest(tmp_path)
    path = manifest.write_json(tmp_path / "qupath_manifest.json")
    loaded = QuPathImageManifest.read_json(path)

    assert loaded == manifest
    assert loaded.cell_geometry.feature_count == 2
    assert loaded.cell_geometry.nucleus_geometry_count == 2
    assert loaded.cell_geometry.object_id_sha256 == hash_identifiers(["a", "b"])
    loaded.validate_current()


def test_qupath_manifest_is_stale_when_any_source_byte_changes(tmp_path):
    manifest = _qupath_manifest(tmp_path)
    Path(manifest.qptiff.path).write_bytes(b"changed qptiff")

    with pytest.raises(StaleArtifactError, match="file (size|content) changed"):
        manifest.validate_current()


def test_fast_qupath_validation_does_not_rehash_unchanged_qptiff(tmp_path, monkeypatch):
    manifest = _qupath_manifest(tmp_path)
    qptiff_path = manifest.qptiff.path
    original = artifacts.sha256_file

    def refuse_qptiff(path):
        if str(path) == qptiff_path:
            raise AssertionError("routine status reread an unchanged qptiff")
        return original(path)

    monkeypatch.setattr(artifacts, "sha256_file", refuse_qptiff)
    manifest.validate_current(mode="fast")
    with pytest.raises(AssertionError, match="reread"):
        manifest.validate_current(mode="content")


def test_qptiff_can_use_an_explicitly_labelled_sampled_fingerprint(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    manifest = QuPathImageManifest.create(
        donor_id="1",
        image_id="image",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        panel_id="panel",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
        segmentation_version="seg-v1",
        qptiff_hash_mode="sampled",
        qptiff_sample_size_bytes=4,
        qptiff_sample_count=3,
    )
    assert manifest.qptiff.hash_mode == "sampled"
    assert manifest.qptiff.sample_size_bytes == 4
    manifest.validate_current(mode="content")


def test_qupath_manifest_fails_on_missing_or_duplicate_geometry(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    _write_geojson(geojson, include_nucleus=False)
    kwargs = dict(
        donor_id="1",
        image_id="image",
        measurement_csv=csv_path,
        qptiff=qptiff,
        cell_geojson=geojson,
        nucleus_geojson=geojson,
        pixel_size_um_x=1.0,
        pixel_size_um_y=1.0,
        panel_id="panel",
        channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
        segmentation_version="seg-v1",
    )
    with pytest.raises(GeometryContractError, match="nucleusGeometry"):
        QuPathImageManifest.create(**kwargs)

    _write_geojson(geojson, object_ids=("same", "same"))
    with pytest.raises(DuplicateObjectError, match="duplicate object_id"):
        QuPathImageManifest.create(**kwargs)

    _write_geojson(geojson, missing_cell_index=1)
    with pytest.raises(GeometryContractError, match="no usable cell geometry"):
        QuPathImageManifest.create(**kwargs)


def _write_partition(root: Path, donor: str, object_ids, values=None):
    partition = root / f"donor_id={donor}"
    partition.mkdir(parents=True, exist_ok=True)
    values = values if values is not None else list(range(len(object_ids)))
    pd.DataFrame(
        {
            "object_id": list(object_ids),
            "donor_id": [str(donor)] * len(object_ids),
            "value": values,
        }
    ).to_parquet(partition / "data_0.parquet", index=False)


def _stage_fixture(tmp_path: Path):
    output = tmp_path / "output"
    _write_partition(output, "1", ["a", "b"], [1.0, 2.0])
    _write_partition(output, "2", ["c", "d"], [3.0, 4.0])
    source = tmp_path / "stage_code.py"
    source.write_text("METHOD = 'v1'\n")
    input_path = tmp_path / "input.txt"
    input_path.write_text("source data\n")
    code = CodeMetadata.capture(tmp_path, source_paths=(source,))
    inputs = (InputArtifact.from_path("raw", input_path),)
    config = {"alpha": 1.0, "emit_compartments": True}
    return output, source, input_path, code, inputs, config


def test_dataset_snapshot_requires_exact_donor_set_schema_and_object_ids(tmp_path):
    output = tmp_path / "output"
    _write_partition(output, "1", ["a", "b"])
    _write_partition(output, "2", ["c", "d"])

    snapshot = snapshot_parquet_dataset(output, expected_donors=("1", "2"))
    assert snapshot.donor_ids == ("1", "2")
    assert snapshot.total_rows == 4

    with pytest.raises(PartialArtifactError, match="extra"):
        snapshot_parquet_dataset(output, expected_donors=("1",))
    with pytest.raises(PartialArtifactError, match="missing"):
        snapshot_parquet_dataset(output, expected_donors=("1", "2", "3"))


def test_dataset_snapshot_rejects_schema_drift_and_duplicate_objects(tmp_path):
    schema_root = tmp_path / "schema"
    _write_partition(schema_root, "1", ["a"], [1])
    _write_partition(schema_root, "2", ["b"], ["one"])
    with pytest.raises(SchemaContractError, match="schema differs"):
        snapshot_parquet_dataset(schema_root, expected_donors=("1", "2"))

    duplicate_root = tmp_path / "duplicates"
    _write_partition(duplicate_root, "1", ["a", "a"])
    with pytest.raises(DuplicateObjectError, match="duplicate object_id"):
        snapshot_parquet_dataset(duplicate_root, expected_donors=("1",))


def test_stage_manifest_is_content_derived_and_validates_current_state(tmp_path):
    output, _, _, code, inputs, config = _stage_fixture(tmp_path)
    first = StageManifest.create(
        stage="redsea",
        method_version="redsea-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )
    second = StageManifest.create(
        stage="redsea",
        method_version="redsea-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )
    # Creation timestamps differ, but the run identity is entirely content-derived.
    assert first.content_id == second.content_id
    output_file = output / "donor_id=1" / "data_0.parquet"
    stat = output_file.stat()
    os.utime(
        output_file,
        ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000),
    )
    touched = StageManifest.create(
        stage="redsea",
        method_version="redsea-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )
    assert touched.content_id == first.content_id

    # Use the freshly captured stat cache after deliberately touching the file.
    path = touched.write_json(tmp_path / "redsea.stage.json")
    loaded = StageManifest.read_json(path)
    loaded.validate_current(config=config)
    assert loaded.output.total_rows == 4
    assert loaded.output.schema_sha256
    assert loaded.output.object_id_sha256


def test_stage_manifest_detects_config_input_code_output_and_partial_drift(tmp_path):
    output, source, input_path, code, inputs, config = _stage_fixture(tmp_path)
    manifest = StageManifest.create(
        stage="geometry_qc",
        method_version="geometry-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )

    with pytest.raises(StaleArtifactError, match="configuration changed"):
        manifest.validate_current(config={**config, "alpha": 0.5})

    input_path.write_text("changed source data\n")
    with pytest.raises(StaleArtifactError, match="input 'raw' changed"):
        manifest.validate_current(config=config)
    input_path.write_text("source data\n")

    source.write_text("METHOD = 'v2'\n")
    with pytest.raises(StaleArtifactError, match="producing code changed"):
        manifest.validate_current(config=config)
    source.write_text("METHOD = 'v1'\n")

    _write_partition(output, "1", ["a", "b"], [100.0, 2.0])
    with pytest.raises(StaleArtifactError, match="content changed"):
        manifest.validate_current(config=config)

    # Restore content, then remove a complete donor partition.
    _write_partition(output, "1", ["a", "b"], [1.0, 2.0])
    (output / "donor_id=2" / "data_0.parquet").unlink()
    (output / "donor_id=2").rmdir()
    with pytest.raises(PartialArtifactError, match="missing"):
        manifest.validate_current(config=config)


def test_stage_manifest_fingerprints_immutable_sidecars(tmp_path):
    output, _, _, code, inputs, config = _stage_fixture(tmp_path)
    sidecar = tmp_path / "model_audit.json"
    sidecar.write_text('{"status":"valid"}\n')
    manifest = StageManifest.create(
        stage="calibrate",
        method_version="calibration-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        sidecar_paths=(sidecar,),
        config=config,
        code=code,
    )
    manifest.validate_current(config=config)
    sidecar.write_text('{"status":"changed"}\n')
    with pytest.raises(StaleArtifactError, match="sidecar changed"):
        manifest.validate_current(config=config)


def test_stage_manifest_refuses_partial_creation(tmp_path):
    output = tmp_path / "output"
    _write_partition(output, "1", ["a"])
    source = tmp_path / "source.py"
    source.write_text("pass\n")
    input_path = tmp_path / "input"
    input_path.write_text("x")
    with pytest.raises(PartialArtifactError, match="missing"):
        StageManifest.create(
            stage="cells",
            method_version="v1",
            expected_donors=("1", "2"),
            inputs=(InputArtifact.from_path("csv", input_path),),
            output_root=output,
            config={},
            code=CodeMetadata.capture(tmp_path, source_paths=(source,)),
        )


def test_run_manifest_references_exact_stage_content(tmp_path):
    output, _, _, code, inputs, config = _stage_fixture(tmp_path)
    stage = StageManifest.create(
        stage="cells",
        method_version="cells-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "cells.stage.json")
    run = RunManifest.create(
        run_name="cohort-run",
        expected_donors=("1", "2"),
        stage_manifest_paths=(stage_path,),
        config=config,
        code=code,
    )
    run_path = run.write_json(tmp_path / "run.json")
    loaded = RunManifest.read_json(run_path)
    loaded.validate_current()
    assert loaded.stages[0].manifest_content_id == stage.content_id


def test_run_manifest_refuses_a_stale_stage_output(tmp_path):
    output, _, _, code, inputs, config = _stage_fixture(tmp_path)
    stage = StageManifest.create(
        stage="cells",
        method_version="cells-v1",
        expected_donors=("1", "2"),
        inputs=inputs,
        output_root=output,
        config=config,
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "cells.stage.json")
    _write_partition(output, "1", ["a", "b"], [100.0, 2.0])

    with pytest.raises(StaleArtifactError, match="content changed"):
        RunManifest.create(
            run_name="cohort-run",
            expected_donors=("1", "2"),
            stage_manifest_paths=(stage_path,),
            config=config,
            code=code,
        )


def test_channel_map_and_measurement_contracts_are_explicit(tmp_path):
    csv_path, qptiff, geojson = _qupath_inputs(tmp_path)
    with pytest.raises(ContractError, match="duplicate markers"):
        QuPathImageManifest.create(
            donor_id="1",
            image_id="image",
            measurement_csv=csv_path,
            qptiff=qptiff,
            cell_geojson=geojson,
            nucleus_geojson=geojson,
            pixel_size_um_x=1.0,
            pixel_size_um_y=1.0,
            panel_id="panel",
            channel_map=(
                ChannelMapping("DAPI", "DAPI", 0),
                ChannelMapping("DAPI", "DAPI-2", 1),
            ),
            segmentation_version="seg-v1",
        )

    bad_csv = tmp_path / "bad.csv"
    bad_csv.write_text("Image,Not Object ID\nimage,a\n")
    with pytest.raises(ContractError, match="Object ID"):
        QuPathImageManifest.create(
            donor_id="1",
            image_id="image",
            measurement_csv=bad_csv,
            qptiff=qptiff,
            cell_geojson=geojson,
            nucleus_geojson=geojson,
            pixel_size_um_x=1.0,
            pixel_size_um_y=1.0,
            panel_id="panel",
            channel_map=(ChannelMapping("DAPI", "DAPI", 0),),
            segmentation_version="seg-v1",
        )


def test_qupath_cohort_manifest_roundtrip_and_pure_ingest_mapping(tmp_path):
    images = (
        _cohort_image(tmp_path, donor_id="1", image_id="donor-1-image"),
        _cohort_image(tmp_path, donor_id="2", image_id="donor-2-image"),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="synthetic-cohort",
        expected_donors=("1", "2"),
        images=images,
    )
    path = cohort.write_json(tmp_path / "qupath_cohort.json")
    loaded = QuPathCohortManifest.read_json(path)

    assert loaded == cohort
    assert loaded.donor_ids == ("1", "2")
    assert loaded.image_ids == ("donor-1-image", "donor-2-image")
    assert len(loaded.panels) == 1
    assert loaded.panels[0].panel_id == "panel-v1"
    expected = {"donor-1-image": "1", "donor-2-image": "2"}
    assert image_id_to_donor_id(loaded) == expected
    assert image_id_to_donor_id(loaded.images) == expected
    assert loaded.image_to_donor() == expected
    loaded.validate_current(mode="fast")

    reference = InputArtifact.from_qupath_cohort_manifest("qupath", path)
    reference.validate_current(mode="fast")


def test_qupath_cohort_requires_exactly_one_unique_image_per_expected_donor(tmp_path):
    first = _cohort_image(tmp_path, donor_id="1", image_id="image-1")
    second = _cohort_image(tmp_path, donor_id="2", image_id="image-2")
    with pytest.raises(PartialArtifactError, match="missing=.*3"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2", "3"),
            images=(first, second),
        )
    with pytest.raises(PartialArtifactError, match="extra=.*2"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1",),
            images=(first, second),
        )

    duplicate_donor = _cohort_image(
        tmp_path, donor_id="1", image_id="another-image"
    )
    with pytest.raises(ContractError, match="duplicate donors"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1",),
            images=(first, duplicate_donor),
        )

    duplicate_image = _cohort_image(
        tmp_path, donor_id="2", image_id="image-1"
    )
    with pytest.raises(ContractError, match="duplicate image_id"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2"),
            images=(first, duplicate_image),
        )


def test_qupath_cohort_enforces_one_exact_channel_map_per_panel(tmp_path):
    first = _cohort_image(tmp_path, donor_id="1", image_id="image-1")
    conflicting = _cohort_image(
        tmp_path,
        donor_id="2",
        image_id="image-2",
        panel_id="panel-v1",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("E_cadherin", "E-cadherin", 2),
        ),
    )
    with pytest.raises(ContractError, match="contradictory channel mappings"):
        QuPathCohortManifest.create(
            cohort_name="cohort",
            expected_donors=("1", "2"),
            images=(first, conflicting),
        )

    # A physical panel change is valid only when it has a distinct panel id.
    second_panel = _cohort_image(
        tmp_path,
        donor_id="2",
        image_id="image-2-panel2",
        panel_id="panel-v2",
        channel_map=(
            ChannelMapping("DAPI", "DAPI", 0),
            ChannelMapping("CD38", "CD38", 1),
        ),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="cohort",
        expected_donors=("1", "2"),
        images=(first, second_panel),
    )
    assert [panel.panel_id for panel in cohort.panels] == ["panel-v1", "panel-v2"]


def test_qupath_cohort_fast_validation_never_rereads_unchanged_qptiffs(
    tmp_path, monkeypatch
):
    images = (
        _cohort_image(tmp_path, donor_id="1", image_id="image-1"),
        _cohort_image(tmp_path, donor_id="2", image_id="image-2"),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="cohort",
        expected_donors=("1", "2"),
        images=images,
    )
    qptiffs = {image.qptiff.path for image in images}
    original = artifacts.sha256_file

    def refuse_qptiffs(path):
        if str(path) in qptiffs:
            raise AssertionError("cohort status reread an unchanged qptiff")
        return original(path)

    monkeypatch.setattr(artifacts, "sha256_file", refuse_qptiffs)
    cohort.validate_current(mode="fast")
    with pytest.raises(AssertionError, match="reread"):
        cohort.validate_current(mode="content")


def test_qupath_cohort_reports_a_changed_embedded_source(tmp_path):
    images = (
        _cohort_image(tmp_path, donor_id="1", image_id="image-1"),
        _cohort_image(tmp_path, donor_id="2", image_id="image-2"),
    )
    cohort = QuPathCohortManifest.create(
        cohort_name="cohort",
        expected_donors=("1", "2"),
        images=images,
    )
    Path(images[1].qptiff.path).write_bytes(b"changed qptiff")
    with pytest.raises(StaleArtifactError, match="file (size|content) changed"):
        cohort.validate_current(mode="fast")
