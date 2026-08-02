"""Donor-scoped receipts and a resource-aware cross-stage task queue."""

from __future__ import annotations

import hashlib
import json
import multiprocessing
import os
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq

from .artifacts import (
    CodeMetadata,
    ConfigSnapshot,
    ContractError,
    DonorDatasetSnapshot,
    FileFingerprint,
    InputArtifact,
    PartialArtifactError,
    SchemaContractError,
    SchemaField,
    StageManifest,
    StaleArtifactError,
    hash_identifiers,
    sha256_json,
)
from .calibration_stage import calibrate_expression
from .expression import (
    build_selected_expression,
    read_single_partition,
    write_donor_partition,
)
from .geometry_qc import GeometryQCConfig
from .geometry_stage import run_geometry_qc_for_image
from .parallel import resolve_jobs, terminate_process_pool_workers
from .redsea import RedseaParams, process_donor
from .reference_controls import (
    build_reference_controls,
    reference_masks_from_wide,
    reference_masks_to_wide,
)
from .state_markers import annotate_proliferation
from .typing_stage import type_calibrated_cells

DONOR_RECEIPT_VERSION = 1
DONOR_STAGES = (
    "geometry",
    "redsea",
    "expression",
    "controls",
    "calibrate",
    "type",
    "states",
)
DONOR_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "geometry": (),
    "redsea": ("geometry",),
    "expression": ("geometry", "redsea"),
    "controls": ("expression",),
    "calibrate": ("expression", "controls"),
    "type": ("calibrate",),
    "states": ("expression",),
}
DONOR_AUDITS = {
    "expression": "expression_availability",
    "controls": "reference_control_models",
    "calibrate": "marker_calibration_models",
    "states": "state_models",
}


def _receipt_identity(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {
            str(key): _receipt_identity(item)
            for key, item in value.items()
            if key not in {"created_at_utc", "mtime_ns"}
        }
    if isinstance(value, (tuple, list)):
        return [_receipt_identity(item) for item in value]
    return value


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    if path.exists():
        raise FileExistsError(f"immutable donor receipt already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def _atomic_write_frame(path: Path, frame: pd.DataFrame) -> Path:
    if path.exists():
        current = pd.read_parquet(path)
        try:
            pd.testing.assert_frame_equal(
                current.reset_index(drop=True),
                frame.reset_index(drop=True),
                check_dtype=False,
                check_categorical=False,
            )
        except AssertionError as exc:
            raise StaleArtifactError(
                f"existing aggregate audit differs: {path}"
            ) from exc
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        frame.to_parquet(temporary, index=False)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def snapshot_donor_partition(
    root: Path | str,
    donor_id: str,
    *,
    object_id_column: str = "object_id",
) -> DonorDatasetSnapshot:
    """Snapshot exactly one donor partition without requiring cohort completion."""

    root = Path(root).expanduser().resolve()
    donor = str(donor_id)
    partition = root / f"donor_id={donor}"
    if not partition.is_dir():
        raise PartialArtifactError(
            f"donor {donor} output partition is missing: {partition}"
        )
    files = sorted(partition.glob("*.parquet"))
    if not files:
        raise PartialArtifactError(f"donor {donor} has no Parquet output")

    ids: list[Any] = []
    row_count = 0
    donor_schema = None
    donor_schema_sha256 = None
    fingerprints = []
    content_entries = []
    for path in files:
        parquet = pq.ParquetFile(path)
        schema = parquet.schema_arrow
        schema_sha256 = hashlib.sha256(
            schema.serialize().to_pybytes()
        ).hexdigest()
        if object_id_column not in schema.names:
            raise SchemaContractError(
                f"{path}: schema has no {object_id_column!r} column"
            )
        if donor_schema is None:
            donor_schema = schema
            donor_schema_sha256 = schema_sha256
        elif not donor_schema.equals(schema, check_metadata=True):
            raise SchemaContractError(
                f"donor {donor} contains inconsistent Parquet schemas"
            )
        columns = [object_id_column]
        if "donor_id" in schema.names:
            columns.append("donor_id")
        table = parquet.read(columns=columns)
        row_count += table.num_rows
        ids.extend(table.column(object_id_column).to_pylist())
        if "donor_id" in columns:
            stored = {
                str(value)
                for value in table.column("donor_id").to_pylist()
                if value is not None
            }
            if stored != {donor}:
                raise ContractError(
                    f"{path}: stored donor_id values {sorted(stored)} "
                    f"do not equal partition {donor}"
                )
        fingerprint = FileFingerprint.capture(path)
        fingerprints.append(fingerprint)
        content_entries.append(
            (
                path.relative_to(partition).as_posix(),
                fingerprint.size_bytes,
                fingerprint.sha256,
            )
        )

    object_hash = hash_identifiers(ids)
    content_hash = sha256_json(
        {
            "donor_id": donor,
            "files": content_entries,
            "rows": row_count,
            "schema": donor_schema_sha256,
            "objects": object_hash,
        }
    )
    return DonorDatasetSnapshot(
        donor_id=donor,
        files=tuple(fingerprints),
        row_count=row_count,
        schema=tuple(
            SchemaField(
                name=field.name,
                type=str(field.type),
                nullable=field.nullable,
            )
            for field in donor_schema
        ),
        schema_sha256=str(donor_schema_sha256),
        object_id_sha256=object_hash,
        content_sha256=content_hash,
    )


def _validate_donor_snapshot(
    root: Path | str,
    snapshot: DonorDatasetSnapshot,
    *,
    mode: str,
) -> None:
    root = Path(root).expanduser().resolve()
    partition = root / f"donor_id={snapshot.donor_id}"
    if not partition.is_dir():
        raise PartialArtifactError(
            f"donor output partition is missing: {partition}"
        )
    actual = {str(path.resolve()) for path in partition.glob("*.parquet")}
    recorded = {item.path for item in snapshot.files}
    if actual != recorded:
        raise PartialArtifactError(
            f"donor {snapshot.donor_id} Parquet inventory changed"
        )
    for fingerprint in snapshot.files:
        fingerprint.validate_current(mode=mode)
    if mode == "content":
        current = snapshot_donor_partition(root, snapshot.donor_id)
        if current.content_sha256 != snapshot.content_sha256:
            raise StaleArtifactError(
                f"donor {snapshot.donor_id} output content changed"
            )


def donor_receipt_path(context, stage: str, donor_id: str) -> Path:
    return (
        context.config.manifests_dir
        / "donors"
        / f"donor_id={donor_id}"
        / f"{stage}.json"
    )


def donor_audit_path(context, stage: str, donor_id: str) -> Path | None:
    name = DONOR_AUDITS.get(stage)
    if name is None:
        return None
    return (
        context.config.audit_dir
        / "donors"
        / name
        / f"donor_id={donor_id}"
        / "data.parquet"
    )


def _donor_stage_config(context, stage: str, donor_id: str) -> dict[str, Any]:
    from .pipeline import _stage_config

    return {
        **_stage_config(context, stage),
        "donor_id": str(donor_id),
        "donor_receipt_version": DONOR_RECEIPT_VERSION,
    }


def _input_from_receipt(
    name: str,
    path: Path,
    *,
    donor_id: str,
) -> InputArtifact:
    receipt = DonorStageReceipt.read_json(path)
    if receipt.stage != str(name) or receipt.donor_id != str(donor_id):
        raise ContractError(
            f"donor receipt at {path} declares "
            f"{receipt.stage}/{receipt.donor_id}, expected {name}/{donor_id}"
        )
    return InputArtifact(
        name=str(name),
        kind="donor_receipt",
        path=str(path.resolve()),
        content_sha256=receipt.content_id,
        files=(),
    )


def _validate_receipt_input(
    artifact: InputArtifact,
    *,
    mode: str,
) -> None:
    if artifact.kind != "donor_receipt":
        artifact.validate_current(mode=mode)
        return
    receipt = DonorStageReceipt.read_json(artifact.path)
    if receipt.content_id != artifact.content_sha256:
        raise StaleArtifactError(
            f"input donor receipt changed: {artifact.name}"
        )
    receipt.validate_current(validation_mode=mode)


@dataclass(frozen=True)
class DonorStageReceipt:
    """Immutable proof that one donor completed one scientific stage."""

    receipt_version: int
    stage: str
    method_version: str
    donor_id: str
    inputs: tuple[InputArtifact, ...]
    output_root: str
    output: DonorDatasetSnapshot
    sidecars: tuple[FileFingerprint, ...]
    config: ConfigSnapshot
    code: CodeMetadata
    created_at_utc: str
    content_id: str

    @classmethod
    def create(
        cls,
        *,
        stage: str,
        method_version: str,
        donor_id: str,
        inputs: Sequence[InputArtifact],
        output_root: Path | str,
        sidecar_paths: Sequence[Path | str],
        config: Any,
        code: CodeMetadata,
    ) -> DonorStageReceipt:
        donor = str(donor_id)
        input_tuple = tuple(sorted(inputs, key=lambda item: item.name))
        if len({item.name for item in input_tuple}) != len(input_tuple):
            raise ContractError("donor receipt input names must be unique")
        for artifact in input_tuple:
            _validate_receipt_input(artifact, mode="fast")
        root = Path(output_root).expanduser().resolve()
        output = snapshot_donor_partition(root, donor)
        sidecars = tuple(
            sorted(
                (
                    FileFingerprint.capture(path)
                    for path in sidecar_paths
                ),
                key=lambda item: item.path,
            )
        )
        draft = cls(
            receipt_version=DONOR_RECEIPT_VERSION,
            stage=str(stage),
            method_version=str(method_version),
            donor_id=donor,
            inputs=input_tuple,
            output_root=str(root),
            output=output,
            sidecars=sidecars,
            config=ConfigSnapshot.capture(config),
            code=code,
            created_at_utc=datetime.now(timezone.utc).isoformat(),
            content_id="",
        )
        return replace(
            draft,
            content_id=sha256_json(
                _receipt_identity(draft.to_dict())
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    def verify_identity(self) -> None:
        if self.receipt_version != DONOR_RECEIPT_VERSION:
            raise ContractError(
                f"unsupported donor receipt version {self.receipt_version}"
            )
        if not self.stage or not self.method_version or not self.donor_id:
            raise ContractError("donor receipt identity fields are required")
        if self.stage not in DONOR_STAGES:
            raise ContractError(f"unknown donor receipt stage {self.stage!r}")
        if self.inputs != tuple(
            sorted(self.inputs, key=lambda item: item.name)
        ) or len({item.name for item in self.inputs}) != len(self.inputs):
            raise ContractError("donor receipt inputs are not canonical")
        if self.sidecars != tuple(
            sorted(self.sidecars, key=lambda item: item.path)
        ):
            raise ContractError("donor receipt sidecars are not canonical")
        if self.output.donor_id != self.donor_id:
            raise ContractError("donor receipt output belongs to another donor")
        payload = self.to_dict()
        payload["content_id"] = ""
        if sha256_json(_receipt_identity(payload)) != self.content_id:
            raise ContractError("donor receipt content_id is invalid")

    def validate_current(
        self,
        *,
        config: Any | None = None,
        validate_code: bool = True,
        validate_inputs: bool = True,
        validation_mode: str = "fast",
    ) -> None:
        self.verify_identity()
        if validate_inputs:
            for artifact in self.inputs:
                _validate_receipt_input(
                    artifact,
                    mode=validation_mode,
                )
        if config is not None and ConfigSnapshot.capture(config) != self.config:
            raise StaleArtifactError(
                f"{self.stage}/{self.donor_id}: configuration changed"
            )
        if validate_code:
            current = self.code.capture_current()
            if current.source_sha256 != self.code.source_sha256:
                raise StaleArtifactError(
                    f"{self.stage}/{self.donor_id}: producing code changed"
                )
        _validate_donor_snapshot(
            self.output_root,
            self.output,
            mode=validation_mode,
        )
        for sidecar in self.sidecars:
            sidecar.validate_current(mode=validation_mode)

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> DonorStageReceipt:
        result = cls(
            receipt_version=int(value["receipt_version"]),
            stage=str(value["stage"]),
            method_version=str(value["method_version"]),
            donor_id=str(value["donor_id"]),
            inputs=tuple(
                InputArtifact.from_dict(item) for item in value["inputs"]
            ),
            output_root=str(value["output_root"]),
            output=DonorDatasetSnapshot.from_dict(value["output"]),
            sidecars=tuple(
                FileFingerprint.from_dict(item)
                for item in value.get("sidecars", ())
            ),
            config=ConfigSnapshot.from_dict(value["config"]),
            code=CodeMetadata.from_dict(value["code"]),
            created_at_utc=str(value["created_at_utc"]),
            content_id=str(value["content_id"]),
        )
        result.verify_identity()
        return result

    @classmethod
    def read_json(cls, path: Path | str) -> DonorStageReceipt:
        with Path(path).open(encoding="utf-8") as handle:
            return cls.from_dict(json.load(handle))


def _validate_task_receipt(
    context,
    stage: str,
    donor: str,
    receipt: DonorStageReceipt,
) -> None:
    from .pipeline import STAGE_METHODS

    if receipt.stage != stage or receipt.donor_id != str(donor):
        raise ContractError(
            f"receipt identity differs from task {stage}/{donor}"
        )
    if receipt.method_version != STAGE_METHODS[stage]:
        raise StaleArtifactError(
            f"{stage}/{donor}: method version changed"
        )
    expected_root = context.output_root(stage).expanduser().resolve()
    if Path(receipt.output_root).expanduser().resolve() != expected_root:
        raise ContractError(
            f"{stage}/{donor}: receipt output root differs from "
            f"{expected_root}"
        )
    receipt.validate_current(
        config=_donor_stage_config(context, stage, donor),
        validation_mode="fast",
    )


def _dependency_input(context, stage: str, donor_id: str) -> InputArtifact:
    receipt_path = donor_receipt_path(context, stage, donor_id)
    if receipt_path.exists():
        return _input_from_receipt(
            stage,
            receipt_path,
            donor_id=donor_id,
        )
    manifest_path = context.stage_manifest_path(stage)
    if manifest_path.exists():
        return InputArtifact.from_stage_manifest(stage, manifest_path)
    raise PartialArtifactError(
        f"donor {donor_id} requires completed {stage!r}"
    )


def donor_stage_inputs(
    context,
    stage: str,
    donor_id: str,
) -> tuple[InputArtifact, ...]:
    inputs: list[InputArtifact] = []
    if stage in {"geometry", "redsea"}:
        inputs.append(
            InputArtifact.from_qupath_cohort_manifest(
                "qupath_cohort",
                context.base_config.qupath_manifest,
            )
        )
    if stage in {"geometry", "redsea", "expression"}:
        ingest_path = context.stage_manifest_path("ingest")
        if not ingest_path.exists():
            raise PartialArtifactError(
                f"donor {donor_id} requires completed ingest"
            )
        inputs.append(InputArtifact.from_stage_manifest("ingest", ingest_path))
    for dependency in DONOR_DEPENDENCIES[stage]:
        inputs.append(_dependency_input(context, dependency, donor_id))
    if stage in {"redsea", "expression", "controls", "calibrate", "type"}:
        inputs.append(
            InputArtifact.from_path(
                "marker_registry",
                context.base_config.marker_registry,
            )
        )
    if stage == "type":
        inputs.append(
            InputArtifact.from_path(
                "typing_rules",
                context.base_config.typing_rules,
            )
        )
        if context.typing_model_bundle is not None:
            from .pipeline import _validate_typing_calibration_compatibility

            _validate_typing_calibration_compatibility(context)
            artifact = context.typing_model_bundle.input_artifact
            if artifact is None:
                raise ContractError(
                    "configured typing model bundle has no captured input artifact"
                )
            inputs.append(artifact)
            promotion = context.typing_model_promotion_manifest
            if promotion is None or promotion.input_artifact is None:
                raise ContractError(
                    "configured typing model bundle requires a captured promotion manifest"
                )
            promotion.validate_for_bundle(context.typing_model_bundle)
            inputs.append(promotion.input_artifact)
    return tuple(inputs)


def _geometry_config(context, image) -> GeometryQCConfig:
    return GeometryQCConfig(
        pixel_size_um_x=image.pixel_size_um_x,
        pixel_size_um_y=image.pixel_size_um_y,
        min_cell_area_um2=context.config.cell_qc_min_cell_area,
        max_cell_area_median_multiple=(
            context.config.cell_qc_max_area_median_multiple
        ),
        min_cell_solidity=context.config.cell_qc_min_cell_solidity,
        min_raster_area_ratio=context.config.cell_qc_min_raster_area_ratio,
        max_raster_area_ratio=context.config.cell_qc_max_raster_area_ratio,
        no_cytoplasm_nc_ratio=context.config.cell_qc_no_cytoplasm_nc_ratio,
        duplicate_radius_um=context.config.cell_qc_duplicate_radius_um,
        duplicate_area_tolerance=context.config.cell_qc_duplicate_area_tol,
        missing_geometry_policy=(
            context.config.cell_qc_missing_geometry_policy
        ),
        duplicate_policy=context.config.cell_qc_duplicate_policy,
    )


def _run_donor_stage_body(context, stage: str, donor_id: str) -> tuple[Path, ...]:
    donor = str(donor_id)
    image_by_donor = {
        image.donor_id: image for image in context.cohort.images
    }
    if stage == "geometry":
        image = image_by_donor[donor]
        result = run_geometry_qc_for_image(
            read_single_partition(context.config.cells_dir, donor),
            image,
            config=_geometry_config(context, image),
        )
        write_donor_partition(
            result.to_frame(),
            context.config.geometry_qc_dir,
            donor_id=donor,
        )
        return ()
    if stage == "redsea":
        process_donor(
            donor,
            context.config,
            RedseaParams.from_config(context.config),
        )
        return ()
    if stage == "expression":
        image = image_by_donor[donor]
        expression, availability = build_selected_expression(
            read_single_partition(context.config.cells_dir, donor),
            read_single_partition(context.config.redsea_dir, donor),
            read_single_partition(context.config.geometry_qc_dir, donor),
            donor_id=donor,
            acquired_markers={
                channel.marker for channel in image.channel_map
            },
            registry=context.registry,
        )
        write_donor_partition(
            expression,
            context.config.selected_expression_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            availability,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "controls":
        expression = read_single_partition(
            context.config.selected_expression_dir,
            donor,
        )
        result = build_reference_controls(
            expression,
            donor_id=donor,
            registry=context.registry,
            value_columns={
                marker.name: marker.name
                for marker in context.registry.markers
            },
        )
        write_donor_partition(
            reference_masks_to_wide(
                result.masks_long,
                registry=context.registry,
            ),
            context.config.reference_controls_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            result.audit,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "calibrate":
        expression = read_single_partition(
            context.config.selected_expression_dir,
            donor,
        )
        masks = read_single_partition(
            context.config.reference_controls_dir,
            donor,
        )
        result = calibrate_expression(
            expression,
            donor_id=donor,
            registry=context.registry,
            reference_controls=reference_masks_from_wide(
                masks,
                registry=context.registry,
            ),
        )
        write_donor_partition(
            result.evidence_wide,
            context.config.marker_evidence_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            result.model_table,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    if stage == "type":
        assignments = type_calibrated_cells(
            read_single_partition(
                context.config.marker_evidence_dir,
                donor,
            ),
            marker_registry=context.registry,
            typing_registry=context.typing_registry,
            model_bundle=context.typing_model_bundle,
            model_promotion_manifest=context.typing_model_promotion_manifest,
        )
        write_donor_partition(
            assignments,
            context.config.assignments_dir,
            donor_id=donor,
        )
        return ()
    if stage == "states":
        states, models = annotate_proliferation(
            read_single_partition(
                context.config.selected_expression_dir,
                donor,
            ),
            donor_id=donor,
        )
        for row in models.itertuples(index=False):
            states[f"{row.marker}__threshold"] = row.threshold
            states[f"{row.marker}__model_status"] = row.status
            states[f"{row.marker}__method"] = row.method
        write_donor_partition(
            states,
            context.config.state_dir,
            donor_id=donor,
        )
        audit_path = donor_audit_path(context, stage, donor)
        write_donor_partition(
            models,
            audit_path.parents[1],
            donor_id=donor,
        )
        return (audit_path,)
    raise ContractError(f"unknown donor stage {stage!r}")


def _partition_has_content(root: Path, donor_id: str) -> bool:
    partition = root / f"donor_id={donor_id}"
    return partition.exists() and any(partition.iterdir())


def execute_donor_stage(context, stage: str, donor_id: str) -> str:
    """Run one ready donor task and atomically publish its receipt."""

    from .pipeline import STAGE_METHODS, _stage_code

    donor = str(donor_id)
    path = donor_receipt_path(context, stage, donor)
    if path.exists():
        receipt = DonorStageReceipt.read_json(path)
        _validate_task_receipt(context, stage, donor, receipt)
        return receipt.content_id
    output_root = context.output_root(stage)
    if _partition_has_content(output_root, donor):
        raise PartialArtifactError(
            f"{stage}/{donor}: output exists without a donor receipt"
        )
    audit_path = donor_audit_path(context, stage, donor)
    if audit_path is not None and _partition_has_content(
        audit_path.parents[1],
        donor,
    ):
        raise PartialArtifactError(
            f"{stage}/{donor}: audit output exists without a donor receipt"
        )

    inputs = donor_stage_inputs(context, stage, donor)
    for artifact in inputs:
        _validate_receipt_input(artifact, mode="fast")
    print(f"[donor-run] {stage}/{donor}", flush=True)
    sidecars = _run_donor_stage_body(context, stage, donor)
    receipt = DonorStageReceipt.create(
        stage=stage,
        method_version=STAGE_METHODS[stage],
        donor_id=donor,
        inputs=inputs,
        output_root=output_root,
        sidecar_paths=sidecars,
        config=_donor_stage_config(context, stage, donor),
        code=_stage_code(context, stage),
    )
    receipt.write_json(path)
    print(
        f"[donor-done] {stage}/{donor}: "
        f"{receipt.output.row_count:,} rows, receipt={receipt.content_id[:12]}",
        flush=True,
    )
    return receipt.content_id


def create_receipt_for_existing_output(
    context,
    stage: str,
    donor_id: str,
    *,
    sidecar_paths: Sequence[Path | str] = (),
) -> DonorStageReceipt:
    """Seal externally verified donor output without recomputing its values."""

    from .pipeline import STAGE_METHODS, _stage_code

    donor = str(donor_id)
    path = donor_receipt_path(context, stage, donor)
    if path.exists():
        receipt = DonorStageReceipt.read_json(path)
        _validate_task_receipt(context, stage, donor, receipt)
        return receipt
    receipt = DonorStageReceipt.create(
        stage=stage,
        method_version=STAGE_METHODS[stage],
        donor_id=donor,
        inputs=donor_stage_inputs(context, stage, donor),
        output_root=context.output_root(stage),
        sidecar_paths=sidecar_paths,
        config=_donor_stage_config(context, stage, donor),
        code=_stage_code(context, stage),
    )
    receipt.write_json(path)
    print(
        f"[donor-receipt] {stage}/{donor}: {receipt.content_id[:12]}",
        flush=True,
    )
    return receipt


@dataclass(frozen=True, order=True)
class DonorTask:
    stage: str
    donor_id: str


@dataclass(frozen=True)
class ResourceRequest:
    group: str
    cpu_slots: int = 1
    gpu_slots: int = 0


@dataclass(frozen=True)
class ResourceLimits:
    total_cpu_slots: int
    gpu_slots: int
    geometry_workers: int
    redsea_workers: int
    downstream_workers: int

    @classmethod
    def from_context(
        cls,
        context,
        *,
        task_count: int,
    ) -> ResourceLimits:
        total = resolve_jobs(context.config.n_jobs, task_count)
        redsea = min(context.config.pipeline_redsea_workers, total)
        if context.config.use_gpu:
            redsea = 1
        return cls(
            total_cpu_slots=total,
            gpu_slots=1 if context.config.use_gpu else 0,
            geometry_workers=min(
                context.config.pipeline_geometry_workers,
                total,
            ),
            redsea_workers=redsea,
            downstream_workers=min(
                context.config.pipeline_downstream_workers,
                total,
            ),
        )

    def group_limit(self, group: str) -> int:
        if group == "geometry":
            return self.geometry_workers
        if group == "redsea":
            return self.redsea_workers
        if group == "downstream":
            return self.downstream_workers
        raise ContractError(f"unknown resource group {group!r}")


@dataclass
class ResourceState:
    cpu_slots: int = 0
    gpu_slots: int = 0
    geometry_workers: int = 0
    redsea_workers: int = 0
    downstream_workers: int = 0

    def _group_value(self, group: str) -> int:
        return int(getattr(self, f"{group}_workers"))

    def can_start(
        self,
        request: ResourceRequest,
        limits: ResourceLimits,
    ) -> bool:
        return (
            self.cpu_slots + request.cpu_slots <= limits.total_cpu_slots
            and self.gpu_slots + request.gpu_slots <= limits.gpu_slots
            and self._group_value(request.group)
            < limits.group_limit(request.group)
        )

    def reserve(self, request: ResourceRequest) -> None:
        self.cpu_slots += request.cpu_slots
        self.gpu_slots += request.gpu_slots
        attribute = f"{request.group}_workers"
        setattr(self, attribute, getattr(self, attribute) + 1)

    def release(self, request: ResourceRequest) -> None:
        self.cpu_slots -= request.cpu_slots
        self.gpu_slots -= request.gpu_slots
        attribute = f"{request.group}_workers"
        setattr(self, attribute, getattr(self, attribute) - 1)


def task_resources(context, task: DonorTask) -> ResourceRequest:
    if task.stage == "geometry":
        return ResourceRequest(group="geometry")
    if task.stage == "redsea":
        return ResourceRequest(
            group="redsea",
            gpu_slots=1 if context.config.use_gpu else 0,
        )
    return ResourceRequest(group="downstream")


def task_ready(task: DonorTask, completed: set[DonorTask]) -> bool:
    return all(
        DonorTask(dependency, task.donor_id) in completed
        for dependency in DONOR_DEPENDENCIES[task.stage]
    )


def _cohort_barrier_ready(context, task: DonorTask) -> bool:
    """Hold probabilistic typing until cohort calibration is sealed.

    Donor-local calibration receipts are sufficient for rules-only typing. A
    fitted model, however, is authorized for one exact cohort calibration
    content hash, which cannot be checked until the aggregate stage manifest
    exists.
    """

    return not (
        task.stage == "type"
        and getattr(context, "typing_model_bundle", None) is not None
        and not context.stage_manifest_path("calibrate").is_file()
    )


def _task_priority(
    task: DonorTask,
    *,
    donor_order: Mapping[str, int],
) -> tuple[int, int, int]:
    stage_index = DONOR_STAGES.index(task.stage)
    if task.stage == "redsea":
        group_priority = 0
    elif task.stage == "geometry":
        group_priority = 1
    else:
        group_priority = 2
    downstream_priority = {
        "expression": 0,
        "controls": 1,
        "calibrate": 2,
        "type": 3,
        "states": 4,
    }
    return (
        group_priority,
        donor_order[task.donor_id],
        (
            downstream_priority[task.stage]
            if group_priority == 2
            else stage_index
        ),
    )


def select_ready_task(
    pending: set[DonorTask],
    completed: set[DonorTask],
    *,
    context,
    resources: ResourceState,
    limits: ResourceLimits,
) -> DonorTask | None:
    donor_order = {
        donor: index for index, donor in enumerate(context.donors)
    }
    candidates = sorted(
        (
            task
            for task in pending
            if task_ready(task, completed)
            and _cohort_barrier_ready(context, task)
            and resources.can_start(task_resources(context, task), limits)
        ),
        key=lambda task: _task_priority(
            task,
            donor_order=donor_order,
        ),
    )
    return candidates[0] if candidates else None


def _aggregate_donor_audit(context, stage: str) -> Path | None:
    from .pipeline import STAGE_SIDECARS

    name = DONOR_AUDITS.get(stage)
    if name is None:
        return None
    frames = []
    for donor in context.donors:
        path = donor_audit_path(context, stage, donor)
        if path is None or not path.is_file():
            raise PartialArtifactError(
                f"{stage}: donor audit is missing for {donor}"
            )
        frames.append(pd.read_parquet(path))
    combined = pd.concat(frames, ignore_index=True)
    relative_paths = STAGE_SIDECARS.get(stage, ())
    if len(relative_paths) != 1:
        raise ContractError(
            f"{stage}: expected exactly one aggregate audit sidecar"
        )
    target = context.run_root / relative_paths[0]
    return _atomic_write_frame(target, combined)


def _manifest_is_current(context, stage: str) -> bool:
    from .pipeline import _stage_config

    path = context.stage_manifest_path(stage)
    if not path.exists():
        return False
    manifest = StageManifest.read_json(path)
    if manifest.stage != stage:
        raise ContractError(
            f"manifest at {path} declares stage {manifest.stage!r}, "
            f"expected {stage!r}"
        )
    manifest.validate_current(
        config=_stage_config(context, stage),
        validation_mode="fast",
    )
    return True


def _receipt_is_current(context, stage: str, donor: str) -> bool:
    path = donor_receipt_path(context, stage, donor)
    if not path.exists():
        return False
    receipt = DonorStageReceipt.read_json(path)
    _validate_task_receipt(context, stage, donor, receipt)
    return True


def _seal_ready_stages(
    context,
    stages: Sequence[str],
    sealed: set[str],
) -> None:
    from .pipeline import (
        STAGE_METHODS,
        _stage_code,
        _stage_config,
        _stage_inputs,
    )

    changed = True
    while changed:
        changed = False
        for stage in stages:
            if stage in sealed:
                continue
            receipt_paths = [
                donor_receipt_path(context, stage, donor)
                for donor in context.donors
            ]
            if not all(path.exists() for path in receipt_paths):
                continue
            for donor, path in zip(
                context.donors,
                receipt_paths,
                strict=True,
            ):
                receipt = DonorStageReceipt.read_json(path)
                _validate_task_receipt(
                    context,
                    stage,
                    donor,
                    receipt,
                )
            try:
                inputs = _stage_inputs(context, stage)
            except PartialArtifactError:
                continue
            for artifact in inputs:
                artifact.validate_current(mode="fast")
            sidecars: list[Path] = list(receipt_paths)
            audit = _aggregate_donor_audit(context, stage)
            if audit is not None:
                sidecars.append(audit)
            if stage == "geometry":
                resume_receipt = (
                    context.config.manifests_dir / "geometry_resume.json"
                )
                if resume_receipt.exists():
                    sidecars.append(resume_receipt)
            manifest = StageManifest.create(
                stage=stage,
                method_version=STAGE_METHODS[stage],
                expected_donors=context.donors,
                inputs=inputs,
                output_root=context.output_root(stage),
                sidecar_paths=sidecars,
                config=_stage_config(context, stage),
                code=_stage_code(context, stage),
            )
            path = context.stage_manifest_path(stage)
            path.parent.mkdir(parents=True, exist_ok=True)
            manifest.write_json(path)
            sealed.add(stage)
            changed = True
            print(
                f"[sealed] {stage}: {manifest.output.total_rows:,} rows, "
                f"manifest={manifest.content_id[:12]}",
                flush=True,
            )


def _initial_task_state(
    context,
    stages: Sequence[str],
) -> tuple[set[DonorTask], set[str]]:
    completed: set[DonorTask] = set()
    sealed: set[str] = set()
    for stage in stages:
        if _manifest_is_current(context, stage):
            sealed.add(stage)
            completed.update(
                DonorTask(stage, donor) for donor in context.donors
            )
            continue
        for donor in context.donors:
            if _receipt_is_current(context, stage, donor):
                completed.add(DonorTask(stage, donor))
                continue
            if _partition_has_content(context.output_root(stage), donor):
                raise PartialArtifactError(
                    f"{stage}/{donor}: output exists without a donor receipt"
                )
            audit = donor_audit_path(context, stage, donor)
            if audit is not None and _partition_has_content(
                audit.parents[1],
                donor,
            ):
                raise PartialArtifactError(
                    f"{stage}/{donor}: audit exists without a donor receipt"
                )
    return completed, sealed


def _run_tasks_serial(
    context,
    pending: set[DonorTask],
    completed: set[DonorTask],
    sealed: set[str],
    stages: Sequence[str],
    limits: ResourceLimits,
) -> None:
    resources = ResourceState()
    while pending:
        _seal_ready_stages(context, stages, sealed)
        task = select_ready_task(
            pending,
            completed,
            context=context,
            resources=resources,
            limits=limits,
        )
        if task is None:
            blocked = sorted(
                f"{item.stage}/{item.donor_id}" for item in pending
            )
            raise PartialArtifactError(
                f"donor task queue is deadlocked: {blocked[:10]}"
            )
        request = task_resources(context, task)
        resources.reserve(request)
        try:
            execute_donor_stage(context, task.stage, task.donor_id)
        finally:
            resources.release(request)
        pending.remove(task)
        completed.add(task)
    _seal_ready_stages(context, stages, sealed)


def _run_tasks_parallel(
    context,
    pending: set[DonorTask],
    completed: set[DonorTask],
    sealed: set[str],
    stages: Sequence[str],
    limits: ResourceLimits,
) -> None:
    resources = ResourceState()
    running = {}
    spawn_context = multiprocessing.get_context("spawn")
    executors = {
        group: ProcessPoolExecutor(
            max_workers=limits.group_limit(group),
            mp_context=spawn_context,
        )
        for group in ("geometry", "redsea", "downstream")
    }
    aborting = False
    try:
        while pending or running:
            _seal_ready_stages(context, stages, sealed)
            while len(running) < limits.total_cpu_slots:
                task = select_ready_task(
                    pending,
                    completed,
                    context=context,
                    resources=resources,
                    limits=limits,
                )
                if task is None:
                    break
                request = task_resources(context, task)
                resources.reserve(request)
                pending.remove(task)
                future = executors[request.group].submit(
                    execute_donor_stage,
                    context,
                    task.stage,
                    task.donor_id,
                )
                running[future] = (task, request)
            if not running:
                blocked = sorted(
                    f"{item.stage}/{item.donor_id}" for item in pending
                )
                raise PartialArtifactError(
                    f"donor task queue is deadlocked: {blocked[:10]}"
                )
            finished, _ = wait(
                tuple(running),
                return_when=FIRST_COMPLETED,
            )
            for future in finished:
                task, request = running.pop(future)
                resources.release(request)
                try:
                    future.result()
                except Exception as exc:
                    raise RuntimeError(
                        f"donor task failed: {task.stage}/{task.donor_id}"
                    ) from exc
                receipt = DonorStageReceipt.read_json(
                    donor_receipt_path(
                        context,
                        task.stage,
                        task.donor_id,
                    )
                )
                _validate_task_receipt(
                    context,
                    task.stage,
                    task.donor_id,
                    receipt,
                )
                completed.add(task)
        _seal_ready_stages(context, stages, sealed)
    except BaseException:
        aborting = True
        for future in running:
            future.cancel()
        _terminate_executor_workers(executors)
        raise
    finally:
        if not aborting:
            for executor in executors.values():
                executor.shutdown(wait=True, cancel_futures=True)


def _terminate_executor_workers(
    executors: Mapping[str, ProcessPoolExecutor],
) -> None:
    """Stop active pool children so an interrupted queue cannot become orphaned.

    Python 3.12 does not expose ``terminate_workers`` on
    :class:`ProcessPoolExecutor`. The executor's process mapping is therefore
    captured before non-waiting shutdown, with TERM followed by a bounded KILL
    fallback. Only children owned by these three queue executors are targeted.
    """

    process_count = sum(
        len(getattr(executor, "_processes", None) or ())
        for executor in executors.values()
    )
    print(
        f"[queue] stopping {process_count} active worker processes",
        flush=True,
    )
    terminate_process_pool_workers(executors)


def run_donor_queue(
    context,
    stages: Sequence[str],
) -> None:
    """Execute a donor-stage DAG under explicit CPU/GPU/group limits."""

    selected = tuple(stages)
    unknown = sorted(set(selected).difference(DONOR_STAGES))
    if unknown:
        raise ContractError(f"unknown donor stages: {unknown}")
    expected_prefix = DONOR_STAGES[: len(selected)]
    if selected != expected_prefix:
        raise ContractError(
            "pipelined donor stages must be a dependency prefix"
        )
    completed, sealed = _initial_task_state(context, selected)
    all_tasks = {
        DonorTask(stage, donor)
        for stage in selected
        for donor in context.donors
    }
    pending = all_tasks.difference(completed)
    limits = ResourceLimits.from_context(
        context,
        task_count=max(1, len(pending)),
    )
    print(
        "[queue] "
        f"cpu={limits.total_cpu_slots}, gpu={limits.gpu_slots}, "
        f"geometry={limits.geometry_workers}, "
        f"redsea={limits.redsea_workers}, "
        f"downstream={limits.downstream_workers}, "
        f"pending={len(pending)}",
        flush=True,
    )
    if not pending:
        _seal_ready_stages(context, selected, sealed)
        return
    if limits.total_cpu_slots == 1:
        _run_tasks_serial(
            context,
            pending,
            completed,
            sealed,
            selected,
            limits,
        )
    else:
        _run_tasks_parallel(
            context,
            pending,
            completed,
            sealed,
            selected,
            limits,
        )


def run_pipelined_pipeline(
    context,
    *,
    stages: Sequence[str],
    export: bool,
):
    """Run ingest as a barrier, then overlap ready donor stages."""

    from .pipeline import (
        STAGE_ORDER,
        export_qupath,
        run_pipeline,
        run_stage,
    )

    selected = tuple(stages)
    expected_prefix = STAGE_ORDER[: len(selected)]
    if selected != expected_prefix:
        raise ContractError(
            "pipelined execution requires a dependency-prefix stage selection"
        )
    context.run_root.mkdir(parents=True, exist_ok=True)
    if not selected:
        return run_pipeline(context, stages=(), export=False)
    run_stage(context, "ingest")
    donor_stages = tuple(
        stage for stage in selected if stage in DONOR_STAGES
    )
    if donor_stages:
        run_donor_queue(context, donor_stages)
    if export and "type" in selected:
        export_qupath(context)
    return run_pipeline(context, stages=(), export=False)


def donor_progress(context, stage: str) -> tuple[int, int]:
    """Return current donor receipt count for a stage."""

    if stage not in DONOR_STAGES:
        return 0, len(context.donors)
    count = sum(
        donor_receipt_path(context, stage, donor).exists()
        for donor in context.donors
    )
    return count, len(context.donors)


__all__ = [
    "DONOR_DEPENDENCIES",
    "DONOR_RECEIPT_VERSION",
    "DONOR_STAGES",
    "DonorStageReceipt",
    "DonorTask",
    "ResourceLimits",
    "ResourceRequest",
    "ResourceState",
    "create_receipt_for_existing_output",
    "donor_progress",
    "donor_receipt_path",
    "donor_stage_inputs",
    "execute_donor_stage",
    "run_donor_queue",
    "run_pipelined_pipeline",
    "select_ready_task",
    "snapshot_donor_partition",
    "task_ready",
    "task_resources",
]
