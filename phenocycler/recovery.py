"""Verified recovery of a stopped pre-receipt run into the donor task queue."""

from __future__ import annotations

import argparse
import hashlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Sequence

import pyarrow as pa
import pyarrow.parquet as pq

from .artifacts import (
    CodeMetadata,
    ContractError,
    DonorDatasetSnapshot,
    PartialArtifactError,
    SchemaContractError,
    StageManifest,
    StaleArtifactError,
    hash_identifiers,
)
from .config import load_config
from .donor_pipeline import (
    create_receipt_for_existing_output,
    run_pipelined_pipeline,
    snapshot_donor_partition,
)
from .pipeline import (
    STAGE_METHODS,
    STAGE_ORDER,
    RunContext,
    _selected_stages,
    resolve_run_context,
)
from .redsea import COMPARTMENTS
from .resume import (
    _canonical_stage_config,
    _content_payload,
    _hardlink_once,
    _reuse_ingest,
    _sampled_fingerprint,
    _source_run_root,
    _validate_ingest_reuse_contract,
    _write_json_once,
)


RECOVERY_CONTRACT_VERSION = 1


@dataclass(frozen=True)
class RecoveredRedseaPartition:
    donor_id: str
    path: str
    row_count: int
    object_id_sha256: str
    schema_sha256: str
    size_bytes: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(frozen=True)
class RecoveryPreparation:
    source_run_id: str
    target_run_id: str
    geometry_donors: tuple[str, ...]
    redsea_donors: tuple[str, ...]


def _expected_redsea_columns(context: RunContext) -> tuple[str, ...]:
    return (
        "object_id",
        "donor_id",
        "cell_area_px",
        "band_area_px",
        "nucleus_area_px",
        "cytoplasm_area_px",
        *(
            f"{compartment}__{marker.name}"
            for compartment in COMPARTMENTS
            for marker in context.registry.markers
        ),
    )


def validate_recovered_redsea_partition(
    path: Path | str,
    *,
    donor_id: str,
    expected: DonorDatasetSnapshot,
    expected_columns: Sequence[str],
    expected_schema: pa.Schema | None = None,
    scan_all_columns: bool = True,
) -> tuple[RecoveredRedseaPartition, pa.Schema]:
    """Prove that one unmanifested REDSEA file is complete and canonical."""

    path = Path(path).expanduser().resolve()
    try:
        parquet = pq.ParquetFile(path)
    except (OSError, pa.ArrowException) as exc:
        raise PartialArtifactError(
            f"{path}: cannot open completed REDSEA Parquet"
        ) from exc
    schema = parquet.schema_arrow
    if tuple(schema.names) != tuple(expected_columns):
        raise SchemaContractError(
            f"{path}: REDSEA columns differ from the current contract"
        )
    if expected_schema is not None and not schema.equals(
        expected_schema,
        check_metadata=True,
    ):
        raise SchemaContractError(
            f"{path}: REDSEA schema differs from preceding recovered donors"
        )
    for name in ("object_id", "donor_id"):
        field = schema.field(name)
        if not (
            pa.types.is_string(field.type)
            or pa.types.is_large_string(field.type)
        ):
            raise SchemaContractError(
                f"{path}: {name} must use an Arrow string type"
            )
    for name in tuple(expected_columns)[2:]:
        if not pa.types.is_floating(schema.field(name).type):
            raise SchemaContractError(
                f"{path}: REDSEA value column {name!r} is not floating point"
            )

    try:
        identity = parquet.read(columns=["object_id", "donor_id"])
        if scan_all_columns and parquet.scan_contents() != parquet.metadata.num_rows:
            raise PartialArtifactError(
                f"{path}: full REDSEA scan returned an incomplete row count"
            )
    except (OSError, pa.ArrowException) as exc:
        raise PartialArtifactError(
            f"{path}: REDSEA data pages are incomplete or unreadable"
        ) from exc
    if identity.num_rows != expected.row_count:
        raise PartialArtifactError(
            f"{path}: REDSEA row count {identity.num_rows:,} differs from "
            f"ingest {expected.row_count:,}"
        )
    if identity.column("donor_id").null_count:
        raise ContractError(f"{path}: donor_id contains null")
    stored_donors = {
        str(value)
        for value in identity.column("donor_id").to_pylist()
        if value is not None
    }
    if stored_donors != {str(donor_id)}:
        raise ContractError(
            f"{path}: stored donor values differ from {donor_id}"
        )
    object_hash = hash_identifiers(
        identity.column("object_id").to_pylist()
    )
    if object_hash != expected.object_id_sha256:
        raise StaleArtifactError(
            f"{path}: REDSEA object universe differs from ingest"
        )
    schema_hash = hashlib.sha256(
        schema.serialize().to_pybytes()
    ).hexdigest()
    return (
        RecoveredRedseaPartition(
            donor_id=str(donor_id),
            path=str(path),
            row_count=identity.num_rows,
            object_id_sha256=object_hash,
            schema_sha256=schema_hash,
            size_bytes=path.stat().st_size,
        ),
        schema,
    )


def _source_geometry_manifest(
    context: RunContext,
    source_root: Path,
    source_ingest: StageManifest,
) -> StageManifest:
    path = source_root / "manifests" / "geometry.json"
    if not path.is_file():
        raise PartialArtifactError(
            f"source run has no completed geometry manifest: {path}"
        )
    manifest = StageManifest.read_json(path)
    if manifest.stage != "geometry":
        raise ContractError("source geometry manifest declares another stage")
    if manifest.method_version != STAGE_METHODS["geometry"]:
        raise StaleArtifactError("source geometry method version differs")
    if manifest.expected_donors != context.donors:
        raise PartialArtifactError("source geometry donor contract differs")
    manifest.validate_current(validation_mode="fast")

    recorded = _canonical_stage_config(manifest)
    if recorded.get("configuration") != context.identity_config:
        raise StaleArtifactError(
            "source geometry scientific configuration differs"
        )
    if recorded.get("registry_fingerprint") != context.registry.fingerprint:
        raise StaleArtifactError("source geometry marker registry differs")
    if (
        recorded.get("typing_rules_fingerprint")
        != context.typing_registry.rules_fingerprint
    ):
        raise StaleArtifactError("source geometry typing registry differs")

    inputs = {artifact.name: artifact for artifact in manifest.inputs}
    if (
        "ingest" not in inputs
        or inputs["ingest"].content_sha256 != source_ingest.content_id
    ):
        raise StaleArtifactError(
            "source geometry does not depend on the source ingest manifest"
        )
    expected_root = (source_root / "geometry_qc").resolve()
    if Path(manifest.output.root).resolve() != expected_root:
        raise ContractError(
            "source geometry manifest records a noncanonical output root"
        )
    if (
        manifest.output.total_rows != source_ingest.output.total_rows
        or manifest.output.object_id_sha256
        != source_ingest.output.object_id_sha256
    ):
        raise StaleArtifactError(
            "source geometry object universe differs from source ingest"
        )
    return manifest


def _link_manifest_dataset(
    context: RunContext,
    *,
    stage: str,
    source: StageManifest,
) -> list[dict[str, Any]]:
    source_root = Path(source.output.root).resolve()
    target_root = context.output_root(stage).resolve()
    expected_targets: set[Path] = set()
    links: list[dict[str, Any]] = []
    for donor in source.output.donors:
        for fingerprint in donor.files:
            source_file = Path(fingerprint.path).resolve()
            try:
                relative = source_file.relative_to(source_root)
            except ValueError as exc:
                raise ContractError(
                    f"{stage}: source file is outside its manifest root"
                ) from exc
            if (
                not relative.parts
                or relative.parts[0] != f"donor_id={donor.donor_id}"
            ):
                raise ContractError(
                    f"{stage}: source file is in the wrong donor partition"
                )
            target_file = target_root / relative
            _hardlink_once(source_file, target_file)
            expected_targets.add(target_file.resolve())
            links.append(
                {
                    "donor_id": donor.donor_id,
                    "source": fingerprint.path,
                    "target": str(target_file.resolve()),
                    "size_bytes": fingerprint.size_bytes,
                    "sha256": fingerprint.sha256,
                }
            )

    actual_targets = (
        {
            path.resolve()
            for path in target_root.rglob("*")
            if path.is_file()
        }
        if target_root.exists()
        else set()
    )
    if actual_targets != expected_targets:
        raise PartialArtifactError(
            f"{stage}: recovery target contains foreign files"
        )
    actual_donors = {
        path.name.split("=", 1)[1]
        for path in target_root.glob("donor_id=*")
        if path.is_dir()
    }
    if actual_donors != set(context.donors):
        raise PartialArtifactError(
            f"{stage}: recovery target donor set differs"
        )
    return links


def _recover_geometry(
    context: RunContext,
    source_root: Path,
    source: StageManifest,
) -> None:
    links = _link_manifest_dataset(
        context,
        stage="geometry",
        source=source,
    )
    receipt_path = context.config.manifests_dir / "geometry_recovery.json"
    payload = _content_payload(
        {
            "recovery_contract_version": RECOVERY_CONTRACT_VERSION,
            "stage": "geometry",
            "transfer": "verified_hard_links",
            "source_run_id": source_root.name,
            "source_manifest_content_id": source.content_id,
            "source_output_content_sha256": source.output.content_sha256,
            "target_run_id": context.run_id,
            "recovery_code": asdict(
                CodeMetadata.capture(
                    context.repository,
                    source_paths=("phenocycler/recovery.py",),
                )
            ),
            "files": links,
        }
    )
    _write_json_once(receipt_path, payload)

    expected = {
        donor.donor_id: donor for donor in source.output.donors
    }
    for donor in context.donors:
        snapshot = snapshot_donor_partition(
            context.config.geometry_qc_dir,
            donor,
        )
        if snapshot.content_sha256 != expected[donor].content_sha256:
            raise StaleArtifactError(
                f"geometry/{donor}: recovered content differs from source"
            )
        create_receipt_for_existing_output(
            context,
            "geometry",
            donor,
            sidecar_paths=(receipt_path,),
        )


def _discover_redsea_sources(
    context: RunContext,
    source_root: Path,
    source_ingest: StageManifest,
) -> tuple[
    tuple[tuple[str, Path, RecoveredRedseaPartition], ...],
    tuple[str, ...],
]:
    root = source_root / "redsea"
    if not root.exists():
        return (), ()
    partitions = {
        path.name.split("=", 1)[1]: path
        for path in root.glob("donor_id=*")
        if path.is_dir()
    }
    extra = sorted(set(partitions).difference(context.donors))
    if extra:
        raise PartialArtifactError(
            f"source REDSEA contains donors outside the cohort: {extra}"
        )
    expected = {
        donor.donor_id: donor for donor in source_ingest.output.donors
    }
    columns = _expected_redsea_columns(context)
    recovered = []
    incomplete = []
    reference_schema = None
    for donor in context.donors:
        partition = partitions.get(donor)
        if partition is None:
            incomplete.append(donor)
            continue
        files = sorted(path for path in partition.iterdir() if path.is_file())
        if not files:
            incomplete.append(donor)
            continue
        if len(files) != 1 or files[0].name != "data_0.parquet":
            raise PartialArtifactError(
                f"source REDSEA donor {donor} has a partial/foreign inventory"
            )
        check, schema = validate_recovered_redsea_partition(
            files[0],
            donor_id=donor,
            expected=expected[donor],
            expected_columns=columns,
            expected_schema=reference_schema,
            scan_all_columns=True,
        )
        reference_schema = schema
        recovered.append((donor, files[0].resolve(), check))
        print(
            f"[recover-check] redsea/{donor}: "
            f"{check.row_count:,} rows, schema={check.schema_sha256[:12]}",
            flush=True,
        )
    return tuple(recovered), tuple(incomplete)


def _recover_redsea(
    context: RunContext,
    source_root: Path,
    recovered: Sequence[
        tuple[str, Path, RecoveredRedseaPartition]
    ],
    incomplete: Sequence[str],
    source_geometry: StageManifest,
) -> None:
    if not recovered:
        return
    target_root = context.config.redsea_dir
    links = []
    expected_targets: set[Path] = set()
    for donor, source_file, check in recovered:
        target = target_root / f"donor_id={donor}" / "data_0.parquet"
        _hardlink_once(source_file, target)
        expected_targets.add(target.resolve())
        links.append(
            {
                "donor_id": donor,
                "source": _sampled_fingerprint(source_file),
                "target": str(target.resolve()),
                "validation": check.to_dict(),
            }
        )
    actual_targets = (
        {
            path.resolve()
            for path in target_root.rglob("*")
            if path.is_file()
        }
        if target_root.exists()
        else set()
    )
    if actual_targets != expected_targets:
        raise PartialArtifactError(
            "redsea: recovery target contains foreign files"
        )

    receipt_path = context.config.manifests_dir / "redsea_recovery.json"
    payload = _content_payload(
        {
            "recovery_contract_version": RECOVERY_CONTRACT_VERSION,
            "stage": "redsea",
            "transfer": "verified_hard_links",
            "source_run_id": source_root.name,
            "source_geometry_manifest_content_id": source_geometry.content_id,
            "target_run_id": context.run_id,
            "recovered_donors": [item[0] for item in recovered],
            "unfinished_donors": list(incomplete),
            "acceptance": (
                "complete Parquet footer and full-column scan; exact current "
                "schema, donor, ingest row count, and UUID universe"
            ),
            "source_publication_note": (
                "source process used the pre-receipt direct writer; only "
                "fully readable final Parquet files are adopted"
            ),
            "recovery_code": asdict(
                CodeMetadata.capture(
                    context.repository,
                    source_paths=("phenocycler/recovery.py",),
                )
            ),
            "files": links,
        }
    )
    _write_json_once(receipt_path, payload)
    for donor, _, _ in recovered:
        create_receipt_for_existing_output(
            context,
            "redsea",
            donor,
            sidecar_paths=(receipt_path,),
        )


def prepare_interrupted_run_recovery(
    context: RunContext,
    *,
    source_run_id: str,
) -> RecoveryPreparation:
    """Adopt only strongly validated completed artifacts from a stopped run."""

    source_root = _source_run_root(context, source_run_id)
    source_ingest_path = source_root / "manifests" / "ingest.json"
    source_ingest = StageManifest.read_json(source_ingest_path)
    source_ingest.validate_current(validation_mode="fast")
    _validate_ingest_reuse_contract(context, source_ingest)
    source_geometry = _source_geometry_manifest(
        context,
        source_root,
        source_ingest,
    )
    recovered, incomplete = _discover_redsea_sources(
        context,
        source_root,
        source_ingest,
    )

    _, target_ingest = _reuse_ingest(context, source_root)
    if (
        target_ingest.output.content_sha256
        != source_ingest.output.content_sha256
    ):
        raise StaleArtifactError(
            "recovered ingest content differs from the stopped run"
        )
    _recover_geometry(context, source_root, source_geometry)
    _recover_redsea(
        context,
        source_root,
        recovered,
        incomplete,
        source_geometry,
    )
    result = RecoveryPreparation(
        source_run_id=source_run_id,
        target_run_id=context.run_id,
        geometry_donors=context.donors,
        redsea_donors=tuple(item[0] for item in recovered),
    )
    print(
        f"[recover-ready] geometry={len(result.geometry_donors)}, "
        f"redsea={len(result.redsea_donors)}, target={result.target_run_id}",
        flush=True,
    )
    return result


def recover(
    *,
    config_path: Path | None,
    source_run_id: str,
    jobs: int | None = None,
    through: str | None = None,
    export: bool = True,
    prepare_only: bool = False,
) -> RunContext:
    overrides = {}
    if jobs is not None:
        overrides["n_jobs"] = jobs
    context = resolve_run_context(load_config(config_path, **overrides))
    prepare_interrupted_run_recovery(
        context,
        source_run_id=source_run_id,
    )
    if not prepare_only:
        run_pipelined_pipeline(
            context,
            stages=_selected_stages(None, through),
            export=export,
        )
    return context


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Recover a stopped pre-receipt run and continue it with the "
            "resource-aware donor queue."
        )
    )
    parser.add_argument("--config", type=Path)
    parser.add_argument("--from-run", required=True)
    parser.add_argument("--jobs", type=int)
    parser.add_argument("--through", choices=STAGE_ORDER[1:])
    parser.add_argument("--no-export", action="store_true")
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args(argv)
    recover(
        config_path=args.config,
        source_run_id=args.from_run,
        jobs=args.jobs,
        through=args.through,
        export=not args.no_export,
        prepare_only=args.prepare_only,
    )
    return 0


__all__ = [
    "RECOVERY_CONTRACT_VERSION",
    "RecoveredRedseaPartition",
    "RecoveryPreparation",
    "main",
    "prepare_interrupted_run_recovery",
    "recover",
    "validate_recovered_redsea_partition",
]
