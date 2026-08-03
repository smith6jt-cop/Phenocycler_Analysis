"""Verified reuse of a completed ingest and duplicate-free geometry prefix."""

from __future__ import annotations

import argparse
import errno
import hashlib
import json
import os
import re
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import pandas as pd
import pyarrow.parquet as pq

from .artifacts import (
    ContractError,
    DonorDatasetSnapshot,
    FileFingerprint,
    PartialArtifactError,
    SchemaContractError,
    StageManifest,
    StaleArtifactError,
    hash_identifiers,
    sha256_json,
)
from .config import load_config
from .geometry_qc import normalize_geometry_output
from .pipeline import (
    STAGE_METHODS,
    STAGE_ORDER,
    STAGE_SIDECARS,
    RunContext,
    _selected_stages,
    _stage_code,
    _stage_config,
    _stage_inputs,
    resolve_run_context,
)


RESUME_CONTRACT_VERSION = 1
RUN_ID_PATTERN = re.compile(r"[0-9a-f]{20}")
GEOMETRY_REQUIRED_COLUMNS = (
    "donor_id",
    "object_id",
    "analysis_eligible",
    "estimation_eligible",
    "spillover_context_eligible",
    "analysis_reasons",
    "estimation_reasons",
    "spillover_context_reasons",
    "duplicate_group",
    "duplicate_survivor",
    "raster_area_ratio",
    "nucleus_cell_area_ratio",
    "analysis_reason",
    "estimation_reason",
    "spillover_context_reason",
)


@dataclass(frozen=True)
class GeometryPartitionCheck:
    donor_id: str
    row_count: int
    object_id_sha256: str
    schema_sha256: str
    duplicate_group_count: int

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _content_payload(payload: Mapping[str, Any]) -> dict[str, Any]:
    result = dict(payload)
    result["content_id"] = sha256_json(result)
    return result


def _write_json_once(path: Path, payload: Mapping[str, Any]) -> Path:
    """Atomically write an immutable JSON receipt, accepting an exact rerun."""

    expected = dict(payload)
    if path.exists():
        with path.open(encoding="utf-8") as handle:
            current = json.load(handle)
        if current != expected:
            raise StaleArtifactError(f"resume receipt differs: {path}")
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(expected, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()
    return path


def _write_parquet_once(frame: pd.DataFrame, path: Path) -> Path:
    """Atomically install one immutable donor Parquet file."""

    if path.exists():
        raise FileExistsError(f"immutable donor artifact already exists: {path}")
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


def _hardlink_once(source: Path, target: Path) -> Path:
    """Install an exact file without allocating another copy of its data."""

    source = source.resolve()
    if not source.is_file():
        raise PartialArtifactError(f"reuse source is missing: {source}")
    if target.exists():
        if target.is_file() and os.path.samefile(source, target):
            return target
        raise PartialArtifactError(
            f"reuse target exists but is not the source hard link: {target}"
        )
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        os.link(source, target)
    except OSError as exc:
        if exc.errno == errno.EXDEV:
            raise ContractError(
                f"{source} and {target.parent} are on different filesystems; "
                "verified cross-run reuse refuses to copy artifact data"
            ) from exc
        raise
    if not os.path.samefile(source, target):
        raise ContractError(f"hard-link verification failed: {target}")
    return target


def _source_run_root(context: RunContext, run_id: str) -> Path:
    if not RUN_ID_PATTERN.fullmatch(run_id):
        raise ContractError("--from-run must be a 20-character lowercase hex run ID")
    runs_root = context.base_config.runs_dir.resolve()
    source = (runs_root / run_id).resolve()
    if source.parent != runs_root or source.name != run_id:
        raise ContractError("--from-run resolves outside the configured runs directory")
    if source == context.run_root.resolve():
        raise ContractError("--from-run must differ from the current resolved run")
    if not source.is_dir():
        raise PartialArtifactError(f"source run is missing: {source}")
    return source


def _canonical_stage_config(manifest: StageManifest) -> dict[str, Any]:
    try:
        value = json.loads(manifest.config.canonical)
    except (TypeError, json.JSONDecodeError) as exc:
        raise ContractError(
            f"{manifest.stage}: configuration snapshot is not valid JSON"
        ) from exc
    if not isinstance(value, dict):
        raise ContractError(f"{manifest.stage}: configuration is not an object")
    return value


def _same_input_identity(left, right) -> bool:
    return (
        left.name,
        left.kind,
        left.content_sha256,
    ) == (
        right.name,
        right.kind,
        right.content_sha256,
    )


def _validate_ingest_reuse_contract(
    context: RunContext,
    source: StageManifest,
) -> tuple:
    if source.stage != "ingest":
        raise ContractError("source ingest manifest declares the wrong stage")
    if source.method_version != STAGE_METHODS["ingest"]:
        raise StaleArtifactError("source ingest method version differs")
    if source.expected_donors != context.donors:
        raise PartialArtifactError("source ingest donor contract differs")
    current_code = _stage_code(context, "ingest")
    if source.code.source_sha256 != current_code.source_sha256:
        raise StaleArtifactError(
            "ingest-producing code differs; source cells cannot be reused"
        )
    current_inputs = _stage_inputs(context, "ingest")
    if len(source.inputs) != len(current_inputs) or not all(
        _same_input_identity(old, new)
        for old, new in zip(source.inputs, current_inputs, strict=True)
    ):
        raise StaleArtifactError("source ingest inputs differ")
    source_config = _canonical_stage_config(source)
    try:
        source_min_area = source_config["configuration"]["cells_min_cell_area"]
    except (KeyError, TypeError) as exc:
        raise ContractError(
            "source ingest configuration lacks cells_min_cell_area"
        ) from exc
    if source_min_area != context.identity_config["cells_min_cell_area"]:
        raise StaleArtifactError("source ingest cell-area filter differs")
    return current_inputs


def _reuse_ingest(
    context: RunContext,
    source_root: Path,
) -> tuple[StageManifest, StageManifest]:
    """Reuse the source ingest with hard links and seal it for the target run."""

    source_path = source_root / "manifests" / "ingest.json"
    source = StageManifest.read_json(source_path)
    source.validate_current(validation_mode="fast")
    current_inputs = _validate_ingest_reuse_contract(context, source)

    target_path = context.stage_manifest_path("ingest")
    if target_path.exists():
        target = StageManifest.read_json(target_path)
        target.validate_current(
            config=_stage_config(context, "ingest"),
            validation_mode="fast",
        )
        if target.output.content_sha256 != source.output.content_sha256:
            raise StaleArtifactError("target ingest does not match the source ingest")
        print(f"[current] ingest: {target.output.total_rows:,} rows")
        return source, target

    source_output = Path(source.output.root).resolve()
    linked_files = []
    for donor in source.output.donors:
        for fingerprint in donor.files:
            source_file = Path(fingerprint.path).resolve()
            try:
                relative = source_file.relative_to(source_output)
            except ValueError as exc:
                raise ContractError(
                    f"source ingest file is outside its recorded root: {source_file}"
                ) from exc
            target_file = context.config.cells_dir / relative
            _hardlink_once(source_file, target_file)
            linked_files.append(
                {
                    "source": str(source_file),
                    "target": str(target_file.resolve()),
                    "size_bytes": source_file.stat().st_size,
                }
            )

    recorded_sidecars = {Path(item.path).resolve() for item in source.sidecars}
    linked_sidecars = []
    for relative in STAGE_SIDECARS["ingest"]:
        source_file = (source_root / relative).resolve()
        if source_file not in recorded_sidecars:
            raise PartialArtifactError(
                f"source ingest manifest does not own sidecar {source_file}"
            )
        target_file = context.run_root / relative
        _hardlink_once(source_file, target_file)
        linked_sidecars.append(
            {
                "source": str(source_file),
                "target": str(target_file.resolve()),
                "size_bytes": source_file.stat().st_size,
            }
        )

    receipt_path = context.config.manifests_dir / "ingest_reuse.json"
    receipt = _content_payload(
        {
            "resume_contract_version": RESUME_CONTRACT_VERSION,
            "stage": "ingest",
            "source_run_id": source_root.name,
            "source_manifest_content_id": source.content_id,
            "source_output_content_sha256": source.output.content_sha256,
            "target_run_id": context.run_id,
            "transfer": "verified_hard_links",
            "files": linked_files,
            "sidecars": linked_sidecars,
        }
    )
    _write_json_once(receipt_path, receipt)

    print(
        f"[reuse] ingest: hard-linked {len(linked_files)} donor files "
        f"from {source_root.name}"
    )
    target = StageManifest.create(
        stage="ingest",
        method_version=STAGE_METHODS["ingest"],
        expected_donors=context.donors,
        inputs=current_inputs,
        output_root=context.config.cells_dir,
        sidecar_paths=(
            *(context.run_root / item for item in STAGE_SIDECARS["ingest"]),
            receipt_path,
        ),
        config=_stage_config(context, "ingest"),
        code=_stage_code(context, "ingest"),
    )
    if target.output.content_sha256 != source.output.content_sha256:
        raise StaleArtifactError(
            "hard-linked ingest content differs from its source snapshot"
        )
    target_path.parent.mkdir(parents=True, exist_ok=True)
    target.write_json(target_path)
    print(
        f"[done] ingest: {target.output.total_rows:,} rows, "
        f"manifest={target.content_id[:12]}"
    )
    return source, target


def _validate_geometry_policy_transition(
    context: RunContext,
    source_ingest: StageManifest,
) -> None:
    source_config = _canonical_stage_config(source_ingest)
    try:
        source_geometry = dict(source_config["configuration"]["geometry_qc"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ContractError(
            "source run does not record a geometry-QC configuration"
        ) from exc
    target_geometry = dict(context.identity_config["geometry_qc"])
    source_policy = source_geometry.pop("duplicate_policy", None)
    target_policy = target_geometry.pop("duplicate_policy", None)
    if source_geometry != target_geometry:
        raise StaleArtifactError(
            "geometry configuration changed beyond duplicate_policy; "
            "the source prefix cannot be reused"
        )
    if (
        source_policy,
        target_policy,
    ) != (
        "fatal",
        "deterministic_keep",
    ):
        raise ContractError(
            "verified prefix reuse requires the policy transition "
            "fatal -> deterministic_keep"
        )


def _single_parquet(root: Path, donor_id: str) -> Path | None:
    partition = root / f"donor_id={donor_id}"
    if not partition.exists():
        return None
    files = sorted(partition.glob("*.parquet"))
    if len(files) != 1:
        raise PartialArtifactError(
            f"expected one Parquet file for donor {donor_id} in {root}; "
            f"found {len(files)}"
        )
    return files[0]


def _schema_sha256(schema) -> str:
    return hashlib.sha256(schema.serialize().to_pybytes()).hexdigest()


def _validate_geometry_partition(
    path: Path,
    *,
    expected: DonorDatasetSnapshot,
    require_no_duplicates: bool,
    expected_schema=None,
) -> tuple[GeometryPartitionCheck, Any]:
    parquet = pq.ParquetFile(path)
    schema = parquet.schema_arrow
    if tuple(schema.names) != GEOMETRY_REQUIRED_COLUMNS:
        raise SchemaContractError(
            f"{path}: geometry columns differ from the current contract"
        )
    if expected_schema is not None and not schema.equals(
        expected_schema, check_metadata=True
    ):
        raise SchemaContractError(
            f"{path}: geometry schema differs from preceding donors"
        )
    columns = [
        "donor_id",
        "object_id",
        "analysis_reasons",
        "spillover_context_reasons",
        "duplicate_group",
        "duplicate_survivor",
    ]
    table = parquet.read(columns=columns)
    if table.num_rows != expected.row_count:
        raise PartialArtifactError(
            f"{path}: row count {table.num_rows:,} differs from ingest "
            f"{expected.row_count:,}"
        )
    values = table.to_pydict()
    stored_donors = {
        str(value) for value in values["donor_id"] if value is not None
    }
    if stored_donors != {expected.donor_id}:
        raise ContractError(
            f"{path}: stored donor values differ from {expected.donor_id}"
        )
    object_hash = hash_identifiers(values["object_id"])
    if object_hash != expected.object_id_sha256:
        raise StaleArtifactError(
            f"{path}: object universe differs from the reused ingest"
        )

    duplicate_groups: dict[str, list[bool]] = {}
    for group, survivor in zip(
        values["duplicate_group"],
        values["duplicate_survivor"],
        strict=True,
    ):
        if (group is None) != (survivor is None):
            raise ContractError(
                f"{path}: duplicate group/survivor nullability is inconsistent"
            )
        if group is not None:
            duplicate_groups.setdefault(str(group), []).append(bool(survivor))
    for group, survivors in duplicate_groups.items():
        if len(survivors) < 2 or sum(survivors) != 1:
            raise ContractError(
                f"{path}: duplicate group {group!r} lacks exactly one survivor"
            )
    if require_no_duplicates and duplicate_groups:
        raise ContractError(
            f"{path}: source prefix contains duplicate geometry groups"
        )

    if require_no_duplicates:
        for analysis, spillover in zip(
            values["analysis_reasons"],
            values["spillover_context_reasons"],
            strict=True,
        ):
            analysis_text = "" if analysis is None else str(analysis)
            spillover_text = "" if spillover is None else str(spillover)
            for reason in (
                "missing_nucleus_geometry",
                "nucleus_raster_empty",
            ):
                if reason in analysis_text and reason not in spillover_text:
                    raise StaleArtifactError(
                        f"{path}: source prefix predates the current nucleus "
                        "exclusion semantics"
                    )

    result = GeometryPartitionCheck(
        donor_id=expected.donor_id,
        row_count=table.num_rows,
        object_id_sha256=object_hash,
        schema_sha256=_schema_sha256(schema),
        duplicate_group_count=len(duplicate_groups),
    )
    return result, schema


def _sampled_fingerprint(path: Path) -> dict[str, Any]:
    return asdict(FileFingerprint.capture(path, hash_mode="sampled"))


def _reuse_geometry_prefix(
    context: RunContext,
    source_root: Path,
    source_ingest: StageManifest,
    target_ingest: StageManifest,
    *,
    from_donor: str,
) -> StageManifest | None:
    manifest_path = context.stage_manifest_path("geometry")
    if manifest_path.exists():
        manifest = StageManifest.read_json(manifest_path)
        manifest.validate_current(
            config=_stage_config(context, "geometry"),
            validation_mode="fast",
        )
        print(f"[current] geometry: {manifest.output.total_rows:,} rows")
        return manifest

    try:
        start = context.donors.index(str(from_donor))
    except ValueError as exc:
        raise ContractError(
            f"--from-donor {from_donor!r} is not in the cohort"
        ) from exc
    prefix = context.donors[:start]
    _validate_geometry_policy_transition(context, source_ingest)

    expected_by_donor = {
        donor.donor_id: donor for donor in target_ingest.output.donors
    }
    output_root = context.config.geometry_qc_dir
    if output_root.exists():
        actual = {
            path.name.split("=", 1)[1]
            for path in output_root.glob("donor_id=*")
            if path.is_dir()
        }
        extra = sorted(actual.difference(context.donors))
        if extra:
            raise PartialArtifactError(
                f"geometry output contains donors outside the cohort: {extra}"
            )

    current_schema = None
    reused_records = []
    for donor in prefix:
        expected = expected_by_donor[donor]
        target_file = _single_parquet(output_root, donor)
        source_file = _single_parquet(
            source_root / "geometry_qc",
            donor,
        )
        if source_file is None:
            raise PartialArtifactError(
                f"source run has no completed geometry for donor {donor}"
            )
        source_check, _ = _validate_geometry_partition(
            source_file,
            expected=expected,
            require_no_duplicates=True,
        )
        if target_file is None:
            print(
                f"[reuse] geometry donor {donor}: "
                "normalizing completed duplicate-free result"
            )
            source_frame = pq.ParquetFile(source_file).read().to_pandas()
            target_file = (
                output_root / f"donor_id={donor}" / "data.parquet"
            )
            _write_parquet_once(
                normalize_geometry_output(source_frame),
                target_file,
            )
        target_check, target_schema = _validate_geometry_partition(
            target_file,
            expected=expected,
            require_no_duplicates=True,
            expected_schema=current_schema,
        )
        current_schema = target_schema
        reused_records.append(
            {
                "donor_id": donor,
                "source": _sampled_fingerprint(source_file),
                "target": _sampled_fingerprint(target_file),
                "source_validation": source_check.to_dict(),
                "target_validation": target_check.to_dict(),
                "transformation": (
                    "values_preserved; nullable duplicate columns "
                    "normalized to string/boolean"
                ),
            }
        )

    receipt_path = context.config.manifests_dir / "geometry_resume.json"
    receipt = _content_payload(
        {
            "resume_contract_version": RESUME_CONTRACT_VERSION,
            "stage": "geometry",
            "source_run_id": source_root.name,
            "source_ingest_manifest_content_id": source_ingest.content_id,
            "target_run_id": context.run_id,
            "target_ingest_manifest_content_id": target_ingest.content_id,
            "from_donor": str(from_donor),
            "reused_prefix_donors": list(prefix),
            "method_version": STAGE_METHODS["geometry"],
            "code_source_sha256": _stage_code(
                context, "geometry"
            ).source_sha256,
            "configuration": context.identity_config["geometry_qc"],
            "reused_partitions": reused_records,
        }
    )
    _write_json_once(receipt_path, receipt)

    from .donor_pipeline import create_receipt_for_existing_output

    for donor in prefix:
        create_receipt_for_existing_output(
            context,
            "geometry",
            donor,
            sidecar_paths=(receipt_path,),
        )
    return None


def resume(
    *,
    config_path: Path | None,
    source_run_id: str,
    from_donor: str,
    jobs: int | None = None,
    through: str | None = None,
    export: bool = True,
) -> RunContext:
    """Prepare verified reuse, finish geometry, and continue selected stages."""

    overrides = {}
    if jobs is not None:
        overrides["n_jobs"] = jobs
    context = resolve_run_context(load_config(config_path, **overrides))
    source_root = _source_run_root(context, source_run_id)
    print(f"[resume] source run {source_run_id}")
    print(f"[resume] target run {context.run_id}")
    source_ingest, target_ingest = _reuse_ingest(context, source_root)
    _reuse_geometry_prefix(
        context,
        source_root,
        source_ingest,
        target_ingest,
        from_donor=from_donor,
    )
    stages = _selected_stages(None, through)
    from .donor_pipeline import run_pipelined_pipeline

    run_pipelined_pipeline(
        context,
        stages=stages,
        export=export,
    )
    return context


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Reuse a verified ingest and duplicate-free geometry prefix, "
            "then continue a corrected content-addressed run."
        )
    )
    parser.add_argument("--config", type=Path)
    parser.add_argument("--from-run", required=True)
    parser.add_argument("--from-donor", required=True)
    parser.add_argument("--jobs", type=int)
    parser.add_argument("--through", choices=STAGE_ORDER[1:])
    parser.add_argument("--no-export", action="store_true")
    args = parser.parse_args(argv)
    resume(
        config_path=args.config,
        source_run_id=args.from_run,
        from_donor=args.from_donor,
        jobs=args.jobs,
        through=args.through,
        export=not args.no_export,
    )
    return 0


__all__ = [
    "GeometryPartitionCheck",
    "main",
    "resume",
]
