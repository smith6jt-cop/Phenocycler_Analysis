"""DAG-aware, provenance-safe adoption of compatible completed stages."""

from __future__ import annotations

import argparse
import hashlib
import json
import shlex
import subprocess
import sys
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Sequence

from .artifacts import (
    CodeMetadata,
    ConfigSnapshot,
    ContractError,
    FileFingerprint,
    PartialArtifactError,
    StageManifest,
    StaleArtifactError,
    sha256_file,
    sha256_json,
    validate_dataset_snapshot_current,
)
from .config import load_config
from .pipeline import (
    STAGE_DEPENDENCIES,
    STAGE_METHODS,
    STAGE_ORDER,
    STAGE_SIDECARS,
    STAGE_SOURCE_PATHS,
    RunContext,
    _stage_code,
    _stage_config,
    _stage_external_inputs,
    _stage_inputs,
    resolve_run_context,
)
from .recovery import _link_manifest_dataset
from .resume import (
    _canonical_stage_config,
    _content_payload,
    _hardlink_once,
    _source_run_root,
    _write_json_once as _write_immutable_json,
)


ADOPTION_CONTRACT_VERSION = 1
ADOPTION_POLICY_PATH = Path(__file__).with_name("adoption_compatibility.json")
ADOPT_ACTIONS = frozenset({"adopt_exact", "adopt_reviewed"})

# Only settings directly consumed by a stage are compared. Dependencies carry
# every upstream setting transitively; operational settings never enter run identity.
STAGE_CONFIGURATION_KEYS: dict[str, tuple[str, ...]] = {
    "ingest": ("run_schema_version", "cells_min_cell_area"),
    "geometry": ("run_schema_version", "geometry_qc"),
    "redsea": ("run_schema_version", "redsea"),
    "expression": ("run_schema_version",),
    "controls": ("run_schema_version", "reference_controls"),
    "calibrate": ("run_schema_version", "marker_calibration"),
    "type": ("run_schema_version", "typing_model"),
    "states": ("run_schema_version",),
}


def _write_adoption_json(path: Path, payload: Mapping[str, Any]) -> Path:
    try:
        return _write_immutable_json(path, payload)
    except StaleArtifactError as exc:
        raise StaleArtifactError(
            f"immutable adoption artifact differs: {path}"
        ) from exc


def _adoption_stage_code(context: RunContext, stage: str) -> CodeMetadata:
    """Bind both the accepted stage contract and the code that derived the target."""

    paths = tuple(
        dict.fromkeys(
            (
                *STAGE_SOURCE_PATHS[stage],
                "phenocycler/adoption.py",
                "phenocycler/adoption_compatibility.json",
            )
        )
    )
    return CodeMetadata.capture(context.repository, source_paths=paths)


def _missing_value() -> dict[str, bool]:
    return {"missing": True}


def project_stage_configuration(
    configuration: Mapping[str, Any],
    stage: str,
) -> dict[str, Any]:
    """Return the scientific configuration slice that can affect ``stage``."""

    result: dict[str, Any] = {}
    for key in STAGE_CONFIGURATION_KEYS[stage]:
        if key in configuration:
            result[key] = configuration[key]
        elif key == "typing_model":
            # Before bundles existed, absence meant the same rules-only lane.
            result[key] = {"mode": "rules_only"}
        else:
            result[key] = _missing_value()
    return result


def _upstream_incompatibilities(
    stage: str,
    decisions: Mapping[str, Mapping[str, Any]],
) -> list[str]:
    return [
        f"upstream_not_adoptable:{dependency}"
        for dependency in STAGE_DEPENDENCIES[stage]
        if decisions.get(dependency, {}).get("action") not in ADOPT_ACTIONS
    ]


def _typing_calibration_plan(
    context: RunContext,
    decisions: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, str | None] | None, str | None]:
    if context.typing_model_bundle is None:
        return None, None
    expected = (
        context.typing_model_bundle.evidence_provenance.source_artifact_content_id
    )
    calibrate = decisions.get("calibrate")
    planned = (
        calibrate.get("source_output_content_sha256")
        if calibrate is not None and calibrate.get("action") in ADOPT_ACTIONS
        else None
    )
    detail = {
        "bundle_evidence_content_sha256": expected,
        "planned_calibrate_content_sha256": planned,
    }
    reason = (
        "typing_bundle_calibration_mismatch"
        if planned is not None and planned != expected
        else None
    )
    return detail, reason


def _load_policy() -> tuple[dict[str, Any], str]:
    try:
        payload = json.loads(ADOPTION_POLICY_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(
            f"cannot load adoption compatibility policy: {ADOPTION_POLICY_PATH}"
        ) from exc
    if payload.get("policy_version") != 1:
        raise ContractError("unsupported adoption compatibility policy version")
    transitions = payload.get("transitions")
    if not isinstance(transitions, list):
        raise ContractError("adoption compatibility policy has no transition list")
    identifiers = [item.get("policy_id") for item in transitions if isinstance(item, dict)]
    if len(identifiers) != len(transitions) or len(set(identifiers)) != len(identifiers):
        raise ContractError("adoption compatibility policy IDs must be unique")
    return payload, sha256_json(payload)


def _git_blob(repository: Path, commit: str, relative: str) -> bytes:
    try:
        return subprocess.run(
            ["git", "-C", str(repository), "show", f"{commit}:{relative}"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as exc:
        raise ContractError(
            f"historical code object is unavailable: {commit}:{relative}"
        ) from exc


def _historical_code_sha256(
    repository: Path,
    commit: str,
    source_paths: Sequence[str],
) -> str:
    entries = []
    for relative in source_paths:
        data = _git_blob(repository, commit, relative)
        entries.append((relative, len(data), hashlib.sha256(data).hexdigest()))
    return sha256_json(tuple(sorted(entries)))


def _changed_source_files(
    repository: Path,
    commit: str,
    source_paths: Sequence[str],
    target_paths: Sequence[str],
) -> dict[str, tuple[str | None, str | None]]:
    source_set = set(source_paths)
    target_set = set(target_paths)
    changed: dict[str, tuple[str | None, str | None]] = {}
    for relative in sorted(source_set | target_set):
        historical = (
            hashlib.sha256(_git_blob(repository, commit, relative)).hexdigest()
            if relative in source_set
            else None
        )
        current_path = repository / relative
        if relative in target_set and not current_path.is_file():
            raise ContractError(f"current stage source is missing: {current_path}")
        current = sha256_file(current_path) if relative in target_set else None
        if historical != current:
            changed[relative] = (historical, current)
    return changed


def _reviewed_transition(
    *,
    context: RunContext,
    stage: str,
    source: StageManifest,
    current_code: CodeMetadata,
    policy: Mapping[str, Any],
) -> tuple[dict[str, Any] | None, str | None]:
    """Match and independently verify one exact checked-in code transition."""

    candidates = [
        item
        for item in policy["transitions"]
        if item.get("stage") == stage
        and item.get("source_method") == source.method_version
        and item.get("target_method") == STAGE_METHODS[stage]
        and item.get("source_code_sha256") == source.code.source_sha256
        and item.get("target_code_sha256") == current_code.source_sha256
    ]
    if not candidates:
        return None, "unreviewed_code_change"
    if len(candidates) != 1:
        raise ContractError(f"{stage}: multiple compatibility policies match")
    record = dict(candidates[0])
    try:
        if tuple(record["configuration_keys"]) != STAGE_CONFIGURATION_KEYS[stage]:
            raise ContractError("reviewed configuration scope differs")
        if source.output.schema_sha256 not in record["allowed_output_schema_sha256"]:
            raise ContractError("source output schema is not reviewed")
        policy_source_paths = tuple(
            map(str, record.get("source_paths", source.code.source_paths))
        )
        policy_target_paths = tuple(
            map(str, record.get("target_paths", current_code.source_paths))
        )
        if tuple(source.code.source_paths) != policy_source_paths:
            raise ContractError("source-path contract differs")
        if tuple(current_code.source_paths) != policy_target_paths:
            raise ContractError("target-path contract differs")
        if source.code.commit is None:
            raise ContractError("source manifest has no reconstructible Git commit")
        historical = _historical_code_sha256(
            context.repository,
            str(source.code.commit),
            policy_source_paths,
        )
        if historical != source.code.source_sha256:
            raise ContractError("source code hash is not reconstructible from its commit")
        changed_files = _changed_source_files(
            context.repository,
            str(source.code.commit),
            policy_source_paths,
            policy_target_paths,
        )
        changed_paths = tuple(changed_files)
        if changed_paths != tuple(sorted(map(str, record["changed_paths"]))):
            raise ContractError("changed source paths differ from the reviewed transition")
        reviewed_files = record.get("changed_files")
        if reviewed_files is None:
            if len(changed_paths) != 1:
                raise ContractError(
                    "legacy transition can review only one changed source file"
                )
            reviewed_files = [
                {
                    "path": changed_paths[0],
                    "source_sha256": record["changed_file_source_sha256"],
                    "target_sha256": record["changed_file_target_sha256"],
                }
            ]
        reviewed_by_path = {
            str(item["path"]): (
                item.get("source_sha256"),
                item.get("target_sha256"),
            )
            for item in reviewed_files
        }
        if reviewed_by_path != changed_files:
            raise ContractError("reviewed changed-file hashes differ")
        diff = subprocess.run(
            [
                "git",
                "-C",
                str(context.repository),
                "diff",
                f"{source.code.commit}..HEAD",
                "--",
                *changed_paths,
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
        if hashlib.sha256(diff).hexdigest() != record["reviewed_diff_sha256"]:
            raise ContractError("reviewed source diff differs")
    except (KeyError, TypeError, subprocess.CalledProcessError, ContractError) as exc:
        return None, f"compatibility_policy_failed:{exc}"
    return record, None


def _assert_within(path: Path, root: Path, *, label: str) -> None:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError as exc:
        raise ContractError(f"{label} resolves outside {root}") from exc


def _validate_archived_manifest(
    context: RunContext,
    source_root: Path,
    stage: str,
    manifest: StageManifest,
    *,
    validation_mode: str,
) -> None:
    """Validate archived bytes without comparing old code to this checkout."""

    if manifest.stage != stage:
        raise ContractError(f"source {stage} manifest declares {manifest.stage!r}")
    if manifest.expected_donors != context.donors:
        raise PartialArtifactError(f"source {stage} donor contract differs")
    target_relative = context.output_root(stage).resolve().relative_to(
        context.run_root.resolve()
    )
    expected_output = (source_root / target_relative).resolve()
    _assert_within(
        expected_output,
        source_root,
        label=f"source {stage} output root",
    )
    if Path(manifest.output.root).resolve() != expected_output:
        raise ContractError(f"source {stage} output root is noncanonical")
    for donor in manifest.output.donors:
        for fingerprint in donor.files:
            _assert_within(
                Path(fingerprint.path),
                expected_output,
                label=f"source {stage} output",
            )
    for sidecar in manifest.sidecars:
        _assert_within(
            Path(sidecar.path),
            source_root,
            label=f"source {stage} sidecar",
        )
    # This verifies identity, output inventory/content cache, and every sidecar.
    # Inputs are checked explicitly below because historical stage inputs would
    # otherwise recurse into current-checkout code validation.
    manifest.validate_current(
        validate_code=False,
        validate_inputs=False,
        validation_mode=validation_mode,
    )
    recorded = _canonical_stage_config(manifest)
    if recorded.get("stage") != stage:
        raise ContractError(f"source {stage} configuration declares another stage")
    if recorded.get("method_version") != manifest.method_version:
        raise ContractError(f"source {stage} configuration method differs")
    if recorded.get("run_id") != source_root.name:
        raise ContractError(f"source {stage} configuration run ID differs")


def _source_input_compatibility(
    *,
    context: RunContext,
    source_root: Path,
    stage: str,
    source: StageManifest,
    source_manifests: Mapping[str, StageManifest],
    validation_mode: str,
) -> tuple[list[dict[str, Any]], list[str]]:
    expected_external = {
        item.name: item for item in _stage_external_inputs(context, stage)
    }
    source_inputs = {item.name: item for item in source.inputs}
    reasons: list[str] = []
    expected_names = set(STAGE_DEPENDENCIES[stage]) | set(expected_external)
    if set(source_inputs) != expected_names:
        reasons.append("input_set_mismatch")

    details = []
    for dependency in STAGE_DEPENDENCIES[stage]:
        item = source_inputs.get(dependency)
        dependency_manifest = source_manifests.get(dependency)
        canonical_path = source_root / "manifests" / f"{dependency}.json"
        if (
            item is None
            or dependency_manifest is None
            or item.kind != "stage_manifest"
            or Path(item.path).resolve() != canonical_path.resolve()
            or item.content_sha256 != dependency_manifest.content_id
        ):
            reasons.append(f"source_dependency_invalid:{dependency}")
            continue
        details.append(
            {
                "name": dependency,
                "kind": "stage_manifest",
                "content_sha256": item.content_sha256,
                "output_content_sha256": dependency_manifest.output.content_sha256,
            }
        )

    for name, current in sorted(expected_external.items()):
        current.validate_current(mode=validation_mode)
        old = source_inputs.get(name)
        matches = old is not None and (
            old.kind,
            old.content_sha256,
        ) == (
            current.kind,
            current.content_sha256,
        )
        if not matches:
            reasons.append(f"external_input_mismatch:{name}")
        details.append(
            {
                "name": name,
                "kind": current.kind,
                "source_content_sha256": old.content_sha256 if old else None,
                "target_content_sha256": current.content_sha256,
                "compatible": matches,
            }
        )
    return details, reasons


def _stage_plan_record(
    *,
    context: RunContext,
    source_root: Path,
    stage: str,
    source: StageManifest | None,
    source_manifests: Mapping[str, StageManifest],
    source_error: str | None,
    decisions: Mapping[str, Mapping[str, Any]],
    policy: Mapping[str, Any],
    policy_content_id: str,
    validation_mode: str,
) -> dict[str, Any]:
    reasons: list[str] = []
    if source is None:
        reasons.append(source_error or "source_manifest_missing")
        return {
            "stage": stage,
            "action": "recompute",
            "reason_codes": reasons,
            "dependencies": list(STAGE_DEPENDENCIES[stage]),
            "source_manifest_content_id": None,
            "source_output_content_sha256": None,
            "source_bytes": 0,
        }

    reasons.extend(_upstream_incompatibilities(stage, decisions))

    input_details, input_reasons = _source_input_compatibility(
        context=context,
        source_root=source_root,
        stage=stage,
        source=source,
        source_manifests=source_manifests,
        validation_mode=validation_mode,
    )
    reasons.extend(input_reasons)

    source_config = _canonical_stage_config(source).get("configuration", {})
    if not isinstance(source_config, dict):
        reasons.append("source_configuration_invalid")
        source_config = {}
    source_projection = project_stage_configuration(source_config, stage)
    target_projection = project_stage_configuration(context.identity_config, stage)
    if source_projection != target_projection:
        reasons.append("configuration_mismatch")
    if source.method_version != STAGE_METHODS[stage]:
        reasons.append("method_version_mismatch")
    typing_calibration = None
    if stage == "type":
        typing_calibration, typing_reason = _typing_calibration_plan(
            context, decisions
        )
        if typing_reason is not None:
            reasons.append(typing_reason)

    current_code = _stage_code(context, stage)
    compatibility_policy = None
    if source.code.source_sha256 == current_code.source_sha256:
        code_basis = "exact_stage_code"
    else:
        compatibility_policy, policy_error = _reviewed_transition(
            context=context,
            stage=stage,
            source=source,
            current_code=current_code,
            policy=policy,
        )
        code_basis = (
            "checked_in_compatibility_policy"
            if compatibility_policy is not None
            else "unreviewed_stage_code"
        )
        if policy_error is not None:
            reasons.append(policy_error)

    source_bytes = sum(
        item.size_bytes for donor in source.output.donors for item in donor.files
    )
    action = "recompute"
    if not reasons:
        action = (
            "adopt_exact"
            if compatibility_policy is None
            else "adopt_reviewed"
        )
    return {
        "stage": stage,
        "action": action,
        "reason_codes": reasons,
        "dependencies": list(STAGE_DEPENDENCIES[stage]),
        "source_manifest_content_id": source.content_id,
        "source_manifest_file_sha256": sha256_file(
            source_root / "manifests" / f"{stage}.json"
        ),
        "source_output_content_sha256": source.output.content_sha256,
        "source_schema_sha256": source.output.schema_sha256,
        "source_object_id_sha256": source.output.object_id_sha256,
        "source_rows": source.output.total_rows,
        "source_bytes": source_bytes,
        "source_method": source.method_version,
        "target_method": STAGE_METHODS[stage],
        "source_code_sha256": source.code.source_sha256,
        "target_code_sha256": current_code.source_sha256,
        "source_code_commit": source.code.commit,
        "code_compatibility_basis": code_basis,
        "compatibility_policy_id": (
            compatibility_policy.get("policy_id")
            if compatibility_policy is not None
            else None
        ),
        "compatibility_policy_content_id": (
            policy_content_id if compatibility_policy is not None else None
        ),
        "configuration_keys": list(STAGE_CONFIGURATION_KEYS[stage]),
        "source_configuration": source_projection,
        "target_configuration": target_projection,
        "inputs": input_details,
        "typing_calibration": typing_calibration,
    }


def build_adoption_plan(
    context: RunContext,
    source_run_id: str,
    *,
    validation_mode: str = "fast",
) -> dict[str, Any]:
    """Plan maximal reuse without modifying either run directory."""

    if validation_mode not in {"fast", "content"}:
        raise ContractError("adoption validation mode must be 'fast' or 'content'")
    source_root = _source_run_root(context, source_run_id)
    policy, policy_content_id = _load_policy()
    source_manifests: dict[str, StageManifest] = {}
    source_errors: dict[str, str] = {}
    for stage in STAGE_ORDER:
        path = source_root / "manifests" / f"{stage}.json"
        if not path.is_file():
            source_errors[stage] = "source_manifest_missing"
            continue
        try:
            _assert_within(
                path,
                source_root,
                label=f"source {stage} manifest",
            )
            manifest = StageManifest.read_json(path)
            _validate_archived_manifest(
                context,
                source_root,
                stage,
                manifest,
                validation_mode=validation_mode,
            )
        except (ContractError, OSError, ValueError) as exc:
            source_errors[stage] = f"source_manifest_invalid:{exc}"
        else:
            source_manifests[stage] = manifest

    decisions: dict[str, dict[str, Any]] = {}
    records = []
    for stage in STAGE_ORDER:
        record = _stage_plan_record(
            context=context,
            source_root=source_root,
            stage=stage,
            source=source_manifests.get(stage),
            source_manifests=source_manifests,
            source_error=source_errors.get(stage),
            decisions=decisions,
            policy=policy,
            policy_content_id=policy_content_id,
            validation_mode=validation_mode,
        )
        decisions[stage] = record
        records.append(record)

    total_bytes = sum(item["source_bytes"] for item in records)
    reusable_bytes = sum(
        item["source_bytes"]
        for item in records
        if item["action"] in ADOPT_ACTIONS
    )
    payload = {
        "adoption_contract_version": ADOPTION_CONTRACT_VERSION,
        "source_run_id": source_run_id,
        "source_run_root": str(source_root),
        "source_has_run_manifest": (source_root / "run.json").is_file(),
        "source_complete_stage_set": len(source_manifests) == len(STAGE_ORDER),
        "target_run_id": context.run_id,
        "target_run_root": str(context.run_root.resolve()),
        "validation_mode": validation_mode,
        "compatibility_policy_path": str(ADOPTION_POLICY_PATH.resolve()),
        "compatibility_policy_content_id": policy_content_id,
        "stage_order": list(STAGE_ORDER),
        "stages": records,
        "summary": {
            "adopt_stages": [
                item["stage"] for item in records if item["action"] in ADOPT_ACTIONS
            ],
            "recompute_stages": [
                item["stage"] for item in records if item["action"] == "recompute"
            ],
            "source_stage_bytes": total_bytes,
            "reusable_stage_bytes": reusable_bytes,
            "reusable_fraction": reusable_bytes / total_bytes if total_bytes else 0.0,
            "continuation_blockers": [
                {
                    "stage": item["stage"],
                    "reason": reason,
                }
                for item in records
                for reason in item["reason_codes"]
                if reason == "typing_bundle_calibration_mismatch"
            ],
        },
    }
    return _content_payload(payload)


def verify_adoption_plan(plan: Mapping[str, Any]) -> None:
    if plan.get("adoption_contract_version") != ADOPTION_CONTRACT_VERSION:
        raise ContractError("unsupported adoption plan version")
    content_id = plan.get("content_id")
    if not isinstance(content_id, str):
        raise ContractError("adoption plan has no content ID")
    payload = dict(plan)
    payload.pop("content_id", None)
    if sha256_json(payload) != content_id:
        raise ContractError("adoption plan content ID is invalid")
    if tuple(plan.get("stage_order", ())) != STAGE_ORDER:
        raise ContractError("adoption plan stage order differs from this workflow")


def write_adoption_plan(plan: Mapping[str, Any], path: Path | str) -> Path:
    verify_adoption_plan(plan)
    return _write_adoption_json(Path(path).expanduser().resolve(), plan)


def read_adoption_plan(path: Path | str) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).expanduser().read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ContractError(f"cannot read adoption plan: {path}") from exc
    if not isinstance(value, dict):
        raise ContractError("adoption plan must be a JSON object")
    verify_adoption_plan(value)
    return value


def _target_sidecars(
    context: RunContext,
    source_root: Path,
    source: StageManifest,
    stage: str,
) -> tuple[list[Path], list[dict[str, Any]]]:
    owned = {Path(item.path).resolve(): item for item in source.sidecars}
    paths = []
    links = []
    for relative in STAGE_SIDECARS.get(stage, ()):
        source_path = (source_root / relative).resolve()
        fingerprint = owned.get(source_path)
        if fingerprint is None:
            raise PartialArtifactError(
                f"source {stage} manifest does not own required sidecar {source_path}"
            )
        target_path = context.run_root / relative
        _assert_within(
            target_path,
            context.run_root,
            label=f"target {stage} sidecar",
        )
        paths.append(target_path)
        links.append(
            {
                "source": str(source_path),
                "target": str(target_path.resolve()),
                "size_bytes": fingerprint.size_bytes,
                "sha256": fingerprint.sha256,
            }
        )
    return paths, links


def _relocated_output_snapshot(
    context: RunContext,
    source: StageManifest,
    stage: str,
):
    """Relocate a verified snapshot onto exact hard links without rereading Parquet."""

    source_root = Path(source.output.root).resolve()
    target_root = context.output_root(stage).resolve()
    donors = []
    for donor in source.output.donors:
        files = []
        for fingerprint in donor.files:
            source_path = Path(fingerprint.path).resolve()
            target_path = target_root / source_path.relative_to(source_root)
            if not target_path.is_file() or not target_path.samefile(source_path):
                raise StaleArtifactError(
                    f"{stage}: target is not the verified source hard link: {target_path}"
                )
            files.append(replace(fingerprint, path=str(target_path.resolve())))
        donors.append(replace(donor, files=tuple(files)))
    output = replace(
        source.output,
        root=str(target_root),
        donors=tuple(donors),
    )
    validate_dataset_snapshot_current(output, mode="fast")
    return output


def _reseal_adopted_stage(
    *,
    context: RunContext,
    stage: str,
    source_run_root: Path,
    source: StageManifest,
    inputs: Sequence[Any],
    scientific_sidecars: Sequence[Path],
    provenance_sidecars: Sequence[Path],
) -> StageManifest:
    """Create a target manifest from hard-link proof, not a full data rescan."""

    input_tuple = tuple(sorted(inputs, key=lambda item: item.name))
    if len({item.name for item in input_tuple}) != len(input_tuple):
        raise ContractError("input artifact names must be unique")
    for item in input_tuple:
        item.validate_current(mode="fast")
    output = _relocated_output_snapshot(context, source, stage)

    source_sidecars = {Path(item.path).resolve(): item for item in source.sidecars}
    relocated_sidecars = []
    for relative, target_path in zip(
        STAGE_SIDECARS.get(stage, ()),
        scientific_sidecars,
        strict=True,
    ):
        source_path = (source_run_root / relative).resolve()
        fingerprint = source_sidecars.get(source_path)
        if fingerprint is None or not target_path.samefile(source_path):
            raise StaleArtifactError(
                f"{stage}: scientific sidecar is not a verified source hard link"
            )
        relocated_sidecars.append(
            replace(fingerprint, path=str(target_path.resolve()))
        )
    sidecars = tuple(
        sorted(
            (
                *relocated_sidecars,
                *(FileFingerprint.capture(path) for path in provenance_sidecars),
            ),
            key=lambda item: item.path,
        )
    )
    if len({item.path for item in sidecars}) != len(sidecars):
        raise ContractError("stage sidecar paths must be unique")
    for sidecar in sidecars:
        sidecar.validate_current(mode="fast")

    draft = StageManifest(
        manifest_version=source.manifest_version,
        stage=stage,
        method_version=STAGE_METHODS[stage],
        expected_donors=tuple(sorted(context.donors)),
        completed_donors=output.donor_ids,
        inputs=input_tuple,
        output=output,
        sidecars=sidecars,
        config=ConfigSnapshot.capture(_stage_config(context, stage)),
        # This is the current compatibility contract used by normal status/run
        # validation. The fingerprinted adoption receipt separately preserves
        # both the original producer and the target derivation code.
        code=_stage_code(context, stage),
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        content_id="",
    )
    target = replace(draft, content_id=sha256_json(draft._identity_payload()))
    target.verify_identity()
    return target


def _preflight_target_stage(
    context: RunContext,
    source_root: Path,
    source: StageManifest,
    stage: str,
) -> None:
    target_root = context.output_root(stage)
    _assert_within(target_root, context.run_root, label=f"target {stage} output")
    expected: dict[Path, Path] = {}
    source_output = Path(source.output.root).resolve()
    for donor in source.output.donors:
        for fingerprint in donor.files:
            source_path = Path(fingerprint.path).resolve()
            relative = source_path.relative_to(source_output)
            expected[(target_root / relative).resolve()] = source_path
    if target_root.exists():
        actual = {path.resolve() for path in target_root.rglob("*") if path.is_file()}
        foreign = actual.difference(expected)
        if foreign:
            raise PartialArtifactError(
                f"{stage}: adoption target contains foreign files: {sorted(map(str, foreign))[:3]}"
            )
        for target in actual:
            if not target.is_file() or not target.samefile(expected[target]):
                raise PartialArtifactError(
                    f"{stage}: existing target is not the reviewed source hard link: {target}"
                )
    sidecars, _ = _target_sidecars(context, source_root, source, stage)
    for relative, target in zip(STAGE_SIDECARS.get(stage, ()), sidecars, strict=True):
        source_path = source_root / relative
        if target.exists() and not (
            target.is_file() and target.samefile(source_path)
        ):
            raise PartialArtifactError(
                f"{stage}: existing adopted sidecar is not the source hard link: {target}"
            )


def _assert_source_matches_plan(
    path: Path,
    manifest: StageManifest,
    record: Mapping[str, Any],
) -> None:
    expected = (
        record.get("source_manifest_content_id"),
        record.get("source_manifest_file_sha256"),
        record.get("source_output_content_sha256"),
        record.get("source_schema_sha256"),
        record.get("source_object_id_sha256"),
    )
    actual = (
        manifest.content_id,
        sha256_file(path),
        manifest.output.content_sha256,
        manifest.output.schema_sha256,
        manifest.output.object_id_sha256,
    )
    if actual != expected:
        raise StaleArtifactError(
            f"source {manifest.stage} manifest differs from the reviewed plan"
        )


def _receipt_payload(
    *,
    context: RunContext,
    plan: Mapping[str, Any],
    record: Mapping[str, Any],
    source: StageManifest,
    dataset_links: Sequence[Mapping[str, Any]],
    sidecar_links: Sequence[Mapping[str, Any]],
    source_snapshot: Path,
    target_inputs: Sequence[Any],
) -> dict[str, Any]:
    adoption_code = CodeMetadata.capture(
        context.repository,
        source_paths=(
            "phenocycler/adoption.py",
            "phenocycler/adoption_compatibility.json",
        ),
    )
    return _content_payload(
        {
            "adoption_contract_version": ADOPTION_CONTRACT_VERSION,
            "plan_content_id": plan["content_id"],
            "stage": record["stage"],
            "decision": record["action"],
            "validation_mode": plan["validation_mode"],
            "source_run_id": plan["source_run_id"],
            "source_manifest_content_id": source.content_id,
            "source_manifest_snapshot": str(source_snapshot.resolve()),
            "source_manifest_file_sha256": record["source_manifest_file_sha256"],
            "source_output_content_sha256": source.output.content_sha256,
            "source_schema_sha256": source.output.schema_sha256,
            "source_object_id_sha256": source.output.object_id_sha256,
            "source_method": source.method_version,
            "source_code": asdict(source.code),
            "target_run_id": context.run_id,
            "target_method": STAGE_METHODS[record["stage"]],
            "target_stage_contract_code": asdict(
                _stage_code(context, record["stage"])
            ),
            "target_manifest_derivation_code": asdict(
                _adoption_stage_code(context, record["stage"])
            ),
            "target_inputs": [
                {
                    "name": item.name,
                    "kind": item.kind,
                    "content_sha256": item.content_sha256,
                }
                for item in target_inputs
            ],
            "configuration_keys": record["configuration_keys"],
            "source_configuration": record["source_configuration"],
            "target_configuration": record["target_configuration"],
            "compatibility_policy_id": record["compatibility_policy_id"],
            "compatibility_policy_content_id": record[
                "compatibility_policy_content_id"
            ],
            "transfer": "verified_hard_links",
            "dataset_files": list(dataset_links),
            "scientific_sidecars": list(sidecar_links),
            "adoption_code": asdict(adoption_code),
        }
    )


def _adopt_stage(
    context: RunContext,
    plan: Mapping[str, Any],
    record: Mapping[str, Any],
) -> StageManifest:
    stage = str(record["stage"])
    source_root = _source_run_root(context, str(plan["source_run_id"]))
    source_path = source_root / "manifests" / f"{stage}.json"
    source = StageManifest.read_json(source_path)
    _assert_source_matches_plan(source_path, source, record)
    target_manifest_path = context.stage_manifest_path(stage)
    receipt_path = context.config.manifests_dir / "adoption" / f"{stage}.json"
    source_snapshot = (
        context.config.manifests_dir
        / "adoption"
        / "source"
        / source_root.name
        / f"{stage}.json"
    )
    plan_snapshot = context.config.manifests_dir / "adoption" / "plan.json"
    for path, label in (
        (target_manifest_path, f"target {stage} manifest"),
        (receipt_path, f"target {stage} adoption receipt"),
        (source_snapshot, f"target {stage} source snapshot"),
        (plan_snapshot, "target adoption plan snapshot"),
    ):
        _assert_within(path, context.run_root, label=label)

    if target_manifest_path.exists():
        current = StageManifest.read_json(target_manifest_path)
        current.validate_current(
            config=_stage_config(context, stage),
            validation_mode="fast",
        )
        if current.output.content_sha256 != source.output.content_sha256:
            raise StaleArtifactError(
                f"existing target {stage} output differs from the adoption source"
            )
        sidecar_paths = {Path(item.path).resolve() for item in current.sidecars}
        if receipt_path.resolve() not in sidecar_paths:
            raise StaleArtifactError(
                f"existing target {stage} is not sealed by its adoption receipt"
            )
        print(f"[adopt-current] {stage}: {current.output.total_rows:,} rows")
        return current

    _hardlink_once(source_path, source_snapshot)
    source = StageManifest.read_json(source_snapshot)
    _assert_source_matches_plan(source_snapshot, source, record)
    dataset_links = _link_manifest_dataset(context, stage=stage, source=source)
    target_sidecars, sidecar_links = _target_sidecars(
        context, source_root, source, stage
    )
    for link in sidecar_links:
        _hardlink_once(Path(link["source"]), Path(link["target"]))
    target_inputs = _stage_inputs(context, stage)
    for item in target_inputs:
        item.validate_current(mode="fast")
    receipt = _receipt_payload(
        context=context,
        plan=plan,
        record=record,
        source=source,
        dataset_links=dataset_links,
        sidecar_links=sidecar_links,
        source_snapshot=source_snapshot,
        target_inputs=target_inputs,
    )
    _write_adoption_json(receipt_path, receipt)
    target = _reseal_adopted_stage(
        context=context,
        stage=stage,
        source_run_root=source_root,
        inputs=target_inputs,
        source=source,
        scientific_sidecars=target_sidecars,
        provenance_sidecars=(
            plan_snapshot,
            source_snapshot,
            receipt_path,
        ),
    )
    if target.output.content_sha256 != source.output.content_sha256:
        raise StaleArtifactError(
            f"adopted {stage} content differs from its source snapshot"
        )
    target_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    target.write_json(target_manifest_path)
    print(
        f"[adopted] {stage}: {target.output.total_rows:,} rows, "
        f"manifest={target.content_id[:12]}"
    )
    return target


def apply_adoption_plan(
    context: RunContext,
    plan: Mapping[str, Any],
) -> tuple[str, ...]:
    """Revalidate and idempotently seal compatible stages one at a time."""

    verify_adoption_plan(plan)
    if plan.get("target_run_id") != context.run_id:
        raise StaleArtifactError(
            "adoption plan target differs from the currently resolved run"
        )
    rebuilt = build_adoption_plan(
        context,
        str(plan["source_run_id"]),
        validation_mode=str(plan["validation_mode"]),
    )
    if rebuilt["content_id"] != plan["content_id"]:
        raise StaleArtifactError(
            "adoption plan is stale; rebuild and review it before applying"
        )

    source_root = _source_run_root(context, str(plan["source_run_id"]))
    _assert_within(
        context.config.manifests_dir,
        context.run_root,
        label="target manifests directory",
    )
    records = [
        item for item in plan["stages"] if item["action"] in ADOPT_ACTIONS
    ]
    # Validate every source and target inventory before the first target write.
    for record in records:
        stage = str(record["stage"])
        source = StageManifest.read_json(
            source_root / "manifests" / f"{stage}.json"
        )
        _assert_source_matches_plan(
            source_root / "manifests" / f"{stage}.json",
            source,
            record,
        )
        _preflight_target_stage(context, source_root, source, stage)

    plan_snapshot = context.config.manifests_dir / "adoption" / "plan.json"
    _write_adoption_json(plan_snapshot, plan)
    adopted_manifests = []
    for record in records:
        adopted_manifests.append(_adopt_stage(context, plan, record))

    aggregate = _content_payload(
        {
            "adoption_contract_version": ADOPTION_CONTRACT_VERSION,
            "plan_content_id": plan["content_id"],
            "source_run_id": plan["source_run_id"],
            "target_run_id": context.run_id,
            "transfer": "verified_hard_links",
            "adopted_stages": [
                {
                    "stage": manifest.stage,
                    "target_manifest_content_id": manifest.content_id,
                    "output_content_sha256": manifest.output.content_sha256,
                }
                for manifest in adopted_manifests
            ],
            "recompute_stages": plan["summary"]["recompute_stages"],
        }
    )
    _write_adoption_json(context.config.manifests_dir / "adoption.json", aggregate)
    return tuple(manifest.stage for manifest in adopted_manifests)


def _print_summary(plan: Mapping[str, Any], *, file=sys.stdout) -> None:
    print(
        f"source {plan['source_run_id']} -> target {plan['target_run_id']}",
        file=file,
    )
    for item in plan["stages"]:
        reason = ", ".join(item["reason_codes"]) or item.get(
            "code_compatibility_basis", "compatible"
        )
        print(f"  {item['stage']:11s} {item['action']:15s} {reason}", file=file)
    summary = plan["summary"]
    gib = summary["reusable_stage_bytes"] / (1024**3)
    print(
        f"reuse {gib:.2f} GiB ({summary['reusable_fraction']:.1%}); "
        f"recompute={summary['recompute_stages']}",
        file=file,
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plan or apply provenance-safe reuse of compatible completed stages."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    plan_parser = subparsers.add_parser(
        "plan",
        help="validate a source run and report compatible stages without changing it",
    )
    plan_parser.add_argument("--config", type=Path, help="pipeline INI file")
    plan_parser.add_argument(
        "--from-run",
        required=True,
        help="20-character source run ID under the configured runs directory",
    )
    plan_parser.add_argument(
        "--out",
        type=Path,
        help="write the immutable reviewed plan here; otherwise emit JSON to stdout",
    )
    plan_parser.add_argument(
        "--validation-mode",
        choices=("fast", "content"),
        default="fast",
        help="use cached fingerprints (fast) or rehash all recorded bytes (content)",
    )
    apply_parser = subparsers.add_parser(
        "apply",
        help="revalidate and hard-link the exact stages authorized by a reviewed plan",
    )
    apply_parser.add_argument("--config", type=Path, help="same pipeline INI used to plan")
    apply_parser.add_argument(
        "--plan", required=True, type=Path, help="immutable plan JSON produced by adopt plan"
    )
    apply_parser.add_argument(
        "--continue",
        dest="continue_run",
        action="store_true",
        help="run incompatible stages and QuPath export after adoption",
    )
    apply_parser.add_argument(
        "--jobs", type=int, help="total donor tasks when --continue is used"
    )
    apply_parser.add_argument(
        "--no-export",
        action="store_true",
        help="with --continue, finish analytical stages but defer QuPath export",
    )
    args = parser.parse_args(argv)

    if args.command == "apply" and not args.continue_run:
        if args.jobs is not None or args.no_export:
            parser.error("--jobs and --no-export require --continue")

    overrides = {}
    if getattr(args, "jobs", None) is not None:
        overrides["n_jobs"] = args.jobs
    context = resolve_run_context(load_config(args.config, **overrides))
    if args.command == "plan":
        plan = build_adoption_plan(
            context,
            args.from_run,
            validation_mode=args.validation_mode,
        )
        if args.out is None:
            _print_summary(plan, file=sys.stderr)
            print(json.dumps(plan, indent=2, sort_keys=True))
        else:
            output_path = args.out.expanduser().resolve()
            for run_root in (
                Path(plan["source_run_root"]),
                Path(plan["target_run_root"]),
            ):
                try:
                    output_path.relative_to(run_root.resolve())
                except ValueError:
                    continue
                parser.error(
                    "--out must be outside the source and target run directories"
                )
            write_adoption_plan(plan, output_path)
            _print_summary(plan)
            print(f"plan {plan['content_id']} written to {output_path}")
        return 0

    plan = read_adoption_plan(args.plan)
    if args.continue_run and plan["summary"].get("continuation_blockers"):
        raise ContractError(
            "adoption plan has continuation blockers: "
            f"{plan['summary']['continuation_blockers']}"
        )
    adopted = apply_adoption_plan(context, plan)
    print(f"[adoption-ready] adopted={list(adopted)}, target={context.run_id}")
    if args.continue_run:
        from .donor_pipeline import run_pipelined_pipeline

        run_pipelined_pipeline(
            context,
            stages=STAGE_ORDER,
            export=not args.no_export,
        )
    else:
        config_option = (
            f" --config {shlex.quote(str(args.config))}" if args.config else ""
        )
        print(
            "continue with: python -m phenocycler.pipeline run"
            f"{config_option} --pipelined"
        )
    return 0


__all__ = [
    "ADOPTION_CONTRACT_VERSION",
    "ADOPTION_POLICY_PATH",
    "STAGE_CONFIGURATION_KEYS",
    "apply_adoption_plan",
    "build_adoption_plan",
    "main",
    "project_stage_configuration",
    "read_adoption_plan",
    "verify_adoption_plan",
    "write_adoption_plan",
]
