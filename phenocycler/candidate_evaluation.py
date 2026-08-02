"""Evaluate an unpromoted probabilistic typing bundle without production use.

This module is the deliberately separate bridge between fitting and promotion.
It may generate assignments for the donor-disjoint calibration or validation
split, or for the locked-test split when that access is explicitly
acknowledged.  It does not accept or create a production promotion manifest.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import pandas as pd

from .artifacts import (
    DatasetSnapshot,
    DonorDatasetSnapshot,
    InputArtifact,
    RunManifest,
    StageManifest,
    StaleArtifactError,
    hash_identifiers,
    sha256_json,
    snapshot_parquet_dataset,
    validate_dataset_snapshot_current,
)
from .hierarchical_typing import AUTHORITATIVE_ANCHOR, INFERRED, TypingRegistry
from .marker_registry import MarkerRegistry, load_registry
from .refinement_contracts import validate_evidence_source
from .typing_model_bundle import (
    TypingModelBundle,
    load_typing_model_bundle,
)
from .typing_stage import (
    CANDIDATE_TYPING_MODE,
    RULE_SCORE_SEMANTICS,
    type_candidate_cells,
)

CANDIDATE_EVALUATION_SCHEMA_VERSION = 1
CANDIDATE_EVALUATION_METHOD_VERSION = "candidate-typing-evaluation-v2"
EVALUATION_SPLITS = ("calibration", "validation", "locked_test")

_REQUIRED_ASSIGNMENT_STRING_COLUMNS = frozenset(
    {
        "donor_id",
        "object_id",
        "broad_label",
        "broad_assignment_status",
        "broad_best_candidate",
        "broad_inference_failure_reasons",
        "broad_reason",
        "broad_authoritative_markers",
        "broad_anchor_classes",
        "ranked_broad_probabilities",
        "broad_type",
        "best_broad_type",
        "assignment_status",
        "ranked_probabilities",
        "subtype_label",
        "subtype_assignment_status",
        "subtype_best_candidate",
        "subtype_inference_failure_reasons",
        "subtype_reason",
        "ranked_subtype_probabilities",
        "cell_type",
        "specific_type",
        "best_specific_type",
        "ranked_specific_probabilities",
        "specific_assignment_status",
        "typing_rules_version",
        "typing_rules_fingerprint",
        "marker_registry_fingerprint",
        "typing_mode",
        "typing_model_bundle_content_id",
        "typing_model_bundle_version",
        "typing_model_method_version",
        "typing_feature_schema_version",
        "typing_score_semantics",
        "typing_broad_score_semantics",
        "typing_subtype_score_semantics",
        "typing_subtype_model_source",
        "candidate_assignment_run_id",
        "candidate_evaluation_split_name",
    }
)
_REQUIRED_ASSIGNMENT_BOOLEAN_COLUMNS = frozenset(
    {
        "broad_unique_winner_pass",
        "broad_probability_pass",
        "broad_margin_pass",
        "broad_stability_pass",
        "broad_valid_evidence_pass",
        "broad_threshold_eligible",
        "broad_all_required_valid",
        "broad_all_required_negative",
        "subtype_unique_winner_pass",
        "subtype_probability_pass",
        "subtype_margin_pass",
        "subtype_stability_pass",
        "subtype_valid_evidence_pass",
        "subtype_threshold_eligible",
        "qc_analysis_eligible",
    }
)
_REQUIRED_ASSIGNMENT_NUMERIC_COLUMNS = frozenset(
    {
        "broad_best_probability",
        "broad_margin",
        "broad_stability",
        "confidence",
        "subtype_best_probability",
        "subtype_margin",
        "subtype_stability",
    }
)
_REQUIRED_ASSIGNMENT_INTEGER_COLUMNS = frozenset(
    {
        "typing_broad_bootstrap_replicates",
        "typing_broad_valid_bootstrap_replicates",
        "typing_subtype_bootstrap_replicates",
        "typing_subtype_valid_bootstrap_replicates",
    }
)
_REQUIRED_ASSIGNMENT_COLUMNS = frozenset(
    _REQUIRED_ASSIGNMENT_STRING_COLUMNS
    | _REQUIRED_ASSIGNMENT_BOOLEAN_COLUMNS
    | _REQUIRED_ASSIGNMENT_NUMERIC_COLUMNS
    | _REQUIRED_ASSIGNMENT_INTEGER_COLUMNS
)

_MANIFEST_FIELDS = {
    "schema_version",
    "method_version",
    "assignment_run_id",
    "model_bundle_content_id",
    "split_manifest_content_id",
    "split_name",
    "source_run_id",
    "source_run_manifest_content_id",
    "source_stage_manifest_content_id",
    "source_artifact_content_id",
    "marker_registry_fingerprint",
    "typing_rules_fingerprint",
    "donors",
    "output",
    "output_content_sha256",
    "content_id",
}
_DATASET_FIELDS = {
    "root",
    "donors",
    "total_rows",
    "schema",
    "schema_sha256",
    "object_id_sha256",
    "content_sha256",
}
_DONOR_SNAPSHOT_FIELDS = {
    "donor_id",
    "files",
    "row_count",
    "schema",
    "schema_sha256",
    "object_id_sha256",
    "content_sha256",
}
_FILE_FINGERPRINT_FIELDS = {
    "path",
    "size_bytes",
    "mtime_ns",
    "sha256",
    "hash_mode",
    "sample_size_bytes",
    "sample_count",
}
_SCHEMA_FIELD_FIELDS = {"name", "type", "nullable"}


def _strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON number is forbidden: {value}")


def _strict_json_float(value: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"non-finite JSON number is forbidden: {value}")
    return result


def _read_strict_json(handle: Any) -> Any:
    return json.load(
        handle,
        object_pairs_hook=_strict_json_object,
        parse_constant=_reject_json_constant,
        parse_float=_strict_json_float,
    )


def _require_mapping(value: Any, *, where: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{where} must be a JSON object")
    return value


def _require_exact_fields(value: Mapping[str, Any], fields: set[str], *, where: str) -> None:
    unknown = sorted(set(value) - fields)
    missing = sorted(fields - set(value))
    if unknown:
        raise ValueError(f"{where} contains unknown fields: {unknown}")
    if missing:
        raise ValueError(f"{where} is missing required fields: {missing}")


def _require_json_int(value: Any, *, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be a JSON integer")
    if value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _require_sha256(value: Any, *, name: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _normalise_donors(values: Sequence[Any]) -> tuple[str, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError("candidate evaluation donors must be a JSON array")
    donors: list[str] = []
    for raw in values:
        if not isinstance(raw, str) or not raw.strip():
            raise TypeError("candidate evaluation donors must be non-empty strings")
        donor = raw.strip()
        if Path(donor).name != donor or "/" in donor or "\\" in donor:
            raise ValueError(f"unsafe donor partition identifier: {donor!r}")
        donors.append(donor)
    if not donors or len(set(donors)) != len(donors):
        raise ValueError("candidate evaluation donors must be non-empty and unique")
    return tuple(sorted(donors))


def _validate_snapshot_json(value: Any) -> Mapping[str, Any]:
    """Reject permissive/unknown nested fields before artifact deserialisation."""

    dataset = _require_mapping(value, where="candidate output snapshot")
    _require_exact_fields(dataset, _DATASET_FIELDS, where="candidate output snapshot")
    _require_json_int(dataset["total_rows"], name="output.total_rows", minimum=1)
    if not isinstance(dataset["root"], str) or not Path(dataset["root"]).is_absolute():
        raise ValueError("output.root must be an absolute path")
    for name in ("schema_sha256", "object_id_sha256", "content_sha256"):
        _require_sha256(dataset[name], name=f"output.{name}")
    donors = dataset["donors"]
    schema = dataset["schema"]
    if isinstance(donors, (str, bytes)) or not isinstance(donors, Sequence):
        raise TypeError("output.donors must be a JSON array")
    if isinstance(schema, (str, bytes)) or not isinstance(schema, Sequence):
        raise TypeError("output.schema must be a JSON array")
    for index, raw_field in enumerate(schema):
        field_value = _require_mapping(raw_field, where=f"output.schema[{index}]")
        _require_exact_fields(
            field_value,
            _SCHEMA_FIELD_FIELDS,
            where=f"output.schema[{index}]",
        )
        if not isinstance(field_value["nullable"], bool):
            raise TypeError(f"output.schema[{index}].nullable must be a JSON boolean")
        if not isinstance(field_value["name"], str) or not isinstance(field_value["type"], str):
            raise TypeError(f"output.schema[{index}] name/type must be strings")
    for donor_index, raw_donor in enumerate(donors):
        donor = _require_mapping(raw_donor, where=f"output.donors[{donor_index}]")
        _require_exact_fields(
            donor,
            _DONOR_SNAPSHOT_FIELDS,
            where=f"output.donors[{donor_index}]",
        )
        _require_json_int(
            donor["row_count"],
            name=f"output.donors[{donor_index}].row_count",
            minimum=1,
        )
        if not isinstance(donor["donor_id"], str) or not donor["donor_id"].strip():
            raise TypeError(f"output.donors[{donor_index}].donor_id must be a string")
        for name in ("schema_sha256", "object_id_sha256", "content_sha256"):
            _require_sha256(donor[name], name=f"output.donors[{donor_index}].{name}")
        donor_schema = donor["schema"]
        files = donor["files"]
        if isinstance(donor_schema, (str, bytes)) or not isinstance(donor_schema, Sequence):
            raise TypeError(f"output.donors[{donor_index}].schema must be an array")
        if isinstance(files, (str, bytes)) or not isinstance(files, Sequence):
            raise TypeError(f"output.donors[{donor_index}].files must be an array")
        if not files:
            raise ValueError(f"output.donors[{donor_index}].files must not be empty")
        for field_index, raw_field in enumerate(donor_schema):
            field_value = _require_mapping(
                raw_field,
                where=f"output.donors[{donor_index}].schema[{field_index}]",
            )
            _require_exact_fields(
                field_value,
                _SCHEMA_FIELD_FIELDS,
                where=f"output.donors[{donor_index}].schema[{field_index}]",
            )
            if not isinstance(field_value["nullable"], bool):
                raise TypeError("output donor schema nullable values must be JSON booleans")
            if not isinstance(field_value["name"], str) or not isinstance(field_value["type"], str):
                raise TypeError("output donor schema name/type must be strings")
        for file_index, raw_file in enumerate(files):
            file_value = _require_mapping(
                raw_file,
                where=f"output.donors[{donor_index}].files[{file_index}]",
            )
            _require_exact_fields(
                file_value,
                _FILE_FINGERPRINT_FIELDS,
                where=f"output.donors[{donor_index}].files[{file_index}]",
            )
            for name in (
                "size_bytes",
                "mtime_ns",
                "sample_size_bytes",
                "sample_count",
            ):
                _require_json_int(
                    file_value[name],
                    name=(f"output.donors[{donor_index}].files[{file_index}].{name}"),
                )
            if (
                not isinstance(file_value["path"], str)
                or not Path(file_value["path"]).is_absolute()
            ):
                raise ValueError("candidate output file paths must be absolute")
            _require_sha256(
                file_value["sha256"],
                name=f"output.donors[{donor_index}].files[{file_index}].sha256",
            )
            if file_value["hash_mode"] not in {"full", "sampled"}:
                raise ValueError("candidate output file hash_mode is unsupported")
    return dataset


def _snapshot_semantic_identity(snapshot: DatasetSnapshot) -> dict[str, Any]:
    """Path-independent fields recomputed by a content validation pass."""

    return {
        "donors": [
            {
                "donor_id": donor.donor_id,
                "row_count": donor.row_count,
                "schema_sha256": donor.schema_sha256,
                "object_id_sha256": donor.object_id_sha256,
                "content_sha256": donor.content_sha256,
                "files": [
                    {
                        "name": Path(item.path).name,
                        "size_bytes": item.size_bytes,
                        "sha256": item.sha256,
                        "hash_mode": item.hash_mode,
                    }
                    for item in donor.files
                ],
            }
            for donor in snapshot.donors
        ],
        "total_rows": snapshot.total_rows,
        "schema_sha256": snapshot.schema_sha256,
        "object_id_sha256": snapshot.object_id_sha256,
        "content_sha256": snapshot.content_sha256,
    }


def _require_constant_assignment_column(
    frame: pd.DataFrame,
    column: str,
    expected: str | int,
    *,
    donor: str,
) -> None:
    values = frame[column]
    if values.isna().any() or not values.astype(str).eq(str(expected)).all():
        raise ValueError(f"candidate assignments for donor {donor} have inconsistent {column}")


def _validate_assignment_frame(
    frame: pd.DataFrame,
    *,
    donor_snapshot: DonorDatasetSnapshot,
    assignment_run_id: str,
    split_name: str,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
) -> None:
    """Prove one assignment partition preserves cells and typing semantics."""

    donor = donor_snapshot.donor_id
    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"candidate assignments for donor {donor} must be a DataFrame")
    if frame.columns.has_duplicates:
        raise ValueError(f"candidate assignments for donor {donor} have duplicate columns")
    broad_probability_columns = {
        f"broad_probability__{class_name}" for class_name in model_bundle.broad.classes
    }
    required = _REQUIRED_ASSIGNMENT_COLUMNS | broad_probability_columns
    missing = sorted(required - set(frame.columns))
    if missing:
        raise KeyError(
            f"candidate assignments for donor {donor} are missing required columns: {missing}"
        )

    for column in _REQUIRED_ASSIGNMENT_STRING_COLUMNS:
        series = frame[column]
        if not (
            pd.api.types.is_object_dtype(series.dtype) or pd.api.types.is_string_dtype(series.dtype)
        ):
            raise TypeError(f"candidate assignment column {column!r} must contain strings")
        if series.isna().any():
            raise ValueError(f"candidate assignment column {column!r} contains nulls")
    for column in _REQUIRED_ASSIGNMENT_BOOLEAN_COLUMNS:
        series = frame[column]
        if not pd.api.types.is_bool_dtype(series.dtype) or series.isna().any():
            raise TypeError(
                f"candidate assignment column {column!r} must contain non-null booleans"
            )
    for column in _REQUIRED_ASSIGNMENT_NUMERIC_COLUMNS | broad_probability_columns:
        series = frame[column]
        if pd.api.types.is_bool_dtype(series.dtype) or not pd.api.types.is_numeric_dtype(
            series.dtype
        ):
            raise TypeError(f"candidate assignment column {column!r} must be numeric")
    for column in _REQUIRED_ASSIGNMENT_INTEGER_COLUMNS:
        series = frame[column]
        if not pd.api.types.is_integer_dtype(series.dtype) or series.isna().any():
            raise TypeError(
                f"candidate assignment column {column!r} must contain non-null integers"
            )
        if (series < 0).any():
            raise ValueError(f"candidate assignment column {column!r} must be non-negative")

    if len(frame) != donor_snapshot.row_count:
        raise ValueError(
            f"candidate assignments for donor {donor} have {len(frame)} rows; "
            f"source evidence has {donor_snapshot.row_count}"
        )
    donor_values = frame["donor_id"].astype(str)
    if set(donor_values) != {donor}:
        raise ValueError(f"candidate assignments for donor {donor} contain another donor")
    object_hash = hash_identifiers(frame["object_id"].astype(str))
    if object_hash != donor_snapshot.object_id_sha256:
        raise ValueError(
            f"candidate assignments for donor {donor} differ from the source object universe"
        )

    expected_constants: dict[str, str | int] = {
        "candidate_assignment_run_id": assignment_run_id,
        "candidate_evaluation_split_name": split_name,
        "typing_model_bundle_content_id": model_bundle.content_id,
        "typing_model_bundle_version": model_bundle.bundle_version,
        "typing_model_method_version": model_bundle.method_version,
        "typing_feature_schema_version": model_bundle.feature_schema_version,
        "typing_score_semantics": model_bundle.score_semantics,
        "typing_broad_score_semantics": model_bundle.score_semantics,
        "typing_mode": CANDIDATE_TYPING_MODE,
        "typing_rules_version": typing_registry.rules_version,
        "typing_rules_fingerprint": typing_registry.rules_fingerprint,
        "marker_registry_fingerprint": marker_registry.fingerprint,
        "typing_broad_bootstrap_replicates": model_bundle.broad.bootstrap_replicates,
        "typing_broad_valid_bootstrap_replicates": (model_bundle.broad.valid_bootstrap_replicates),
    }
    for column, expected in expected_constants.items():
        _require_constant_assignment_column(
            frame,
            column,
            expected,
            donor=donor,
        )

    configured_parents = {str(rule.parent) for rule in typing_registry.subtype_rules}
    bundled_parents = set(model_bundle.subtype_by_parent)
    accepted_parent = frame["broad_assignment_status"].isin({AUTHORITATIVE_ANCHOR, INFERRED})
    subtype_applicable = accepted_parent & frame["broad_label"].isin(configured_parents)
    expected_source = pd.Series("not_applicable", index=frame.index, dtype="object")
    expected_source.loc[subtype_applicable & frame["broad_label"].isin(bundled_parents)] = (
        "probabilistic_bundle"
    )
    expected_source.loc[subtype_applicable & ~frame["broad_label"].isin(bundled_parents)] = (
        "rules_fallback"
    )
    observed_source = frame["typing_subtype_model_source"].astype(str)
    if not observed_source.eq(expected_source).all():
        raise ValueError(f"candidate assignments for donor {donor} misstate subtype model source")
    expected_semantics = expected_source.map(
        {
            "not_applicable": "",
            "probabilistic_bundle": model_bundle.score_semantics,
            "rules_fallback": RULE_SCORE_SEMANTICS,
        }
    )
    if not frame["typing_subtype_score_semantics"].astype(str).eq(expected_semantics).all():
        raise ValueError(
            f"candidate assignments for donor {donor} misstate subtype score semantics"
        )

    subtype_replicates = {
        str(ensemble.parent): ensemble.bootstrap_replicates for ensemble in model_bundle.subtypes
    }
    valid_subtype_replicates = {
        str(ensemble.parent): ensemble.valid_bootstrap_replicates
        for ensemble in model_bundle.subtypes
    }
    expected_replicates = frame["broad_label"].map(subtype_replicates).fillna(0).astype(int)
    expected_valid_replicates = (
        frame["broad_label"].map(valid_subtype_replicates).fillna(0).astype(int)
    )
    if not frame["typing_subtype_bootstrap_replicates"].eq(expected_replicates).all():
        raise ValueError(
            f"candidate assignments for donor {donor} misstate subtype bootstrap count"
        )
    if not frame["typing_subtype_valid_bootstrap_replicates"].eq(expected_valid_replicates).all():
        raise ValueError(
            f"candidate assignments for donor {donor} misstate valid subtype bootstrap count"
        )


def candidate_assignment_run_id(
    *,
    model_bundle_content_id: str,
    split_manifest_content_id: str,
    split_name: str,
    source_run_id: str,
    source_run_manifest_content_id: str,
    source_stage_manifest_content_id: str,
    source_artifact_content_id: str,
    marker_registry_fingerprint: str,
    typing_rules_fingerprint: str,
    donors: Sequence[str],
) -> str:
    """Return the deterministic scientific identity of an evaluation run."""

    normalised_donors = _normalise_donors(donors)
    if split_name not in EVALUATION_SPLITS:
        raise ValueError(f"candidate evaluation split must be one of {EVALUATION_SPLITS}")
    payload = {
        "schema_version": CANDIDATE_EVALUATION_SCHEMA_VERSION,
        "method_version": CANDIDATE_EVALUATION_METHOD_VERSION,
        "model_bundle_content_id": _require_sha256(
            model_bundle_content_id, name="model_bundle_content_id"
        ),
        "split_manifest_content_id": _require_sha256(
            split_manifest_content_id, name="split_manifest_content_id"
        ),
        "split_name": split_name,
        "source_run_id": str(source_run_id),
        "source_run_manifest_content_id": _require_sha256(
            source_run_manifest_content_id,
            name="source_run_manifest_content_id",
        ),
        "source_stage_manifest_content_id": _require_sha256(
            source_stage_manifest_content_id,
            name="source_stage_manifest_content_id",
        ),
        "source_artifact_content_id": _require_sha256(
            source_artifact_content_id, name="source_artifact_content_id"
        ),
        "marker_registry_fingerprint": _require_sha256(
            marker_registry_fingerprint, name="marker_registry_fingerprint"
        ),
        "typing_rules_fingerprint": _require_sha256(
            typing_rules_fingerprint, name="typing_rules_fingerprint"
        ),
        "donors": list(normalised_donors),
    }
    if not payload["source_run_id"].strip():
        raise ValueError("source_run_id must be a non-empty string")
    return sha256_json(payload)


@dataclass(frozen=True)
class CandidateEvaluationManifest:
    """Strict manifest for non-production assignments from one frozen bundle."""

    assignment_run_id: str
    model_bundle_content_id: str
    split_manifest_content_id: str
    split_name: str
    source_run_id: str
    source_run_manifest_content_id: str
    source_stage_manifest_content_id: str
    source_artifact_content_id: str
    marker_registry_fingerprint: str
    typing_rules_fingerprint: str
    donors: tuple[str, ...]
    output: DatasetSnapshot
    output_content_sha256: str
    schema_version: int = CANDIDATE_EVALUATION_SCHEMA_VERSION
    method_version: str = CANDIDATE_EVALUATION_METHOD_VERSION
    content_id: str = ""
    input_artifact: InputArtifact | None = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != CANDIDATE_EVALUATION_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported candidate-evaluation schema {self.schema_version!r}")
        if self.method_version != CANDIDATE_EVALUATION_METHOD_VERSION:
            raise ValueError(f"unsupported candidate-evaluation method {self.method_version!r}")
        if self.split_name not in EVALUATION_SPLITS:
            raise ValueError(f"candidate evaluation split must be one of {EVALUATION_SPLITS}")
        if not isinstance(self.source_run_id, str) or not self.source_run_id.strip():
            raise ValueError("source_run_id must be a non-empty string")
        object.__setattr__(self, "source_run_id", self.source_run_id.strip())
        for name in (
            "assignment_run_id",
            "model_bundle_content_id",
            "split_manifest_content_id",
            "source_run_manifest_content_id",
            "source_stage_manifest_content_id",
            "source_artifact_content_id",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
            "output_content_sha256",
        ):
            object.__setattr__(self, name, _require_sha256(getattr(self, name), name=name))
        donors = _normalise_donors(self.donors)
        object.__setattr__(self, "donors", donors)
        if not isinstance(self.output, DatasetSnapshot):
            raise TypeError("candidate output must be a DatasetSnapshot")
        if self.output.donor_ids != donors:
            raise ValueError("candidate output donor set/order differs from the evaluation split")
        if self.output.total_rows < 1 or any(donor.row_count < 1 for donor in self.output.donors):
            raise ValueError("candidate output requires at least one row per donor")
        if self.output.total_rows != sum(donor.row_count for donor in self.output.donors):
            raise ValueError("candidate output total row count is inconsistent")
        if self.output.content_sha256 != self.output_content_sha256:
            raise ValueError("candidate output content identity is inconsistent")
        expected_run_id = candidate_assignment_run_id(
            model_bundle_content_id=self.model_bundle_content_id,
            split_manifest_content_id=self.split_manifest_content_id,
            split_name=self.split_name,
            source_run_id=self.source_run_id,
            source_run_manifest_content_id=self.source_run_manifest_content_id,
            source_stage_manifest_content_id=self.source_stage_manifest_content_id,
            source_artifact_content_id=self.source_artifact_content_id,
            marker_registry_fingerprint=self.marker_registry_fingerprint,
            typing_rules_fingerprint=self.typing_rules_fingerprint,
            donors=self.donors,
        )
        if self.assignment_run_id != expected_run_id:
            raise ValueError("candidate assignment_run_id is invalid")
        expected_content_id = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected_content_id:
            raise ValueError("candidate evaluation manifest content_id is invalid")
        object.__setattr__(self, "content_id", expected_content_id)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "method_version": self.method_version,
            "assignment_run_id": self.assignment_run_id,
            "model_bundle_content_id": self.model_bundle_content_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "split_name": self.split_name,
            "source_run_id": self.source_run_id,
            "source_run_manifest_content_id": (self.source_run_manifest_content_id),
            "source_stage_manifest_content_id": (self.source_stage_manifest_content_id),
            "source_artifact_content_id": self.source_artifact_content_id,
            "marker_registry_fingerprint": self.marker_registry_fingerprint,
            "typing_rules_fingerprint": self.typing_rules_fingerprint,
            "donors": list(self.donors),
            # Absolute output roots and file paths are validation caches, not
            # scientific identity.  Keeping them out of the manifest ID makes
            # relocation/rerun of byte-identical assignments unable to alter
            # bound review draws or bootstrap seeds.
            "output": _snapshot_semantic_identity(self.output),
            "output_content_sha256": self.output_content_sha256,
        }

    def to_dict(self) -> dict[str, Any]:
        value = self.identity_payload()
        # Persist mtimes as validation caches even though they do not determine
        # scientific identity.
        value["output"] = self.output.to_dict()
        value["content_id"] = self.content_id
        return value

    @classmethod
    def create(
        cls,
        *,
        model_bundle: TypingModelBundle,
        split_name: str,
        source_run_manifest: RunManifest,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
        output: DatasetSnapshot,
    ) -> CandidateEvaluationManifest:
        provenance = model_bundle.evidence_provenance
        donors = tuple(getattr(model_bundle.split_manifest, split_name))
        run_id = candidate_assignment_run_id(
            model_bundle_content_id=model_bundle.content_id,
            split_manifest_content_id=model_bundle.split_manifest.content_id,
            split_name=split_name,
            source_run_id=model_bundle.split_manifest.source_run_id,
            source_run_manifest_content_id=provenance.source_run_manifest_content_id,
            source_stage_manifest_content_id=provenance.source_stage_manifest_content_id,
            source_artifact_content_id=provenance.source_artifact_content_id,
            marker_registry_fingerprint=marker_registry.fingerprint,
            typing_rules_fingerprint=typing_registry.rules_fingerprint,
            donors=donors,
        )
        return cls(
            assignment_run_id=run_id,
            model_bundle_content_id=model_bundle.content_id,
            split_manifest_content_id=model_bundle.split_manifest.content_id,
            split_name=split_name,
            source_run_id=model_bundle.split_manifest.source_run_id,
            source_run_manifest_content_id=source_run_manifest.content_id,
            source_stage_manifest_content_id=(provenance.source_stage_manifest_content_id),
            source_artifact_content_id=provenance.source_artifact_content_id,
            marker_registry_fingerprint=marker_registry.fingerprint,
            typing_rules_fingerprint=typing_registry.rules_fingerprint,
            donors=donors,
            output=output,
            output_content_sha256=output.content_sha256,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> CandidateEvaluationManifest:
        manifest = _require_mapping(value, where="candidate evaluation manifest")
        _require_exact_fields(manifest, _MANIFEST_FIELDS, where="candidate evaluation manifest")
        schema_version = _require_json_int(
            manifest["schema_version"], name="schema_version", minimum=1
        )
        output_value = _validate_snapshot_json(manifest["output"])
        raw_donors = manifest["donors"]
        if not isinstance(raw_donors, list):
            raise TypeError("candidate evaluation donors must be a JSON array")
        donors = _normalise_donors(raw_donors)
        return cls(
            schema_version=schema_version,
            method_version=manifest["method_version"],
            assignment_run_id=manifest["assignment_run_id"],
            model_bundle_content_id=manifest["model_bundle_content_id"],
            split_manifest_content_id=manifest["split_manifest_content_id"],
            split_name=manifest["split_name"],
            source_run_id=manifest["source_run_id"],
            source_run_manifest_content_id=(manifest["source_run_manifest_content_id"]),
            source_stage_manifest_content_id=(manifest["source_stage_manifest_content_id"]),
            source_artifact_content_id=manifest["source_artifact_content_id"],
            marker_registry_fingerprint=manifest["marker_registry_fingerprint"],
            typing_rules_fingerprint=manifest["typing_rules_fingerprint"],
            donors=donors,
            output=DatasetSnapshot.from_dict(output_value),
            output_content_sha256=manifest["output_content_sha256"],
            content_id=manifest["content_id"],
        )

    @classmethod
    def read_json(cls, path: Path | str) -> CandidateEvaluationManifest:
        resolved = Path(path).expanduser().resolve()
        artifact = InputArtifact.from_path("candidate_evaluation_manifest", resolved)
        with resolved.open(encoding="utf-8") as handle:
            value = _read_strict_json(handle)
        result = cls.from_dict(value)
        artifact.validate_current(mode="fast")
        return replace(result, input_artifact=artifact)

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    def validate_for_bundle(
        self,
        *,
        model_bundle: TypingModelBundle,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
    ) -> None:
        model_bundle.validate_current(
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        if model_bundle.input_artifact is not None:
            model_bundle.input_artifact.validate_current(mode="content")
        provenance = model_bundle.evidence_provenance
        expected = {
            "model_bundle_content_id": model_bundle.content_id,
            "split_manifest_content_id": model_bundle.split_manifest.content_id,
            "source_run_id": model_bundle.split_manifest.source_run_id,
            "source_run_manifest_content_id": (provenance.source_run_manifest_content_id),
            "source_stage_manifest_content_id": (provenance.source_stage_manifest_content_id),
            "source_artifact_content_id": provenance.source_artifact_content_id,
            "marker_registry_fingerprint": marker_registry.fingerprint,
            "typing_rules_fingerprint": typing_registry.rules_fingerprint,
        }
        for name, wanted in expected.items():
            if getattr(self, name) != wanted:
                raise ValueError(f"candidate evaluation {name} differs from its bundle")
        split_donors = tuple(getattr(model_bundle.split_manifest, self.split_name))
        if self.donors != split_donors:
            raise ValueError("candidate evaluation donor set differs from its bundle split")

    def validate_current(
        self,
        *,
        model_bundle: TypingModelBundle,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
        source_run_manifest: RunManifest,
        allow_locked_test: bool = False,
    ) -> None:
        _authorize_split(self.split_name, allow_locked_test=allow_locked_test)
        if self.input_artifact is not None:
            self.input_artifact.validate_current(mode="content")
        # Re-running construction validates both semantic identifiers.
        self.__post_init__()
        self.validate_for_bundle(
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        _evidence_root, source_stage = _validate_source_lineage(
            model_bundle=model_bundle,
            source_run_manifest=source_run_manifest,
        )
        validate_dataset_snapshot_current(self.output, mode="content")
        current = snapshot_parquet_dataset(self.output.root, expected_donors=self.donors)
        if _snapshot_semantic_identity(current) != _snapshot_semantic_identity(self.output):
            raise StaleArtifactError(
                "candidate assignment DatasetSnapshot differs from current output"
            )
        _validate_assignment_dataset(
            output_root=Path(self.output.root),
            output_snapshot=current,
            source_snapshot=source_stage.output,
            donors=self.donors,
            assignment_run_id=self.assignment_run_id,
            split_name=self.split_name,
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )


def _atomic_write_json(path: Path, value: Mapping[str, Any]) -> Path:
    path = path.expanduser().resolve()
    if path.exists():
        raise FileExistsError(f"immutable candidate manifest already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _authorize_split(split_name: str, *, allow_locked_test: bool) -> tuple[str, ...]:
    if split_name not in EVALUATION_SPLITS:
        raise ValueError(
            "candidate assignments may only read calibration, validation, or locked_test donors"
        )
    if split_name == "locked_test" and not allow_locked_test:
        raise PermissionError(
            "locked_test candidate assignments require explicit allow_locked_test=True"
        )
    return EVALUATION_SPLITS


def _validate_source_lineage(
    *,
    model_bundle: TypingModelBundle,
    source_run_manifest: RunManifest,
    evidence_path: Path | str | None = None,
) -> tuple[Path, StageManifest]:
    references = [
        reference for reference in source_run_manifest.stages if reference.stage == "calibrate"
    ]
    if len(references) != 1:
        raise ValueError("source run must contain exactly one calibrate stage")
    stage = StageManifest.read_json(references[0].manifest_path)
    root = Path(stage.output.root).expanduser().resolve()
    requested = root if evidence_path is None else Path(evidence_path).expanduser().resolve()
    validate_evidence_source(
        requested,
        provenance=model_bundle.evidence_provenance,
        splits=model_bundle.split_manifest,
        source_run_manifest=source_run_manifest,
    )
    return requested, stage


def _partition_path(root: Path, donor: str) -> Path:
    partition = (root / f"donor_id={donor}").resolve()
    if partition.parent != root.resolve():
        raise ValueError(f"unsafe donor partition identifier: {donor!r}")
    return partition


def _read_evidence_partition(root: Path, donor_snapshot: DonorDatasetSnapshot) -> pd.DataFrame:
    """Read exactly one authorized donor and prove its complete cell universe."""

    donor = donor_snapshot.donor_id
    partition = _partition_path(root, donor)
    files = sorted(partition.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no evidence Parquet files for donor {donor}")
    frames = [pd.read_parquet(path) for path in files]
    frame = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    if "object_id" not in frame:
        raise KeyError(f"evidence for donor {donor} is missing object_id")
    if "donor_id" not in frame:
        frame.insert(0, "donor_id", donor)
    frame["donor_id"] = frame["donor_id"].astype(str)
    frame["object_id"] = frame["object_id"].astype(str)
    if set(frame["donor_id"]) != {donor}:
        raise ValueError(f"evidence partition {donor} contains another donor")
    if len(frame) != donor_snapshot.row_count:
        raise ValueError(f"evidence partition {donor} row count differs from its snapshot")
    if hash_identifiers(frame["object_id"]) != donor_snapshot.object_id_sha256:
        raise ValueError(f"evidence partition {donor} object universe differs from its snapshot")
    return frame


def _read_assignment_partition(root: Path, *, donor: str) -> pd.DataFrame:
    partition = _partition_path(root, donor)
    files = sorted(partition.glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no candidate assignment Parquet files for donor {donor}")
    frames = [pd.read_parquet(path) for path in files]
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]


def _candidate_assignments_from_evidence(
    evidence: pd.DataFrame,
    *,
    assignment_run_id: str,
    split_name: str,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
) -> pd.DataFrame:
    """Return the deterministic candidate frame for one complete donor."""

    assignments = type_candidate_cells(
        evidence,
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )
    assignments["candidate_assignment_run_id"] = assignment_run_id
    assignments["candidate_evaluation_split_name"] = split_name
    return assignments


def _validate_deterministic_assignment_replay(
    observed: pd.DataFrame,
    expected: pd.DataFrame,
    *,
    donor: str,
) -> None:
    """Reject a re-snapshotted output that differs from deterministic inference."""

    if tuple(observed.columns) != tuple(expected.columns):
        raise StaleArtifactError(
            f"candidate assignments for donor {donor} differ from deterministic "
            "replay: output columns/order changed"
        )

    try:
        pd.testing.assert_frame_equal(
            observed.reset_index(drop=True),
            expected.reset_index(drop=True),
            check_dtype=False,
            check_exact=True,
        )
    except AssertionError as error:
        raise StaleArtifactError(
            f"candidate assignments for donor {donor} differ from deterministic "
            "replay of the exact source evidence"
        ) from error


def _validate_assignment_dataset(
    *,
    output_root: Path,
    output_snapshot: DatasetSnapshot,
    source_snapshot: DatasetSnapshot,
    donors: Sequence[str],
    assignment_run_id: str,
    split_name: str,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
) -> None:
    """Replay and validate every persisted donor against calibrated source."""

    expected_donors = _normalise_donors(donors)
    if output_snapshot.donor_ids != expected_donors:
        raise ValueError("candidate output donors differ from the authorized split")
    source_by_donor = {snapshot.donor_id: snapshot for snapshot in source_snapshot.donors}
    output_by_donor = {snapshot.donor_id: snapshot for snapshot in output_snapshot.donors}
    missing_source = sorted(set(expected_donors) - set(source_by_donor))
    if missing_source:
        raise ValueError(f"source evidence is missing evaluation donors: {missing_source}")
    for donor in expected_donors:
        source_donor = source_by_donor[donor]
        output_donor = output_by_donor[donor]
        if output_donor.row_count != source_donor.row_count:
            raise ValueError(
                f"candidate assignments for donor {donor} have a different row count "
                "from source evidence"
            )
        if output_donor.object_id_sha256 != source_donor.object_id_sha256:
            raise ValueError(
                f"candidate assignments for donor {donor} have a different object "
                "universe from source evidence"
            )
        frame = _read_assignment_partition(output_root, donor=donor)
        _validate_assignment_frame(
            frame,
            donor_snapshot=source_donor,
            assignment_run_id=assignment_run_id,
            split_name=split_name,
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        # A fresh snapshot can truthfully attest to edited bytes.  It cannot
        # attest that those bytes came from this model.  Re-run inference from
        # the exact source donor and compare all generated fields, keeping peak
        # memory bounded to one source/output/replay frame at a time.
        evidence = _read_evidence_partition(Path(source_snapshot.root), source_donor)
        expected = _candidate_assignments_from_evidence(
            evidence,
            assignment_run_id=assignment_run_id,
            split_name=split_name,
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        _validate_deterministic_assignment_replay(frame, expected, donor=donor)


def _write_assignment_partition(frame: pd.DataFrame, root: Path, *, donor: str) -> Path:
    partition = _partition_path(root, donor)
    partition.mkdir(parents=False, exist_ok=False)
    destination = partition / "assignments.parquet"
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=".assignments.", suffix=".tmp", dir=partition
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        frame.to_parquet(temporary, index=False)
        os.link(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def evaluate_candidate_bundle(
    *,
    model_bundle: TypingModelBundle,
    evidence_path: Path | str,
    source_run_manifest: RunManifest,
    output_root: Path | str,
    manifest_path: Path | str,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    split_name: str = "validation",
    allow_locked_test: bool = False,
) -> CandidateEvaluationManifest:
    """Generate immutable, donor-partitioned assignments for one review split."""

    _authorize_split(split_name, allow_locked_test=allow_locked_test)
    model_bundle.validate_current(
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )
    donors = tuple(getattr(model_bundle.split_manifest, split_name))
    if not donors:
        raise ValueError(f"candidate bundle has no donors in split {split_name!r}")
    donors = _normalise_donors(donors)
    evidence_root, source_stage = _validate_source_lineage(
        model_bundle=model_bundle,
        source_run_manifest=source_run_manifest,
        evidence_path=evidence_path,
    )

    output = Path(output_root).expanduser().resolve()
    manifest_target = Path(manifest_path).expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"immutable candidate output root already exists: {output}")
    if manifest_target.exists():
        raise FileExistsError(f"immutable candidate manifest already exists: {manifest_target}")
    if output == manifest_target:
        raise ValueError("candidate output root and manifest path must differ")
    if output.is_relative_to(evidence_root) or manifest_target.is_relative_to(evidence_root):
        raise ValueError("candidate assignments and manifest must be outside source evidence")

    source_snapshots = {donor.donor_id: donor for donor in source_stage.output.donors}
    missing = sorted(set(donors) - set(source_snapshots))
    if missing:
        raise ValueError(f"source evidence is missing evaluation donors: {missing}")

    run_id = candidate_assignment_run_id(
        model_bundle_content_id=model_bundle.content_id,
        split_manifest_content_id=model_bundle.split_manifest.content_id,
        split_name=split_name,
        source_run_id=model_bundle.split_manifest.source_run_id,
        source_run_manifest_content_id=(
            model_bundle.evidence_provenance.source_run_manifest_content_id
        ),
        source_stage_manifest_content_id=(
            model_bundle.evidence_provenance.source_stage_manifest_content_id
        ),
        source_artifact_content_id=(model_bundle.evidence_provenance.source_artifact_content_id),
        marker_registry_fingerprint=marker_registry.fingerprint,
        typing_rules_fingerprint=typing_registry.rules_fingerprint,
        donors=donors,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(exist_ok=False)
    for donor in donors:
        evidence = _read_evidence_partition(evidence_root, source_snapshots[donor])
        assignments = _candidate_assignments_from_evidence(
            evidence,
            assignment_run_id=run_id,
            split_name=split_name,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
            model_bundle=model_bundle,
        )
        if set(assignments["donor_id"].astype(str)) != {donor}:
            raise AssertionError("candidate typing changed the authorized donor")
        _validate_assignment_frame(
            assignments,
            donor_snapshot=source_snapshots[donor],
            assignment_run_id=run_id,
            split_name=split_name,
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        _write_assignment_partition(assignments, output, donor=donor)

    output_snapshot = snapshot_parquet_dataset(output, expected_donors=donors)
    _validate_assignment_dataset(
        output_root=output,
        output_snapshot=output_snapshot,
        source_snapshot=source_stage.output,
        donors=donors,
        assignment_run_id=run_id,
        split_name=split_name,
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )
    manifest = CandidateEvaluationManifest.create(
        model_bundle=model_bundle,
        split_name=split_name,
        source_run_manifest=source_run_manifest,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        output=output_snapshot,
    )
    if manifest.assignment_run_id != run_id:
        raise AssertionError("candidate assignment identity changed after writing")
    manifest.write_json(manifest_target)
    return manifest


def run_candidate_evaluation(
    *,
    model_bundle_path: Path | str,
    evidence_path: Path | str,
    source_run_manifest_path: Path | str,
    output_root: Path | str,
    manifest_path: Path | str,
    marker_registry_path: Path | str | None = None,
    typing_rules_path: Path | str | None = None,
    split_name: str = "validation",
    allow_locked_test: bool = False,
) -> CandidateEvaluationManifest:
    """Load an unpromoted bundle and execute the isolated evaluation path."""

    markers = load_registry(marker_registry_path)
    rules = TypingRegistry.from_marker_registry(markers, typing_rules=typing_rules_path)
    bundle = load_typing_model_bundle(
        model_bundle_path,
        marker_registry=markers,
        typing_registry=rules,
    )
    source_run = RunManifest.read_json(source_run_manifest_path)
    return evaluate_candidate_bundle(
        model_bundle=bundle,
        evidence_path=evidence_path,
        source_run_manifest=source_run,
        output_root=output_root,
        manifest_path=manifest_path,
        marker_registry=markers,
        typing_registry=rules,
        split_name=split_name,
        allow_locked_test=allow_locked_test,
    )


def load_candidate_evaluation_manifest(
    path: Path | str,
    *,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    source_run_manifest: RunManifest,
    allow_locked_test: bool = False,
) -> CandidateEvaluationManifest:
    """Load and content-validate an exact candidate assignment artifact."""

    manifest = CandidateEvaluationManifest.read_json(path)
    manifest.validate_current(
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=allow_locked_test,
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-bundle", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--source-run-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--marker-registry", type=Path)
    parser.add_argument("--typing-rules", type=Path)
    parser.add_argument("--split", choices=EVALUATION_SPLITS, default="validation")
    parser.add_argument(
        "--allow-locked-test",
        action="store_true",
        help="explicitly authorize one locked-test assignment run",
    )
    args = parser.parse_args(argv)
    run_candidate_evaluation(
        model_bundle_path=args.model_bundle,
        evidence_path=args.evidence,
        source_run_manifest_path=args.source_run_manifest,
        output_root=args.output_root,
        manifest_path=args.manifest,
        marker_registry_path=args.marker_registry,
        typing_rules_path=args.typing_rules,
        split_name=args.split,
        allow_locked_test=args.allow_locked_test,
    )
    return 0


__all__ = [
    "CANDIDATE_EVALUATION_METHOD_VERSION",
    "CANDIDATE_EVALUATION_SCHEMA_VERSION",
    "EVALUATION_SPLITS",
    "CandidateEvaluationManifest",
    "candidate_assignment_run_id",
    "evaluate_candidate_bundle",
    "load_candidate_evaluation_manifest",
    "main",
    "run_candidate_evaluation",
]


if __name__ == "__main__":  # pragma: no cover - exercised through the CLI
    raise SystemExit(main())
