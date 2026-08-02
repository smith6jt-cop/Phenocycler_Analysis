"""Fail-closed calibration provenance for probabilistic typing thresholds.

The public builder deliberately does not accept an operating-point table or
caller-authored metrics.  It joins blinded calibration samples and adjudicated
review ledgers to one content-valid candidate-assignment artifact, derives the
production-aware operating grid, and freezes the deterministic selection.  The
normalized scientific rows and the complete derived grid are embedded so a
strict JSON load can repeat the calculation without reopening mutable files.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import InitVar, dataclass, field
from numbers import Real
from pathlib import Path
from types import MappingProxyType
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from .artifacts import (
    dataframe_content_sha256,
    sha256_json,
    snapshot_parquet_dataset,
    validate_dataset_snapshot_current,
)
from .refinement_contracts import (
    REVIEW_LABEL_COLUMNS,
    REVIEW_SAMPLE_AUDIT_COLUMNS,
    REVIEW_SAMPLE_BLINDED_COLUMNS,
    DonorSplitManifest,
    ReviewLabelProvenance,
    blinded_review_sample_from_dataset,
    review_label_frame_fingerprint,
    review_sample_frame_fingerprint,
    reviewer_frame_context_artifact_content_id,
    reviewer_frame_fingerprint,
    validate_review_labels,
    validate_review_sampling_against_assignments,
)

if TYPE_CHECKING:  # pragma: no cover - import cycle documentation only
    from .artifacts import RunManifest
    from .candidate_evaluation import CandidateEvaluationManifest
    from .hierarchical_typing import TypingRegistry
    from .marker_registry import MarkerRegistry
    from .typing_model_bundle import TypingModelBundle


POLICY_SCHEMA_VERSION = 1
PROVENANCE_SCHEMA_VERSION = 3
DECISION_INPUT_SCHEMA_VERSION = 2
OPERATING_POINT_SCHEMA_VERSION = 1

DEFAULT_PROBABILITY_THRESHOLDS = (0.70, 0.75, 0.80, 0.85, 0.90)
DEFAULT_MARGIN_THRESHOLDS = (0.10, 0.20, 0.25, 0.30)
DEFAULT_STABILITY_THRESHOLDS = (0.80, 0.85, 0.90, 0.95)

CALIBRATION_FRAME_KEY_COLUMNS = ("donor_id", "object_id", "level", "parent")
OPERATING_POINT_COLUMNS = (
    "probability_threshold",
    "margin_threshold",
    "stability_threshold",
    "coverage",
    "selective_error",
    "correct_call_rate",
    "false_positive_weight",
    "false_negative_weight",
    "rescue_eligible_weight",
    "rescue_called_weight",
    "rescue_coverage",
    "rescue_selective_error",
    "error_cost",
)
SELECTED_METRIC_COLUMNS = OPERATING_POINT_COLUMNS[3:]
_THRESHOLD_COLUMNS = OPERATING_POINT_COLUMNS[:3]
_ASSIGNMENT_THRESHOLD_FIELDS = (
    "inferred_probability",
    "inferred_margin",
    "inferred_stability",
    "anchor_probability",
    "negative_probability",
)
_OPTIONAL_RATE_COLUMNS = (
    "selective_error",
    "rescue_coverage",
    "rescue_selective_error",
)
_RATE_COLUMNS = (
    "coverage",
    "selective_error",
    "correct_call_rate",
    "rescue_coverage",
    "rescue_selective_error",
)
_WEIGHT_COLUMNS = (
    "false_positive_weight",
    "false_negative_weight",
    "rescue_eligible_weight",
    "rescue_called_weight",
)
_RESCUE_REASONS = {
    "nonunique_model_winner",
    "valid_positive_support_required",
    "posterior_margin_or_stability_failed",
}
_ALLOWED_STATUSES = {"anchor", "inferred", "other", "ambiguous", "unavailable"}
_DISPOSITION_LABEL = {
    "other": "Other",
    "ambiguous": "Ambiguous",
    "unavailable": "Unavailable",
}


def _reject_unknown(value: Mapping[str, Any], allowed: set[str], *, where: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(f"{where} contains unknown fields: {unknown}")


def _require_exact_fields(value: Mapping[str, Any], required: set[str], *, where: str) -> None:
    _reject_unknown(value, required, where=where)
    missing = sorted(required - set(value))
    if missing:
        raise ValueError(f"{where} is missing required fields: {missing}")


def _strict_string(value: Any, *, name: str, allow_blank: bool = False) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{name} must be a JSON string")
    result = value.strip()
    if not allow_blank and not result:
        raise ValueError(f"{name} must be non-empty")
    return result


def _strict_sha256(value: Any, *, name: str) -> str:
    text = _strict_string(value, name=name)
    if re.fullmatch(r"[0-9a-f]{64}", text) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return text


def _strict_positive_int(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be a JSON integer")
    if value < 1:
        raise ValueError(f"{name} must be positive")
    return value


def _strict_bool(value: Any, *, name: str) -> bool:
    if not isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a JSON boolean")
    return bool(value)


def _strict_finite_float(value: Any, *, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite JSON number")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _unit_interval(value: Any, *, name: str) -> float:
    result = _strict_finite_float(value, name=name)
    if result < 0.0 or result > 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return result


def _strict_thresholds(values: Sequence[Any], *, name: str) -> tuple[float, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence):
        raise TypeError(f"{name} must be a JSON array/sequence")
    result = tuple(_unit_interval(value, name=name) for value in values)
    if not result:
        raise ValueError(f"{name} must not be empty")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} contains duplicate thresholds")
    return tuple(sorted(result))


def _strict_assignment_thresholds(value: Mapping[str, Any]) -> Mapping[str, float]:
    if not isinstance(value, Mapping):
        raise TypeError("preliminary_assignment_thresholds must be a mapping")
    _require_exact_fields(
        value,
        set(_ASSIGNMENT_THRESHOLD_FIELDS),
        where="preliminary_assignment_thresholds",
    )
    return MappingProxyType(
        {
            name: _unit_interval(
                value[name], name=f"preliminary_assignment_thresholds.{name}"
            )
            for name in _ASSIGNMENT_THRESHOLD_FIELDS
        }
    )


def threshold_grid_content_id(
    *,
    probability_thresholds: Sequence[float],
    margin_thresholds: Sequence[float],
    stability_thresholds: Sequence[float],
) -> str:
    """Return the preregisterable identity of a threshold search grid."""

    return sha256_json(
        {
            "schema_version": 1,
            "probability_thresholds": list(
                _strict_thresholds(
                    probability_thresholds, name="probability_thresholds"
                )
            ),
            "margin_thresholds": list(
                _strict_thresholds(margin_thresholds, name="margin_thresholds")
            ),
            "stability_thresholds": list(
                _strict_thresholds(
                    stability_thresholds, name="stability_thresholds"
                )
            ),
        }
    )


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path = path.expanduser().resolve()
    if path.exists():
        raise FileExistsError(f"immutable threshold artifact already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)
    return path


def _strict_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite number {value!r}")


def _read_strict_json(path: Path | str) -> Mapping[str, Any]:
    with Path(path).expanduser().open(encoding="utf-8") as handle:
        value = json.load(
            handle,
            object_pairs_hook=_strict_object,
            parse_constant=_reject_json_constant,
        )
    if not isinstance(value, Mapping):
        raise TypeError("threshold artifact must contain one JSON object")
    return value


@dataclass(frozen=True)
class ThresholdSelectionPolicy:
    """Predeclared FP/FN costs and total plus rescue acceptance gates."""

    policy_version: str
    false_positive_cost: float
    false_negative_cost: float
    maximum_selective_error: float
    minimum_coverage: float
    maximum_rescue_selective_error: float
    minimum_rescue_coverage: float
    schema_version: int = POLICY_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != POLICY_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported threshold-policy schema_version {self.schema_version!r}")
        object.__setattr__(
            self, "policy_version", _strict_string(self.policy_version, name="policy_version")
        )
        for name in ("false_positive_cost", "false_negative_cost"):
            number = _strict_finite_float(getattr(self, name), name=name)
            if number < 0:
                raise ValueError(f"{name} must be >= 0")
            object.__setattr__(self, name, number)
        if self.false_positive_cost == 0 and self.false_negative_cost == 0:
            raise ValueError("at least one false-positive/false-negative cost must be positive")
        for name in (
            "maximum_selective_error",
            "minimum_coverage",
            "maximum_rescue_selective_error",
            "minimum_rescue_coverage",
        ):
            object.__setattr__(self, name, _unit_interval(getattr(self, name), name=name))
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("threshold-selection policy content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "policy_version": self.policy_version,
            "false_positive_cost": self.false_positive_cost,
            "false_negative_cost": self.false_negative_cost,
            "maximum_selective_error": self.maximum_selective_error,
            "minimum_coverage": self.minimum_coverage,
            "maximum_rescue_selective_error": self.maximum_rescue_selective_error,
            "minimum_rescue_coverage": self.minimum_rescue_coverage,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> ThresholdSelectionPolicy:
        if not isinstance(value, Mapping):
            raise TypeError("threshold-selection policy must be a mapping")
        fields = {
            "schema_version",
            "policy_version",
            "false_positive_cost",
            "false_negative_cost",
            "maximum_selective_error",
            "minimum_coverage",
            "maximum_rescue_selective_error",
            "minimum_rescue_coverage",
            "content_id",
        }
        _require_exact_fields(value, fields, where="threshold-selection policy")
        if isinstance(value["schema_version"], bool) or not isinstance(
            value["schema_version"], int
        ):
            raise TypeError("threshold-selection policy schema_version must be an integer")
        return cls(
            schema_version=value["schema_version"],
            policy_version=value["policy_version"],
            false_positive_cost=value["false_positive_cost"],
            false_negative_cost=value["false_negative_cost"],
            maximum_selective_error=value["maximum_selective_error"],
            minimum_coverage=value["minimum_coverage"],
            maximum_rescue_selective_error=value["maximum_rescue_selective_error"],
            minimum_rescue_coverage=value["minimum_rescue_coverage"],
            content_id=_strict_sha256(value["content_id"], name="content_id"),
        )

    @classmethod
    def read_json(cls, path: Path | str) -> ThresholdSelectionPolicy:
        return cls.from_dict(_read_strict_json(path))


@dataclass(frozen=True)
class CalibrationDecisionRow:
    """One reviewed cell-level decision, normalized for deterministic replay."""

    donor_id: str
    object_id: str
    reference_label: str
    reviewer_1: str
    reviewer_2: str
    adjudicated: bool
    root_cause: str
    evidence_sufficient: bool
    sampling_stratum: str
    sampling_method: str
    sampling_salt: str
    stratum_population: int
    stratum_sampled: int
    sampling_weight: float
    selection_probability: float
    current_prediction: str
    current_assignment_status: str
    current_reason: str
    candidate_prediction: str
    threshold_eligible: bool
    valid_evidence: bool
    stability: float | None
    probabilities: tuple[float, ...]
    schema_version: int = DECISION_INPUT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != DECISION_INPUT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported calibration-row schema {self.schema_version!r}")
        for name in (
            "donor_id",
            "object_id",
            "reference_label",
            "reviewer_1",
            "reviewer_2",
            "root_cause",
            "sampling_stratum",
            "sampling_method",
            "sampling_salt",
            "current_prediction",
            "current_reason",
        ):
            object.__setattr__(self, name, _strict_string(getattr(self, name), name=name))
        status = _strict_string(
            self.current_assignment_status, name="current_assignment_status"
        ).lower()
        if status not in _ALLOWED_STATUSES:
            raise ValueError(f"unknown current_assignment_status {status!r}")
        object.__setattr__(self, "current_assignment_status", status)
        object.__setattr__(
            self, "adjudicated", _strict_bool(self.adjudicated, name="adjudicated")
        )
        object.__setattr__(
            self,
            "evidence_sufficient",
            _strict_bool(self.evidence_sufficient, name="evidence_sufficient"),
        )
        if not self.adjudicated or not self.evidence_sufficient:
            raise ValueError(
                "threshold calibration rows require adjudicated, evidence-sufficient review"
            )
        for name in ("stratum_population", "stratum_sampled"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
                raise TypeError(f"{name} must be an integer")
            if int(value) < 1:
                raise ValueError(f"{name} must be positive")
            object.__setattr__(self, name, int(value))
        if self.stratum_sampled > self.stratum_population:
            raise ValueError("stratum_sampled cannot exceed stratum_population")
        object.__setattr__(
            self,
            "candidate_prediction",
            _strict_string(
                self.candidate_prediction, name="candidate_prediction", allow_blank=True
            ),
        )
        weight = _strict_finite_float(self.sampling_weight, name="sampling_weight")
        if weight <= 0:
            raise ValueError("sampling_weight must be positive")
        object.__setattr__(self, "sampling_weight", weight)
        selection_probability = _unit_interval(
            self.selection_probability, name="selection_probability"
        )
        if selection_probability <= 0:
            raise ValueError("selection_probability must be positive")
        expected_weight = self.stratum_population / self.stratum_sampled
        if not np.isclose(weight, expected_weight, rtol=1e-12, atol=1e-12):
            raise ValueError("sampling_weight differs from stratum population/sample count")
        if not np.isclose(
            selection_probability, 1.0 / weight, rtol=1e-12, atol=1e-12
        ):
            raise ValueError("selection_probability differs from inverse sampling weight")
        object.__setattr__(self, "selection_probability", selection_probability)
        object.__setattr__(
            self,
            "threshold_eligible",
            _strict_bool(self.threshold_eligible, name="threshold_eligible"),
        )
        object.__setattr__(
            self, "valid_evidence", _strict_bool(self.valid_evidence, name="valid_evidence")
        )
        if self.stability is not None:
            object.__setattr__(self, "stability", _unit_interval(self.stability, name="stability"))
        probabilities = tuple(
            _unit_interval(value, name="probability") for value in self.probabilities
        )
        if len(probabilities) < 2 or not np.isclose(sum(probabilities), 1.0, atol=1e-6, rtol=0.0):
            raise ValueError("calibration probabilities require >=2 classes and must sum to one")
        object.__setattr__(self, "probabilities", probabilities)

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "donor_id": self.donor_id,
            "object_id": self.object_id,
            "reference_label": self.reference_label,
            "reviewer_1": self.reviewer_1,
            "reviewer_2": self.reviewer_2,
            "adjudicated": self.adjudicated,
            "root_cause": self.root_cause,
            "evidence_sufficient": self.evidence_sufficient,
            "sampling_stratum": self.sampling_stratum,
            "sampling_method": self.sampling_method,
            "sampling_salt": self.sampling_salt,
            "stratum_population": self.stratum_population,
            "stratum_sampled": self.stratum_sampled,
            "sampling_weight": self.sampling_weight,
            "selection_probability": self.selection_probability,
            "current_prediction": self.current_prediction,
            "current_assignment_status": self.current_assignment_status,
            "current_reason": self.current_reason,
            "candidate_prediction": self.candidate_prediction,
            "threshold_eligible": self.threshold_eligible,
            "valid_evidence": self.valid_evidence,
            "stability": self.stability,
            "probabilities": list(self.probabilities),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> CalibrationDecisionRow:
        fields = {
            "schema_version",
            "donor_id",
            "object_id",
            "reference_label",
            "reviewer_1",
            "reviewer_2",
            "adjudicated",
            "root_cause",
            "evidence_sufficient",
            "sampling_stratum",
            "sampling_method",
            "sampling_salt",
            "stratum_population",
            "stratum_sampled",
            "sampling_weight",
            "selection_probability",
            "current_prediction",
            "current_assignment_status",
            "current_reason",
            "candidate_prediction",
            "threshold_eligible",
            "valid_evidence",
            "stability",
            "probabilities",
        }
        if not isinstance(value, Mapping):
            raise TypeError("calibration decision row must be a mapping")
        _require_exact_fields(value, fields, where="calibration decision row")
        if not isinstance(value["probabilities"], list):
            raise TypeError("calibration probabilities must be a JSON array")
        if value["stability"] is not None and isinstance(value["stability"], bool):
            raise TypeError("stability must be null or a JSON number")
        return cls(
            schema_version=value["schema_version"],
            donor_id=value["donor_id"],
            object_id=value["object_id"],
            reference_label=value["reference_label"],
            reviewer_1=value["reviewer_1"],
            reviewer_2=value["reviewer_2"],
            adjudicated=value["adjudicated"],
            root_cause=value["root_cause"],
            evidence_sufficient=value["evidence_sufficient"],
            sampling_stratum=value["sampling_stratum"],
            sampling_method=value["sampling_method"],
            sampling_salt=value["sampling_salt"],
            stratum_population=value["stratum_population"],
            stratum_sampled=value["stratum_sampled"],
            sampling_weight=value["sampling_weight"],
            selection_probability=value["selection_probability"],
            current_prediction=value["current_prediction"],
            current_assignment_status=value["current_assignment_status"],
            current_reason=value["current_reason"],
            candidate_prediction=value["candidate_prediction"],
            threshold_eligible=value["threshold_eligible"],
            valid_evidence=value["valid_evidence"],
            stability=value["stability"],
            probabilities=tuple(value["probabilities"]),
        )


@dataclass(frozen=True)
class CalibrationTargetDecisionInput:
    """All reviewed inputs for one broad or probabilistic subtype target."""

    level: str
    parent: str
    classes: tuple[str, ...]
    rows: tuple[CalibrationDecisionRow, ...]
    schema_version: int = DECISION_INPUT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != DECISION_INPUT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported calibration-target schema {self.schema_version!r}")
        level = _strict_string(self.level, name="level").lower()
        parent = _strict_string(self.parent, name="parent", allow_blank=True)
        if (
            level not in {"broad", "subtype"}
            or (level == "broad" and parent)
            or (level == "subtype" and not parent)
        ):
            raise ValueError("broad target requires blank parent; subtype target requires parent")
        object.__setattr__(self, "level", level)
        object.__setattr__(self, "parent", parent)
        if isinstance(self.classes, (str, bytes)) or not isinstance(self.classes, Sequence):
            raise TypeError("target classes must be a sequence")
        classes = tuple(_strict_string(value, name="class") for value in self.classes)
        if len(classes) < 2 or len(set(classes)) != len(classes):
            raise ValueError("probabilistic target requires at least two unique classes")
        object.__setattr__(self, "classes", classes)
        if isinstance(self.rows, (str, bytes)) or not isinstance(self.rows, Sequence):
            raise TypeError("target rows must be a sequence")
        rows = tuple(self.rows)
        if not rows or any(not isinstance(row, CalibrationDecisionRow) for row in rows):
            raise ValueError("target rows must contain CalibrationDecisionRow values")
        rows = tuple(sorted(rows, key=lambda row: (row.donor_id, row.object_id)))
        keys = [(row.donor_id, row.object_id) for row in rows]
        if len(set(keys)) != len(keys):
            raise ValueError("calibration target contains duplicate cell keys")
        for row in rows:
            if len(row.probabilities) != len(classes):
                raise ValueError("calibration row probability count differs from target classes")
            if row.reference_label not in set(classes) | {"Other"}:
                raise ValueError("reference label is outside target classes plus Other")
            status = row.current_assignment_status
            if status in _DISPOSITION_LABEL:
                if row.current_prediction != _DISPOSITION_LABEL[status]:
                    raise ValueError("current disposition label differs from assignment status")
            elif row.current_prediction not in classes:
                raise ValueError("anchor/inferred current prediction is outside target classes")
            expected_eligible = status == "inferred" or (
                status == "ambiguous" and row.current_reason in _RESCUE_REASONS
            )
            if row.threshold_eligible != expected_eligible:
                raise ValueError(
                    "threshold_eligible differs from assignment status/reason contract"
                )
            best = classes[int(np.argmax(np.asarray(row.probabilities, dtype=float)))]
            if row.threshold_eligible and row.candidate_prediction != best:
                raise ValueError("eligible candidate differs from probability argmax")
            if row.threshold_eligible and row.stability is None:
                raise ValueError("threshold-eligible rows require finite stability")
        references = [row.reference_label for row in rows]
        for class_name in classes:
            positive = any(value == class_name for value in references)
            negative = any(value != class_name for value in references)
            if not positive or not negative:
                raise ValueError(
                    f"target class {class_name!r} lacks nonzero reviewed positive "
                    "or negative support"
                )
        object.__setattr__(self, "rows", rows)

    @property
    def target(self) -> tuple[str, str]:
        return self.level, self.parent

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "level": self.level,
            "parent": self.parent,
            "classes": list(self.classes),
            "rows": [row.to_dict() for row in self.rows],
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> CalibrationTargetDecisionInput:
        fields = {"schema_version", "level", "parent", "classes", "rows"}
        if not isinstance(value, Mapping):
            raise TypeError("calibration target input must be a mapping")
        _require_exact_fields(value, fields, where="calibration target input")
        if not isinstance(value["classes"], list) or not isinstance(value["rows"], list):
            raise TypeError("calibration target classes/rows must be JSON arrays")
        return cls(
            schema_version=value["schema_version"],
            level=value["level"],
            parent=value["parent"],
            classes=tuple(value["classes"]),
            rows=tuple(CalibrationDecisionRow.from_dict(row) for row in value["rows"]),
        )


def _embedded_sampling_audit(target: CalibrationTargetDecisionInput) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "donor_id": row.donor_id,
                "object_id": row.object_id,
                "sampling_stratum": row.sampling_stratum,
                "sampling_method": row.sampling_method,
                "sampling_salt": row.sampling_salt,
                "stratum_population": row.stratum_population,
                "stratum_sampled": row.stratum_sampled,
                "sampling_weight": row.sampling_weight,
                "selection_probability": row.selection_probability,
            }
            for row in target.rows
        ],
        columns=list(REVIEW_SAMPLE_AUDIT_COLUMNS),
    )


def _embedded_review_ledger(target: CalibrationTargetDecisionInput) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "donor_id": row.donor_id,
                "object_id": row.object_id,
                "reference_label": row.reference_label,
                "reviewer_1": row.reviewer_1,
                "reviewer_2": row.reviewer_2,
                "adjudicated": row.adjudicated,
                "root_cause": row.root_cause,
                "evidence_sufficient": row.evidence_sufficient,
            }
            for row in target.rows
        ],
        columns=list(REVIEW_LABEL_COLUMNS),
    )


def _normalise_operating_grid(table: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(table, pd.DataFrame) or table.empty:
        raise ValueError("operating-point grid must be a nonempty pandas DataFrame")
    if table.columns.has_duplicates:
        raise ValueError("operating-point grid columns must be unique")
    observed = set(map(str, table.columns))
    expected = set(OPERATING_POINT_COLUMNS)
    if observed != expected:
        raise ValueError(
            "operating-point grid schema differs from the acceptance API: "
            f"missing={sorted(expected - observed)}, unknown={sorted(observed - expected)}"
        )
    work = table.loc[:, list(OPERATING_POINT_COLUMNS)].copy()
    for column in OPERATING_POINT_COLUMNS:
        if is_bool_dtype(work[column].dtype) or not is_numeric_dtype(work[column].dtype):
            raise TypeError(f"operating-point column {column!r} must be numeric")
        work[column] = work[column].astype(float)
    for column in _THRESHOLD_COLUMNS:
        if not np.isfinite(work[column]).all() or not work[column].between(0, 1).all():
            raise ValueError(f"operating-point {column} must lie in [0, 1]")
    for column in _RATE_COLUMNS:
        finite = np.isfinite(work[column])
        if np.isinf(work[column]).any() or not work.loc[finite, column].between(0, 1).all():
            raise ValueError(f"operating-point {column} must be null or lie in [0, 1]")
    for column in (*_WEIGHT_COLUMNS, "error_cost"):
        if not np.isfinite(work[column]).all() or (work[column] < 0).any():
            raise ValueError(f"operating-point {column} must be finite and >= 0")
    if work.duplicated(list(_THRESHOLD_COLUMNS)).any():
        raise ValueError("operating-point grid contains duplicate threshold rows")
    expected_size = int(np.prod([work[column].nunique() for column in _THRESHOLD_COLUMNS]))
    if expected_size != len(work):
        raise ValueError("operating-point table is not a complete Cartesian grid")
    return work.sort_values(list(_THRESHOLD_COLUMNS), kind="mergesort").reset_index(drop=True)


def operating_grid_fingerprint(table: pd.DataFrame) -> str:
    return dataframe_content_sha256(
        _normalise_operating_grid(table), sort_by=list(_THRESHOLD_COLUMNS)
    )


@dataclass(frozen=True)
class FrozenOperatingPoint:
    probability_threshold: float
    margin_threshold: float
    stability_threshold: float
    coverage: float
    selective_error: float | None
    correct_call_rate: float
    false_positive_weight: float
    false_negative_weight: float
    rescue_eligible_weight: float
    rescue_called_weight: float
    rescue_coverage: float | None
    rescue_selective_error: float | None
    error_cost: float
    schema_version: int = OPERATING_POINT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != OPERATING_POINT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported operating-point schema {self.schema_version!r}")
        for name in _THRESHOLD_COLUMNS:
            object.__setattr__(self, name, _unit_interval(getattr(self, name), name=name))
        for name in _RATE_COLUMNS:
            value = getattr(self, name)
            if value is None:
                if name not in _OPTIONAL_RATE_COLUMNS:
                    raise ValueError(f"{name} cannot be null")
            else:
                object.__setattr__(self, name, _unit_interval(value, name=name))
        for name in (*_WEIGHT_COLUMNS, "error_cost"):
            number = _strict_finite_float(getattr(self, name), name=name)
            if number < 0:
                raise ValueError(f"{name} must be >= 0")
            object.__setattr__(self, name, number)
        if self.correct_call_rate > self.coverage + 1e-12:
            raise ValueError("correct_call_rate cannot exceed coverage")
        if self.rescue_called_weight > self.rescue_eligible_weight + 1e-12:
            raise ValueError("rescue calls cannot exceed rescue eligibility")

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            **{name: getattr(self, name) for name in OPERATING_POINT_COLUMNS},
        }

    @classmethod
    def from_series(cls, row: pd.Series) -> FrozenOperatingPoint:
        values: dict[str, Any] = {}
        for name in OPERATING_POINT_COLUMNS:
            value = float(row[name])
            values[name] = None if name in _OPTIONAL_RATE_COLUMNS and np.isnan(value) else value
        return cls(**values)

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> FrozenOperatingPoint:
        fields = {"schema_version", *OPERATING_POINT_COLUMNS}
        if not isinstance(value, Mapping):
            raise TypeError("operating point must be a mapping")
        _require_exact_fields(value, fields, where="operating point")
        return cls(**dict(value))


def _grid_from_points(points: Sequence[FrozenOperatingPoint]) -> pd.DataFrame:
    return _normalise_operating_grid(
        pd.DataFrame(
            [{name: getattr(point, name) for name in OPERATING_POINT_COLUMNS} for point in points]
        )
    )


def _decision_inputs_fingerprint(inputs: Sequence[CalibrationTargetDecisionInput]) -> str:
    return sha256_json(
        {
            "schema_version": DECISION_INPUT_SCHEMA_VERSION,
            "targets": [item.to_dict() for item in inputs],
        }
    )


def calibration_frame_fingerprint(frame: pd.DataFrame) -> str:
    """Fingerprint every column of an externally supplied calibration frame.

    This utility remains useful to notebook callers.  Threshold construction
    itself fingerprints its stricter embedded :class:`CalibrationTargetDecisionInput` rows.
    """

    if not isinstance(frame, pd.DataFrame) or frame.empty:
        raise ValueError("calibration frame must be a nonempty pandas DataFrame")
    if frame.columns.has_duplicates:
        raise ValueError("calibration frame columns must be unique")
    missing = sorted(set(CALIBRATION_FRAME_KEY_COLUMNS) - set(frame.columns))
    if missing:
        raise KeyError(f"calibration frame is missing key columns: {missing}")
    work = frame.copy()
    for column in CALIBRATION_FRAME_KEY_COLUMNS:
        work[column] = work[column].astype("string").fillna("").str.strip()
    if (work[["donor_id", "object_id", "level"]] == "").any().any():
        raise ValueError("calibration frame contains blank keys")
    if work.duplicated(list(CALIBRATION_FRAME_KEY_COLUMNS)).any():
        raise ValueError("calibration frame contains duplicate target/cell keys")
    ordered = [
        *CALIBRATION_FRAME_KEY_COLUMNS,
        *sorted(set(work.columns) - set(CALIBRATION_FRAME_KEY_COLUMNS)),
    ]
    return dataframe_content_sha256(
        work.loc[:, ordered], sort_by=list(CALIBRATION_FRAME_KEY_COLUMNS)
    )


def _derive_operating_grid(
    inputs: Sequence[CalibrationTargetDecisionInput],
    *,
    probability_thresholds: Sequence[float],
    margin_thresholds: Sequence[float],
    stability_thresholds: Sequence[float],
    policy: ThresholdSelectionPolicy,
) -> pd.DataFrame:
    from .probabilistic_training import acceptance_operating_points

    targets = tuple(inputs)
    global_classes: list[str] = []
    class_namespaces: list[dict[str, str]] = []
    for index, target in enumerate(targets):
        namespace = {class_name: f"T{index}::{class_name}" for class_name in target.classes}
        class_namespaces.append(namespace)
        global_classes.extend(namespace.values())
    references: list[str] = []
    probabilities: list[list[float]] = []
    stability: list[float] = []
    weights: list[float] = []
    current: list[str] = []
    statuses: list[str] = []
    reasons: list[str] = []
    candidates: list[str] = []
    eligible: list[bool] = []
    evidence: list[bool] = []
    offsets = np.cumsum((0, *(len(target.classes) for target in targets)))
    for target_index, target in enumerate(targets):
        namespace = class_namespaces[target_index]
        start, stop = int(offsets[target_index]), int(offsets[target_index + 1])
        for row in target.rows:
            vector = [0.0] * len(global_classes)
            vector[start:stop] = row.probabilities
            probabilities.append(vector)
            references.append(
                "Other" if row.reference_label == "Other" else namespace[row.reference_label]
            )
            current.append(
                row.current_prediction
                if row.current_prediction in _DISPOSITION_LABEL.values()
                else namespace[row.current_prediction]
            )
            candidates.append(
                namespace.get(
                    row.candidate_prediction, row.candidate_prediction or global_classes[0]
                )
            )
            statuses.append(row.current_assignment_status)
            reasons.append(row.current_reason)
            eligible.append(row.threshold_eligible)
            evidence.append(row.valid_evidence)
            stability.append(np.nan if row.stability is None else row.stability)
            weights.append(row.sampling_weight)
    return _normalise_operating_grid(
        acceptance_operating_points(
            reference=references,
            probabilities=np.asarray(probabilities, dtype=float),
            classes=tuple(global_classes),
            stability=stability,
            probability_thresholds=probability_thresholds,
            margin_thresholds=margin_thresholds,
            stability_thresholds=stability_thresholds,
            weights=weights,
            false_positive_cost=policy.false_positive_cost,
            false_negative_cost=policy.false_negative_cost,
            current_prediction=current,
            current_assignment_status=statuses,
            current_reason=reasons,
            candidate_prediction=candidates,
            threshold_eligible=eligible,
            valid_evidence=evidence,
        )
    )


def _select_targetwise_operating_point(
    inputs: Sequence[CalibrationTargetDecisionInput],
    aggregate_grid: pd.DataFrame,
    *,
    probability_thresholds: Sequence[float],
    margin_thresholds: Sequence[float],
    stability_thresholds: Sequence[float],
    policy: ThresholdSelectionPolicy,
) -> pd.Series:
    """Select only among threshold tuples feasible for every modeled target."""

    threshold_columns = list(_THRESHOLD_COLUMNS)
    common: set[tuple[float, float, float]] | None = None
    failed_targets: list[tuple[str, str]] = []
    for target in inputs:
        target_grid = _derive_operating_grid(
            (target,),
            probability_thresholds=probability_thresholds,
            margin_thresholds=margin_thresholds,
            stability_thresholds=stability_thresholds,
            policy=policy,
        )
        feasible = target_grid.loc[
            target_grid["selective_error"].le(policy.maximum_selective_error)
            & target_grid["coverage"].ge(policy.minimum_coverage)
            & target_grid["rescue_selective_error"].le(policy.maximum_rescue_selective_error)
            & target_grid["rescue_coverage"].ge(policy.minimum_rescue_coverage)
        ]
        keys = {
            tuple(float(row[column]) for column in threshold_columns)
            for _, row in feasible.iterrows()
        }
        if not keys:
            failed_targets.append(target.target)
        common = keys if common is None else common & keys
    if failed_targets or not common:
        detail = failed_targets if failed_targets else [item.target for item in inputs]
        raise ValueError(
            "no common operating point meets all total/rescue gates for every "
            f"probabilistic target; failing_targets={detail}"
        )
    aggregate_keys = pd.MultiIndex.from_frame(aggregate_grid.loc[:, threshold_columns])
    allowed = pd.MultiIndex.from_tuples(sorted(common), names=threshold_columns)
    restricted = aggregate_grid.loc[aggregate_keys.isin(allowed)].copy()
    from .probabilistic_training import choose_operating_point

    return choose_operating_point(
        restricted,
        maximum_selective_error=policy.maximum_selective_error,
        minimum_coverage=policy.minimum_coverage,
        maximum_rescue_selective_error=policy.maximum_rescue_selective_error,
        minimum_rescue_coverage=policy.minimum_rescue_coverage,
    )


def _selected_metrics(value: Mapping[str, Any]) -> Mapping[str, float]:
    if not isinstance(value, Mapping):
        raise TypeError("selected_metrics must be a mapping")
    _require_exact_fields(value, set(SELECTED_METRIC_COLUMNS), where="selected_metrics")
    result = {
        name: _strict_finite_float(value[name], name=f"selected_metrics.{name}")
        for name in SELECTED_METRIC_COLUMNS
    }
    for name in _RATE_COLUMNS:
        if not 0 <= result[name] <= 1:
            raise ValueError(f"selected_metrics.{name} must lie in [0, 1]")
    for name in (*_WEIGHT_COLUMNS, "error_cost"):
        if result[name] < 0:
            raise ValueError(f"selected_metrics.{name} must be >= 0")
    return MappingProxyType(result)


@dataclass(frozen=True)
class ThresholdSelectionProvenance:
    """Replayable calibration lineage and selector-derived threshold decision."""

    selection_version: str
    preliminary_bundle_content_id: str
    model_fit_content_id: str
    preliminary_assignment_thresholds: Mapping[str, float]
    split_manifest: DonorSplitManifest
    calibration_candidate_manifest_content_id: str
    calibration_candidate_output_content_id: str
    calibration_review_provenance: tuple[ReviewLabelProvenance, ...]
    assignment_run_id: str
    calibration_donors: tuple[str, ...]
    calibration_inputs: tuple[CalibrationTargetDecisionInput, ...]
    calibration_row_count: int
    calibration_frame_fingerprint: str
    probability_thresholds: tuple[float, ...]
    margin_thresholds: tuple[float, ...]
    stability_thresholds: tuple[float, ...]
    operating_grid: tuple[FrozenOperatingPoint, ...]
    operating_grid_row_count: int
    operating_grid_fingerprint: str
    policy: ThresholdSelectionPolicy
    selected_probability_threshold: float
    selected_margin_threshold: float
    selected_stability_threshold: float
    selected_metrics: Mapping[str, float]
    schema_version: int = PROVENANCE_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or not isinstance(self.schema_version, int)
            or self.schema_version != PROVENANCE_SCHEMA_VERSION
        ):
            raise ValueError(
                f"unsupported threshold-provenance schema_version {self.schema_version!r}"
            )
        object.__setattr__(
            self,
            "selection_version",
            _strict_string(self.selection_version, name="selection_version"),
        )
        for name in (
            "preliminary_bundle_content_id",
            "model_fit_content_id",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
        ):
            object.__setattr__(self, name, _strict_sha256(getattr(self, name), name=name))
        preliminary_thresholds = _strict_assignment_thresholds(
            self.preliminary_assignment_thresholds
        )
        object.__setattr__(
            self, "preliminary_assignment_thresholds", preliminary_thresholds
        )
        if not isinstance(self.split_manifest, DonorSplitManifest):
            raise TypeError("split_manifest must be DonorSplitManifest")
        if not isinstance(self.policy, ThresholdSelectionPolicy):
            raise TypeError("policy must be ThresholdSelectionPolicy")
        object.__setattr__(
            self,
            "assignment_run_id",
            _strict_sha256(self.assignment_run_id, name="assignment_run_id"),
        )
        donors = tuple(
            sorted(
                _strict_string(value, name="calibration donor") for value in self.calibration_donors
            )
        )
        if len(set(donors)) != len(donors) or donors != self.split_manifest.calibration:
            raise ValueError("calibration_donors must exactly equal the calibration split")
        object.__setattr__(self, "calibration_donors", donors)

        inputs = tuple(self.calibration_inputs)
        if not inputs or any(
            not isinstance(item, CalibrationTargetDecisionInput) for item in inputs
        ):
            raise ValueError("calibration_inputs must contain target decision inputs")
        inputs = tuple(sorted(inputs, key=lambda item: (item.level != "broad", item.parent)))
        targets = [item.target for item in inputs]
        if len(set(targets)) != len(targets):
            raise ValueError("calibration_inputs contains duplicate targets")
        object.__setattr__(self, "calibration_inputs", inputs)
        observed_donors = tuple(sorted({row.donor_id for item in inputs for row in item.rows}))
        if observed_donors != donors:
            raise ValueError("embedded calibration inputs must cover every calibration donor")
        row_count = _strict_positive_int(self.calibration_row_count, name="calibration_row_count")
        if row_count != sum(len(item.rows) for item in inputs):
            raise ValueError("calibration_row_count differs from embedded scientific rows")
        object.__setattr__(self, "calibration_row_count", row_count)
        expected_frame_fingerprint = _decision_inputs_fingerprint(inputs)
        if (
            _strict_sha256(self.calibration_frame_fingerprint, name="calibration_frame_fingerprint")
            != expected_frame_fingerprint
        ):
            raise ValueError("calibration_frame_fingerprint differs from embedded scientific rows")
        object.__setattr__(self, "calibration_frame_fingerprint", expected_frame_fingerprint)

        reviews = tuple(self.calibration_review_provenance)
        if not reviews or any(not isinstance(review, ReviewLabelProvenance) for review in reviews):
            raise ValueError("calibration_review_provenance must contain review manifests")
        reviews = tuple(
            sorted(reviews, key=lambda item: (item.level != "broad", str(item.parent or "")))
        )
        review_targets = [(review.level, str(review.parent or "")) for review in reviews]
        if review_targets != targets:
            raise ValueError("calibration review target set differs from embedded target inputs")
        for review, target_input in zip(reviews, inputs, strict=True):
            if review.split_name != "calibration" or not review.complete_evidence_sufficient:
                raise ValueError("threshold selection requires complete calibration reviews")
            if (
                review.source_run_id != self.split_manifest.source_run_id
                or review.split_manifest_content_id != self.split_manifest.content_id
            ):
                raise ValueError("calibration review refers to a different source/split")
            if review.assignment_run_id != self.assignment_run_id:
                raise ValueError("calibration review refers to a different assignment run")
            if (
                review.candidate_evaluation_manifest_content_id
                != self.calibration_candidate_manifest_content_id
                or review.candidate_assignment_output_content_id
                != self.calibration_candidate_output_content_id
            ):
                raise ValueError(
                    "calibration review refers to a different candidate assignment artifact"
                )
            if review.typing_model_bundle_content_id != self.preliminary_bundle_content_id:
                raise ValueError(
                    "calibration review must bind the preliminary full bundle content_id"
                )
            if review.sample_row_count != len(target_input.rows):
                raise ValueError("calibration review row count differs from embedded target rows")
            embedded_sample_fingerprint = review_sample_frame_fingerprint(
                _embedded_sampling_audit(target_input)
            )
            if embedded_sample_fingerprint != review.sample_frame_fingerprint:
                raise ValueError(
                    "embedded calibration sampling audit differs from review provenance"
                )
            embedded_label_fingerprint = review_label_frame_fingerprint(
                _embedded_review_ledger(target_input)
            )
            if embedded_label_fingerprint != review.fingerprint:
                raise ValueError(
                    "embedded calibration ledger differs from review provenance"
                )
        object.__setattr__(self, "calibration_review_provenance", reviews)

        probability_thresholds = _strict_thresholds(
            self.probability_thresholds, name="probability_thresholds"
        )
        margin_thresholds = _strict_thresholds(self.margin_thresholds, name="margin_thresholds")
        stability_thresholds = _strict_thresholds(
            self.stability_thresholds, name="stability_thresholds"
        )
        object.__setattr__(self, "probability_thresholds", probability_thresholds)
        object.__setattr__(self, "margin_thresholds", margin_thresholds)
        object.__setattr__(self, "stability_thresholds", stability_thresholds)
        if (
            self.split_manifest.threshold_selection_policy_content_id
            != self.policy.content_id
        ):
            raise ValueError(
                "threshold-selection policy differs from donor-split precommitment"
            )
        expected_grid_content_id = threshold_grid_content_id(
            probability_thresholds=probability_thresholds,
            margin_thresholds=margin_thresholds,
            stability_thresholds=stability_thresholds,
        )
        if self.split_manifest.threshold_grid_content_id != expected_grid_content_id:
            raise ValueError("threshold grid differs from donor-split precommitment")
        preliminary_floors = (
            (
                "probability",
                min(probability_thresholds),
                preliminary_thresholds["inferred_probability"],
            ),
            (
                "margin",
                min(margin_thresholds),
                preliminary_thresholds["inferred_margin"],
            ),
            (
                "stability",
                min(stability_thresholds),
                preliminary_thresholds["inferred_stability"],
            ),
        )
        more_permissive = [
            name
            for name, grid_minimum, preliminary_gate in preliminary_floors
            if grid_minimum < preliminary_gate
        ]
        if more_permissive:
            raise ValueError(
                "calibration grid cannot be more permissive than the preliminary "
                "assignment gates because subtype scores are materialized only after "
                "an accepted broad parent; violations="
                f"{more_permissive}"
            )
        derived_grid = _derive_operating_grid(
            inputs,
            probability_thresholds=probability_thresholds,
            margin_thresholds=margin_thresholds,
            stability_thresholds=stability_thresholds,
            policy=self.policy,
        )
        points = tuple(self.operating_grid)
        if not points or any(not isinstance(point, FrozenOperatingPoint) for point in points):
            raise ValueError("operating_grid must contain frozen operating points")
        points = tuple(
            sorted(
                points,
                key=lambda item: (
                    item.probability_threshold,
                    item.margin_threshold,
                    item.stability_threshold,
                ),
            )
        )
        supplied_grid = _grid_from_points(points)
        derived_points = tuple(
            FrozenOperatingPoint.from_series(row) for _, row in derived_grid.iterrows()
        )
        if [point.to_dict() for point in points] != [point.to_dict() for point in derived_points]:
            raise ValueError("embedded operating grid differs from calibration scientific inputs")
        object.__setattr__(self, "operating_grid", points)
        grid_count = _strict_positive_int(
            self.operating_grid_row_count, name="operating_grid_row_count"
        )
        if grid_count != len(points):
            raise ValueError("operating_grid_row_count differs from embedded grid")
        object.__setattr__(self, "operating_grid_row_count", grid_count)
        expected_grid_fingerprint = operating_grid_fingerprint(supplied_grid)
        if (
            _strict_sha256(self.operating_grid_fingerprint, name="operating_grid_fingerprint")
            != expected_grid_fingerprint
        ):
            raise ValueError("operating_grid_fingerprint differs from embedded grid")
        object.__setattr__(self, "operating_grid_fingerprint", expected_grid_fingerprint)

        selected = _select_targetwise_operating_point(
            inputs,
            derived_grid,
            probability_thresholds=probability_thresholds,
            margin_thresholds=margin_thresholds,
            stability_thresholds=stability_thresholds,
            policy=self.policy,
        )
        for field_name, column in (
            ("selected_probability_threshold", "probability_threshold"),
            ("selected_margin_threshold", "margin_threshold"),
            ("selected_stability_threshold", "stability_threshold"),
        ):
            observed = _unit_interval(getattr(self, field_name), name=field_name)
            expected = float(selected[column])
            if observed != expected:
                raise ValueError(f"{field_name} is not the deterministic selected grid member")
            object.__setattr__(self, field_name, observed)
        metrics = _selected_metrics(self.selected_metrics)
        expected_metrics = {name: float(selected[name]) for name in SELECTED_METRIC_COLUMNS}
        if dict(metrics) != expected_metrics:
            raise ValueError("selected_metrics differs from the deterministic selected grid row")
        object.__setattr__(self, "selected_metrics", metrics)
        expected_content_id = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected_content_id:
            raise ValueError("threshold-selection provenance content_id is invalid")
        object.__setattr__(self, "content_id", expected_content_id)

    @property
    def evaluated_targets(self) -> tuple[tuple[str, str], ...]:
        return tuple(item.target for item in self.calibration_inputs)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "selection_version": self.selection_version,
            "preliminary_bundle_content_id": self.preliminary_bundle_content_id,
            "model_fit_content_id": self.model_fit_content_id,
            "preliminary_assignment_thresholds": dict(
                self.preliminary_assignment_thresholds
            ),
            "split_manifest": self.split_manifest.to_dict(),
            "calibration_candidate_manifest_content_id": self.calibration_candidate_manifest_content_id,
            "calibration_candidate_output_content_id": self.calibration_candidate_output_content_id,
            "calibration_review_provenance": [
                review.to_dict() for review in self.calibration_review_provenance
            ],
            "assignment_run_id": self.assignment_run_id,
            "calibration_donors": list(self.calibration_donors),
            "calibration_inputs": [item.to_dict() for item in self.calibration_inputs],
            "calibration_row_count": self.calibration_row_count,
            "calibration_frame_fingerprint": self.calibration_frame_fingerprint,
            "probability_thresholds": list(self.probability_thresholds),
            "margin_thresholds": list(self.margin_thresholds),
            "stability_thresholds": list(self.stability_thresholds),
            "operating_grid": [point.to_dict() for point in self.operating_grid],
            "operating_grid_row_count": self.operating_grid_row_count,
            "operating_grid_fingerprint": self.operating_grid_fingerprint,
            "policy": self.policy.to_dict(),
            "selected_probability_threshold": self.selected_probability_threshold,
            "selected_margin_threshold": self.selected_margin_threshold,
            "selected_stability_threshold": self.selected_stability_threshold,
            "selected_metrics": dict(self.selected_metrics),
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> ThresholdSelectionProvenance:
        fields = {
            "schema_version",
            "selection_version",
            "preliminary_bundle_content_id",
            "model_fit_content_id",
            "preliminary_assignment_thresholds",
            "split_manifest",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
            "calibration_review_provenance",
            "assignment_run_id",
            "calibration_donors",
            "calibration_inputs",
            "calibration_row_count",
            "calibration_frame_fingerprint",
            "probability_thresholds",
            "margin_thresholds",
            "stability_thresholds",
            "operating_grid",
            "operating_grid_row_count",
            "operating_grid_fingerprint",
            "policy",
            "selected_probability_threshold",
            "selected_margin_threshold",
            "selected_stability_threshold",
            "selected_metrics",
            "content_id",
        }
        if not isinstance(value, Mapping):
            raise TypeError("threshold-selection provenance must be a mapping")
        _require_exact_fields(value, fields, where="threshold-selection provenance")
        for name in (
            "calibration_review_provenance",
            "calibration_donors",
            "calibration_inputs",
            "probability_thresholds",
            "margin_thresholds",
            "stability_thresholds",
            "operating_grid",
        ):
            if not isinstance(value[name], list):
                raise TypeError(f"{name} must be a JSON array")
        return cls(
            schema_version=value["schema_version"],
            selection_version=value["selection_version"],
            preliminary_bundle_content_id=value["preliminary_bundle_content_id"],
            model_fit_content_id=value["model_fit_content_id"],
            preliminary_assignment_thresholds=value[
                "preliminary_assignment_thresholds"
            ],
            split_manifest=DonorSplitManifest.from_dict(value["split_manifest"]),
            calibration_candidate_manifest_content_id=value[
                "calibration_candidate_manifest_content_id"
            ],
            calibration_candidate_output_content_id=value[
                "calibration_candidate_output_content_id"
            ],
            calibration_review_provenance=tuple(
                ReviewLabelProvenance.from_dict(item)
                for item in value["calibration_review_provenance"]
            ),
            assignment_run_id=value["assignment_run_id"],
            calibration_donors=tuple(value["calibration_donors"]),
            calibration_inputs=tuple(
                CalibrationTargetDecisionInput.from_dict(item)
                for item in value["calibration_inputs"]
            ),
            calibration_row_count=value["calibration_row_count"],
            calibration_frame_fingerprint=value["calibration_frame_fingerprint"],
            probability_thresholds=tuple(value["probability_thresholds"]),
            margin_thresholds=tuple(value["margin_thresholds"]),
            stability_thresholds=tuple(value["stability_thresholds"]),
            operating_grid=tuple(
                FrozenOperatingPoint.from_dict(item) for item in value["operating_grid"]
            ),
            operating_grid_row_count=value["operating_grid_row_count"],
            operating_grid_fingerprint=value["operating_grid_fingerprint"],
            policy=ThresholdSelectionPolicy.from_dict(value["policy"]),
            selected_probability_threshold=value["selected_probability_threshold"],
            selected_margin_threshold=value["selected_margin_threshold"],
            selected_stability_threshold=value["selected_stability_threshold"],
            selected_metrics=value["selected_metrics"],
            content_id=value["content_id"],
        )

    @classmethod
    def read_json(cls, path: Path | str) -> ThresholdSelectionProvenance:
        return cls.from_dict(_read_strict_json(path))


def _normalise_target_mapping(
    value: Mapping[tuple[str, str], pd.DataFrame], *, name: str
) -> dict[tuple[str, str], pd.DataFrame]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a target-keyed mapping")
    result: dict[tuple[str, str], pd.DataFrame] = {}
    for raw_target, frame in value.items():
        if not isinstance(raw_target, tuple) or len(raw_target) != 2:
            raise TypeError(f"{name} keys must be (level, parent) tuples")
        level = _strict_string(raw_target[0], name=f"{name} level").lower()
        parent = _strict_string(raw_target[1], name=f"{name} parent", allow_blank=True)
        target = (level, parent)
        if target in result:
            raise ValueError(f"{name} contains duplicate normalized target {target}")
        if not isinstance(frame, pd.DataFrame):
            raise TypeError(f"{name}[{target}] must be a pandas DataFrame")
        result[target] = frame
    return result


def _load_current_assignments(manifest: CandidateEvaluationManifest) -> pd.DataFrame:
    if manifest.input_artifact is not None:
        manifest.input_artifact.validate_current(mode="content")
    validate_dataset_snapshot_current(manifest.output, mode="content")
    current = snapshot_parquet_dataset(manifest.output.root, expected_donors=manifest.donors)
    expected_semantics = (
        manifest.output.total_rows,
        manifest.output.schema_sha256,
        manifest.output.object_id_sha256,
        manifest.output.content_sha256,
        tuple(
            (item.donor_id, item.row_count, item.object_id_sha256, item.content_sha256)
            for item in manifest.output.donors
        ),
    )
    observed_semantics = (
        current.total_rows,
        current.schema_sha256,
        current.object_id_sha256,
        current.content_sha256,
        tuple(
            (item.donor_id, item.row_count, item.object_id_sha256, item.content_sha256)
            for item in current.donors
        ),
    )
    if observed_semantics != expected_semantics:
        raise ValueError("candidate assignment output differs from its immutable snapshot")
    frames: list[pd.DataFrame] = []
    for donor in manifest.output.donors:
        for fingerprint in donor.files:
            frame = pd.read_parquet(fingerprint.path)
            if "donor_id" not in frame:
                frame["donor_id"] = donor.donor_id
            elif not frame["donor_id"].astype(str).eq(donor.donor_id).all():
                raise ValueError("candidate assignment file escapes its donor partition")
            frames.append(frame)
    assignments = pd.concat(frames, ignore_index=True)
    for column in ("donor_id", "object_id"):
        if column not in assignments:
            raise KeyError(f"candidate assignments are missing {column!r}")
        assignments[column] = assignments[column].astype("string").fillna("").str.strip()
    if (
        len(assignments) != manifest.output.total_rows
        or assignments[["donor_id", "object_id"]].eq("").any().any()
    ):
        raise ValueError("candidate assignment row/key universe differs from its snapshot")
    if assignments.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("candidate assignments contain duplicate cell keys")
    return assignments


def _explicit_bool(value: Any, *, name: str) -> bool:
    if pd.isna(value) or value not in (True, False, 0, 1):
        raise ValueError(f"{name} must be an explicit boolean")
    return bool(value)


def _ranking_probabilities(
    value: Any, classes: tuple[str, ...], *, allow_empty: bool
) -> tuple[float, ...]:
    if not isinstance(value, str):
        raise TypeError("ranked probabilities must be JSON strings")
    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as exc:
        raise ValueError("ranked probabilities contain invalid JSON") from exc
    if parsed == [] and allow_empty:
        return tuple(1.0 / len(classes) for _ in classes)
    if not isinstance(parsed, list) or len(parsed) != len(classes):
        raise ValueError("ranked probabilities do not cover every target class")
    observed: dict[str, float] = {}
    for entry in parsed:
        if not isinstance(entry, list) or len(entry) != 2 or not isinstance(entry[0], str):
            raise ValueError("ranked probability entries must be [class, probability]")
        name = entry[0]
        if name in observed:
            raise ValueError("ranked probabilities contain duplicate classes")
        observed[name] = _unit_interval(entry[1], name="ranked probability")
    if set(observed) != set(classes):
        raise ValueError("ranked probability classes differ from the preliminary bundle")
    probabilities = tuple(observed[name] for name in classes)
    if not np.isclose(sum(probabilities), 1.0, atol=1e-6, rtol=0.0):
        raise ValueError("ranked probabilities must sum to one")
    return probabilities


def _expected_target_classes(
    model_bundle: TypingModelBundle,
) -> dict[tuple[str, str], tuple[str, ...]]:
    result = {("broad", ""): tuple(model_bundle.broad.classes)}
    result.update(
        {("subtype", str(model.parent)): tuple(model.classes) for model in model_bundle.subtypes}
    )
    return result


def build_threshold_selection_provenance(
    *,
    selection_version: str,
    model_bundle: TypingModelBundle,
    calibration_candidate_manifest: CandidateEvaluationManifest,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    source_run_manifest: RunManifest,
    calibration_sample_frames: Mapping[tuple[str, str], pd.DataFrame],
    calibration_reviewer_frames: Mapping[tuple[str, str], pd.DataFrame],
    calibration_review_ledgers: Mapping[tuple[str, str], pd.DataFrame],
    calibration_review_provenance: Sequence[ReviewLabelProvenance],
    review_blinding_key: str | bytes,
    policy: ThresholdSelectionPolicy,
    probability_thresholds: Sequence[float] = DEFAULT_PROBABILITY_THRESHOLDS,
    margin_thresholds: Sequence[float] = DEFAULT_MARGIN_THRESHOLDS,
    stability_thresholds: Sequence[float] = DEFAULT_STABILITY_THRESHOLDS,
) -> ThresholdSelectionProvenance:
    """Derive and freeze thresholds from one exact calibration candidate run.

    There is intentionally no ``operating_points`` or ``selected_metrics``
    argument: all metrics are recomputed from the candidate assignments and
    exact blinded/adjudicated review inputs.
    """

    from .artifacts import RunManifest
    from .candidate_evaluation import CandidateEvaluationManifest
    from .hierarchical_typing import TypingRegistry
    from .marker_registry import MarkerRegistry
    from .typing_model_bundle import TypingModelBundle

    if not isinstance(model_bundle, TypingModelBundle):
        raise TypeError("model_bundle must be a preliminary TypingModelBundle")
    if not isinstance(calibration_candidate_manifest, CandidateEvaluationManifest):
        raise TypeError("calibration_candidate_manifest must be CandidateEvaluationManifest")
    if not isinstance(marker_registry, MarkerRegistry):
        raise TypeError("marker_registry must be MarkerRegistry")
    if not isinstance(typing_registry, TypingRegistry):
        raise TypeError("typing_registry must be TypingRegistry")
    if not isinstance(source_run_manifest, RunManifest):
        raise TypeError("source_run_manifest must be RunManifest")
    if not isinstance(policy, ThresholdSelectionPolicy):
        raise TypeError("policy must be ThresholdSelectionPolicy")
    probability_values = _strict_thresholds(
        probability_thresholds, name="probability_thresholds"
    )
    margin_values = _strict_thresholds(margin_thresholds, name="margin_thresholds")
    stability_values = _strict_thresholds(
        stability_thresholds, name="stability_thresholds"
    )
    splits = model_bundle.split_manifest
    if not splits.threshold_selection_policy_content_id:
        raise ValueError(
            "donor split must preregister threshold_selection_policy_content_id "
            "before calibration labels are opened"
        )
    if splits.threshold_selection_policy_content_id != policy.content_id:
        raise ValueError("threshold-selection policy differs from donor-split precommitment")
    grid_content_id = threshold_grid_content_id(
        probability_thresholds=probability_values,
        margin_thresholds=margin_values,
        stability_thresholds=stability_values,
    )
    if not splits.threshold_grid_content_id:
        raise ValueError(
            "donor split must preregister threshold_grid_content_id before "
            "calibration labels are opened"
        )
    if splits.threshold_grid_content_id != grid_content_id:
        raise ValueError("threshold grid differs from donor-split precommitment")
    if model_bundle.threshold_selection_provenance is not None:
        raise ValueError("threshold selection requires a preliminary/unfrozen model bundle")
    manifest = calibration_candidate_manifest
    if manifest.split_name != "calibration":
        raise ValueError("threshold selection requires a calibration candidate manifest")
    expected_manifest = {
        "model_bundle_content_id": model_bundle.content_id,
        "split_manifest_content_id": model_bundle.split_manifest.content_id,
        "source_run_id": model_bundle.split_manifest.source_run_id,
        "source_run_manifest_content_id": model_bundle.evidence_provenance.source_run_manifest_content_id,
        "source_stage_manifest_content_id": model_bundle.evidence_provenance.source_stage_manifest_content_id,
        "source_artifact_content_id": model_bundle.evidence_provenance.source_artifact_content_id,
        "marker_registry_fingerprint": model_bundle.marker_registry_fingerprint,
        "typing_rules_fingerprint": model_bundle.typing_rules_fingerprint,
    }
    for name, expected in expected_manifest.items():
        if getattr(manifest, name) != expected:
            raise ValueError(f"calibration candidate {name} differs from the preliminary bundle")
    if manifest.donors != model_bundle.split_manifest.calibration:
        raise ValueError("calibration candidate donors differ from the calibration split")
    manifest.validate_current(
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=False,
    )
    assignments = _load_current_assignments(manifest)
    metadata = {
        "typing_model_bundle_content_id": model_bundle.content_id,
        "typing_mode": "probabilistic_candidate",
        "typing_rules_fingerprint": model_bundle.typing_rules_fingerprint,
        "marker_registry_fingerprint": model_bundle.marker_registry_fingerprint,
        "typing_feature_schema_version": model_bundle.feature_schema_version,
        "typing_score_semantics": model_bundle.score_semantics,
    }
    for name, expected in metadata.items():
        if name not in assignments:
            raise KeyError(f"candidate assignments are missing metadata column {name!r}")
        observed = assignments[name].astype("string").fillna("").str.strip()
        if not observed.eq(str(expected)).all():
            raise ValueError(f"candidate assignment {name} differs from the preliminary bundle")

    expected_targets = _expected_target_classes(model_bundle)
    samples = _normalise_target_mapping(calibration_sample_frames, name="calibration_sample_frames")
    reviewer_frames = _normalise_target_mapping(
        calibration_reviewer_frames, name="calibration_reviewer_frames"
    )
    ledgers = _normalise_target_mapping(
        calibration_review_ledgers, name="calibration_review_ledgers"
    )
    if (
        set(samples) != set(expected_targets)
        or set(reviewer_frames) != set(expected_targets)
        or set(ledgers) != set(expected_targets)
    ):
        raise ValueError(
            "calibration sample/reviewer/ledger target set must equal broad plus "
            "every probabilistic subtype: "
            f"expected={sorted(expected_targets)}, samples={sorted(samples)}, "
            f"reviewer_frames={sorted(reviewer_frames)}, ledgers={sorted(ledgers)}"
        )
    reviews: dict[tuple[str, str], ReviewLabelProvenance] = {}
    for review in calibration_review_provenance:
        if not isinstance(review, ReviewLabelProvenance):
            raise TypeError("calibration_review_provenance must contain ReviewLabelProvenance")
        target = (review.level, str(review.parent or ""))
        if target in reviews:
            raise ValueError(f"duplicate calibration review target: {target}")
        reviews[target] = review
    if set(reviews) != set(expected_targets):
        raise ValueError("calibration review provenance target set is incomplete or unrelated")

    indexed_assignments = assignments.set_index(["donor_id", "object_id"])
    if not indexed_assignments.index.is_unique:  # defensive; loader already rejects this
        raise ValueError("candidate assignments contain duplicate cell keys")
    target_inputs: list[CalibrationTargetDecisionInput] = []
    observed_sample_donors: set[str] = set()
    ingest_references = tuple(
        reference
        for reference in source_run_manifest.stages
        if reference.stage == "ingest"
    )
    if len(ingest_references) != 1:
        raise ValueError(
            "threshold selection requires exactly one ingest stage in the bound "
            "source run"
        )
    authorized_review_context_id = ingest_references[0].output_content_sha256
    for target in sorted(expected_targets, key=lambda item: (item[0] != "broad", item[1])):
        level, parent = target
        review = reviews[target]
        if review.level != level or str(review.parent or "") != parent:
            raise ValueError("review provenance target metadata differs from its mapping key")
        if review.split_name != "calibration" or not review.complete_evidence_sufficient:
            raise ValueError("threshold selection requires complete calibration reviews")
        if review.typing_model_bundle_content_id != model_bundle.content_id:
            raise ValueError("calibration review must bind the preliminary full bundle content_id")
        if (
            review.candidate_evaluation_manifest_content_id != manifest.content_id
            or review.candidate_assignment_output_content_id
            != manifest.output_content_sha256
        ):
            raise ValueError(
                "calibration review must bind the exact candidate manifest and output"
            )
        if review.review_context_artifact_content_id != authorized_review_context_id:
            raise ValueError(
                "calibration review context is not the sole ingest output of the "
                "bound source run"
            )
        if set(samples[target].columns) != set(REVIEW_SAMPLE_AUDIT_COLUMNS):
            raise ValueError(
                "threshold selection requires the sealed sampling-audit columns "
                "only; reviewer locations belong in the separate blinded frame"
            )
        if (
            reviewer_frame_fingerprint(reviewer_frames[target])
            != review.reviewer_frame_fingerprint
        ):
            raise ValueError(
                "calibration reviewer frame differs from review provenance"
            )
        validate_review_sampling_against_assignments(
            samples[target],
            assignments,
            splits=model_bundle.split_manifest,
            split_name="calibration",
            candidate_assignment_run_id=manifest.assignment_run_id,
            level=level,
            parent=(None if level == "broad" else parent),
            blinding_key=review_blinding_key,
        )
        canonical_labels = validate_review_labels(
            ledgers[target],
            provenance=review,
            splits=model_bundle.split_manifest,
            sample_frame=samples[target],
            reviewer_frame=reviewer_frames[target],
            assignment_run_id=manifest.assignment_run_id,
            typing_model_bundle_content_id=model_bundle.content_id,
            candidate_evaluation_manifest_content_id=manifest.content_id,
            candidate_assignment_output_content_id=manifest.output_content_sha256,
        )
        sample = samples[target].copy()
        sample["donor_id"] = sample["donor_id"].astype(str)
        sample["object_id"] = sample["object_id"].astype(str)
        donors = set(sample["donor_id"])
        if not donors <= set(manifest.donors):
            raise ValueError("calibration sample contains donors outside the candidate manifest")
        observed_sample_donors.update(donors)
        keys = list(map(tuple, sample[["donor_id", "object_id"]].to_numpy()))
        missing_keys = sorted(set(keys) - set(indexed_assignments.index.tolist()))
        if missing_keys:
            raise ValueError(
                f"calibration sample contains unrelated assignment keys: {missing_keys[:5]}"
            )
        selected_assignments = indexed_assignments.loc[keys].reset_index()
        if list(map(tuple, selected_assignments[["donor_id", "object_id"]].to_numpy())) != keys:
            raise AssertionError("candidate assignment join changed the blinded sample order")
        aligned = (
            sample.loc[:, list(REVIEW_SAMPLE_AUDIT_COLUMNS)]
            .merge(
                canonical_labels.loc[:, list(REVIEW_LABEL_COLUMNS)],
                on=["donor_id", "object_id"],
                how="inner",
                validate="1:1",
            )
            .merge(selected_assignments, on=["donor_id", "object_id"], how="inner", validate="1:1")
        )
        if len(aligned) != len(sample):
            raise ValueError("sample, review ledger, and assignments do not have exact cell keys")
        prefix = "broad" if level == "broad" else "subtype"
        required_columns = {
            f"{prefix}_label",
            f"{prefix}_assignment_status",
            f"{prefix}_reason",
            f"{prefix}_best_candidate",
            f"{prefix}_stability",
            f"{prefix}_threshold_eligible",
            f"{prefix}_valid_evidence_pass",
            f"ranked_{prefix}_probabilities",
        }
        missing_columns = sorted(required_columns - set(aligned.columns))
        if missing_columns:
            raise KeyError(f"candidate assignments are missing target columns: {missing_columns}")
        if level == "subtype":
            for column, expected in (
                ("broad_label", parent),
                ("typing_subtype_model_source", "probabilistic_bundle"),
                ("typing_subtype_score_semantics", model_bundle.score_semantics),
            ):
                if column not in aligned or not aligned[column].astype(str).eq(expected).all():
                    raise ValueError(
                        f"subtype calibration rows are unrelated to probabilistic parent {parent!r}"
                    )

        classes = expected_targets[target]
        rows: list[CalibrationDecisionRow] = []
        for record in aligned.to_dict(orient="records"):
            status = _strict_string(
                record[f"{prefix}_assignment_status"], name="assignment status"
            ).lower()
            raw_label = _strict_string(record[f"{prefix}_label"], name="assignment label")
            current = _DISPOSITION_LABEL.get(status, raw_label)
            eligible = _explicit_bool(
                record[f"{prefix}_threshold_eligible"], name="threshold_eligible"
            )
            probabilities = _ranking_probabilities(
                record[f"ranked_{prefix}_probabilities"],
                classes,
                allow_empty=(not eligible and status == "unavailable"),
            )
            stability_value = pd.to_numeric(
                pd.Series([record[f"{prefix}_stability"]]), errors="coerce"
            ).iloc[0]
            stability = None if pd.isna(stability_value) else float(stability_value)
            rows.append(
                CalibrationDecisionRow(
                    donor_id=str(record["donor_id"]),
                    object_id=str(record["object_id"]),
                    reference_label=str(record["reference_label"]),
                    reviewer_1=str(record["reviewer_1"]),
                    reviewer_2=str(record["reviewer_2"]),
                    adjudicated=_explicit_bool(
                        record["adjudicated"], name="adjudicated"
                    ),
                    root_cause=str(record["root_cause"]),
                    evidence_sufficient=_explicit_bool(
                        record["evidence_sufficient"], name="evidence_sufficient"
                    ),
                    sampling_stratum=str(record["sampling_stratum"]),
                    sampling_method=str(record["sampling_method"]),
                    sampling_salt=str(record["sampling_salt"]),
                    stratum_population=int(record["stratum_population"]),
                    stratum_sampled=int(record["stratum_sampled"]),
                    sampling_weight=record["sampling_weight"],
                    selection_probability=record["selection_probability"],
                    current_prediction=current,
                    current_assignment_status=status,
                    current_reason=str(record[f"{prefix}_reason"]),
                    candidate_prediction=str(record[f"{prefix}_best_candidate"] or ""),
                    threshold_eligible=eligible,
                    valid_evidence=_explicit_bool(
                        record[f"{prefix}_valid_evidence_pass"], name="valid_evidence_pass"
                    ),
                    stability=stability,
                    probabilities=probabilities,
                )
            )
        target_inputs.append(
            CalibrationTargetDecisionInput(
                level=level, parent=parent, classes=classes, rows=tuple(rows)
            )
        )
    if tuple(sorted(observed_sample_donors)) != manifest.donors:
        raise ValueError(
            "calibration review samples must collectively cover every calibration donor"
        )

    preliminary_floors = (
        (
            "probability",
            min(probability_values),
            float(model_bundle.thresholds.inferred_probability),
        ),
        ("margin", min(margin_values), float(model_bundle.thresholds.inferred_margin)),
        (
            "stability",
            min(stability_values),
            float(model_bundle.thresholds.inferred_stability),
        ),
    )
    more_permissive = [
        name
        for name, grid_minimum, preliminary_gate in preliminary_floors
        if grid_minimum < preliminary_gate
    ]
    if more_permissive:
        raise ValueError(
            "calibration grid cannot be more permissive than the preliminary "
            "assignment gates because subtype scores are materialized only after "
            "an accepted broad parent; refit the preliminary bundle with the grid "
            f"minima (violations={more_permissive})"
        )
    grid = _derive_operating_grid(
        target_inputs,
        probability_thresholds=probability_values,
        margin_thresholds=margin_values,
        stability_thresholds=stability_values,
        policy=policy,
    )
    selected = _select_targetwise_operating_point(
        target_inputs,
        grid,
        probability_thresholds=probability_values,
        margin_thresholds=margin_values,
        stability_thresholds=stability_values,
        policy=policy,
    )
    # Close the check/use window before sealing provenance.  This repeats the
    # exact source-stage universe, assignment schema, and content validation.
    manifest.validate_current(
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=False,
    )
    frozen_grid = tuple(FrozenOperatingPoint.from_series(row) for _, row in grid.iterrows())
    return ThresholdSelectionProvenance(
        selection_version=selection_version,
        preliminary_bundle_content_id=model_bundle.content_id,
        model_fit_content_id=model_bundle.model_fit_content_id,
        preliminary_assignment_thresholds={
            name: float(getattr(model_bundle.thresholds, name))
            for name in _ASSIGNMENT_THRESHOLD_FIELDS
        },
        split_manifest=model_bundle.split_manifest,
        calibration_candidate_manifest_content_id=manifest.content_id,
        calibration_candidate_output_content_id=manifest.output_content_sha256,
        calibration_review_provenance=tuple(reviews.values()),
        assignment_run_id=manifest.assignment_run_id,
        calibration_donors=manifest.donors,
        calibration_inputs=tuple(target_inputs),
        calibration_row_count=sum(len(item.rows) for item in target_inputs),
        calibration_frame_fingerprint=_decision_inputs_fingerprint(target_inputs),
        probability_thresholds=probability_values,
        margin_thresholds=margin_values,
        stability_thresholds=stability_values,
        operating_grid=frozen_grid,
        operating_grid_row_count=len(frozen_grid),
        operating_grid_fingerprint=operating_grid_fingerprint(grid),
        policy=policy,
        selected_probability_threshold=float(selected["probability_threshold"]),
        selected_margin_threshold=float(selected["margin_threshold"]),
        selected_stability_threshold=float(selected["stability_threshold"]),
        selected_metrics={name: float(selected[name]) for name in SELECTED_METRIC_COLUMNS},
    )


@dataclass(frozen=True)
class ThresholdSelectionReplayInputs:
    """Protected, in-memory inputs for full calibration source replay.

    This object deliberately has no JSON representation.  In particular, the
    review blinding key must remain outside threshold, model, evaluation, and
    release artifacts.  The mappings retain the exact sealed audit frames,
    reviewer frames, and adjudicated ledgers used for the original decision;
    mutating a retained frame makes the subsequent replay fail closed.
    """

    preliminary_bundle: TypingModelBundle = field(repr=False)
    calibration_candidate_manifest: CandidateEvaluationManifest = field(repr=False)
    marker_registry: MarkerRegistry = field(repr=False)
    typing_registry: TypingRegistry = field(repr=False)
    source_run_manifest: RunManifest = field(repr=False)
    calibration_sample_frames: Mapping[tuple[str, str], pd.DataFrame] = field(
        repr=False
    )
    calibration_reviewer_frames: Mapping[tuple[str, str], pd.DataFrame] = field(
        repr=False
    )
    calibration_review_ledgers: Mapping[tuple[str, str], pd.DataFrame] = field(
        repr=False
    )
    review_blinding_key: str | bytes = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        from .artifacts import RunManifest
        from .candidate_evaluation import CandidateEvaluationManifest
        from .hierarchical_typing import TypingRegistry
        from .marker_registry import MarkerRegistry
        from .typing_model_bundle import TypingModelBundle

        expected_types = (
            ("preliminary_bundle", self.preliminary_bundle, TypingModelBundle),
            (
                "calibration_candidate_manifest",
                self.calibration_candidate_manifest,
                CandidateEvaluationManifest,
            ),
            ("marker_registry", self.marker_registry, MarkerRegistry),
            ("typing_registry", self.typing_registry, TypingRegistry),
            ("source_run_manifest", self.source_run_manifest, RunManifest),
        )
        for name, value, expected_type in expected_types:
            if not isinstance(value, expected_type):
                raise TypeError(f"{name} must be {expected_type.__name__}")
        if self.preliminary_bundle.threshold_selection_provenance is not None:
            raise ValueError(
                "threshold replay inputs require the original preliminary/unfrozen bundle"
            )
        for name in (
            "calibration_sample_frames",
            "calibration_reviewer_frames",
            "calibration_review_ledgers",
        ):
            normalized = _normalise_target_mapping(getattr(self, name), name=name)
            object.__setattr__(self, name, MappingProxyType(normalized))
        key = self.review_blinding_key
        if isinstance(key, str):
            protected_key: str | bytes = key
        elif isinstance(key, bytes):
            protected_key = bytes(key)
        else:
            raise TypeError("review_blinding_key must be str or bytes")
        object.__setattr__(self, "review_blinding_key", protected_key)


_THRESHOLD_REPLAY_PROOF_TOKEN = object()


@dataclass(frozen=True)
class ThresholdSelectionSourceReplayProof:
    """In-memory proof emitted only after exact threshold source replay.

    The proof deliberately has no serialization API.  Its guarded constructor
    prevents ordinary callers from manufacturing a release-authorization token
    merely by copying content identifiers out of persisted artifacts.
    """

    threshold_selection_provenance_content_id: str
    preliminary_bundle_content_id: str
    model_fit_content_id: str
    calibration_candidate_manifest_content_id: str
    calibration_candidate_output_content_id: str
    assignment_run_id: str
    source_run_manifest_content_id: str
    source_stage_manifest_content_id: str
    source_artifact_content_id: str
    review_context_artifact_content_id: str
    _verification_token: InitVar[object] = field(repr=False)

    def __post_init__(self, _verification_token: object) -> None:
        if _verification_token is not _THRESHOLD_REPLAY_PROOF_TOKEN:
            raise TypeError(
                "ThresholdSelectionSourceReplayProof is created only by successful "
                "threshold source replay"
            )
        for name in (
            "threshold_selection_provenance_content_id",
            "preliminary_bundle_content_id",
            "model_fit_content_id",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
            "assignment_run_id",
            "source_run_manifest_content_id",
            "source_stage_manifest_content_id",
            "source_artifact_content_id",
            "review_context_artifact_content_id",
        ):
            object.__setattr__(
                self, name, _strict_sha256(getattr(self, name), name=name)
            )


def _validate_replayed_reviewer_context(
    replay_inputs: ThresholdSelectionReplayInputs,
) -> str:
    """Rebuild every reviewer frame from the sole immutable ingest output."""

    samples = replay_inputs.calibration_sample_frames
    reviewer_frames = replay_inputs.calibration_reviewer_frames
    if set(samples) != set(reviewer_frames):
        raise ValueError(
            "calibration sample and reviewer-frame target sets differ during source replay"
        )
    if not samples:
        raise ValueError("calibration source replay requires reviewer targets")
    context_ids: set[str] = set()
    for target in sorted(samples):
        observed = reviewer_frames[target]
        location_columns = tuple(
            column
            for column in observed.columns
            if column not in REVIEW_SAMPLE_BLINDED_COLUMNS
        )
        replayed = blinded_review_sample_from_dataset(
            samples[target],
            source_run_manifest=replay_inputs.source_run_manifest,
            splits=replay_inputs.preliminary_bundle.split_manifest,
            location_columns=location_columns,
        )
        observed_context = reviewer_frame_context_artifact_content_id(observed)
        replayed_context = reviewer_frame_context_artifact_content_id(replayed)
        context_ids.add(replayed_context)
        if observed_context != replayed_context:
            raise ValueError(
                f"calibration reviewer context for target {target} differs from "
                "the sole ingest artifact"
            )
        if reviewer_frame_fingerprint(observed) != reviewer_frame_fingerprint(replayed):
            raise ValueError(
                f"calibration reviewer frame for target {target} differs from "
                "source replay of the sole ingest artifact"
            )
    if len(context_ids) != 1:  # defensive: one source run has one authorized ingest
        raise ValueError("calibration reviewer targets use different ingest artifacts")
    return next(iter(context_ids))


def validate_threshold_selection_source_replay(
    frozen_bundle: TypingModelBundle,
    replay_inputs: ThresholdSelectionReplayInputs,
) -> ThresholdSelectionSourceReplayProof:
    """Rebuild a frozen threshold decision from its protected source inputs.

    Replay always runs candidate inference and calibration against the original
    preliminary bundle.  It never refits the model and never applies the frozen
    operating thresholds while reconstructing the calibration candidate.
    Exact provenance content identity is required before the frozen bundle's
    selected thresholds and scientific-fit linkage are accepted.
    """

    from .typing_model_bundle import TypingModelBundle

    if not isinstance(frozen_bundle, TypingModelBundle):
        raise TypeError("frozen_bundle must be TypingModelBundle")
    if not isinstance(replay_inputs, ThresholdSelectionReplayInputs):
        raise TypeError("replay_inputs must be ThresholdSelectionReplayInputs")
    provenance = frozen_bundle.threshold_selection_provenance
    if provenance is None:
        raise ValueError("threshold source replay requires a frozen model bundle")
    if (
        not isinstance(provenance, ThresholdSelectionProvenance)
        or provenance.schema_version != PROVENANCE_SCHEMA_VERSION
    ):
        raise ValueError(
            f"threshold source replay requires schema-v{PROVENANCE_SCHEMA_VERSION} provenance"
        )

    preliminary = replay_inputs.preliminary_bundle
    candidate = replay_inputs.calibration_candidate_manifest
    # Reconstructing both bundles from their public payloads rechecks their
    # content identities and all internal linkage without fitting a model.
    TypingModelBundle.from_dict(preliminary.to_dict())
    TypingModelBundle.from_dict(frozen_bundle.to_dict())

    if provenance.preliminary_bundle_content_id != preliminary.content_id:
        raise ValueError(
            "frozen threshold provenance refers to a different preliminary bundle"
        )
    if not (
        preliminary.model_fit_content_id
        == provenance.model_fit_content_id
        == frozen_bundle.model_fit_content_id
    ):
        raise ValueError("preliminary, provenance, and frozen model-fit identities differ")
    if not (
        preliminary.split_manifest.content_id
        == provenance.split_manifest.content_id
        == frozen_bundle.split_manifest.content_id
    ):
        raise ValueError("preliminary, provenance, and frozen donor-split identities differ")
    if preliminary.evidence_provenance.content_id != (
        frozen_bundle.evidence_provenance.content_id
    ):
        raise ValueError("preliminary and frozen evidence identities differ")
    for name in ("marker_registry_fingerprint", "typing_rules_fingerprint"):
        if getattr(preliminary, name) != getattr(frozen_bundle, name):
            raise ValueError(f"preliminary and frozen {name} values differ")
    if replay_inputs.source_run_manifest.content_id != (
        preliminary.evidence_provenance.source_run_manifest_content_id
    ):
        raise ValueError("replay source run manifest differs from preliminary evidence")
    if (
        candidate.content_id
        != provenance.calibration_candidate_manifest_content_id
        or candidate.output_content_sha256
        != provenance.calibration_candidate_output_content_id
        or candidate.assignment_run_id != provenance.assignment_run_id
    ):
        raise ValueError(
            "replay candidate manifest, output, or assignment run differs from "
            "frozen threshold provenance"
        )

    preliminary_thresholds = {
        name: float(getattr(preliminary.thresholds, name))
        for name in _ASSIGNMENT_THRESHOLD_FIELDS
    }
    if preliminary_thresholds != dict(provenance.preliminary_assignment_thresholds):
        raise ValueError(
            "replay preliminary thresholds differ from frozen threshold provenance"
        )

    review_context_artifact_content_id = _validate_replayed_reviewer_context(
        replay_inputs
    )
    replayed = build_threshold_selection_provenance(
        selection_version=provenance.selection_version,
        model_bundle=preliminary,
        calibration_candidate_manifest=candidate,
        marker_registry=replay_inputs.marker_registry,
        typing_registry=replay_inputs.typing_registry,
        source_run_manifest=replay_inputs.source_run_manifest,
        calibration_sample_frames=replay_inputs.calibration_sample_frames,
        calibration_reviewer_frames=replay_inputs.calibration_reviewer_frames,
        calibration_review_ledgers=replay_inputs.calibration_review_ledgers,
        calibration_review_provenance=provenance.calibration_review_provenance,
        review_blinding_key=replay_inputs.review_blinding_key,
        policy=provenance.policy,
        probability_thresholds=provenance.probability_thresholds,
        margin_thresholds=provenance.margin_thresholds,
        stability_thresholds=provenance.stability_thresholds,
    )
    if replayed.content_id != provenance.content_id:
        raise ValueError(
            "threshold selection source replay differs from frozen provenance"
        )

    selected_thresholds = {
        "inferred_probability": replayed.selected_probability_threshold,
        "inferred_margin": replayed.selected_margin_threshold,
        "inferred_stability": replayed.selected_stability_threshold,
        "anchor_probability": preliminary_thresholds["anchor_probability"],
        "negative_probability": preliminary_thresholds["negative_probability"],
    }
    observed_thresholds = {
        name: float(getattr(frozen_bundle.thresholds, name))
        for name in _ASSIGNMENT_THRESHOLD_FIELDS
    }
    if observed_thresholds != selected_thresholds:
        raise ValueError(
            "frozen bundle thresholds differ from the source-replayed decision"
        )
    if frozen_bundle.threshold_selection_provenance.content_id != replayed.content_id:
        raise ValueError("frozen bundle is linked to a different threshold decision")
    return ThresholdSelectionSourceReplayProof(
        threshold_selection_provenance_content_id=replayed.content_id,
        preliminary_bundle_content_id=preliminary.content_id,
        model_fit_content_id=preliminary.model_fit_content_id,
        calibration_candidate_manifest_content_id=candidate.content_id,
        calibration_candidate_output_content_id=candidate.output_content_sha256,
        assignment_run_id=candidate.assignment_run_id,
        source_run_manifest_content_id=(
            preliminary.evidence_provenance.source_run_manifest_content_id
        ),
        source_stage_manifest_content_id=(
            preliminary.evidence_provenance.source_stage_manifest_content_id
        ),
        source_artifact_content_id=(
            preliminary.evidence_provenance.source_artifact_content_id
        ),
        review_context_artifact_content_id=review_context_artifact_content_id,
        _verification_token=_THRESHOLD_REPLAY_PROOF_TOKEN,
    )


__all__ = [
    "CALIBRATION_FRAME_KEY_COLUMNS",
    "DEFAULT_MARGIN_THRESHOLDS",
    "DEFAULT_PROBABILITY_THRESHOLDS",
    "DEFAULT_STABILITY_THRESHOLDS",
    "OPERATING_POINT_COLUMNS",
    "POLICY_SCHEMA_VERSION",
    "PROVENANCE_SCHEMA_VERSION",
    "SELECTED_METRIC_COLUMNS",
    "CalibrationDecisionRow",
    "CalibrationTargetDecisionInput",
    "FrozenOperatingPoint",
    "ThresholdSelectionPolicy",
    "ThresholdSelectionProvenance",
    "ThresholdSelectionReplayInputs",
    "ThresholdSelectionSourceReplayProof",
    "build_threshold_selection_provenance",
    "calibration_frame_fingerprint",
    "operating_grid_fingerprint",
    "threshold_grid_content_id",
    "validate_threshold_selection_source_replay",
]
