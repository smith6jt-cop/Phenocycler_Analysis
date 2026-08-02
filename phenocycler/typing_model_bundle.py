"""Strict JSON model bundles for the optional probabilistic rescue lane.

Production typing remains rules-only unless one of these bundles is explicitly
configured.  A configured bundle is an immutable type-stage input: schema,
content identity, registry compatibility, coefficient constraints, and every
bootstrap member are validated before any assignment output is written.
"""

from __future__ import annotations

import hashlib
import hmac
import importlib.metadata
import json
import os
import platform
import re
import sys
import sysconfig
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field, replace
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .artifacts import (
    CodeMetadata,
    InputArtifact,
    RunManifest,
    hash_identifiers,
    sha256_file,
    sha256_json,
    validate_dataset_snapshot_current,
)
from .hierarchical_typing import (
    AMBIGUOUS,
    AUTHORITATIVE_ANCHOR,
    FEATURE_SCHEMA_VERSION,
    INFERRED,
    OTHER,
    UNAVAILABLE,
    AdditiveMultinomialModel,
    AssignmentThresholds,
    TypingRegistry,
    classifier_from_registry,
    evidence_feature_matrix,
)
from .marker_registry import MarkerRegistry
from .refinement_contracts import (
    REVIEW_LABEL_COLUMNS,
    REVIEW_SAMPLING_METHOD,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
    ReviewLabelProvenance,
    blinded_review_sample_from_dataset,
    review_label_frame_fingerprint,
    review_sample_frame_fingerprint,
    reviewer_frame_fingerprint,
)
from .threshold_selection import (
    ThresholdSelectionProvenance,
    ThresholdSelectionReplayInputs,
    ThresholdSelectionSourceReplayProof,
    validate_threshold_selection_source_replay,
)

BUNDLE_SCHEMA_VERSION = 2
BUNDLE_METHOD_VERSION = "probabilistic-rescue-bundle-v1"
PROMOTION_SCHEMA_VERSION = 3
PROMOTION_COMPONENT_SCHEMA_VERSION = 3
EVALUATION_REPORT_METHOD_VERSION = "probabilistic-typing-heldout-evaluation-v3"
EVALUATION_BOOTSTRAP_UNIT = "donor_then_stratified_cell"
EVALUATION_MULTIPLICITY_METHOD = (
    "studywise_two_report_bonferroni_v1"
    "+target_supnorm_half_alpha_bonferroni_targets_finite_b_v1"
    "+cluster_kish_cp_boundary_guard_v1"
)
MINIMUM_SIMULTANEOUS_TAIL_DRAWS = 20.0
PROMOTION_REPORT_FAMILY_COUNT = 2
PROMOTION_REPLAY_RECEIPT_METHOD_VERSION = "probabilistic-promotion-replay-receipt-v2"
PROMOTION_DECISION = "approved_for_production"
PRODUCTION_MAX_ABS_WEIGHT = 8.0
MINIMUM_VALID_BOOTSTRAP_FRACTION = 0.80
TEMPERATURE_MIN = 0.05
TEMPERATURE_MAX = 20.0
SCORE_SEMANTICS = "donor_class_balanced_temperature_scaled_softmax_probability"

_FIT_PARAMETER_FIELDS = {
    "n_bootstrap",
    "random_state",
    "l2",
    "max_abs_weight",
    "maxiter",
    "require_all_subtype_parents",
    "structural_rules_only_parents",
    "missing_probabilistic_subtype_parents",
    "final_specific_classes",
    "minimum_valid_bootstrap_fraction",
    "class_priors",
    "donor_weighting",
}

_PROMOTION_REPLAY_SOURCE_PATHS = (
    "phenocycler",
    "pyproject.toml",
    "environment.yml",
    "environment-integration.yml",
)
_PROMOTION_RUNTIME_DISTRIBUTIONS = ("numpy", "pandas", "pyarrow", "scipy")


def promotion_replay_code_sha256() -> str:
    """Hash the complete package and declared environment specifications."""

    repository = Path(__file__).resolve().parents[1]
    return CodeMetadata.capture(
        repository,
        source_paths=_PROMOTION_REPLAY_SOURCE_PATHS,
    ).source_sha256


@lru_cache(maxsize=1)
def promotion_replay_environment_sha256() -> str:
    """Hash the numerical runtime and installed files used by release replay.

    Distribution versions alone cannot distinguish rebuilt wheels or changed
    compiled kernels.  Promotion is rare, so this deliberately hashes the
    installed bytes for the numerical/dataframe/Parquet stack and the Python
    executable.  The result is cached for the process lifetime.  Deployments
    may additionally bind an immutable image or solved-environment digest via
    ``PHENOCYCLER_RUNTIME_IMAGE_DIGEST``.
    """

    distributions: list[dict[str, Any]] = []
    for name in _PROMOTION_RUNTIME_DISTRIBUTIONS:
        try:
            distribution = importlib.metadata.distribution(name)
        except importlib.metadata.PackageNotFoundError as exc:
            raise RuntimeError(
                f"promotion replay runtime is missing required distribution {name!r}"
            ) from exc
        files = distribution.files
        if not files:
            raise RuntimeError(
                f"installed distribution {name!r} has no auditable file inventory"
            )
        entries: list[tuple[str, int, str]] = []
        for item in sorted(files, key=lambda value: str(value)):
            path = Path(distribution.locate_file(item)).resolve()
            if "__pycache__" in path.parts:
                continue
            if not path.is_file():
                raise RuntimeError(
                    f"installed distribution {name!r} is missing declared file {item!s}"
                )
            entries.append((str(item), path.stat().st_size, sha256_file(path)))
        distributions.append(
            {
                "name": name,
                "version": distribution.version,
                "files": entries,
            }
        )

    executable = Path(sys.executable).resolve()
    image_digest = os.environ.get("PHENOCYCLER_RUNTIME_IMAGE_DIGEST", "").strip()
    if image_digest and re.fullmatch(r"[0-9a-f]{64}", image_digest) is None:
        raise ValueError(
            "PHENOCYCLER_RUNTIME_IMAGE_DIGEST must be a lowercase SHA-256 digest"
        )
    library_dir = sysconfig.get_config_var("LIBDIR")
    library_name = sysconfig.get_config_var("LDLIBRARY")
    runtime_library = (
        Path(library_dir, library_name).resolve()
        if library_dir and library_name
        else None
    )
    payload = {
        "python": {
            "implementation": sys.implementation.name,
            "cache_tag": sys.implementation.cache_tag,
            "version": list(sys.version_info),
            "soabi": sysconfig.get_config_var("SOABI"),
            "executable_sha256": sha256_file(executable),
            "runtime_library_sha256": (
                sha256_file(runtime_library)
                if runtime_library is not None and runtime_library.is_file()
                else ""
            ),
        },
        "platform": {
            "system": platform.system(),
            "release": platform.release(),
            "machine": platform.machine(),
            "libc": list(platform.libc_ver()),
        },
        "runtime_image_digest": image_digest,
        "distributions": distributions,
    }
    return sha256_json(payload)


def _reject_unknown(value: Mapping[str, Any], allowed: set[str], *, where: str) -> None:
    unknown = sorted(set(value) - allowed)
    if unknown:
        raise ValueError(f"{where} contains unknown fields: {unknown}")


def _strict_bool(value: Any, *, name: str) -> bool:
    if not isinstance(value, bool):
        raise TypeError(f"{name} must be a JSON boolean")
    return value


def _promotion_report_confidence_level(policy_confidence_level: float) -> float:
    """Allocate study-wise alpha across validation and locked-test reports."""

    confidence = _strict_unit_float(
        policy_confidence_level, name="policy.confidence_level"
    )
    return 1.0 - (1.0 - confidence) / PROMOTION_REPORT_FAMILY_COUNT


def _require_sha256(value: Any, *, name: str) -> str:
    text = str(value)
    if re.fullmatch(r"[0-9a-f]{64}", text) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return text


def _require_exact_fields(value: Mapping[str, Any], fields: set[str], *, where: str) -> None:
    _reject_unknown(value, fields, where=where)
    missing = sorted(fields - set(value))
    if missing:
        raise ValueError(f"{where} is missing required fields: {missing}")


def _strict_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"JSON contains duplicate key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"JSON contains non-finite number {value!r}")


def _read_strict_json(path: Path) -> Mapping[str, Any]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(
            handle,
            object_pairs_hook=_strict_json_object,
            parse_constant=_reject_json_constant,
        )
    if not isinstance(value, Mapping):
        raise TypeError("JSON artifact must contain one object")
    return value


def _strict_finite_float(value: Any, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError(f"{name} must be a JSON number")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _strict_unit_float(value: Any, *, name: str) -> float:
    result = _strict_finite_float(value, name=name)
    if result < 0 or result > 1:
        raise ValueError(f"{name} must lie in [0, 1]")
    return result


def _strict_nonnegative_float(value: Any, *, name: str) -> float:
    result = _strict_finite_float(value, name=name)
    if result < 0:
        raise ValueError(f"{name} must be >= 0")
    return result


def _strict_positive_int(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise TypeError(f"{name} must be a positive JSON integer")
    return value


def _validated_fit_parameters(value: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError("fit_parameters must be a mapping")
    _reject_unknown(value, _FIT_PARAMETER_FIELDS, where="fit_parameters")
    missing = sorted(_FIT_PARAMETER_FIELDS - set(value))
    if missing:
        raise ValueError(f"fit_parameters is missing required fields: {missing}")

    n_bootstrap = _strict_positive_int(value["n_bootstrap"], name="fit_parameters.n_bootstrap")
    random_state = value["random_state"]
    if random_state is not None and (
        isinstance(random_state, bool) or not isinstance(random_state, int)
    ):
        raise TypeError("fit_parameters.random_state must be a JSON integer or null")
    l2 = float(value["l2"])
    if not np.isfinite(l2) or l2 < 0:
        raise ValueError("fit_parameters.l2 must be finite and >= 0")
    max_abs_weight = float(value["max_abs_weight"])
    if (
        not np.isfinite(max_abs_weight)
        or max_abs_weight <= 0
        or max_abs_weight > PRODUCTION_MAX_ABS_WEIGHT
    ):
        raise ValueError(
            "fit_parameters.max_abs_weight must be finite, > 0, and <= "
            f"{PRODUCTION_MAX_ABS_WEIGHT:g}"
        )
    maxiter = _strict_positive_int(value["maxiter"], name="fit_parameters.maxiter")
    require_all = _strict_bool(
        value["require_all_subtype_parents"],
        name="fit_parameters.require_all_subtype_parents",
    )
    def _parent_array(name: str) -> list[str]:
        raw = value[name]
        if isinstance(raw, (str, bytes)) or not isinstance(raw, Sequence):
            raise TypeError(f"fit_parameters.{name} must be a JSON array")
        result: list[str] = []
        for item in raw:
            if not isinstance(item, str) or not item.strip():
                raise TypeError(
                    f"fit_parameters.{name} must contain non-empty strings"
                )
            result.append(item.strip())
        if len(set(result)) != len(result):
            raise ValueError(f"fit_parameters.{name} contains duplicates")
        return sorted(result)

    structural = _parent_array("structural_rules_only_parents")
    missing_probabilistic = _parent_array(
        "missing_probabilistic_subtype_parents"
    )
    final_specific_classes = _parent_array("final_specific_classes")
    if len(final_specific_classes) < 2:
        raise ValueError(
            "fit_parameters.final_specific_classes must contain at least two classes"
        )
    overlap = sorted(set(structural) & set(missing_probabilistic))
    if overlap:
        raise ValueError(
            "fit_parameters structural and missing probabilistic subtype parents "
            f"overlap: {overlap}"
        )
    valid_fraction = float(value["minimum_valid_bootstrap_fraction"])
    if valid_fraction != MINIMUM_VALID_BOOTSTRAP_FRACTION:
        raise ValueError(
            "fit_parameters.minimum_valid_bootstrap_fraction must equal "
            f"{MINIMUM_VALID_BOOTSTRAP_FRACTION:g}"
        )
    expected_strings = {
        "class_priors": "uniform_no_intercept",
        "donor_weighting": "class_then_donor_balanced",
    }
    for name, expected in expected_strings.items():
        if value[name] != expected:
            raise ValueError(f"fit_parameters.{name} must equal {expected!r}")
    return {
        "n_bootstrap": n_bootstrap,
        "random_state": random_state,
        "l2": l2,
        "max_abs_weight": max_abs_weight,
        "maxiter": maxiter,
        "require_all_subtype_parents": require_all,
        "structural_rules_only_parents": structural,
        "missing_probabilistic_subtype_parents": missing_probabilistic,
        "final_specific_classes": final_specific_classes,
        "minimum_valid_bootstrap_fraction": valid_fraction,
        **expected_strings,
    }


def _finite_matrix(
    value: Any,
    *,
    rows: int,
    columns: int,
    name: str,
) -> tuple[tuple[float, ...], ...]:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be numeric") from exc
    if array.shape != (rows, columns):
        raise ValueError(f"{name} shape {array.shape} does not match {(rows, columns)}")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains non-finite coefficients")
    if np.max(np.abs(array), initial=0.0) > PRODUCTION_MAX_ABS_WEIGHT:
        raise ValueError(
            f"{name} exceeds the production coefficient cap {PRODUCTION_MAX_ABS_WEIGHT:g}"
        )
    return tuple(tuple(float(item) for item in row) for row in array)


def _atomic_write_json(path: Path, payload: Mapping[str, Any]) -> Path:
    path = path.expanduser().resolve()
    if path.exists():
        raise FileExistsError(f"immutable JSON artifact already exists: {path}")
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
        if temporary.exists():
            temporary.unlink()
    return path


@dataclass(frozen=True)
class FrozenModelEnsemble:
    """One calibrated base model plus all whole-donor bootstrap refits."""

    level: str
    parent: str | None
    classes: tuple[str, ...]
    markers: tuple[str, ...]
    coefficients: tuple[tuple[float, ...], ...]
    bootstrap_coefficients: tuple[tuple[tuple[float, ...], ...], ...]
    training_donors: tuple[str, ...]
    calibration_donors: tuple[str, ...]
    bootstrap_draws: tuple[tuple[str, ...], ...]
    bootstrap_valid: tuple[bool, ...]
    bootstrap_failures: tuple[str, ...]
    temperature: float = 1.0
    n_development_labels: int = 0
    n_calibration_labels: int = 0
    fit_converged: bool = True

    def __post_init__(self) -> None:
        level = str(self.level).strip().lower()
        object.__setattr__(self, "level", level)
        if level not in {"broad", "subtype"}:
            raise ValueError(f"unknown model level {self.level!r}")
        parent = None if self.parent in (None, "") else str(self.parent)
        object.__setattr__(self, "parent", parent)
        if (level == "broad") != (parent is None):
            raise ValueError("broad model requires no parent; subtype model requires one")
        classes = tuple(str(value) for value in self.classes)
        markers = tuple(str(value) for value in self.markers)
        if len(classes) < 2 or len(set(classes)) != len(classes):
            raise ValueError("probabilistic rescue models require at least two unique classes")
        if not markers or len(set(markers)) != len(markers):
            raise ValueError("model markers must be non-empty and unique")
        object.__setattr__(self, "classes", classes)
        object.__setattr__(self, "markers", markers)
        matrix = _finite_matrix(
            self.coefficients,
            rows=len(classes),
            columns=len(markers),
            name=f"{level} base coefficients",
        )
        object.__setattr__(self, "coefficients", matrix)
        if not self.bootstrap_coefficients:
            raise ValueError("at least one complete bootstrap model is required")
        bootstrap = tuple(
            _finite_matrix(
                value,
                rows=len(classes),
                columns=len(markers),
                name=f"{level} bootstrap coefficients {index}",
            )
            for index, value in enumerate(self.bootstrap_coefficients)
        )
        object.__setattr__(self, "bootstrap_coefficients", bootstrap)
        training = tuple(sorted({str(value) for value in self.training_donors}))
        calibration = tuple(sorted({str(value) for value in self.calibration_donors}))
        if len(training) < 2:
            raise ValueError("bootstrap model requires at least two training donors")
        if not calibration:
            raise ValueError("probability calibration requires held-out calibration donors")
        if set(training) & set(calibration):
            raise ValueError("training and calibration donors overlap")
        object.__setattr__(self, "training_donors", training)
        object.__setattr__(self, "calibration_donors", calibration)
        draws = tuple(tuple(str(value) for value in draw) for draw in self.bootstrap_draws)
        if len(draws) != len(bootstrap):
            raise ValueError("one donor draw is required for every bootstrap model")
        for index, draw in enumerate(draws):
            if len(draw) != len(training):
                raise ValueError(
                    f"bootstrap draw {index} has {len(draw)} donors; expected {len(training)}"
                )
            unknown = sorted(set(draw) - set(training))
            if unknown:
                raise ValueError(f"bootstrap draw {index} contains non-training donors: {unknown}")
        object.__setattr__(self, "bootstrap_draws", draws)
        if len(self.bootstrap_valid) != len(bootstrap):
            raise ValueError("one validity flag is required for every bootstrap draw")
        validity = tuple(
            _strict_bool(value, name=f"bootstrap_valid[{index}]")
            for index, value in enumerate(self.bootstrap_valid)
        )
        if len(self.bootstrap_failures) != len(bootstrap):
            raise ValueError("one failure reason is required for every bootstrap draw")
        failures: list[str] = []
        for index, value in enumerate(self.bootstrap_failures):
            if not isinstance(value, str):
                raise TypeError(f"bootstrap_failures[{index}] must be a JSON string")
            failure = value.strip()
            if validity[index] and failure:
                raise ValueError(f"valid bootstrap draw {index} cannot have a failure reason")
            if not validity[index] and not failure:
                raise ValueError(f"invalid bootstrap draw {index} requires a failure reason")
            failures.append(failure)
        object.__setattr__(self, "bootstrap_valid", validity)
        object.__setattr__(self, "bootstrap_failures", tuple(failures))
        temperature = float(self.temperature)
        if (
            not np.isfinite(temperature)
            or temperature < TEMPERATURE_MIN
            or temperature > TEMPERATURE_MAX
        ):
            raise ValueError(
                f"temperature must be finite and lie in [{TEMPERATURE_MIN:g}, {TEMPERATURE_MAX:g}]"
            )
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            scaled = np.asarray(matrix, dtype=float) / temperature
        if not np.isfinite(scaled).all():
            raise ValueError("temperature produces non-finite calibrated coefficients")
        object.__setattr__(self, "temperature", temperature)
        if int(self.n_development_labels) < 1 or int(self.n_calibration_labels) < 1:
            raise ValueError("development and calibration label counts must be positive")
        if not _strict_bool(self.fit_converged, name="fit_converged"):
            raise ValueError("a non-converged base fit cannot be promoted")

    @property
    def bootstrap_replicates(self) -> int:
        return len(self.bootstrap_coefficients)

    @property
    def valid_bootstrap_replicates(self) -> int:
        return sum(self.bootstrap_valid)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "level": self.level,
            "parent": self.parent,
            "classes": list(self.classes),
            "markers": list(self.markers),
            "coefficients": [list(row) for row in self.coefficients],
            "bootstrap_coefficients": [
                [list(row) for row in matrix] for matrix in self.bootstrap_coefficients
            ],
            "training_donors": list(self.training_donors),
            "calibration_donors": list(self.calibration_donors),
            "bootstrap_draws": [list(draw) for draw in self.bootstrap_draws],
            "bootstrap_valid": list(self.bootstrap_valid),
            "bootstrap_failures": list(self.bootstrap_failures),
            "temperature": self.temperature,
            "n_development_labels": int(self.n_development_labels),
            "n_calibration_labels": int(self.n_calibration_labels),
            "fit_converged": bool(self.fit_converged),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> FrozenModelEnsemble:
        if not isinstance(value, Mapping):
            raise TypeError("frozen model ensemble must be a mapping")
        allowed = {
            "level",
            "parent",
            "classes",
            "markers",
            "coefficients",
            "bootstrap_coefficients",
            "training_donors",
            "calibration_donors",
            "bootstrap_draws",
            "bootstrap_valid",
            "bootstrap_failures",
            "temperature",
            "n_development_labels",
            "n_calibration_labels",
            "fit_converged",
        }
        _reject_unknown(value, allowed, where="model ensemble")
        return cls(
            level=str(value.get("level", "")),
            parent=value.get("parent"),
            classes=tuple(value.get("classes", ())),
            markers=tuple(value.get("markers", ())),
            coefficients=tuple(value.get("coefficients", ())),
            bootstrap_coefficients=tuple(value.get("bootstrap_coefficients", ())),
            training_donors=tuple(value.get("training_donors", ())),
            calibration_donors=tuple(value.get("calibration_donors", ())),
            bootstrap_draws=tuple(value.get("bootstrap_draws", ())),
            bootstrap_valid=tuple(value.get("bootstrap_valid", ())),
            bootstrap_failures=tuple(value.get("bootstrap_failures", ())),
            temperature=float(value.get("temperature", np.nan)),
            n_development_labels=int(value.get("n_development_labels", 0)),
            n_calibration_labels=int(value.get("n_calibration_labels", 0)),
            fit_converged=_strict_bool(value.get("fit_converged"), name="fit_converged"),
        )

    def validate_registry(self, registry: TypingRegistry) -> None:
        expected = classifier_from_registry(registry, level=self.level, parent=self.parent)
        if self.classes != expected.classes:
            raise ValueError(
                f"{self.level} bundle classes differ from typing rules: "
                f"{self.classes} != {expected.classes}"
            )
        if self.markers != expected.markers:
            raise ValueError(
                f"{self.level} bundle markers differ from typing rules: "
                f"{self.markers} != {expected.markers}"
            )
        matrices = (self.coefficients, *self.bootstrap_coefficients)
        tolerance = 1e-12
        for model_index, value in enumerate(matrices):
            coefficients = np.asarray(value, dtype=float)
            for class_index, rule in enumerate(expected.rules):
                positive = set(rule.anchor_markers + rule.support_markers)
                negative = set(rule.exclusion_markers)
                for marker_index, marker in enumerate(expected.markers):
                    coefficient = coefficients[class_index, marker_index]
                    if marker in positive and coefficient < -tolerance:
                        raise ValueError(
                            f"model {model_index}: {rule.name}/{marker} violates positive sign"
                        )
                    if marker in negative and coefficient > tolerance:
                        raise ValueError(
                            f"model {model_index}: {rule.name}/{marker} violates exclusion sign"
                        )
                    if marker not in positive | negative and abs(coefficient) > tolerance:
                        raise ValueError(
                            f"model {model_index}: unregistered {rule.name}/{marker} coefficient is nonzero"
                        )

    def materialize_model(self, registry: TypingRegistry) -> AdditiveMultinomialModel:
        self.validate_registry(registry)
        template = classifier_from_registry(registry, level=self.level, parent=self.parent)
        return replace(
            template,
            coefficients=np.asarray(self.coefficients, dtype=float) / self.temperature,
            converged=True,
            message=(
                f"{BUNDLE_METHOD_VERSION}; temperature={self.temperature:.8g}; "
                f"bootstrap={self.bootstrap_replicates}"
            ),
            n_seeds=int(self.n_development_labels),
        )


@dataclass(frozen=True)
class TypingModelBundle:
    """A semantically content-addressed model and stability artifact."""

    bundle_version: str
    marker_registry_fingerprint: str
    typing_rules_fingerprint: str
    split_manifest: DonorSplitManifest
    label_provenance: DevelopmentLabelProvenance
    evidence_provenance: DevelopmentEvidenceProvenance
    thresholds: AssignmentThresholds
    broad: FrozenModelEnsemble
    subtypes: tuple[FrozenModelEnsemble, ...] = ()
    fit_parameters: Mapping[str, Any] = field(default_factory=dict)
    threshold_selection_provenance: ThresholdSelectionProvenance | None = None
    created_at_utc: str = ""
    schema_version: int = BUNDLE_SCHEMA_VERSION
    method_version: str = BUNDLE_METHOD_VERSION
    feature_schema_version: str = FEATURE_SCHEMA_VERSION
    score_semantics: str = SCORE_SEMANTICS
    model_fit_content_id: str = ""
    content_id: str = ""
    input_artifact: InputArtifact | None = field(default=None, repr=False, compare=False)

    def __post_init__(self) -> None:
        if int(self.schema_version) != BUNDLE_SCHEMA_VERSION:
            raise ValueError(f"unsupported typing bundle schema_version {self.schema_version!r}")
        if self.method_version != BUNDLE_METHOD_VERSION:
            raise ValueError(f"unsupported typing bundle method {self.method_version!r}")
        if self.feature_schema_version != FEATURE_SCHEMA_VERSION:
            raise ValueError(
                "typing bundle feature contract differs from production: "
                f"{self.feature_schema_version!r} != {FEATURE_SCHEMA_VERSION!r}"
            )
        if self.score_semantics != SCORE_SEMANTICS:
            raise ValueError(f"unsupported score semantics {self.score_semantics!r}")
        for name in (
            "bundle_version",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
        ):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} must be a non-empty string")
        if len(self.marker_registry_fingerprint) != 64:
            raise ValueError("marker registry fingerprint must be SHA-256")
        if len(self.typing_rules_fingerprint) != 64:
            raise ValueError("typing rules fingerprint must be SHA-256")
        if self.split_manifest.source_run_id != self.label_provenance.source_run_id:
            raise ValueError("split manifest and label bundle refer to different source runs")
        if self.label_provenance.split_manifest_content_id != self.split_manifest.content_id:
            raise ValueError("label bundle refers to a different donor split")
        if self.split_manifest.source_run_id != self.evidence_provenance.source_run_id:
            raise ValueError("split manifest and training evidence refer to different runs")
        if self.evidence_provenance.split_manifest_content_id != self.split_manifest.content_id:
            raise ValueError("training evidence refers to a different donor split")
        if self.evidence_provenance.feature_schema_version != FEATURE_SCHEMA_VERSION:
            raise ValueError("training evidence uses a different feature schema")
        if self.broad.level != "broad":
            raise ValueError("bundle broad model has the wrong level")
        parents = [model.parent for model in self.subtypes]
        if any(model.level != "subtype" for model in self.subtypes):
            raise ValueError("all subtype bundle models must have subtype level")
        if len(set(parents)) != len(parents):
            raise ValueError("bundle contains duplicate subtype parents")
        fit_parameters = _validated_fit_parameters(self.fit_parameters)
        object.__setattr__(self, "fit_parameters", fit_parameters)
        expected_bootstrap = fit_parameters["n_bootstrap"]
        declared_weight_cap = fit_parameters["max_abs_weight"]
        minimum_valid = fit_parameters["minimum_valid_bootstrap_fraction"]
        expected_train = set(self.split_manifest.development)
        expected_calibration = set(self.split_manifest.calibration)
        for model in (self.broad, *self.subtypes):
            if set(model.training_donors) - expected_train:
                raise ValueError("model training donors escape the development split")
            if set(model.calibration_donors) - expected_calibration:
                raise ValueError("model calibration donors escape the calibration split")
            if model.bootstrap_replicates != expected_bootstrap:
                raise ValueError(f"{model.level} bootstrap count differs from fit_parameters")
            valid_fraction = model.valid_bootstrap_replicates / model.bootstrap_replicates
            if valid_fraction < minimum_valid:
                raise ValueError(
                    f"{model.level} valid bootstrap fraction {valid_fraction:.3f} "
                    f"is below {minimum_valid:.3f}"
                )
            for matrix in (model.coefficients, *model.bootstrap_coefficients):
                if (
                    np.max(np.abs(np.asarray(matrix, dtype=float)), initial=0.0)
                    > declared_weight_cap
                ):
                    raise ValueError(f"{model.level} coefficients exceed declared max_abs_weight")
        if self.created_at_utc:
            try:
                datetime.fromisoformat(self.created_at_utc.replace("Z", "+00:00"))
            except ValueError as exc:
                raise ValueError("created_at_utc is not an ISO timestamp") from exc
        expected_fit_id = sha256_json(self.model_fit_identity_payload())
        if self.model_fit_content_id and self.model_fit_content_id != expected_fit_id:
            raise ValueError("typing model bundle model_fit_content_id is invalid")
        object.__setattr__(self, "model_fit_content_id", expected_fit_id)
        provenance = self.threshold_selection_provenance
        if provenance is not None:
            if not isinstance(provenance, ThresholdSelectionProvenance):
                raise TypeError(
                    "threshold_selection_provenance must be ThresholdSelectionProvenance or null"
                )
            if provenance.model_fit_content_id != expected_fit_id:
                raise ValueError("threshold selection refers to a different fitted model")
            if provenance.split_manifest.source_run_id != self.split_manifest.source_run_id:
                raise ValueError("threshold selection refers to a different source run")
            if provenance.split_manifest.content_id != self.split_manifest.content_id:
                raise ValueError("threshold selection refers to a different donor split")
            if provenance.calibration_donors != self.split_manifest.calibration:
                raise ValueError("threshold selection donors differ from the calibration split")
            expected_targets = {
                ("broad", ""),
                *{("subtype", str(model.parent)) for model in self.subtypes},
            }
            observed_targets = set(provenance.evaluated_targets)
            if observed_targets != expected_targets:
                missing = sorted(expected_targets - observed_targets)
                extra = sorted(observed_targets - expected_targets)
                raise ValueError(
                    "threshold selection target coverage is incomplete: "
                    f"missing={missing}, extra={extra}"
                )
            selected_thresholds = (
                provenance.selected_probability_threshold,
                provenance.selected_margin_threshold,
                provenance.selected_stability_threshold,
            )
            configured_thresholds = (
                float(self.thresholds.inferred_probability),
                float(self.thresholds.inferred_margin),
                float(self.thresholds.inferred_stability),
            )
            if selected_thresholds != configured_thresholds:
                raise ValueError(
                    "bundle inferred thresholds differ from the calibration "
                    "threshold-selection decision"
                )
            preliminary_thresholds = provenance.preliminary_assignment_thresholds
            for name in ("anchor_probability", "negative_probability"):
                if float(getattr(self.thresholds, name)) != float(
                    preliminary_thresholds[name]
                ):
                    raise ValueError(
                        f"bundle {name} differs from the calibration preliminary bundle"
                    )
        payload = self.identity_payload()
        expected = sha256_json(payload)
        if self.content_id and self.content_id != expected:
            raise ValueError("typing model bundle content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    @property
    def subtype_by_parent(self) -> dict[str, FrozenModelEnsemble]:
        return {str(model.parent): model for model in self.subtypes}

    def model_fit_identity_payload(self) -> dict[str, Any]:
        """Return the scientific fit payload independent of release decisions.

        The fit identity intentionally excludes the bundle schema/release version,
        creation timestamp, assignment thresholds, and threshold-selection
        provenance.  Re-freezing the same scientific fit at a new operating point
        therefore retains one stable calibration lineage identifier.
        """

        return {
            "method_version": self.method_version,
            "feature_schema_version": self.feature_schema_version,
            "score_semantics": self.score_semantics,
            "marker_registry_fingerprint": self.marker_registry_fingerprint,
            "typing_rules_fingerprint": self.typing_rules_fingerprint,
            "split_manifest": self.split_manifest.to_dict(),
            "label_provenance": self.label_provenance.to_dict(),
            "evidence_provenance": self.evidence_provenance.to_dict(),
            "broad": self.broad.identity_payload(),
            "subtypes": [
                model.identity_payload()
                for model in sorted(self.subtypes, key=lambda item: str(item.parent))
            ],
            "fit_parameters": dict(self.fit_parameters),
        }

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": int(self.schema_version),
            "method_version": self.method_version,
            "feature_schema_version": self.feature_schema_version,
            "score_semantics": self.score_semantics,
            "bundle_version": self.bundle_version,
            "model_fit_content_id": self.model_fit_content_id,
            "marker_registry_fingerprint": self.marker_registry_fingerprint,
            "typing_rules_fingerprint": self.typing_rules_fingerprint,
            "split_manifest": self.split_manifest.to_dict(),
            "label_provenance": self.label_provenance.to_dict(),
            "evidence_provenance": self.evidence_provenance.to_dict(),
            "thresholds": {
                name: float(getattr(self.thresholds, name))
                for name in (
                    "inferred_probability",
                    "inferred_margin",
                    "inferred_stability",
                    "anchor_probability",
                    "negative_probability",
                )
            },
            "broad": self.broad.identity_payload(),
            "subtypes": [
                model.identity_payload()
                for model in sorted(self.subtypes, key=lambda item: str(item.parent))
            ],
            "fit_parameters": dict(self.fit_parameters),
            "threshold_selection_provenance": (
                None
                if self.threshold_selection_provenance is None
                else self.threshold_selection_provenance.to_dict()
            ),
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.identity_payload(),
            "created_at_utc": self.created_at_utc,
            "content_id": self.content_id,
        }

    @classmethod
    def create(
        cls,
        *,
        bundle_version: str,
        marker_registry_fingerprint: str,
        typing_rules_fingerprint: str,
        split_manifest: DonorSplitManifest,
        label_provenance: DevelopmentLabelProvenance,
        evidence_provenance: DevelopmentEvidenceProvenance,
        thresholds: AssignmentThresholds,
        broad: FrozenModelEnsemble,
        subtypes: Sequence[FrozenModelEnsemble] = (),
        fit_parameters: Mapping[str, Any] | None = None,
        threshold_selection_provenance: ThresholdSelectionProvenance | None = None,
        created_at_utc: str | None = None,
    ) -> TypingModelBundle:
        return cls(
            bundle_version=str(bundle_version),
            marker_registry_fingerprint=str(marker_registry_fingerprint),
            typing_rules_fingerprint=str(typing_rules_fingerprint),
            split_manifest=split_manifest,
            label_provenance=label_provenance,
            evidence_provenance=evidence_provenance,
            thresholds=thresholds,
            broad=broad,
            subtypes=tuple(subtypes),
            fit_parameters=dict(fit_parameters or {}),
            threshold_selection_provenance=threshold_selection_provenance,
            created_at_utc=("" if created_at_utc is None else str(created_at_utc)),
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> TypingModelBundle:
        if not isinstance(value, Mapping):
            raise TypeError("typing model bundle must be a mapping")
        allowed = {
            "schema_version",
            "method_version",
            "feature_schema_version",
            "score_semantics",
            "bundle_version",
            "model_fit_content_id",
            "created_at_utc",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
            "split_manifest",
            "label_provenance",
            "evidence_provenance",
            "thresholds",
            "broad",
            "subtypes",
            "fit_parameters",
            "threshold_selection_provenance",
            "content_id",
        }
        if not str(value.get("content_id", "")).strip():
            raise ValueError("typing model bundle is missing content_id")
        _require_exact_fields(value, allowed, where="typing model bundle")
        provenance_value = value["threshold_selection_provenance"]
        if provenance_value is not None and not isinstance(provenance_value, Mapping):
            raise TypeError(
                "typing model bundle threshold_selection_provenance must be an object or null"
            )
        threshold_value = value.get("thresholds")
        if not isinstance(threshold_value, Mapping):
            raise ValueError("typing model bundle is missing thresholds")
        threshold_fields = {
            "inferred_probability",
            "inferred_margin",
            "inferred_stability",
            "anchor_probability",
            "negative_probability",
        }
        _reject_unknown(threshold_value, threshold_fields, where="thresholds")
        if set(threshold_value) != threshold_fields:
            raise ValueError("typing model bundle thresholds are incomplete")
        return cls(
            schema_version=int(value.get("schema_version", -1)),
            method_version=str(value.get("method_version", "")),
            feature_schema_version=str(value.get("feature_schema_version", "")),
            score_semantics=str(value.get("score_semantics", "")),
            bundle_version=str(value.get("bundle_version", "")),
            created_at_utc=str(value.get("created_at_utc", "")),
            marker_registry_fingerprint=str(value.get("marker_registry_fingerprint", "")),
            typing_rules_fingerprint=str(value.get("typing_rules_fingerprint", "")),
            split_manifest=DonorSplitManifest.from_dict(value.get("split_manifest", {})),
            label_provenance=DevelopmentLabelProvenance.from_dict(
                value.get("label_provenance", {})
            ),
            evidence_provenance=DevelopmentEvidenceProvenance.from_dict(
                value.get("evidence_provenance", {})
            ),
            thresholds=AssignmentThresholds(
                **{name: float(threshold_value[name]) for name in threshold_fields}
            ),
            broad=FrozenModelEnsemble.from_dict(value.get("broad", {})),
            subtypes=tuple(
                FrozenModelEnsemble.from_dict(item) for item in value.get("subtypes", ())
            ),
            fit_parameters=dict(value.get("fit_parameters", {})),
            threshold_selection_provenance=(
                None
                if provenance_value is None
                else ThresholdSelectionProvenance.from_dict(provenance_value)
            ),
            model_fit_content_id=str(value["model_fit_content_id"]),
            content_id=str(value.get("content_id", "")),
        )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def read_json(cls, path: Path | str) -> TypingModelBundle:
        path = Path(path).expanduser().resolve()
        artifact = InputArtifact.from_path("typing_model_bundle", path)
        result = cls.from_dict(_read_strict_json(path))
        artifact.validate_current(mode="fast")
        return replace(result, input_artifact=artifact)

    def validate_current(
        self,
        *,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
    ) -> None:
        if self.input_artifact is not None:
            self.input_artifact.validate_current(mode="fast")
        if self.marker_registry_fingerprint != marker_registry.fingerprint:
            raise ValueError("typing bundle marker registry fingerprint is stale")
        if self.typing_rules_fingerprint != typing_registry.rules_fingerprint:
            raise ValueError("typing bundle rules fingerprint is stale")
        self.broad.validate_registry(typing_registry)
        configured_parents = {str(rule.parent) for rule in typing_registry.subtype_rules}
        unknown = sorted(set(self.subtype_by_parent) - configured_parents)
        if unknown:
            raise ValueError(f"typing bundle contains unknown subtype parents: {unknown}")
        structural = set(self.fit_parameters["structural_rules_only_parents"])
        unknown_structural = sorted(structural - configured_parents)
        if unknown_structural:
            raise ValueError(
                f"structural rules-only parents are not configured: {unknown_structural}"
            )
        overlapping = sorted(structural & set(self.subtype_by_parent))
        if overlapping:
            raise ValueError(
                f"structural rules-only parents cannot have probabilistic models: {overlapping}"
            )
        one_class_parents = {
            parent
            for parent in configured_parents
            if len(typing_registry.rules_for(level="subtype", parent=parent)) == 1
        }
        if structural != one_class_parents:
            raise ValueError(
                "structural rules-only parents must exactly match configured "
                f"one-class subtype parents: expected={sorted(one_class_parents)}, "
                f"observed={sorted(structural)}"
            )
        expected_missing = configured_parents - structural - set(self.subtype_by_parent)
        declared_missing = set(
            self.fit_parameters["missing_probabilistic_subtype_parents"]
        )
        if declared_missing != expected_missing:
            raise ValueError(
                "declared missing probabilistic subtype parents differ from the "
                f"registry/model inventory: expected={sorted(expected_missing)}, "
                f"observed={sorted(declared_missing)}"
            )
        if self.fit_parameters["require_all_subtype_parents"]:
            if declared_missing:
                raise ValueError(
                    "required probabilistic subtype models are missing for parents: "
                    f"{sorted(declared_missing)}"
                )
        expected_specific_classes: set[str] = set()
        for broad_class in self.broad.classes:
            subtype_rules = typing_registry.rules_for(
                level="subtype", parent=broad_class
            )
            if subtype_rules:
                expected_specific_classes.update(rule.name for rule in subtype_rules)
            else:
                expected_specific_classes.add(broad_class)
        declared_specific_classes = set(
            self.fit_parameters["final_specific_classes"]
        )
        if declared_specific_classes != expected_specific_classes:
            raise ValueError(
                "declared final specific classes differ from the typing ontology: "
                f"expected={sorted(expected_specific_classes)}, "
                f"observed={sorted(declared_specific_classes)}"
            )
        for model in self.subtypes:
            model.validate_registry(typing_registry)

    def materialize_models(
        self,
        registry: TypingRegistry,
    ) -> tuple[AdditiveMultinomialModel, dict[str, AdditiveMultinomialModel]]:
        broad = self.broad.materialize_model(registry)
        subtype = {str(model.parent): model.materialize_model(registry) for model in self.subtypes}
        return broad, subtype

    def freeze_threshold_selection(
        self,
        provenance: ThresholdSelectionProvenance,
        *,
        bundle_version: str | None = None,
        created_at_utc: str | None = None,
    ) -> TypingModelBundle:
        """Attach one reviewed calibration decision without refitting models."""

        if self.threshold_selection_provenance is not None:
            raise ValueError("typing model bundle already has frozen thresholds")
        if not isinstance(provenance, ThresholdSelectionProvenance):
            raise TypeError("provenance must be ThresholdSelectionProvenance")
        if provenance.preliminary_bundle_content_id != self.content_id:
            raise ValueError(
                "threshold selection refers to a different preliminary bundle"
            )
        observed_preliminary = {
            name: float(getattr(self.thresholds, name))
            for name in (
                "inferred_probability",
                "inferred_margin",
                "inferred_stability",
                "anchor_probability",
                "negative_probability",
            )
        }
        if observed_preliminary != dict(provenance.preliminary_assignment_thresholds):
            raise ValueError(
                "threshold selection preliminary thresholds differ from the "
                "preliminary bundle"
            )
        thresholds = replace(
            self.thresholds,
            inferred_probability=provenance.selected_probability_threshold,
            inferred_margin=provenance.selected_margin_threshold,
            inferred_stability=provenance.selected_stability_threshold,
        )
        frozen = replace(
            self,
            bundle_version=(
                self.bundle_version if bundle_version is None else str(bundle_version)
            ),
            thresholds=thresholds,
            threshold_selection_provenance=provenance,
            created_at_utc=(
                self.created_at_utc if created_at_utc is None else str(created_at_utc)
            ),
            model_fit_content_id="",
            content_id="",
            input_artifact=None,
        )
        if frozen.model_fit_content_id != self.model_fit_content_id:
            raise AssertionError("freezing thresholds changed the fitted-model identity")
        return frozen


def _normalized_json_scalar(value: Any, *, name: str) -> str | int | float | bool | None:
    """Return one losslessly JSON-serializable scalar or fail closed."""

    if value is None or value is pd.NA:
        return None
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        result = float(value)
        if not np.isfinite(result):
            raise ValueError(f"{name} must be finite")
        return result
    if isinstance(value, str):
        return value
    raise TypeError(f"{name} must be a JSON scalar, not {type(value).__name__}")


@dataclass(frozen=True)
class EvaluatedDecisionRow:
    """Normalized scientific evidence for one reviewed assignment decision.

    The complete blinded sample record and adjudication fields are embedded so
    a report reader can replay both review fingerprints.  Operational flags are
    retained for transparency, but their relationships to assignment status are
    validated again before any metric is derived.
    """

    donor_id: str
    object_id: str
    sample_record: Mapping[str, Any] | tuple[tuple[str, Any], ...]
    reviewer_record: Mapping[str, Any] | tuple[tuple[str, Any], ...]
    reference_label: str
    reviewer_1: str
    reviewer_2: str
    adjudicated: bool
    root_cause: str
    evidence_sufficient: bool
    prediction: str
    assignment_status: str
    assignment_reason: str
    target_scope: bool
    called: bool
    rescue_eligible: bool
    rescue_called: bool
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(
                f"unsupported evaluated-row schema_version {self.schema_version!r}"
            )
        donor_id = str(self.donor_id).strip()
        object_id = str(self.object_id).strip()
        if not donor_id or not object_id:
            raise ValueError("evaluated rows require non-empty donor/object keys")
        object.__setattr__(self, "donor_id", donor_id)
        object.__setattr__(self, "object_id", object_id)

        raw_sample = self.sample_record
        if isinstance(raw_sample, Mapping):
            sample_items = tuple(raw_sample.items())
        elif isinstance(raw_sample, tuple):
            sample_items = raw_sample
        else:
            raise TypeError("evaluated-row sample_record must be a JSON object")
        sample: dict[str, str | int | float | bool | None] = {}
        for raw_key, raw_value in sample_items:
            if not isinstance(raw_key, str) or not raw_key:
                raise TypeError("evaluated-row sample_record keys must be non-empty strings")
            if raw_key in sample:
                raise ValueError(f"evaluated-row sample_record repeats column {raw_key!r}")
            sample[raw_key] = _normalized_json_scalar(
                raw_value,
                name=f"evaluated_row.sample_record[{raw_key!r}]",
            )
        required_sample = {
            "donor_id",
            "object_id",
            "sampling_stratum",
            "sampling_method",
            "sampling_salt",
            "stratum_population",
            "stratum_sampled",
            "sampling_weight",
            "selection_probability",
        }
        missing_sample = sorted(required_sample - set(sample))
        if missing_sample:
            raise ValueError(
                f"evaluated-row sample_record is missing required fields: {missing_sample}"
            )
        if str(sample["donor_id"]).strip() != donor_id:
            raise ValueError("evaluated-row donor_id differs from its sample record")
        if str(sample["object_id"]).strip() != object_id:
            raise ValueError("evaluated-row object_id differs from its sample record")
        weight = _strict_nonnegative_float(
            sample["sampling_weight"],
            name="evaluated_row.sample_record.sampling_weight",
        )
        if weight <= 0:
            raise ValueError("evaluated-row sampling_weight must be positive")
        sample["donor_id"] = donor_id
        sample["object_id"] = object_id
        sample["sampling_weight"] = weight
        object.__setattr__(self, "sample_record", tuple(sorted(sample.items())))

        raw_reviewer = self.reviewer_record
        if isinstance(raw_reviewer, Mapping):
            reviewer_items = tuple(raw_reviewer.items())
        elif isinstance(raw_reviewer, tuple):
            reviewer_items = raw_reviewer
        else:
            raise TypeError("evaluated-row reviewer_record must be a JSON object")
        reviewer: dict[str, str | int | float | bool | None] = {}
        for raw_key, raw_value in reviewer_items:
            if not isinstance(raw_key, str) or not raw_key:
                raise TypeError(
                    "evaluated-row reviewer_record keys must be non-empty strings"
                )
            if raw_key in reviewer:
                raise ValueError(
                    f"evaluated-row reviewer_record repeats column {raw_key!r}"
                )
            reviewer[raw_key] = _normalized_json_scalar(
                raw_value,
                name=f"evaluated_row.reviewer_record[{raw_key!r}]",
            )
        required_reviewer = {"donor_id", "object_id", "sampling_stratum"}
        missing_reviewer = sorted(required_reviewer - set(reviewer))
        if missing_reviewer:
            raise ValueError(
                "evaluated-row reviewer_record is missing required fields: "
                f"{missing_reviewer}"
            )
        if str(reviewer["donor_id"]).strip() != donor_id or str(
            reviewer["object_id"]
        ).strip() != object_id:
            raise ValueError(
                "evaluated-row reviewer keys differ from its donor/object keys"
            )
        if str(reviewer["sampling_stratum"]).strip() != str(
            sample["sampling_stratum"]
        ).strip():
            raise ValueError(
                "evaluated-row reviewer stratum differs from its sample record"
            )
        reviewer["donor_id"] = donor_id
        reviewer["object_id"] = object_id
        reviewer["sampling_stratum"] = str(
            reviewer["sampling_stratum"]
        ).strip()
        object.__setattr__(
            self, "reviewer_record", tuple(sorted(reviewer.items()))
        )

        for name in (
            "reference_label",
            "reviewer_1",
            "reviewer_2",
            "root_cause",
            "assignment_status",
            "assignment_reason",
        ):
            value = str(getattr(self, name)).strip()
            if not value:
                raise ValueError(f"evaluated-row {name} must be non-empty")
            object.__setattr__(self, name, value)
        object.__setattr__(self, "prediction", str(self.prediction).strip())
        for name in (
            "adjudicated",
            "evidence_sufficient",
            "target_scope",
            "called",
            "rescue_eligible",
            "rescue_called",
        ):
            object.__setattr__(
                self,
                name,
                _strict_bool(getattr(self, name), name=f"evaluated_row.{name}"),
            )
        if not self.adjudicated or not self.evidence_sufficient:
            raise ValueError("evaluated rows require adjudicated evidence-sufficient reviews")
        if self.reviewer_1 == self.reviewer_2:
            raise ValueError("evaluated rows require two distinct reviewers")

        known_statuses = {
            AUTHORITATIVE_ANCHOR,
            INFERRED,
            AMBIGUOUS,
            OTHER,
            UNAVAILABLE,
            "not_applicable",
        }
        if self.assignment_status not in known_statuses:
            raise ValueError(
                f"evaluated row has unknown assignment status {self.assignment_status!r}"
            )
        expected_called = self.target_scope and self.assignment_status in {
            AUTHORITATIVE_ANCHOR,
            INFERRED,
        }
        if self.called != expected_called:
            raise ValueError("evaluated-row called flag conflicts with scope/status")
        if self.target_scope and self.assignment_status == INFERRED and not self.rescue_called:
            raise ValueError("an inferred evaluated row must be recorded as a rescue call")
        if self.rescue_eligible and not self.target_scope:
            raise ValueError("evaluated-row rescue eligibility conflicts with scope/status")
        if self.rescue_called and (not self.rescue_eligible or not self.called):
            raise ValueError(
                "evaluated-row rescue call is outside the eligible accepted lane"
            )

        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("evaluated-row content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    @property
    def sampling_weight(self) -> float:
        return float(dict(self.sample_record)["sampling_weight"])

    def sample_dict(self) -> dict[str, Any]:
        return dict(self.sample_record)

    def reviewer_dict(self) -> dict[str, Any]:
        return dict(self.reviewer_record)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "donor_id": self.donor_id,
            "object_id": self.object_id,
            "sample_record": self.sample_dict(),
            "reviewer_record": self.reviewer_dict(),
            "reference_label": self.reference_label,
            "reviewer_1": self.reviewer_1,
            "reviewer_2": self.reviewer_2,
            "adjudicated": self.adjudicated,
            "root_cause": self.root_cause,
            "evidence_sufficient": self.evidence_sufficient,
            "prediction": self.prediction,
            "assignment_status": self.assignment_status,
            "assignment_reason": self.assignment_reason,
            "target_scope": self.target_scope,
            "called": self.called,
            "rescue_eligible": self.rescue_eligible,
            "rescue_called": self.rescue_called,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> EvaluatedDecisionRow:
        if not isinstance(value, Mapping):
            raise TypeError("evaluated decision row must be a mapping")
        fields = {
            "schema_version",
            "donor_id",
            "object_id",
            "sample_record",
            "reviewer_record",
            "reference_label",
            "reviewer_1",
            "reviewer_2",
            "adjudicated",
            "root_cause",
            "evidence_sufficient",
            "prediction",
            "assignment_status",
            "assignment_reason",
            "target_scope",
            "called",
            "rescue_eligible",
            "rescue_called",
            "content_id",
        }
        _require_exact_fields(value, fields, where="evaluated decision row")
        return cls(**{name: value[name] for name in fields})


@dataclass(frozen=True)
class MetricInterval:
    """One held-out metric and its donor-bootstrap confidence interval."""

    estimate: float
    lower: float
    upper: float
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported metric interval schema_version {self.schema_version!r}")
        for name in ("estimate", "lower", "upper"):
            object.__setattr__(
                self,
                name,
                _strict_unit_float(getattr(self, name), name=f"metric.{name}"),
            )
        if not self.lower <= self.estimate <= self.upper:
            raise ValueError("metric interval must contain its estimate")
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("metric interval content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "estimate": self.estimate,
            "lower": self.lower,
            "upper": self.upper,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> MetricInterval:
        if not isinstance(value, Mapping):
            raise TypeError("metric interval must be a mapping")
        fields = {"schema_version", "estimate", "lower", "upper", "content_id"}
        _require_exact_fields(value, fields, where="metric interval")
        return cls(
            schema_version=value["schema_version"],
            estimate=value["estimate"],
            lower=value["lower"],
            upper=value["upper"],
            content_id=value["content_id"],
        )


@dataclass(frozen=True)
class ClassErrorResult:
    """False-positive and false-negative evidence for one modeled class."""

    label: str
    reference_positive_rows: int
    reference_negative_rows: int
    reference_positive_donors: tuple[str, ...]
    reference_negative_donors: tuple[str, ...]
    reference_positive_weight: float
    reference_negative_weight: float
    false_positive_weight: float
    false_negative_weight: float
    false_positive_rate: MetricInterval
    false_negative_rate: MetricInterval
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported class error schema_version {self.schema_version!r}")
        label = str(self.label).strip()
        if not label:
            raise ValueError("class error label must be non-empty")
        object.__setattr__(self, "label", label)
        for name in ("reference_positive_rows", "reference_negative_rows"):
            object.__setattr__(
                self,
                name,
                _strict_positive_int(getattr(self, name), name=f"class_error.{name}"),
            )
        for name in ("reference_positive_donors", "reference_negative_donors"):
            raw = getattr(self, name)
            if isinstance(raw, (str, bytes)):
                raise TypeError(f"class_error.{name} must be a JSON array")
            donors = tuple(sorted(str(value).strip() for value in raw))
            if not donors or any(not donor for donor in donors) or len(set(donors)) != len(donors):
                raise ValueError(f"class_error.{name} must contain unique non-empty donors")
            object.__setattr__(self, name, donors)
        for name in (
            "reference_positive_weight",
            "reference_negative_weight",
            "false_positive_weight",
            "false_negative_weight",
        ):
            object.__setattr__(
                self,
                name,
                _strict_nonnegative_float(getattr(self, name), name=f"class_error.{name}"),
            )
        if self.reference_positive_weight <= 0 or self.reference_negative_weight <= 0:
            raise ValueError("each modeled class requires nonzero positive and negative support")
        if self.false_positive_weight > self.reference_negative_weight:
            raise ValueError("false-positive weight exceeds negative support")
        if self.false_negative_weight > self.reference_positive_weight:
            raise ValueError("false-negative weight exceeds positive support")
        if not isinstance(self.false_positive_rate, MetricInterval) or not isinstance(
            self.false_negative_rate, MetricInterval
        ):
            raise TypeError("class error rates must be MetricInterval objects")
        expected_rates = {
            "false_positive_rate": (self.false_positive_weight / self.reference_negative_weight),
            "false_negative_rate": (self.false_negative_weight / self.reference_positive_weight),
        }
        for name, expected_rate in expected_rates.items():
            observed = getattr(self, name).estimate
            if not np.isclose(observed, expected_rate, rtol=1e-12, atol=1e-12):
                raise ValueError(f"{name} estimate differs from its support weights")
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("class error content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "label": self.label,
            "reference_positive_rows": self.reference_positive_rows,
            "reference_negative_rows": self.reference_negative_rows,
            "reference_positive_donors": list(self.reference_positive_donors),
            "reference_negative_donors": list(self.reference_negative_donors),
            "reference_positive_weight": self.reference_positive_weight,
            "reference_negative_weight": self.reference_negative_weight,
            "false_positive_weight": self.false_positive_weight,
            "false_negative_weight": self.false_negative_weight,
            "false_positive_rate": self.false_positive_rate.to_dict(),
            "false_negative_rate": self.false_negative_rate.to_dict(),
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> ClassErrorResult:
        if not isinstance(value, Mapping):
            raise TypeError("class error result must be a mapping")
        fields = {
            "schema_version",
            "label",
            "reference_positive_rows",
            "reference_negative_rows",
            "reference_positive_donors",
            "reference_negative_donors",
            "reference_positive_weight",
            "reference_negative_weight",
            "false_positive_weight",
            "false_negative_weight",
            "false_positive_rate",
            "false_negative_rate",
            "content_id",
        }
        _require_exact_fields(value, fields, where="class error result")
        return cls(
            schema_version=value["schema_version"],
            label=value["label"],
            reference_positive_rows=value["reference_positive_rows"],
            reference_negative_rows=value["reference_negative_rows"],
            reference_positive_donors=tuple(value["reference_positive_donors"]),
            reference_negative_donors=tuple(value["reference_negative_donors"]),
            reference_positive_weight=value["reference_positive_weight"],
            reference_negative_weight=value["reference_negative_weight"],
            false_positive_weight=value["false_positive_weight"],
            false_negative_weight=value["false_negative_weight"],
            false_positive_rate=MetricInterval.from_dict(value["false_positive_rate"]),
            false_negative_rate=MetricInterval.from_dict(value["false_negative_rate"]),
            content_id=value["content_id"],
        )


def _fail_closed_ratio(
    numerator: float,
    denominator: float,
    *,
    zero_denominator: float,
) -> float:
    if denominator <= 0:
        return zero_denominator
    return float(np.clip(numerator / denominator, 0.0, 1.0))


def _target_simultaneous_intervals(
    estimates: Mapping[str, float],
    bootstrap_values: Mapping[str, np.ndarray],
    denominator_contributions: Mapping[str, np.ndarray],
    *,
    donor_ids: np.ndarray,
    confidence_level: float,
    multiplicity_target_count: int,
) -> dict[str, MetricInterval]:
    if (
        set(estimates) != set(bootstrap_values)
        or set(estimates) != set(denominator_contributions)
        or not estimates
    ):
        raise ValueError("simultaneous interval metric keys are incomplete")
    target_count = _strict_positive_int(
        multiplicity_target_count, name="multiplicity_target_count"
    )
    ordered = tuple(sorted(estimates))
    vectors = [np.asarray(bootstrap_values[name], dtype=float) for name in ordered]
    replicate_counts = {len(vector) for vector in vectors if vector.ndim == 1}
    if len(replicate_counts) != 1 or any(vector.ndim != 1 for vector in vectors):
        raise ValueError("simultaneous bootstrap vectors must have one equal length")
    n_replicates = next(iter(replicate_counts))
    if n_replicates < 1 or any(not np.isfinite(vector).all() for vector in vectors):
        raise ValueError("simultaneous bootstrap vectors must be finite and nonempty")
    alpha = 1.0 - confidence_level
    alpha_sup = alpha / (2.0 * target_count)
    order_statistic = int(np.ceil((n_replicates + 1) * (1.0 - alpha_sup)))
    if (
        n_replicates * alpha_sup < MINIMUM_SIMULTANEOUS_TAIL_DRAWS
        or order_statistic > n_replicates
    ):
        return {
            name: MetricInterval(estimate=estimates[name], lower=0.0, upper=1.0)
            for name in ordered
        }
    matrix = np.column_stack(vectors)
    point = np.asarray([estimates[name] for name in ordered], dtype=float)
    max_deviation = np.max(np.abs(matrix - point[None, :]), axis=1)
    # Finite-B conservative order statistic: k=ceil((B+1)(1-alpha_sup)),
    # indexed from one. Every metric for this target uses the same critical
    # value. Half of report alpha is reserved for this target-wise sup-norm
    # band and half for boundary guards below; Bonferroni spans targets.
    critical = float(np.sort(max_deviation)[order_statistic - 1])
    intervals = {
        name: MetricInterval(
            estimate=estimates[name],
            lower=max(0.0, estimates[name] - critical),
            upper=min(1.0, estimates[name] + critical),
        )
        for name in ordered
    }

    favorable_boundary = {
        name
        for name in ordered
        if (
            (name in {"coverage", "rescue_coverage"} and estimates[name] == 1.0)
            or (
                name not in {"coverage", "rescue_coverage"}
                and estimates[name] == 0.0
            )
        )
    }
    if not favorable_boundary:
        return intervals
    donors_array = np.asarray(donor_ids, dtype=str)
    if donors_array.ndim != 1:
        raise ValueError("boundary-guard donor IDs must be one-dimensional")
    alpha_boundary = alpha / (
        2.0 * target_count * len(favorable_boundary)
    )
    for name in favorable_boundary:
        estimate = estimates[name]
        contributions = np.asarray(denominator_contributions[name], dtype=float)
        if (
            contributions.ndim != 1
            or len(contributions) != len(donors_array)
            or not np.isfinite(contributions).all()
            or np.any(contributions < 0)
        ):
            raise ValueError(
                f"boundary denominator contributions are invalid for {name!r}"
            )
        total = float(contributions.sum())
        cell_square_sum = float(np.square(contributions).sum())
        donor_totals = np.asarray(
            [
                contributions[donors_array == donor].sum()
                for donor in sorted(set(donors_array))
            ],
            dtype=float,
        )
        donor_square_sum = float(np.square(donor_totals).sum())
        if total <= 0 or cell_square_sum <= 0 or donor_square_sum <= 0:
            intervals[name] = MetricInterval(
                estimate=estimate, lower=0.0, upper=1.0
            )
            continue
        n_cell = total * total / cell_square_sum
        n_cluster = total * total / donor_square_sum
        n_effective = int(np.floor(min(n_cell, n_cluster)))
        if n_effective < 2:
            intervals[name] = MetricInterval(
                estimate=estimate, lower=0.0, upper=1.0
            )
            continue
        zero_adverse_upper = 1.0 - alpha_boundary ** (1.0 / n_effective)
        current = intervals[name]
        if name in {"coverage", "rescue_coverage"}:
            intervals[name] = MetricInterval(
                estimate=estimate,
                lower=min(current.lower, 1.0 - zero_adverse_upper),
                upper=current.upper,
            )
        else:
            intervals[name] = MetricInterval(
                estimate=estimate,
                lower=current.lower,
                upper=max(current.upper, zero_adverse_upper),
            )
    return intervals


def _evaluation_target_seed(random_state: int, key: tuple[str, str | None]) -> int:
    digest = sha256_json(
        {"random_state": random_state, "level": key[0], "parent": key[1]}
    )
    return int(digest[:16], 16)


def _two_stage_survey_bootstrap(
    *,
    rows: tuple[EvaluatedDecisionRow, ...],
    values: np.ndarray,
    donors: tuple[str, ...],
    n_bootstrap: int,
    random_state: int,
    target_key: tuple[str, str | None],
) -> np.ndarray:
    """Rao-Wu-style cell resampling inside a donor cluster bootstrap.

    Donors are first sampled with replacement.  Within each selected donor,
    every non-census review stratum is independently resampled using nonnegative
    Rao-Wu rescaled multipliers and the SRSWOR finite-population correction.
    Census strata contribute their exact totals and therefore add no artificial
    second-stage variance.
    """

    groups: dict[tuple[str, str], list[int]] = {}
    metadata: dict[tuple[str, str], tuple[int, int]] = {}
    for index, row in enumerate(rows):
        sample = row.sample_dict()
        stratum = str(sample["sampling_stratum"]).strip()
        method = str(sample["sampling_method"]).strip()
        salt = str(sample["sampling_salt"]).strip()
        if not stratum or not salt:
            raise ValueError("evaluation sample design has a blank stratum or salt")
        if method != REVIEW_SAMPLING_METHOD:
            raise ValueError(
                "evaluation sample design uses an unsupported sampling method"
            )
        population = _strict_positive_int(
            sample["stratum_population"],
            name="evaluated_row.sample_record.stratum_population",
        )
        sampled = _strict_positive_int(
            sample["stratum_sampled"],
            name="evaluated_row.sample_record.stratum_sampled",
        )
        if sampled > population:
            raise ValueError("evaluation stratum sample exceeds its population")
        expected_weight = population / sampled
        expected_probability = sampled / population
        observed_probability = _strict_unit_float(
            sample["selection_probability"],
            name="evaluated_row.sample_record.selection_probability",
        )
        if not np.isclose(
            row.sampling_weight, expected_weight, rtol=1e-12, atol=1e-12
        ) or not np.isclose(
            observed_probability,
            expected_probability,
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError(
                "evaluation sample weight/probability differs from its stratum design"
            )
        key = row.donor_id, stratum
        previous = metadata.setdefault(key, (population, sampled))
        if previous != (population, sampled):
            raise ValueError("evaluation stratum metadata is inconsistent")
        groups.setdefault(key, []).append(index)

    for key, indices in groups.items():
        population, sampled = metadata[key]
        if len(indices) != sampled:
            raise ValueError(
                f"evaluation stratum {key} embeds {len(indices)} rows but declares "
                f"stratum_sampled={sampled}"
            )
        if population > sampled and sampled < 2:
            raise ValueError(
                f"evaluation stratum {key} is non-census but has fewer than two "
                "reviewed cells; second-stage sampling variance is not estimable"
            )

    groups_by_donor: dict[str, tuple[tuple[np.ndarray, int, int], ...]] = {}
    for donor in donors:
        donor_groups: list[tuple[np.ndarray, int, int]] = []
        for key in sorted(group for group in groups if group[0] == donor):
            population, sampled = metadata[key]
            donor_groups.append(
                (np.asarray(groups[key], dtype=int), population, sampled)
            )
        if not donor_groups:
            raise ValueError(f"evaluation donor {donor!r} has no review strata")
        groups_by_donor[donor] = tuple(donor_groups)

    rng = np.random.default_rng(_evaluation_target_seed(random_state, target_key))
    result = np.zeros((n_bootstrap, values.shape[1]), dtype=float)
    for replicate in range(n_bootstrap):
        sampled_donor_indices = rng.integers(0, len(donors), size=len(donors))
        for donor_index in sampled_donor_indices:
            donor = donors[int(donor_index)]
            for indices, population, sampled in groups_by_donor[donor]:
                stratum_values = values[indices]
                if population == sampled:
                    result[replicate] += stratum_values.sum(axis=0)
                    continue
                counts = rng.multinomial(
                    sampled - 1,
                    np.full(sampled, 1.0 / sampled, dtype=float),
                ).astype(float)
                rescaled = (sampled / (sampled - 1.0)) * counts
                fpc = np.sqrt((population - sampled) / (population - 1.0))
                multipliers = 1.0 + fpc * (rescaled - 1.0)
                result[replicate] += (
                    stratum_values * multipliers[:, None]
                ).sum(axis=0)
    return result


def _replay_target_evaluation(
    *,
    rows: tuple[EvaluatedDecisionRow, ...],
    classes: tuple[str, ...],
    target_key: tuple[str, str | None],
    confidence_level: float,
    n_bootstrap: int,
    random_state: int,
    multiplicity_target_count: int,
) -> dict[str, Any]:
    """Replay all promotion-relevant statistics from normalized cell rows."""

    weight = np.asarray([row.sampling_weight for row in rows], dtype=float)
    reference = np.asarray([row.reference_label for row in rows], dtype=str)
    prediction = np.asarray([row.prediction for row in rows], dtype=str)
    called = np.asarray([row.called for row in rows], dtype=bool)
    eligible = np.asarray([row.rescue_eligible for row in rows], dtype=bool)
    rescue_called = np.asarray([row.rescue_called for row in rows], dtype=bool)
    wrong = called & (prediction != reference)
    rescue_wrong = rescue_called & (prediction != reference)

    component_names = [
        "population",
        "called",
        "wrong",
        "rescue_eligible",
        "rescue_called",
        "rescue_wrong",
    ]
    component_values: list[np.ndarray] = [
        weight,
        weight * called,
        weight * wrong,
        weight * eligible,
        weight * rescue_called,
        weight * rescue_wrong,
    ]
    for label in classes:
        positive = reference == label
        negative = ~positive
        component_names.extend(
            [f"{label}:positive", f"{label}:negative", f"{label}:fp", f"{label}:fn"]
        )
        component_values.extend(
            [
                weight * positive,
                weight * negative,
                weight * (called & (prediction == label) & negative),
                weight * (positive & (~called | (prediction != label))),
            ]
        )
    values = np.column_stack(component_values).astype(float, copy=False)
    donors = tuple(sorted({row.donor_id for row in rows}))
    donor_index = {donor: index for index, donor in enumerate(donors)}
    donor_values = np.zeros((len(donors), values.shape[1]), dtype=float)
    for row_index, row in enumerate(rows):
        donor_values[donor_index[row.donor_id]] += values[row_index]
    point = donor_values.sum(axis=0)
    position = {name: index for index, name in enumerate(component_names)}

    if point[position["called"]] <= 0:
        raise ValueError(f"target {target_key} has no accepted calls")
    if point[position["rescue_eligible"]] <= 0:
        raise ValueError(f"target {target_key} has no rescue-eligible cells")
    if point[position["rescue_called"]] <= 0:
        raise ValueError(f"target {target_key} has no rescue calls")
    for label in classes:
        if point[position[f"{label}:positive"]] <= 0:
            raise ValueError(f"target {target_key} class {label!r} has no positive support")
        if point[position[f"{label}:negative"]] <= 0:
            raise ValueError(f"target {target_key} class {label!r} has no negative support")

    bootstrap = _two_stage_survey_bootstrap(
        rows=rows,
        values=values,
        donors=donors,
        n_bootstrap=n_bootstrap,
        random_state=random_state,
        target_key=target_key,
    )

    metric_specs = {
        "coverage": ("called", "population", 0.0),
        "selective_error": ("wrong", "called", 1.0),
        "rescue_coverage": ("rescue_called", "rescue_eligible", 0.0),
        "rescue_selective_error": ("rescue_wrong", "rescue_called", 1.0),
    }
    metric_estimates: dict[str, float] = {}
    metric_bootstrap: dict[str, np.ndarray] = {}
    metric_denominators: dict[str, np.ndarray] = {}
    for name, (numerator, denominator, fallback) in metric_specs.items():
        numerator_index = position[numerator]
        denominator_index = position[denominator]
        estimate = _fail_closed_ratio(
            point[numerator_index],
            point[denominator_index],
            zero_denominator=fallback,
        )
        bootstrap_metric = np.asarray(
            [
                _fail_closed_ratio(
                    numerator_value,
                    denominator_value,
                    zero_denominator=fallback,
                )
                for numerator_value, denominator_value in zip(
                    bootstrap[:, numerator_index],
                    bootstrap[:, denominator_index],
                    strict=True,
                )
            ],
            dtype=float,
        )
        metric_estimates[name] = estimate
        metric_bootstrap[name] = bootstrap_metric
        metric_denominators[name] = values[:, denominator_index]

    class_support: list[dict[str, Any]] = []
    for label in classes:
        positive_index = position[f"{label}:positive"]
        negative_index = position[f"{label}:negative"]
        fp_index = position[f"{label}:fp"]
        fn_index = position[f"{label}:fn"]
        fp_estimate = _fail_closed_ratio(
            point[fp_index], point[negative_index], zero_denominator=1.0
        )
        fn_estimate = _fail_closed_ratio(
            point[fn_index], point[positive_index], zero_denominator=1.0
        )
        fp_bootstrap = np.asarray(
            [
                _fail_closed_ratio(numerator, denominator, zero_denominator=1.0)
                for numerator, denominator in zip(
                    bootstrap[:, fp_index],
                    bootstrap[:, negative_index],
                    strict=True,
                )
            ],
            dtype=float,
        )
        fn_bootstrap = np.asarray(
            [
                _fail_closed_ratio(numerator, denominator, zero_denominator=1.0)
                for numerator, denominator in zip(
                    bootstrap[:, fn_index],
                    bootstrap[:, positive_index],
                    strict=True,
                )
            ],
            dtype=float,
        )
        fp_name = f"{label}:false_positive_rate"
        fn_name = f"{label}:false_negative_rate"
        metric_estimates[fp_name] = fp_estimate
        metric_estimates[fn_name] = fn_estimate
        metric_bootstrap[fp_name] = fp_bootstrap
        metric_bootstrap[fn_name] = fn_bootstrap
        metric_denominators[fp_name] = values[:, negative_index]
        metric_denominators[fn_name] = values[:, positive_index]
        class_support.append(
            {
                "label": label,
                "positive_index": positive_index,
                "negative_index": negative_index,
                "fp_index": fp_index,
                "fn_index": fn_index,
                "fp_name": fp_name,
                "fn_name": fn_name,
            }
        )

    intervals = _target_simultaneous_intervals(
        metric_estimates,
        metric_bootstrap,
        metric_denominators,
        donor_ids=np.asarray([row.donor_id for row in rows], dtype=str),
        confidence_level=confidence_level,
        multiplicity_target_count=multiplicity_target_count,
    )
    class_errors: list[ClassErrorResult] = []
    for support in class_support:
        label = support["label"]
        positive_index = support["positive_index"]
        negative_index = support["negative_index"]
        fp_index = support["fp_index"]
        fn_index = support["fn_index"]
        class_errors.append(
            ClassErrorResult(
                label=label,
                reference_positive_rows=int(np.count_nonzero(reference == label)),
                reference_negative_rows=int(np.count_nonzero(reference != label)),
                reference_positive_donors=tuple(
                    sorted(
                        {
                            row.donor_id
                            for row in rows
                            if row.reference_label == label
                        }
                    )
                ),
                reference_negative_donors=tuple(
                    sorted(
                        {
                            row.donor_id
                            for row in rows
                            if row.reference_label != label
                        }
                    )
                ),
                reference_positive_weight=float(point[positive_index]),
                reference_negative_weight=float(point[negative_index]),
                false_positive_weight=float(point[fp_index]),
                false_negative_weight=float(point[fn_index]),
                false_positive_rate=intervals[support["fp_name"]],
                false_negative_rate=intervals[support["fn_name"]],
            )
        )

    return {
        "row_count": len(rows),
        "donors": donors,
        "population_weight": float(point[position["population"]]),
        "called_weight": float(point[position["called"]]),
        "wrong_call_weight": float(point[position["wrong"]]),
        "rescue_eligible_weight": float(point[position["rescue_eligible"]]),
        "rescue_called_weight": float(point[position["rescue_called"]]),
        "rescue_wrong_call_weight": float(point[position["rescue_wrong"]]),
        "coverage": intervals["coverage"],
        "selective_error": intervals["selective_error"],
        "rescue_coverage": intervals["rescue_coverage"],
        "rescue_selective_error": intervals["rescue_selective_error"],
        "class_errors": tuple(class_errors),
    }


@dataclass(frozen=True)
class TargetEvaluationResult:
    """Replayable held-out results for one fitted model.

    Summary weights and intervals remain serialized for readability, but they
    are never trusted.  Construction and deserialization deterministically
    rederive every value from ``evaluated_rows`` and reject any difference.
    """

    level: str
    parent: str | None
    classes: tuple[str, ...]
    evaluated_rows: tuple[EvaluatedDecisionRow, ...]
    confidence_level: float
    n_bootstrap: int
    random_state: int
    multiplicity_target_count: int
    row_count: int
    donors: tuple[str, ...]
    population_weight: float
    called_weight: float
    wrong_call_weight: float
    rescue_eligible_weight: float
    rescue_called_weight: float
    rescue_wrong_call_weight: float
    coverage: MetricInterval
    selective_error: MetricInterval
    rescue_coverage: MetricInterval
    rescue_selective_error: MetricInterval
    class_errors: tuple[ClassErrorResult, ...]
    review_provenance: ReviewLabelProvenance
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported target result schema_version {self.schema_version!r}")
        level = str(self.level).strip().lower()
        parent = None if self.parent in (None, "") else str(self.parent).strip()
        if level not in {"broad", "subtype", "specific"}:
            raise ValueError("target level must be broad, subtype, or specific")
        if (level in {"broad", "specific"} and parent is not None) or (
            level == "subtype" and parent is None
        ):
            raise ValueError(
                "broad/specific targets have no parent; subtype target requires one"
            )
        object.__setattr__(self, "level", level)
        object.__setattr__(self, "parent", parent)
        classes = tuple(str(value).strip() for value in self.classes)
        if len(classes) < 2 or any(not value for value in classes):
            raise ValueError("evaluated probabilistic target requires at least two classes")
        if len(set(classes)) != len(classes):
            raise ValueError("target classes must be unique")
        object.__setattr__(self, "classes", classes)
        confidence = _strict_unit_float(
            self.confidence_level,
            name="target.confidence_level",
        )
        if confidence <= 0 or confidence >= 1:
            raise ValueError("target confidence_level must lie in (0, 1)")
        object.__setattr__(self, "confidence_level", confidence)
        object.__setattr__(
            self,
            "n_bootstrap",
            _strict_positive_int(self.n_bootstrap, name="target.n_bootstrap"),
        )
        if (
            isinstance(self.random_state, bool)
            or not isinstance(self.random_state, int)
            or self.random_state < 0
        ):
            raise TypeError("target random_state must be a nonnegative JSON integer")
        object.__setattr__(
            self,
            "multiplicity_target_count",
            _strict_positive_int(
                self.multiplicity_target_count,
                name="target.multiplicity_target_count",
            ),
        )

        rows = tuple(self.evaluated_rows)
        if not rows or any(not isinstance(value, EvaluatedDecisionRow) for value in rows):
            raise TypeError("target evaluated_rows must contain EvaluatedDecisionRow objects")
        rows = tuple(sorted(rows, key=lambda row: (row.donor_id, row.object_id)))
        keys = [(row.donor_id, row.object_id) for row in rows]
        if len(set(keys)) != len(keys):
            raise ValueError("target evaluated rows contain duplicate donor/object keys")
        sample_columns = set(dict(rows[0].sample_record))
        if any(set(dict(row.sample_record)) != sample_columns for row in rows):
            raise ValueError("target evaluated rows have inconsistent sample-record columns")
        reviewer_columns = set(dict(rows[0].reviewer_record))
        if any(set(dict(row.reviewer_record)) != reviewer_columns for row in rows):
            raise ValueError(
                "target evaluated rows have inconsistent reviewer-record columns"
            )
        supported_references = set(classes) | {"Other"}
        unsupported = sorted({row.reference_label for row in rows} - supported_references)
        if unsupported:
            raise ValueError(f"target review contains unsupported labels: {unsupported}")
        invalid_calls = sorted(
            {
                row.prediction
                for row in rows
                if row.called and row.prediction not in set(classes)
            }
        )
        if invalid_calls:
            raise ValueError(
                "target accepted predictions fall outside modeled classes: "
                f"{invalid_calls}"
            )
        object.__setattr__(self, "evaluated_rows", rows)

        row_count = _strict_positive_int(self.row_count, name="target.row_count")
        object.__setattr__(self, "row_count", row_count)
        if isinstance(self.donors, (str, bytes)):
            raise TypeError("target donors must be a JSON array")
        donors = tuple(sorted(str(value).strip() for value in self.donors))
        if not donors or any(not donor for donor in donors):
            raise ValueError("target donors must contain non-empty donor IDs")
        if len(set(donors)) != len(donors):
            raise ValueError("target donors contain duplicates")
        object.__setattr__(self, "donors", donors)
        for name in (
            "population_weight",
            "called_weight",
            "wrong_call_weight",
            "rescue_eligible_weight",
            "rescue_called_weight",
            "rescue_wrong_call_weight",
        ):
            object.__setattr__(
                self,
                name,
                _strict_nonnegative_float(getattr(self, name), name=f"target.{name}"),
            )
        if self.population_weight <= 0:
            raise ValueError("target population weight must be positive")
        if self.called_weight <= 0 or self.called_weight > self.population_weight:
            raise ValueError("target called weight must lie in (0, population]")
        if self.wrong_call_weight > self.called_weight:
            raise ValueError("target wrong-call weight exceeds called weight")
        if self.rescue_eligible_weight <= 0 or self.rescue_eligible_weight > self.population_weight:
            raise ValueError("target rescue-eligible weight must lie in (0, population]")
        if (
            self.rescue_called_weight <= 0
            or self.rescue_called_weight > self.rescue_eligible_weight
            or self.rescue_called_weight > self.called_weight
        ):
            raise ValueError(
                "target rescue-called weight must be nonzero and no larger than "
                "eligible/total called weight"
            )
        if self.rescue_wrong_call_weight > self.rescue_called_weight:
            raise ValueError("target rescue wrong-call weight exceeds rescue calls")
        intervals = (
            self.coverage,
            self.selective_error,
            self.rescue_coverage,
            self.rescue_selective_error,
        )
        if any(not isinstance(value, MetricInterval) for value in intervals):
            raise TypeError("target performance metrics must be MetricInterval objects")
        expected_metrics = {
            "coverage": self.called_weight / self.population_weight,
            "selective_error": self.wrong_call_weight / self.called_weight,
            "rescue_coverage": (self.rescue_called_weight / self.rescue_eligible_weight),
            "rescue_selective_error": (self.rescue_wrong_call_weight / self.rescue_called_weight),
        }
        for name, expected_value in expected_metrics.items():
            if not np.isclose(
                getattr(self, name).estimate,
                expected_value,
                rtol=1e-12,
                atol=1e-12,
            ):
                raise ValueError(f"target {name} estimate differs from its weights")
        errors = tuple(self.class_errors)
        if any(not isinstance(value, ClassErrorResult) for value in errors):
            raise TypeError("target class_errors must contain ClassErrorResult objects")
        if tuple(value.label for value in errors) != classes:
            raise ValueError("target class-error labels must exactly match modeled class order")
        if not np.isclose(
            sum(value.false_positive_weight for value in errors),
            self.wrong_call_weight,
            rtol=1e-12,
            atol=1e-12,
        ):
            raise ValueError("class false-positive weights must sum to target wrong-call weight")
        object.__setattr__(self, "class_errors", errors)
        if not isinstance(self.review_provenance, ReviewLabelProvenance):
            raise TypeError("target review_provenance must be ReviewLabelProvenance")
        expected_parent = "" if parent is None else parent
        provenance = self.review_provenance
        if (provenance.level, provenance.parent) != (level, expected_parent):
            raise ValueError("target review provenance level/parent differs from the model target")
        if not provenance.complete_evidence_sufficient:
            raise ValueError(
                "target review provenance must certify complete evidence-sufficient review"
            )
        if provenance.sample_row_count != row_count:
            raise ValueError("target evaluation row_count differs from its complete review sample")

        sample_frame = pd.DataFrame([row.sample_dict() for row in rows])
        if review_sample_frame_fingerprint(sample_frame) != provenance.sample_frame_fingerprint:
            raise ValueError(
                "target embedded sample rows differ from review sample provenance"
            )
        reviewer_frame = pd.DataFrame(
            [row.reviewer_dict() for row in rows]
        )
        if (
            reviewer_frame_fingerprint(reviewer_frame)
            != provenance.reviewer_frame_fingerprint
        ):
            raise ValueError(
                "target embedded reviewer rows differ from reviewer-frame provenance"
            )
        review_frame = pd.DataFrame(
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
                for row in rows
            ],
            columns=list(REVIEW_LABEL_COLUMNS),
        )
        if review_label_frame_fingerprint(review_frame) != provenance.fingerprint:
            raise ValueError(
                "target embedded review rows differ from review-label provenance"
            )

        replayed = _replay_target_evaluation(
            rows=rows,
            classes=classes,
            target_key=(level, parent),
            confidence_level=confidence,
            n_bootstrap=self.n_bootstrap,
            random_state=self.random_state,
            multiplicity_target_count=self.multiplicity_target_count,
        )
        scalar_fields = (
            "row_count",
            "donors",
            "population_weight",
            "called_weight",
            "wrong_call_weight",
            "rescue_eligible_weight",
            "rescue_called_weight",
            "rescue_wrong_call_weight",
        )
        mismatched_scalars = [
            name for name in scalar_fields if getattr(self, name) != replayed[name]
        ]
        if mismatched_scalars:
            raise ValueError(
                "target summary differs from replayed evaluated rows: "
                f"{sorted(mismatched_scalars)}"
            )
        metric_fields = (
            "coverage",
            "selective_error",
            "rescue_coverage",
            "rescue_selective_error",
        )
        mismatched_metrics = [
            name
            for name in metric_fields
            if getattr(self, name).to_dict() != replayed[name].to_dict()
        ]
        if mismatched_metrics:
            raise ValueError(
                "target metric intervals differ from deterministic donor-bootstrap replay: "
                f"{sorted(mismatched_metrics)}"
            )
        if tuple(value.to_dict() for value in self.class_errors) != tuple(
            value.to_dict() for value in replayed["class_errors"]
        ):
            raise ValueError(
                "target class errors differ from deterministic donor-bootstrap replay"
            )
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("target evaluation content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    @property
    def target_key(self) -> tuple[str, str | None]:
        return self.level, self.parent

    @property
    def evaluated_rows_fingerprint(self) -> str:
        return sha256_json(
            {
                "schema_version": PROMOTION_COMPONENT_SCHEMA_VERSION,
                "target": [self.level, self.parent],
                "rows": [row.to_dict() for row in self.evaluated_rows],
            }
        )

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "level": self.level,
            "parent": self.parent,
            "classes": list(self.classes),
            "evaluated_rows": [row.to_dict() for row in self.evaluated_rows],
            "confidence_level": self.confidence_level,
            "n_bootstrap": self.n_bootstrap,
            "random_state": self.random_state,
            "multiplicity_target_count": self.multiplicity_target_count,
            "row_count": self.row_count,
            "donors": list(self.donors),
            "population_weight": self.population_weight,
            "called_weight": self.called_weight,
            "wrong_call_weight": self.wrong_call_weight,
            "rescue_eligible_weight": self.rescue_eligible_weight,
            "rescue_called_weight": self.rescue_called_weight,
            "rescue_wrong_call_weight": self.rescue_wrong_call_weight,
            "coverage": self.coverage.to_dict(),
            "selective_error": self.selective_error.to_dict(),
            "rescue_coverage": self.rescue_coverage.to_dict(),
            "rescue_selective_error": self.rescue_selective_error.to_dict(),
            "class_errors": [value.to_dict() for value in self.class_errors],
            "review_provenance": self.review_provenance.to_dict(),
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    @classmethod
    def from_evaluated_rows(
        cls,
        *,
        level: str,
        parent: str | None,
        classes: Sequence[str],
        evaluated_rows: Sequence[EvaluatedDecisionRow],
        confidence_level: float,
        n_bootstrap: int,
        random_state: int,
        multiplicity_target_count: int,
        review_provenance: ReviewLabelProvenance,
    ) -> TargetEvaluationResult:
        """Construct a target report without accepting caller-authored metrics."""

        normalized_level = str(level).strip().lower()
        normalized_parent = None if parent in (None, "") else str(parent).strip()
        normalized_classes = tuple(str(value).strip() for value in classes)
        rows = tuple(sorted(evaluated_rows, key=lambda row: (row.donor_id, row.object_id)))
        replayed = _replay_target_evaluation(
            rows=rows,
            classes=normalized_classes,
            target_key=(normalized_level, normalized_parent),
            confidence_level=float(confidence_level),
            n_bootstrap=int(n_bootstrap),
            random_state=int(random_state),
            multiplicity_target_count=int(multiplicity_target_count),
        )
        return cls(
            level=normalized_level,
            parent=normalized_parent,
            classes=normalized_classes,
            evaluated_rows=rows,
            confidence_level=confidence_level,
            n_bootstrap=n_bootstrap,
            random_state=random_state,
            multiplicity_target_count=multiplicity_target_count,
            review_provenance=review_provenance,
            **replayed,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> TargetEvaluationResult:
        if not isinstance(value, Mapping):
            raise TypeError("target evaluation result must be a mapping")
        fields = {
            "schema_version",
            "level",
            "parent",
            "classes",
            "evaluated_rows",
            "confidence_level",
            "n_bootstrap",
            "random_state",
            "multiplicity_target_count",
            "row_count",
            "donors",
            "population_weight",
            "called_weight",
            "wrong_call_weight",
            "rescue_eligible_weight",
            "rescue_called_weight",
            "rescue_wrong_call_weight",
            "coverage",
            "selective_error",
            "rescue_coverage",
            "rescue_selective_error",
            "class_errors",
            "review_provenance",
            "content_id",
        }
        _require_exact_fields(value, fields, where="target evaluation result")
        return cls(
            schema_version=value["schema_version"],
            level=value["level"],
            parent=value["parent"],
            classes=tuple(value["classes"]),
            evaluated_rows=tuple(
                EvaluatedDecisionRow.from_dict(item) for item in value["evaluated_rows"]
            ),
            confidence_level=value["confidence_level"],
            n_bootstrap=value["n_bootstrap"],
            random_state=value["random_state"],
            multiplicity_target_count=value["multiplicity_target_count"],
            row_count=value["row_count"],
            donors=tuple(value["donors"]),
            population_weight=value["population_weight"],
            called_weight=value["called_weight"],
            wrong_call_weight=value["wrong_call_weight"],
            rescue_eligible_weight=value["rescue_eligible_weight"],
            rescue_called_weight=value["rescue_called_weight"],
            rescue_wrong_call_weight=value["rescue_wrong_call_weight"],
            coverage=MetricInterval.from_dict(value["coverage"]),
            selective_error=MetricInterval.from_dict(value["selective_error"]),
            rescue_coverage=MetricInterval.from_dict(value["rescue_coverage"]),
            rescue_selective_error=MetricInterval.from_dict(value["rescue_selective_error"]),
            class_errors=tuple(ClassErrorResult.from_dict(item) for item in value["class_errors"]),
            review_provenance=ReviewLabelProvenance.from_dict(value["review_provenance"]),
            content_id=value["content_id"],
        )


@dataclass(frozen=True)
class PromotionGatePolicy:
    """Prespecified confidence-bound gates applied to every fitted model."""

    minimum_donors: int
    minimum_rows: int
    minimum_coverage_lower: float
    maximum_selective_error_upper: float
    minimum_rescue_coverage_lower: float
    maximum_rescue_selective_error_upper: float
    maximum_class_false_positive_rate_upper: float
    maximum_class_false_negative_rate_upper: float
    minimum_class_positive_rows: int
    minimum_class_negative_rows: int
    minimum_class_positive_donors: int
    minimum_class_negative_donors: int
    confidence_level: float
    bootstrap_replicates: int
    bootstrap_unit: str = EVALUATION_BOOTSTRAP_UNIT
    multiplicity_method: str = EVALUATION_MULTIPLICITY_METHOD
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(f"unsupported promotion policy schema_version {self.schema_version!r}")
        for name in (
            "minimum_donors",
            "minimum_rows",
            "minimum_class_positive_rows",
            "minimum_class_negative_rows",
            "minimum_class_positive_donors",
            "minimum_class_negative_donors",
            "bootstrap_replicates",
        ):
            object.__setattr__(
                self,
                name,
                _strict_positive_int(getattr(self, name), name=f"policy.{name}"),
            )
        for name in (
            "minimum_coverage_lower",
            "maximum_selective_error_upper",
            "minimum_rescue_coverage_lower",
            "maximum_rescue_selective_error_upper",
            "maximum_class_false_positive_rate_upper",
            "maximum_class_false_negative_rate_upper",
        ):
            object.__setattr__(
                self,
                name,
                _strict_unit_float(getattr(self, name), name=f"policy.{name}"),
            )
        confidence = _strict_unit_float(self.confidence_level, name="policy.confidence_level")
        if confidence <= 0 or confidence >= 1:
            raise ValueError("policy confidence_level must lie in (0, 1)")
        object.__setattr__(self, "confidence_level", confidence)
        absolute_minima = {
            "minimum_donors": 5,
            "minimum_rows": 20,
            "minimum_class_positive_rows": 3,
            "minimum_class_negative_rows": 3,
            "minimum_class_positive_donors": 3,
            "minimum_class_negative_donors": 3,
            "bootstrap_replicates": 1000,
        }
        too_weak = {
            name: (getattr(self, name), floor)
            for name, floor in absolute_minima.items()
            if getattr(self, name) < floor
        }
        if too_weak:
            raise ValueError(
                "promotion policy is weaker than absolute production floors: "
                f"{too_weak}"
            )
        if confidence < 0.90:
            raise ValueError("promotion policy confidence_level must be >= 0.90")
        if self.minimum_coverage_lower <= 0 or self.minimum_rescue_coverage_lower <= 0:
            raise ValueError("promotion coverage lower-bound gates must be > 0")
        error_gates = (
            self.maximum_selective_error_upper,
            self.maximum_rescue_selective_error_upper,
            self.maximum_class_false_positive_rate_upper,
            self.maximum_class_false_negative_rate_upper,
        )
        if any(value >= 1 for value in error_gates):
            raise ValueError("promotion error upper-bound gates must be < 1")
        if self.bootstrap_unit != EVALUATION_BOOTSTRAP_UNIT:
            raise ValueError(
                "promotion policy bootstrap_unit must equal "
                f"{EVALUATION_BOOTSTRAP_UNIT!r}"
            )
        if self.multiplicity_method != EVALUATION_MULTIPLICITY_METHOD:
            raise ValueError(
                "promotion multiplicity_method must equal "
                f"{EVALUATION_MULTIPLICITY_METHOD!r}"
            )
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("promotion policy content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "minimum_donors": self.minimum_donors,
            "minimum_rows": self.minimum_rows,
            "minimum_coverage_lower": self.minimum_coverage_lower,
            "maximum_selective_error_upper": self.maximum_selective_error_upper,
            "minimum_rescue_coverage_lower": self.minimum_rescue_coverage_lower,
            "maximum_rescue_selective_error_upper": (self.maximum_rescue_selective_error_upper),
            "maximum_class_false_positive_rate_upper": (
                self.maximum_class_false_positive_rate_upper
            ),
            "maximum_class_false_negative_rate_upper": (
                self.maximum_class_false_negative_rate_upper
            ),
            "minimum_class_positive_rows": self.minimum_class_positive_rows,
            "minimum_class_negative_rows": self.minimum_class_negative_rows,
            "minimum_class_positive_donors": self.minimum_class_positive_donors,
            "minimum_class_negative_donors": self.minimum_class_negative_donors,
            "confidence_level": self.confidence_level,
            "bootstrap_replicates": self.bootstrap_replicates,
            "bootstrap_unit": self.bootstrap_unit,
            "multiplicity_method": self.multiplicity_method,
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        """Persist the exact preregistered policy as strict JSON."""

        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> PromotionGatePolicy:
        if not isinstance(value, Mapping):
            raise TypeError("promotion policy must be a mapping")
        fields = {
            "schema_version",
            "minimum_donors",
            "minimum_rows",
            "minimum_coverage_lower",
            "maximum_selective_error_upper",
            "minimum_rescue_coverage_lower",
            "maximum_rescue_selective_error_upper",
            "maximum_class_false_positive_rate_upper",
            "maximum_class_false_negative_rate_upper",
            "minimum_class_positive_rows",
            "minimum_class_negative_rows",
            "minimum_class_positive_donors",
            "minimum_class_negative_donors",
            "confidence_level",
            "bootstrap_replicates",
            "bootstrap_unit",
            "multiplicity_method",
            "content_id",
        }
        _require_exact_fields(value, fields, where="promotion policy")
        return cls(**{name: value[name] for name in fields})

    @classmethod
    def read_json(cls, path: Path | str) -> PromotionGatePolicy:
        """Load and content-validate an exact preregistered policy."""

        return cls.from_dict(_read_strict_json(path))


def evaluation_bootstrap_random_state(
    *,
    split_name: str,
    model_bundle_content_id: str,
    assignment_run_id: str,
) -> int:
    """Derive the only valid report bootstrap seed from immutable lineage."""

    digest = sha256_json(
        {
            "method_version": EVALUATION_REPORT_METHOD_VERSION,
            "split_name": str(split_name),
            "model_bundle_content_id": str(model_bundle_content_id),
            "assignment_run_id": str(assignment_run_id),
        }
    )
    return int(digest[:8], 16)


def _evaluation_frame_fingerprint(
    *,
    method_version: str,
    split_name: str,
    model_bundle_content_id: str,
    split_manifest_content_id: str,
    source_run_id: str,
    candidate_evaluation_manifest_content_id: str,
    assignment_run_id: str,
    assignment_artifact_content_id: str,
    evaluation_donors: tuple[str, ...],
    targets: tuple[TargetEvaluationResult, ...],
) -> str:
    return sha256_json(
        {
            "schema_version": PROMOTION_COMPONENT_SCHEMA_VERSION,
            "method_version": method_version,
            "split_name": split_name,
            "model_bundle_content_id": model_bundle_content_id,
            "split_manifest_content_id": split_manifest_content_id,
            "source_run_id": source_run_id,
            "candidate_evaluation_manifest_content_id": (
                candidate_evaluation_manifest_content_id
            ),
            "assignment_run_id": assignment_run_id,
            "assignment_artifact_content_id": assignment_artifact_content_id,
            "evaluation_donors": list(evaluation_donors),
            "targets": [
                {
                    "level": target.level,
                    "parent": target.parent,
                    "review_provenance_content_id": target.review_provenance.content_id,
                    "sample_frame_fingerprint": (
                        target.review_provenance.sample_frame_fingerprint
                    ),
                    "review_label_fingerprint": target.review_provenance.fingerprint,
                    "evaluated_rows_fingerprint": target.evaluated_rows_fingerprint,
                }
                for target in targets
            ],
        }
    )


def _candidate_manifest_for_report(candidate_value: Mapping[str, Any]) -> Any:
    """Parse embedded lineage without touching large candidate partitions."""

    # Local import avoids the candidate module's intentional dependency on the
    # model-bundle contract during module initialization.
    from .candidate_evaluation import CandidateEvaluationManifest

    return CandidateEvaluationManifest.from_dict(candidate_value)


def _candidate_assignment_partition_for_report(
    donor_snapshot: Any,
    *,
    columns: Sequence[str],
) -> pd.DataFrame:
    """Read one projected donor partition and recheck its exact cell universe."""

    projected = tuple(dict.fromkeys(("donor_id", "object_id", *columns)))
    donor_frames = [
        pd.read_parquet(Path(fingerprint.path), columns=list(projected))
        for fingerprint in donor_snapshot.files
    ]
    frame = (
        pd.concat(donor_frames, ignore_index=True)
        if len(donor_frames) > 1
        else donor_frames[0]
    )
    for column in ("donor_id", "object_id"):
        frame[column] = frame[column].astype("string").fillna("").str.strip()
    if (frame[["donor_id", "object_id"]] == "").any().any():
        raise ValueError("embedded candidate output contains blank donor/object keys")
    if set(frame["donor_id"]) != {donor_snapshot.donor_id}:
        raise ValueError("embedded candidate partition contains a different donor")
    if len(frame) != donor_snapshot.row_count:
        raise ValueError("embedded candidate partition row count differs from its snapshot")
    if frame.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("embedded candidate partition contains duplicate donor/object keys")
    if hash_identifiers(frame["object_id"]) != donor_snapshot.object_id_sha256:
        raise ValueError("embedded candidate object universe differs from its snapshot")
    return frame


def _verify_target_evidence_against_candidate(
    target: TargetEvaluationResult,
    assignments: pd.DataFrame,
) -> bool:
    """Bind embedded reviewed decisions to exact candidate assignment bytes."""

    assignment_donors = set(assignments["donor_id"].astype(str))
    rows = tuple(row for row in target.evaluated_rows if row.donor_id in assignment_donors)
    if target.level in {"broad", "specific"}:
        population_scope = pd.Series(True, index=assignments.index, dtype=bool)
    else:
        population_scope = (
            assignments["broad_assignment_status"]
            .astype("string")
            .fillna("")
            .str.strip()
            .isin({AUTHORITATIVE_ANCHOR, INFERRED})
            & assignments["broad_label"]
            .astype("string")
            .fillna("")
            .str.strip()
            .eq(str(target.parent))
            & assignments["typing_subtype_model_source"]
            .astype("string")
            .fillna("")
            .str.strip()
            .eq("probabilistic_bundle")
        )
    population_nonempty = bool(population_scope.any())
    if not population_nonempty:
        if rows:
            raise ValueError(
                f"target {target.target_key} has reviewed rows for an empty donor scope"
            )
        return False
    if not rows:
        raise ValueError(
            f"target {target.target_key} omits a donor with a nonempty candidate population"
        )

    prefix = "broad" if target.level == "broad" else "subtype"
    required = {"donor_id", "object_id"}
    if target.level == "specific":
        required.update(
            {
                "specific_type",
                "specific_assignment_status",
                "broad_assignment_status",
                "broad_reason",
                "broad_threshold_eligible",
                "subtype_assignment_status",
                "subtype_reason",
                "subtype_threshold_eligible",
            }
        )
    else:
        required.update(
            {
                f"{prefix}_label",
                f"{prefix}_assignment_status",
                f"{prefix}_reason",
                f"{prefix}_threshold_eligible",
            }
        )
    if target.level == "subtype":
        required.update(
            {"broad_label", "broad_assignment_status", "typing_subtype_model_source"}
        )
    missing = sorted(required - set(assignments.columns))
    if missing:
        raise KeyError(
            f"candidate assignments for target {target.target_key} are missing columns: {missing}"
        )
    indexed = assignments.loc[:, sorted(required)].copy()
    reviewed = pd.DataFrame(
        [(row.donor_id, row.object_id) for row in rows],
        columns=["donor_id", "object_id"],
    )
    selected = reviewed.merge(
        indexed,
        on=["donor_id", "object_id"],
        how="left",
        validate="1:1",
        indicator=True,
    )
    if selected["_merge"].ne("both").any():
        raise ValueError("embedded evaluated rows are missing from candidate assignments")
    selected = selected.drop(columns="_merge")
    if target.level == "specific":
        prediction = selected["specific_type"].astype("string").fillna("").str.strip()
        status = (
            selected["specific_assignment_status"]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        broad_status = (
            selected["broad_assignment_status"]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        subtype_status = (
            selected["subtype_assignment_status"]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        use_subtype = broad_status.isin({AUTHORITATIVE_ANCHOR, INFERRED}) & (
            subtype_status != "not_applicable"
        )
        broad_reason = selected["broad_reason"].astype("string").fillna("").str.strip()
        subtype_reason = (
            selected["subtype_reason"].astype("string").fillna("").str.strip()
        )
        reason = broad_reason.where(~use_subtype, subtype_reason)
        broad_eligible = selected["broad_threshold_eligible"]
        subtype_eligible = selected["subtype_threshold_eligible"]
        for name, raw in (
            ("broad_threshold_eligible", broad_eligible),
            ("subtype_threshold_eligible", subtype_eligible),
        ):
            if raw.isna().any() or not raw.isin([True, False, 0, 1]).all():
                raise ValueError(f"candidate {name} must contain explicit booleans")
        eligible = broad_eligible.astype(bool) | (
            use_subtype & subtype_eligible.astype(bool)
        )
    else:
        prediction = selected[f"{prefix}_label"].astype("string").fillna("").str.strip()
        status = (
            selected[f"{prefix}_assignment_status"]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        reason = selected[f"{prefix}_reason"].astype("string").fillna("").str.strip()
        raw_eligible = selected[f"{prefix}_threshold_eligible"]
        if raw_eligible.isna().any() or not raw_eligible.isin([True, False, 0, 1]).all():
            raise ValueError("candidate threshold eligibility must contain explicit booleans")
        eligible = raw_eligible.astype(bool)
    if target.level in {"broad", "specific"}:
        scope = pd.Series(True, index=selected.index, dtype=bool)
    else:
        broad_status = (
            selected["broad_assignment_status"].astype("string").fillna("").str.strip()
        )
        broad_label = selected["broad_label"].astype("string").fillna("").str.strip()
        subtype_source = (
            selected["typing_subtype_model_source"]
            .astype("string")
            .fillna("")
            .str.strip()
        )
        scope = (
            broad_status.isin({AUTHORITATIVE_ANCHOR, INFERRED})
            & broad_label.eq(str(target.parent))
            & subtype_source.eq("probabilistic_bundle")
        )
    called = scope & status.isin({AUTHORITATIVE_ANCHOR, INFERRED})
    rescue_eligible = scope & eligible
    if target.level == "specific":
        rescue_called = called & (
            broad_status.eq(INFERRED) | (use_subtype & subtype_status.eq(INFERRED))
        )
    else:
        rescue_called = scope & status.eq(INFERRED)
    for index, row in enumerate(rows):
        observed = {
            "prediction": str(prediction.iloc[index]),
            "assignment_status": str(status.iloc[index]),
            "assignment_reason": str(reason.iloc[index]),
            "target_scope": bool(scope.iloc[index]),
            "called": bool(called.iloc[index]),
            "rescue_eligible": bool(rescue_eligible.iloc[index]),
            "rescue_called": bool(rescue_called.iloc[index]),
        }
        expected = {name: getattr(row, name) for name in observed}
        mismatched = sorted(name for name in observed if observed[name] != expected[name])
        if mismatched:
            raise ValueError(
                f"target {target.target_key} evaluated row "
                f"{(row.donor_id, row.object_id)} differs from exact candidate output: "
                f"{mismatched}"
            )
    return True


@dataclass(frozen=True)
class SplitEvaluationReport:
    """Content-addressed metrics for one exact held-out assignment artifact."""

    split_name: str
    model_bundle_content_id: str
    split_manifest_content_id: str
    source_run_id: str
    source_run_manifest: Mapping[str, Any]
    candidate_evaluation_manifest_content_id: str
    candidate_evaluation_manifest: Mapping[str, Any]
    assignment_run_id: str
    assignment_artifact_content_id: str
    evaluation_frame_fingerprint: str
    evaluation_donors: tuple[str, ...]
    promotion_policy: PromotionGatePolicy
    confidence_level: float
    n_bootstrap: int
    random_state: int
    targets: tuple[TargetEvaluationResult, ...]
    bootstrap_unit: str = EVALUATION_BOOTSTRAP_UNIT
    method_version: str = EVALUATION_REPORT_METHOD_VERSION
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(
                f"unsupported evaluation report schema_version {self.schema_version!r}"
            )
        if self.method_version != EVALUATION_REPORT_METHOD_VERSION:
            raise ValueError(f"unsupported evaluation report method {self.method_version!r}")
        split_name = str(self.split_name).strip()
        if split_name not in {"validation", "locked_test"}:
            raise ValueError("evaluation report split must be validation or locked_test")
        object.__setattr__(self, "split_name", split_name)
        source_run_id = str(self.source_run_id).strip()
        assignment_run_id = str(self.assignment_run_id).strip()
        if not source_run_id or not assignment_run_id:
            raise ValueError("evaluation source and assignment run IDs must be non-empty")
        object.__setattr__(self, "source_run_id", source_run_id)
        object.__setattr__(self, "assignment_run_id", assignment_run_id)
        for name in (
            "model_bundle_content_id",
            "split_manifest_content_id",
            "candidate_evaluation_manifest_content_id",
            "assignment_artifact_content_id",
            "evaluation_frame_fingerprint",
        ):
            object.__setattr__(self, name, _require_sha256(getattr(self, name), name=name))
        if isinstance(self.evaluation_donors, (str, bytes)):
            raise TypeError("evaluation_donors must be a JSON array")
        donors = tuple(sorted(str(value).strip() for value in self.evaluation_donors))
        if not donors or any(not donor for donor in donors):
            raise ValueError("evaluation_donors must contain non-empty donor IDs")
        if len(set(donors)) != len(donors):
            raise ValueError("evaluation_donors contains duplicates")
        object.__setattr__(self, "evaluation_donors", donors)
        if not isinstance(self.promotion_policy, PromotionGatePolicy):
            raise TypeError("evaluation promotion_policy must be PromotionGatePolicy")
        confidence = _strict_unit_float(self.confidence_level, name="evaluation.confidence_level")
        if confidence <= 0 or confidence >= 1:
            raise ValueError("evaluation confidence_level must lie in (0, 1)")
        object.__setattr__(self, "confidence_level", confidence)
        object.__setattr__(
            self,
            "n_bootstrap",
            _strict_positive_int(self.n_bootstrap, name="evaluation.n_bootstrap"),
        )
        expected_confidence = _promotion_report_confidence_level(
            self.promotion_policy.confidence_level
        )
        if confidence != expected_confidence:
            raise ValueError("evaluation confidence differs from its precommitted policy")
        if self.n_bootstrap != self.promotion_policy.bootstrap_replicates:
            raise ValueError("evaluation bootstrap count differs from its precommitted policy")
        if (
            isinstance(self.random_state, bool)
            or not isinstance(self.random_state, int)
            or self.random_state < 0
        ):
            raise TypeError("evaluation random_state must be a nonnegative JSON integer")
        expected_random_state = evaluation_bootstrap_random_state(
            split_name=split_name,
            model_bundle_content_id=self.model_bundle_content_id,
            assignment_run_id=self.assignment_run_id,
        )
        if self.random_state != expected_random_state:
            raise ValueError(
                "evaluation random_state differs from its immutable-lineage-derived seed"
            )
        if self.bootstrap_unit != EVALUATION_BOOTSTRAP_UNIT:
            raise ValueError(
                "evaluation bootstrap_unit must equal "
                f"{EVALUATION_BOOTSTRAP_UNIT!r}"
            )
        targets = tuple(self.targets)
        if not targets or any(not isinstance(value, TargetEvaluationResult) for value in targets):
            raise TypeError("evaluation targets must contain TargetEvaluationResult objects")
        keys = [value.target_key for value in targets]
        if len(set(keys)) != len(keys):
            raise ValueError("evaluation report contains duplicate model targets")
        expected_target_count = len(targets)
        for target in targets:
            unknown_target_donors = sorted(set(target.donors) - set(donors))
            if unknown_target_donors:
                raise ValueError(
                    f"target {target.target_key} contains donors outside the "
                    f"evaluation split: {unknown_target_donors}"
                )
            if (
                target.confidence_level != confidence
                or target.n_bootstrap != self.n_bootstrap
                or target.random_state != self.random_state
                or target.multiplicity_target_count != expected_target_count
            ):
                raise ValueError(
                    f"target {target.target_key} bootstrap configuration differs "
                    "from the evaluation report"
                )
            provenance = target.review_provenance
            lineage = {
                "split_name": (provenance.split_name, self.split_name),
                "source_run_id": (provenance.source_run_id, self.source_run_id),
                "split_manifest_content_id": (
                    provenance.split_manifest_content_id,
                    self.split_manifest_content_id,
                ),
                "model_bundle_content_id": (
                    provenance.typing_model_bundle_content_id,
                    self.model_bundle_content_id,
                ),
                "assignment_run_id": (
                    provenance.assignment_run_id,
                    self.assignment_run_id,
                ),
                "candidate_evaluation_manifest_content_id": (
                    provenance.candidate_evaluation_manifest_content_id,
                    self.candidate_evaluation_manifest_content_id,
                ),
                "candidate_assignment_output_content_id": (
                    provenance.candidate_assignment_output_content_id,
                    self.assignment_artifact_content_id,
                ),
            }
            mismatched = [name for name, pair in lineage.items() if pair[0] != pair[1]]
            if mismatched:
                raise ValueError(
                    f"target {target.target_key} review provenance lineage "
                    f"mismatch: {sorted(mismatched)}"
                )
        broad = [target for target in targets if target.level == "broad"]
        if len(broad) != 1:
            raise ValueError("evaluation report requires exactly one broad target")
        specific = [target for target in targets if target.level == "specific"]
        if len(specific) != 1:
            raise ValueError(
                "evaluation report requires exactly one all-cell specific target"
            )
        if not isinstance(self.candidate_evaluation_manifest, Mapping):
            raise TypeError("evaluation candidate_evaluation_manifest must be a JSON object")
        if not isinstance(self.source_run_manifest, Mapping):
            raise TypeError("evaluation source_run_manifest must be a JSON object")
        source_manifest = RunManifest.from_dict(self.source_run_manifest)
        candidate = _candidate_manifest_for_report(self.candidate_evaluation_manifest)
        candidate_lineage = {
            "content_id": (
                candidate.content_id,
                self.candidate_evaluation_manifest_content_id,
            ),
            "model_bundle_content_id": (
                candidate.model_bundle_content_id,
                self.model_bundle_content_id,
            ),
            "split_manifest_content_id": (
                candidate.split_manifest_content_id,
                self.split_manifest_content_id,
            ),
            "split_name": (candidate.split_name, self.split_name),
            "source_run_id": (candidate.source_run_id, self.source_run_id),
            "source_run_manifest_content_id": (
                candidate.source_run_manifest_content_id,
                source_manifest.content_id,
            ),
            "assignment_run_id": (candidate.assignment_run_id, self.assignment_run_id),
            "assignment_artifact_content_id": (
                candidate.output_content_sha256,
                self.assignment_artifact_content_id,
            ),
            "evaluation_donors": (candidate.donors, self.evaluation_donors),
        }
        candidate_mismatches = sorted(
            name for name, pair in candidate_lineage.items() if pair[0] != pair[1]
        )
        if candidate_mismatches:
            raise ValueError(
                "embedded candidate manifest lineage differs from evaluation report: "
                f"{candidate_mismatches}"
            )
        ingest_references = tuple(
            reference
            for reference in source_manifest.stages
            if reference.stage == "ingest"
        )
        if len(ingest_references) != 1:
            raise ValueError(
                "evaluation source run must contain exactly one ingest review context"
            )
        review_context_ids = {ingest_references[0].output_content_sha256}
        for target in targets:
            context_id = target.review_provenance.review_context_artifact_content_id
            if context_id not in review_context_ids:
                raise ValueError(
                    f"target {target.target_key} review context is not an output "
                    "artifact in the bound source run"
                )
        object.__setattr__(self, "candidate_evaluation_manifest", candidate.to_dict())
        object.__setattr__(self, "source_run_manifest", source_manifest.to_dict())
        object.__setattr__(self, "targets", targets)
        replayed_fingerprint = _evaluation_frame_fingerprint(
            method_version=self.method_version,
            split_name=self.split_name,
            model_bundle_content_id=self.model_bundle_content_id,
            split_manifest_content_id=self.split_manifest_content_id,
            source_run_id=self.source_run_id,
            candidate_evaluation_manifest_content_id=(
                self.candidate_evaluation_manifest_content_id
            ),
            assignment_run_id=self.assignment_run_id,
            assignment_artifact_content_id=self.assignment_artifact_content_id,
            evaluation_donors=self.evaluation_donors,
            targets=targets,
        )
        if self.evaluation_frame_fingerprint != replayed_fingerprint:
            raise ValueError(
                "evaluation_frame_fingerprint differs from replayable target evidence"
            )
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("split evaluation report content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def validate_candidate_artifact_binding(self) -> None:
        """Audit reviewed rows against the current exact candidate artifact.

        This donor-streaming audit is intentionally separate from routine
        manifest loading so production startup does not scan a large image
        cohort.  Metric replay from the embedded reviewed rows remains part of
        ordinary construction.
        """

        candidate = _candidate_manifest_for_report(
            self.candidate_evaluation_manifest
        )
        validate_dataset_snapshot_current(candidate.output, mode="content")
        projection_columns = {
            "candidate_assignment_run_id",
            "typing_model_bundle_content_id",
            "broad_label",
            "broad_assignment_status",
            "broad_reason",
            "broad_threshold_eligible",
        }
        if any(target.level in {"subtype", "specific"} for target in self.targets):
            projection_columns.update(
                {
                    "subtype_label",
                    "subtype_assignment_status",
                    "subtype_reason",
                    "subtype_threshold_eligible",
                    "typing_subtype_model_source",
                    "specific_type",
                    "specific_assignment_status",
                }
            )
        population_donors: dict[tuple[str, str | None], list[str]] = {
            target.target_key: [] for target in self.targets
        }
        for donor_snapshot in candidate.output.donors:
            assignments = _candidate_assignment_partition_for_report(
                donor_snapshot,
                columns=tuple(sorted(projection_columns)),
            )
            constants = {
                "candidate_assignment_run_id": candidate.assignment_run_id,
                "typing_model_bundle_content_id": candidate.model_bundle_content_id,
            }
            for column, expected_value in constants.items():
                values = assignments[column].astype("string").fillna("").str.strip()
                if not values.eq(expected_value).all():
                    raise ValueError(
                        f"embedded candidate partition has inconsistent {column}"
                    )
            for target in self.targets:
                if _verify_target_evidence_against_candidate(target, assignments):
                    population_donors[target.target_key].append(donor_snapshot.donor_id)
        for target in self.targets:
            expected_donors = tuple(sorted(population_donors[target.target_key]))
            if target.donors != expected_donors:
                raise ValueError(
                    f"target {target.target_key} donors differ from donors with a "
                    "nonempty candidate target population"
                )

    def validate_candidate_source_replay(
        self,
        *,
        model_bundle: TypingModelBundle,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
        review_blinding_key: str | bytes | None = None,
    ) -> None:
        """Replay the exact candidate output from its immutable source evidence."""

        from .candidate_evaluation import CandidateEvaluationManifest

        candidate = CandidateEvaluationManifest.from_dict(
            self.candidate_evaluation_manifest
        )
        source_run = RunManifest.from_dict(self.source_run_manifest)
        self.validate_candidate_artifact_binding()
        candidate.validate_current(
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
            source_run_manifest=source_run,
            allow_locked_test=self.split_name == "locked_test",
        )
        if review_blinding_key is None:
            return
        from .refinement_contracts import validate_review_sampling_against_assignments

        sampling_columns = {
            "broad_label",
            "broad_assignment_status",
            "subtype_label",
            "subtype_assignment_status",
            "typing_subtype_model_source",
            "specific_type",
            "specific_assignment_status",
        }
        for target in self.targets:
            sample_frame = pd.DataFrame(
                [row.sample_dict() for row in target.evaluated_rows]
            )
            embedded_reviewer = pd.DataFrame(
                [row.reviewer_dict() for row in target.evaluated_rows]
            )
            context_columns = tuple(
                column
                for column in embedded_reviewer.columns
                if column not in {"donor_id", "object_id", "sampling_stratum"}
            )
            replayed_reviewer = blinded_review_sample_from_dataset(
                sample_frame,
                source_run_manifest=source_run,
                splits=model_bundle.split_manifest,
                location_columns=context_columns,
            )
            if reviewer_frame_fingerprint(replayed_reviewer) != (
                target.review_provenance.reviewer_frame_fingerprint
            ):
                raise ValueError(
                    f"target {target.target_key} reviewer context differs from "
                    "the exact ingest artifact"
                )
            expected_records = embedded_reviewer.sort_values(
                ["donor_id", "object_id"], kind="mergesort", ignore_index=True
            ).to_dict(orient="records")
            replayed_records = replayed_reviewer.sort_values(
                ["donor_id", "object_id"], kind="mergesort", ignore_index=True
            ).to_dict(orient="records")
            if replayed_records != expected_records:
                raise ValueError(
                    f"target {target.target_key} embedded reviewer records differ "
                    "from exact ingest replay"
                )
        for donor_snapshot in candidate.output.donors:
            assignments = _candidate_assignment_partition_for_report(
                donor_snapshot,
                columns=tuple(sorted(sampling_columns)),
            )
            for target in self.targets:
                if donor_snapshot.donor_id not in target.donors:
                    continue
                sample_frame = pd.DataFrame(
                    [
                        row.sample_dict()
                        for row in target.evaluated_rows
                        if row.donor_id == donor_snapshot.donor_id
                    ]
                )
                validate_review_sampling_against_assignments(
                    sample_frame,
                    assignments,
                    splits=model_bundle.split_manifest,
                    split_name=self.split_name,
                    candidate_assignment_run_id=candidate.assignment_run_id,
                    level=target.level,
                    parent=target.parent,
                    blinding_key=review_blinding_key,
                    expected_donors=(donor_snapshot.donor_id,),
                )

    def validate_for_promotion(
        self,
        *,
        model_bundle: TypingModelBundle,
        policy: PromotionGatePolicy,
    ) -> None:
        """Fail closed unless this report passes the exact production policy.

        This public check is the supported validation-stage gate before an
        operator is allowed to open the locked-test split.  It performs no
        source reads; use :meth:`validate_candidate_source_replay` for the
        separate donor-streaming artifact audit.
        """

        if not isinstance(model_bundle, TypingModelBundle):
            raise TypeError("model_bundle must be TypingModelBundle")
        if not isinstance(policy, PromotionGatePolicy):
            raise TypeError("policy must be PromotionGatePolicy")
        split = model_bundle.split_manifest
        if not split.promotion_policy_content_id or (
            split.promotion_policy_content_id != policy.content_id
        ):
            raise ValueError(
                "promotion policy differs from donor-split preregistration"
            )
        if self.promotion_policy.content_id != policy.content_id:
            raise ValueError("evaluation report uses a different promotion policy")
        if self.model_bundle_content_id != model_bundle.content_id:
            raise ValueError("evaluation report refers to a different model bundle")
        if self.split_manifest_content_id != split.content_id:
            raise ValueError("evaluation report refers to a different donor split")
        if self.source_run_id != split.source_run_id:
            raise ValueError("evaluation report refers to a different source run")
        expected_donors = tuple(getattr(split, self.split_name))
        if self.evaluation_donors != expected_donors:
            raise ValueError(
                f"{self.split_name} evaluation donors do not exactly cover the split"
            )
        expected_models = {
            ("broad", None): model_bundle.broad.classes,
            ("specific", None): tuple(
                model_bundle.fit_parameters["final_specific_classes"]
            ),
            **{
                ("subtype", model.parent): model.classes
                for model in model_bundle.subtypes
            },
        }
        observed = {target.target_key: target for target in self.targets}
        missing = sorted(
            set(expected_models) - set(observed), key=lambda item: str(item)
        )
        extra = sorted(
            set(observed) - set(expected_models), key=lambda item: str(item)
        )
        if missing or extra:
            raise ValueError(
                f"{self.split_name} model-target coverage is incomplete: "
                f"missing={missing}, extra={extra}"
            )
        if self.bootstrap_unit != policy.bootstrap_unit:
            raise ValueError("evaluation bootstrap unit differs from promotion policy")
        if not np.isclose(
            self.confidence_level,
            _promotion_report_confidence_level(policy.confidence_level),
            rtol=0.0,
            atol=1e-12,
        ):
            raise ValueError("evaluation confidence level differs from promotion policy")
        if self.n_bootstrap != policy.bootstrap_replicates:
            raise ValueError(
                "evaluation bootstrap count differs from the precommitted policy"
            )
        for key, expected_classes in expected_models.items():
            target = observed[key]
            if target.classes != expected_classes:
                raise ValueError(
                    f"{self.split_name} target {key} classes differ from the model"
                )
            unknown_donors = sorted(set(target.donors) - set(expected_donors))
            if unknown_donors:
                raise ValueError(
                    f"{self.split_name} target {key} contains wrong-split donors: "
                    f"{unknown_donors}"
                )
            if len(target.donors) < policy.minimum_donors:
                raise ValueError(
                    f"{self.split_name} target {key} has too few reviewed donors"
                )
            if target.row_count < policy.minimum_rows:
                raise ValueError(
                    f"{self.split_name} target {key} has too few reviewed rows"
                )
            failed: list[str] = []
            if target.coverage.lower < policy.minimum_coverage_lower:
                failed.append("coverage")
            if target.selective_error.upper > policy.maximum_selective_error_upper:
                failed.append("selective_error")
            if target.rescue_coverage.lower < policy.minimum_rescue_coverage_lower:
                failed.append("rescue_coverage")
            if (
                target.rescue_selective_error.upper
                > policy.maximum_rescue_selective_error_upper
            ):
                failed.append("rescue_selective_error")
            for class_error in target.class_errors:
                if (
                    class_error.reference_positive_rows
                    < policy.minimum_class_positive_rows
                ):
                    failed.append(f"{class_error.label}:positive_rows")
                if (
                    class_error.reference_negative_rows
                    < policy.minimum_class_negative_rows
                ):
                    failed.append(f"{class_error.label}:negative_rows")
                if (
                    len(class_error.reference_positive_donors)
                    < policy.minimum_class_positive_donors
                ):
                    failed.append(f"{class_error.label}:positive_donors")
                if (
                    len(class_error.reference_negative_donors)
                    < policy.minimum_class_negative_donors
                ):
                    failed.append(f"{class_error.label}:negative_donors")
                if (
                    class_error.false_positive_rate.upper
                    > policy.maximum_class_false_positive_rate_upper
                ):
                    failed.append(f"{class_error.label}:false_positive_rate")
                if (
                    class_error.false_negative_rate.upper
                    > policy.maximum_class_false_negative_rate_upper
                ):
                    failed.append(f"{class_error.label}:false_negative_rate")
            if failed:
                raise ValueError(
                    f"{self.split_name} target {key} fails promotion policy: "
                    f"{sorted(failed)}"
                )

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "method_version": self.method_version,
            "split_name": self.split_name,
            "model_bundle_content_id": self.model_bundle_content_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "source_run_id": self.source_run_id,
            "source_run_manifest": dict(self.source_run_manifest),
            "candidate_evaluation_manifest_content_id": (
                self.candidate_evaluation_manifest_content_id
            ),
            "candidate_evaluation_manifest": dict(self.candidate_evaluation_manifest),
            "assignment_run_id": self.assignment_run_id,
            "assignment_artifact_content_id": self.assignment_artifact_content_id,
            "evaluation_frame_fingerprint": self.evaluation_frame_fingerprint,
            "evaluation_donors": list(self.evaluation_donors),
            "promotion_policy": self.promotion_policy.to_dict(),
            "confidence_level": self.confidence_level,
            "n_bootstrap": self.n_bootstrap,
            "random_state": self.random_state,
            "bootstrap_unit": self.bootstrap_unit,
            "targets": [value.to_dict() for value in self.targets],
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    def write_json(self, path: Path | str) -> Path:
        """Persist this frozen held-out report as strict JSON."""

        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def from_target_results(
        cls,
        *,
        split_name: str,
        model_bundle_content_id: str,
        split_manifest_content_id: str,
        source_run_id: str,
        source_run_manifest: Mapping[str, Any],
        candidate_evaluation_manifest_content_id: str,
        candidate_evaluation_manifest: Mapping[str, Any],
        assignment_run_id: str,
        assignment_artifact_content_id: str,
        evaluation_donors: Sequence[str],
        promotion_policy: PromotionGatePolicy,
        targets: Sequence[TargetEvaluationResult],
    ) -> SplitEvaluationReport:
        """Construct a report while deriving all replay and seed identities."""

        donors = tuple(sorted(str(value).strip() for value in evaluation_donors))
        target_values = tuple(targets)
        random_state = evaluation_bootstrap_random_state(
            split_name=split_name,
            model_bundle_content_id=model_bundle_content_id,
            assignment_run_id=assignment_run_id,
        )
        fingerprint = _evaluation_frame_fingerprint(
            method_version=EVALUATION_REPORT_METHOD_VERSION,
            split_name=split_name,
            model_bundle_content_id=model_bundle_content_id,
            split_manifest_content_id=split_manifest_content_id,
            source_run_id=source_run_id,
            candidate_evaluation_manifest_content_id=(
                candidate_evaluation_manifest_content_id
            ),
            assignment_run_id=assignment_run_id,
            assignment_artifact_content_id=assignment_artifact_content_id,
            evaluation_donors=donors,
            targets=target_values,
        )
        return cls(
            split_name=split_name,
            model_bundle_content_id=model_bundle_content_id,
            split_manifest_content_id=split_manifest_content_id,
            source_run_id=source_run_id,
            source_run_manifest=source_run_manifest,
            candidate_evaluation_manifest_content_id=(
                candidate_evaluation_manifest_content_id
            ),
            candidate_evaluation_manifest=candidate_evaluation_manifest,
            assignment_run_id=assignment_run_id,
            assignment_artifact_content_id=assignment_artifact_content_id,
            evaluation_frame_fingerprint=fingerprint,
            evaluation_donors=donors,
            promotion_policy=promotion_policy,
            confidence_level=_promotion_report_confidence_level(
                promotion_policy.confidence_level
            ),
            n_bootstrap=promotion_policy.bootstrap_replicates,
            random_state=random_state,
            targets=target_values,
        )

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> SplitEvaluationReport:
        if not isinstance(value, Mapping):
            raise TypeError("split evaluation report must be a mapping")
        fields = {
            "schema_version",
            "method_version",
            "split_name",
            "model_bundle_content_id",
            "split_manifest_content_id",
            "source_run_id",
            "source_run_manifest",
            "candidate_evaluation_manifest_content_id",
            "candidate_evaluation_manifest",
            "assignment_run_id",
            "assignment_artifact_content_id",
            "evaluation_frame_fingerprint",
            "evaluation_donors",
            "promotion_policy",
            "confidence_level",
            "n_bootstrap",
            "random_state",
            "bootstrap_unit",
            "targets",
            "content_id",
        }
        _require_exact_fields(value, fields, where="split evaluation report")
        return cls(
            schema_version=value["schema_version"],
            method_version=value["method_version"],
            split_name=value["split_name"],
            model_bundle_content_id=value["model_bundle_content_id"],
            split_manifest_content_id=value["split_manifest_content_id"],
            source_run_id=value["source_run_id"],
            source_run_manifest=value["source_run_manifest"],
            candidate_evaluation_manifest_content_id=value[
                "candidate_evaluation_manifest_content_id"
            ],
            candidate_evaluation_manifest=value["candidate_evaluation_manifest"],
            assignment_run_id=value["assignment_run_id"],
            assignment_artifact_content_id=value["assignment_artifact_content_id"],
            evaluation_frame_fingerprint=value["evaluation_frame_fingerprint"],
            evaluation_donors=tuple(value["evaluation_donors"]),
            promotion_policy=PromotionGatePolicy.from_dict(value["promotion_policy"]),
            confidence_level=value["confidence_level"],
            n_bootstrap=value["n_bootstrap"],
            random_state=value["random_state"],
            bootstrap_unit=value["bootstrap_unit"],
            targets=tuple(TargetEvaluationResult.from_dict(item) for item in value["targets"]),
            content_id=value["content_id"],
        )

    @classmethod
    def read_json(cls, path: Path | str) -> SplitEvaluationReport:
        """Load and fully replay a frozen held-out report from strict JSON."""

        return cls.from_dict(_read_strict_json(path))


def _release_key_bytes(value: str | bytes) -> bytes:
    if isinstance(value, str):
        key = value.encode("utf-8")
    elif isinstance(value, bytes):
        key = value
    else:
        raise TypeError("release signing key must be text or bytes")
    if len(key) < 16:
        raise ValueError("release signing key must contain at least 16 bytes")
    return key


def _canonical_hmac_payload(value: Mapping[str, Any]) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


@dataclass(frozen=True)
class PromotionReplayReceipt:
    """Deployment-key attestation that the full release replay succeeded."""

    typing_model_bundle_content_id: str
    model_fit_content_id: str
    split_manifest_content_id: str
    threshold_selection_provenance_content_id: str
    preliminary_bundle_content_id: str
    calibration_candidate_manifest_content_id: str
    calibration_candidate_output_content_id: str
    calibration_assignment_run_id: str
    calibration_source_run_manifest_content_id: str
    calibration_source_stage_manifest_content_id: str
    calibration_source_artifact_content_id: str
    calibration_review_context_artifact_content_id: str
    promotion_policy_content_id: str
    validation_report_content_id: str
    locked_test_report_content_id: str
    marker_registry_fingerprint: str
    typing_rules_fingerprint: str
    replay_code_sha256: str
    replay_environment_sha256: str
    issuer: str
    issued_at_utc: str
    replay_method_version: str = EVALUATION_REPORT_METHOD_VERSION
    method_version: str = PROMOTION_REPLAY_RECEIPT_METHOD_VERSION
    schema_version: int = PROMOTION_COMPONENT_SCHEMA_VERSION
    signature_hmac_sha256: str = ""
    content_id: str = ""

    def __post_init__(self) -> None:
        if (
            isinstance(self.schema_version, bool)
            or self.schema_version != PROMOTION_COMPONENT_SCHEMA_VERSION
        ):
            raise ValueError(
                f"unsupported promotion replay receipt schema {self.schema_version!r}"
            )
        if self.method_version != PROMOTION_REPLAY_RECEIPT_METHOD_VERSION:
            raise ValueError("unsupported promotion replay receipt method")
        if self.replay_method_version != EVALUATION_REPORT_METHOD_VERSION:
            raise ValueError("promotion receipt uses a stale replay implementation")
        for name in (
            "typing_model_bundle_content_id",
            "model_fit_content_id",
            "split_manifest_content_id",
            "threshold_selection_provenance_content_id",
            "preliminary_bundle_content_id",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
            "calibration_assignment_run_id",
            "calibration_source_run_manifest_content_id",
            "calibration_source_stage_manifest_content_id",
            "calibration_source_artifact_content_id",
            "calibration_review_context_artifact_content_id",
            "promotion_policy_content_id",
            "validation_report_content_id",
            "locked_test_report_content_id",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
            "replay_code_sha256",
            "replay_environment_sha256",
        ):
            object.__setattr__(
                self, name, _require_sha256(getattr(self, name), name=name)
            )
        issuer = str(self.issuer).strip()
        if not issuer:
            raise ValueError("promotion replay receipt issuer must be non-empty")
        object.__setattr__(self, "issuer", issuer)
        issued_at = str(self.issued_at_utc).strip()
        try:
            parsed = datetime.fromisoformat(issued_at.replace("Z", "+00:00"))
        except ValueError as exc:
            raise ValueError("receipt issued_at_utc must be an ISO timestamp") from exc
        if parsed.tzinfo is None or parsed.utcoffset() is None:
            raise ValueError("receipt issued_at_utc must include a timezone")
        if parsed.utcoffset().total_seconds() != 0:
            raise ValueError("receipt issued_at_utc must use UTC")
        if issued_at != parsed.isoformat():
            raise ValueError("receipt issued_at_utc must use canonical ISO UTC form")
        object.__setattr__(self, "issued_at_utc", issued_at)
        signature = _require_sha256(
            self.signature_hmac_sha256,
            name="receipt.signature_hmac_sha256",
        )
        object.__setattr__(self, "signature_hmac_sha256", signature)
        expected = sha256_json(
            {
                "identity": self.identity_payload(),
                "signature_hmac_sha256": signature,
            }
        )
        if self.content_id and self.content_id != expected:
            raise ValueError("promotion replay receipt content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "method_version": self.method_version,
            "replay_method_version": self.replay_method_version,
            "typing_model_bundle_content_id": self.typing_model_bundle_content_id,
            "model_fit_content_id": self.model_fit_content_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "threshold_selection_provenance_content_id": (
                self.threshold_selection_provenance_content_id
            ),
            "preliminary_bundle_content_id": self.preliminary_bundle_content_id,
            "calibration_candidate_manifest_content_id": (
                self.calibration_candidate_manifest_content_id
            ),
            "calibration_candidate_output_content_id": (
                self.calibration_candidate_output_content_id
            ),
            "calibration_assignment_run_id": self.calibration_assignment_run_id,
            "calibration_source_run_manifest_content_id": (
                self.calibration_source_run_manifest_content_id
            ),
            "calibration_source_stage_manifest_content_id": (
                self.calibration_source_stage_manifest_content_id
            ),
            "calibration_source_artifact_content_id": (
                self.calibration_source_artifact_content_id
            ),
            "calibration_review_context_artifact_content_id": (
                self.calibration_review_context_artifact_content_id
            ),
            "promotion_policy_content_id": self.promotion_policy_content_id,
            "validation_report_content_id": self.validation_report_content_id,
            "locked_test_report_content_id": self.locked_test_report_content_id,
            "marker_registry_fingerprint": self.marker_registry_fingerprint,
            "typing_rules_fingerprint": self.typing_rules_fingerprint,
            "replay_code_sha256": self.replay_code_sha256,
            "replay_environment_sha256": self.replay_environment_sha256,
            "issuer": self.issuer,
            "issued_at_utc": self.issued_at_utc,
        }

    def to_dict(self) -> dict[str, Any]:
        return {
            **self.identity_payload(),
            "signature_hmac_sha256": self.signature_hmac_sha256,
            "content_id": self.content_id,
        }

    @classmethod
    def create(
        cls,
        *,
        model_bundle: TypingModelBundle,
        policy: PromotionGatePolicy,
        validation_report: SplitEvaluationReport,
        locked_test_report: SplitEvaluationReport,
        threshold_source_replay_proof: ThresholdSelectionSourceReplayProof,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
        release_signing_key: str | bytes,
        issuer: str,
        issued_at_utc: str,
    ) -> PromotionReplayReceipt:
        provenance = model_bundle.threshold_selection_provenance
        if provenance is None:
            raise ValueError("receipt requires frozen threshold-selection provenance")
        if not isinstance(
            threshold_source_replay_proof, ThresholdSelectionSourceReplayProof
        ):
            raise TypeError(
                "receipt requires a successful ThresholdSelectionSourceReplayProof"
            )
        proof = threshold_source_replay_proof
        expected_proof = {
            "threshold_selection_provenance_content_id": provenance.content_id,
            "preliminary_bundle_content_id": provenance.preliminary_bundle_content_id,
            "model_fit_content_id": model_bundle.model_fit_content_id,
            "calibration_candidate_manifest_content_id": (
                provenance.calibration_candidate_manifest_content_id
            ),
            "calibration_candidate_output_content_id": (
                provenance.calibration_candidate_output_content_id
            ),
            "assignment_run_id": provenance.assignment_run_id,
            "source_run_manifest_content_id": (
                model_bundle.evidence_provenance.source_run_manifest_content_id
            ),
            "source_stage_manifest_content_id": (
                model_bundle.evidence_provenance.source_stage_manifest_content_id
            ),
            "source_artifact_content_id": (
                model_bundle.evidence_provenance.source_artifact_content_id
            ),
        }
        proof_mismatch = sorted(
            name
            for name, expected_value in expected_proof.items()
            if getattr(proof, name) != expected_value
        )
        review_context_ids = {
            item.review_context_artifact_content_id
            for item in provenance.calibration_review_provenance
        }
        if proof_mismatch or review_context_ids != {
            proof.review_context_artifact_content_id
        }:
            raise ValueError(
                "threshold source replay proof differs from the frozen calibration "
                f"decision: fields={proof_mismatch}"
            )
        key = _release_key_bytes(release_signing_key)
        commitment = hashlib.sha256(key).hexdigest()
        expected_commitment = model_bundle.split_manifest.release_signing_key_sha256
        if not expected_commitment or not hmac.compare_digest(
            commitment, expected_commitment
        ):
            raise ValueError(
                "release signing key differs from the preregistered split commitment"
            )
        unsigned = {
            "schema_version": PROMOTION_COMPONENT_SCHEMA_VERSION,
            "method_version": PROMOTION_REPLAY_RECEIPT_METHOD_VERSION,
            "replay_method_version": EVALUATION_REPORT_METHOD_VERSION,
            "typing_model_bundle_content_id": model_bundle.content_id,
            "model_fit_content_id": model_bundle.model_fit_content_id,
            "split_manifest_content_id": model_bundle.split_manifest.content_id,
            "threshold_selection_provenance_content_id": provenance.content_id,
            "preliminary_bundle_content_id": proof.preliminary_bundle_content_id,
            "calibration_candidate_manifest_content_id": (
                proof.calibration_candidate_manifest_content_id
            ),
            "calibration_candidate_output_content_id": (
                proof.calibration_candidate_output_content_id
            ),
            "calibration_assignment_run_id": proof.assignment_run_id,
            "calibration_source_run_manifest_content_id": (
                proof.source_run_manifest_content_id
            ),
            "calibration_source_stage_manifest_content_id": (
                proof.source_stage_manifest_content_id
            ),
            "calibration_source_artifact_content_id": proof.source_artifact_content_id,
            "calibration_review_context_artifact_content_id": (
                proof.review_context_artifact_content_id
            ),
            "promotion_policy_content_id": policy.content_id,
            "validation_report_content_id": validation_report.content_id,
            "locked_test_report_content_id": locked_test_report.content_id,
            "marker_registry_fingerprint": marker_registry.fingerprint,
            "typing_rules_fingerprint": typing_registry.rules_fingerprint,
            "replay_code_sha256": promotion_replay_code_sha256(),
            "replay_environment_sha256": promotion_replay_environment_sha256(),
            "issuer": str(issuer).strip(),
            "issued_at_utc": str(issued_at_utc).strip(),
        }
        signature = hmac.new(
            key,
            _canonical_hmac_payload(unsigned),
            hashlib.sha256,
        ).hexdigest()
        return cls(**unsigned, signature_hmac_sha256=signature)

    def verify(self, release_signing_key: str | bytes) -> None:
        key = _release_key_bytes(release_signing_key)
        expected = hmac.new(
            key,
            _canonical_hmac_payload(self.identity_payload()),
            hashlib.sha256,
        ).hexdigest()
        if not hmac.compare_digest(expected, self.signature_hmac_sha256):
            raise ValueError("promotion replay receipt signature is invalid")

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> PromotionReplayReceipt:
        if not isinstance(value, Mapping):
            raise TypeError("promotion replay receipt must be a mapping")
        fields = {
            "schema_version",
            "method_version",
            "replay_method_version",
            "typing_model_bundle_content_id",
            "model_fit_content_id",
            "split_manifest_content_id",
            "threshold_selection_provenance_content_id",
            "preliminary_bundle_content_id",
            "calibration_candidate_manifest_content_id",
            "calibration_candidate_output_content_id",
            "calibration_assignment_run_id",
            "calibration_source_run_manifest_content_id",
            "calibration_source_stage_manifest_content_id",
            "calibration_source_artifact_content_id",
            "calibration_review_context_artifact_content_id",
            "promotion_policy_content_id",
            "validation_report_content_id",
            "locked_test_report_content_id",
            "marker_registry_fingerprint",
            "typing_rules_fingerprint",
            "replay_code_sha256",
            "replay_environment_sha256",
            "issuer",
            "issued_at_utc",
            "signature_hmac_sha256",
            "content_id",
        }
        _require_exact_fields(value, fields, where="promotion replay receipt")
        return cls(**{name: value[name] for name in fields})


@dataclass(frozen=True)
class TypingModelPromotionManifest:
    """Verified release approval for one exact fitted model bundle."""

    promotion_version: str
    model_bundle_content_id: str
    split_manifest_content_id: str
    policy: PromotionGatePolicy
    validation_report: SplitEvaluationReport
    locked_test_report: SplitEvaluationReport
    replay_receipt: PromotionReplayReceipt
    decision: str = PROMOTION_DECISION
    schema_version: int = PROMOTION_SCHEMA_VERSION
    content_id: str = ""
    input_artifact: InputArtifact | None = field(default=None, repr=False, compare=False)
    _signature_verified: bool = field(
        default=False, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if isinstance(self.schema_version, bool) or self.schema_version != PROMOTION_SCHEMA_VERSION:
            raise ValueError(
                f"unsupported promotion schema_version {self.schema_version!r}; "
                "schema v1 summary hashes are unverifiable and must be regenerated"
            )
        if not isinstance(self.promotion_version, str) or not self.promotion_version.strip():
            raise ValueError("promotion_version must be a non-empty string")
        object.__setattr__(self, "promotion_version", self.promotion_version.strip())
        if self.decision != PROMOTION_DECISION:
            raise ValueError(f"promotion decision must be {PROMOTION_DECISION!r}")
        for name in ("model_bundle_content_id", "split_manifest_content_id"):
            object.__setattr__(self, name, _require_sha256(getattr(self, name), name=name))
        if not isinstance(self.policy, PromotionGatePolicy):
            raise TypeError("policy must be PromotionGatePolicy")
        if not isinstance(self.validation_report, SplitEvaluationReport) or not isinstance(
            self.locked_test_report, SplitEvaluationReport
        ):
            raise TypeError("promotion reports must be SplitEvaluationReport objects")
        if not isinstance(self.replay_receipt, PromotionReplayReceipt):
            raise TypeError("promotion replay_receipt must be PromotionReplayReceipt")
        if self.validation_report.split_name != "validation":
            raise ValueError("validation report has the wrong split")
        if self.locked_test_report.split_name != "locked_test":
            raise ValueError("locked-test report has the wrong split")
        for report in (self.validation_report, self.locked_test_report):
            if report.promotion_policy.content_id != self.policy.content_id:
                raise ValueError(
                    "promotion report uses a different precommitted gate policy"
                )
            if report.model_bundle_content_id != self.model_bundle_content_id:
                raise ValueError("promotion report refers to a different model bundle")
            if report.split_manifest_content_id != self.split_manifest_content_id:
                raise ValueError("promotion report refers to a different donor split")
        receipt_lineage = {
            "typing_model_bundle_content_id": (
                self.replay_receipt.typing_model_bundle_content_id,
                self.model_bundle_content_id,
            ),
            "split_manifest_content_id": (
                self.replay_receipt.split_manifest_content_id,
                self.split_manifest_content_id,
            ),
            "promotion_policy_content_id": (
                self.replay_receipt.promotion_policy_content_id,
                self.policy.content_id,
            ),
            "validation_report_content_id": (
                self.replay_receipt.validation_report_content_id,
                self.validation_report.content_id,
            ),
            "locked_test_report_content_id": (
                self.replay_receipt.locked_test_report_content_id,
                self.locked_test_report.content_id,
            ),
        }
        mismatched_receipt = sorted(
            name for name, pair in receipt_lineage.items() if pair[0] != pair[1]
        )
        if mismatched_receipt:
            raise ValueError(
                "promotion replay receipt lineage mismatch: "
                f"{mismatched_receipt}"
            )
        expected = sha256_json(self.identity_payload())
        if self.content_id and self.content_id != expected:
            raise ValueError("typing promotion manifest content_id is invalid")
        object.__setattr__(self, "content_id", expected)

    def identity_payload(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "promotion_version": self.promotion_version,
            "decision": self.decision,
            "model_bundle_content_id": self.model_bundle_content_id,
            "split_manifest_content_id": self.split_manifest_content_id,
            "policy": self.policy.to_dict(),
            "validation_report": self.validation_report.to_dict(),
            "locked_test_report": self.locked_test_report.to_dict(),
            "replay_receipt": self.replay_receipt.to_dict(),
        }

    def to_dict(self) -> dict[str, Any]:
        return {**self.identity_payload(), "content_id": self.content_id}

    @classmethod
    def create(
        cls,
        *,
        promotion_version: str,
        model_bundle: TypingModelBundle,
        policy: PromotionGatePolicy,
        validation_report: SplitEvaluationReport,
        locked_test_report: SplitEvaluationReport,
        threshold_source_replay_inputs: ThresholdSelectionReplayInputs,
        marker_registry: MarkerRegistry,
        typing_registry: TypingRegistry,
        review_blinding_key: str | bytes,
        release_signing_key: str | bytes,
        issuer: str,
        issued_at_utc: str,
    ) -> TypingModelPromotionManifest:
        model_bundle.validate_current(
            marker_registry=marker_registry,
            typing_registry=typing_registry,
        )
        if model_bundle.split_manifest.promotion_policy_content_id != policy.content_id:
            raise ValueError("promotion policy was not preregistered in the donor split")
        threshold_source_replay_proof = validate_threshold_selection_source_replay(
            model_bundle,
            threshold_source_replay_inputs,
        )
        for report in (validation_report, locked_test_report):
            report.validate_candidate_source_replay(
                model_bundle=model_bundle,
                marker_registry=marker_registry,
                typing_registry=typing_registry,
                review_blinding_key=review_blinding_key,
            )
        receipt = PromotionReplayReceipt.create(
            model_bundle=model_bundle,
            policy=policy,
            validation_report=validation_report,
            locked_test_report=locked_test_report,
            threshold_source_replay_proof=threshold_source_replay_proof,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
            release_signing_key=release_signing_key,
            issuer=issuer,
            issued_at_utc=issued_at_utc,
        )
        result = cls(
            promotion_version=promotion_version,
            model_bundle_content_id=model_bundle.content_id,
            split_manifest_content_id=model_bundle.split_manifest.content_id,
            policy=policy,
            validation_report=validation_report,
            locked_test_report=locked_test_report,
            replay_receipt=receipt,
        )
        result.validate_for_bundle(
            model_bundle,
            release_signing_key=release_signing_key,
        )
        object.__setattr__(result, "_signature_verified", True)
        return result

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> TypingModelPromotionManifest:
        if not isinstance(value, Mapping):
            raise TypeError("typing promotion manifest must be a mapping")
        schema_version = value.get("schema_version")
        if schema_version == 1:
            raise ValueError(
                "promotion schema v1 contains unverifiable summary hashes; "
                "regenerate a schema v3 promotion"
            )
        allowed = {
            "schema_version",
            "promotion_version",
            "decision",
            "model_bundle_content_id",
            "split_manifest_content_id",
            "policy",
            "validation_report",
            "locked_test_report",
            "replay_receipt",
            "content_id",
        }
        _require_exact_fields(value, allowed, where="typing promotion manifest")
        if not str(value.get("content_id", "")).strip():
            raise ValueError("typing promotion manifest is missing content_id")
        if isinstance(schema_version, bool) or not isinstance(schema_version, int):
            raise TypeError("promotion schema_version must be a JSON integer")
        return cls(
            schema_version=schema_version,
            promotion_version=value["promotion_version"],
            decision=value["decision"],
            model_bundle_content_id=value["model_bundle_content_id"],
            split_manifest_content_id=value["split_manifest_content_id"],
            policy=PromotionGatePolicy.from_dict(value["policy"]),
            validation_report=SplitEvaluationReport.from_dict(value["validation_report"]),
            locked_test_report=SplitEvaluationReport.from_dict(value["locked_test_report"]),
            replay_receipt=PromotionReplayReceipt.from_dict(value["replay_receipt"]),
            content_id=value["content_id"],
        )

    def write_json(self, path: Path | str) -> Path:
        return _atomic_write_json(Path(path), self.to_dict())

    @classmethod
    def read_json(cls, path: Path | str) -> TypingModelPromotionManifest:
        path = Path(path).expanduser().resolve()
        artifact = InputArtifact.from_path("typing_model_promotion", path)
        result = cls.from_dict(_read_strict_json(path))
        artifact.validate_current(mode="fast")
        return replace(result, input_artifact=artifact)

    def validate_for_bundle(
        self,
        model_bundle: TypingModelBundle,
        *,
        marker_registry: MarkerRegistry | None = None,
        typing_registry: TypingRegistry | None = None,
        review_blinding_key: str | bytes | None = None,
        release_signing_key: str | bytes | None = None,
        require_source_replay: bool = False,
    ) -> None:
        if self.input_artifact is not None:
            self.input_artifact.validate_current(mode="content")
        if self.model_bundle_content_id != model_bundle.content_id:
            raise ValueError("typing promotion refers to a different model bundle")
        if model_bundle.threshold_selection_provenance is None:
            raise ValueError(
                "production promotion requires complete threshold-selection "
                "provenance from the calibration split"
            )
        missing_subtypes = model_bundle.fit_parameters[
            "missing_probabilistic_subtype_parents"
        ]
        if missing_subtypes:
            raise ValueError(
                "production promotion requires probabilistic models for every "
                "multi-class subtype parent; missing="
                f"{missing_subtypes}"
            )
        split = model_bundle.split_manifest
        if self.split_manifest_content_id != split.content_id:
            raise ValueError("typing promotion refers to a different donor split")
        if not split.validation or not split.locked_test:
            raise ValueError(
                "production promotion requires non-empty validation and locked-test splits"
            )
        if not split.promotion_policy_content_id or (
            split.promotion_policy_content_id != self.policy.content_id
        ):
            raise ValueError(
                "production promotion policy differs from donor-split preregistration"
            )
        if not split.release_signing_key_sha256:
            raise ValueError("production split has no release signing-key commitment")
        key: bytes | None = None
        if release_signing_key is None:
            if not self._signature_verified:
                raise ValueError(
                    "production promotion validation requires either the external "
                    "release key or a loader-verified in-memory capability"
                )
        else:
            key = _release_key_bytes(release_signing_key)
            if not hmac.compare_digest(
                hashlib.sha256(key).hexdigest(), split.release_signing_key_sha256
            ):
                raise ValueError("release signing key differs from donor-split commitment")
        provenance = model_bundle.threshold_selection_provenance
        assert provenance is not None
        receipt_expected = {
            "typing_model_bundle_content_id": model_bundle.content_id,
            "model_fit_content_id": model_bundle.model_fit_content_id,
            "split_manifest_content_id": split.content_id,
            "threshold_selection_provenance_content_id": provenance.content_id,
            "preliminary_bundle_content_id": provenance.preliminary_bundle_content_id,
            "calibration_candidate_manifest_content_id": (
                provenance.calibration_candidate_manifest_content_id
            ),
            "calibration_candidate_output_content_id": (
                provenance.calibration_candidate_output_content_id
            ),
            "calibration_assignment_run_id": provenance.assignment_run_id,
            "calibration_source_run_manifest_content_id": (
                model_bundle.evidence_provenance.source_run_manifest_content_id
            ),
            "calibration_source_stage_manifest_content_id": (
                model_bundle.evidence_provenance.source_stage_manifest_content_id
            ),
            "calibration_source_artifact_content_id": (
                model_bundle.evidence_provenance.source_artifact_content_id
            ),
            "promotion_policy_content_id": self.policy.content_id,
            "validation_report_content_id": self.validation_report.content_id,
            "locked_test_report_content_id": self.locked_test_report.content_id,
            "marker_registry_fingerprint": model_bundle.marker_registry_fingerprint,
            "typing_rules_fingerprint": model_bundle.typing_rules_fingerprint,
            "replay_code_sha256": promotion_replay_code_sha256(),
            "replay_environment_sha256": promotion_replay_environment_sha256(),
        }
        calibration_review_context_ids = {
            item.review_context_artifact_content_id
            for item in provenance.calibration_review_provenance
        }
        if len(calibration_review_context_ids) != 1:
            raise ValueError(
                "calibration reviews do not share one source-ingest context"
            )
        receipt_expected["calibration_review_context_artifact_content_id"] = next(
            iter(calibration_review_context_ids)
        )
        receipt_mismatch = sorted(
            name
            for name, expected_value in receipt_expected.items()
            if getattr(self.replay_receipt, name) != expected_value
        )
        if receipt_mismatch:
            raise ValueError(
                "promotion replay receipt differs from the frozen release: "
                f"{receipt_mismatch}"
            )
        if key is not None:
            self.replay_receipt.verify(key)
        if (marker_registry is None) != (typing_registry is None):
            raise TypeError("source replay requires both marker and typing registries")
        if require_source_replay and marker_registry is None:
            raise ValueError("production promotion validation requires source inference replay")
        for report in (self.validation_report, self.locked_test_report):
            self._validate_report(
                report,
                model_bundle=model_bundle,
            )
            if require_source_replay:
                assert marker_registry is not None and typing_registry is not None
                report.validate_candidate_source_replay(
                    model_bundle=model_bundle,
                    marker_registry=marker_registry,
                    typing_registry=typing_registry,
                    review_blinding_key=review_blinding_key,
                )

    def _validate_report(
        self,
        report: SplitEvaluationReport,
        *,
        model_bundle: TypingModelBundle,
    ) -> None:
        report.validate_for_promotion(
            model_bundle=model_bundle,
            policy=self.policy,
        )


def load_typing_model_promotion_manifest(
    path: Path | str,
    *,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry | None = None,
    typing_registry: TypingRegistry | None = None,
    review_blinding_key: str | bytes | None = None,
    release_signing_key: str | bytes,
    audit_source_replay: bool = False,
) -> TypingModelPromotionManifest:
    """Load a release authorization, with optional donor-streaming source audit."""

    promotion = TypingModelPromotionManifest.read_json(path)
    promotion.validate_for_bundle(
        model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        review_blinding_key=review_blinding_key,
        release_signing_key=release_signing_key,
        require_source_replay=audit_source_replay,
    )
    object.__setattr__(promotion, "_signature_verified", True)
    return promotion


def load_typing_model_bundle(
    path: Path | str,
    *,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
) -> TypingModelBundle:
    """Load a bundle, capture its bytes, and validate scientific compatibility."""

    bundle = TypingModelBundle.read_json(path)
    bundle.validate_current(
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )
    return bundle


def ensemble_stability_wide(
    evidence: pd.DataFrame,
    *,
    ensemble: FrozenModelEnsemble,
    chunk_size: int = 250_000,
) -> pd.DataFrame:
    """Compute bootstrap agreement without a cells × bootstraps tensor.

    Production marker evidence is one row per cell. Each chunk is transformed
    once, then integer agreement counts are accumulated across complete
    bootstrap coefficient matrices.
    """

    if not isinstance(evidence, pd.DataFrame):
        raise TypeError("evidence must be a pandas DataFrame")
    if "marker" in evidence:
        raise ValueError("chunked bundle stability requires wide production evidence")
    if isinstance(chunk_size, bool) or int(chunk_size) != chunk_size or chunk_size < 1:
        raise ValueError("chunk_size must be an integer >= 1")
    pieces: list[pd.DataFrame] = []
    base_coefficients = np.asarray(ensemble.coefficients, dtype=float)
    bootstrap = tuple(np.asarray(value, dtype=float) for value in ensemble.bootstrap_coefficients)
    bootstrap_valid = ensemble.bootstrap_valid
    classes = np.asarray(ensemble.classes, dtype=object)
    for start in range(0, len(evidence), int(chunk_size)):
        chunk = evidence.iloc[start : start + int(chunk_size)]
        keys, features, _state, _valid = evidence_feature_matrix(chunk, ensemble.markers)
        base_logits = features @ base_coefficients.T
        if not np.isfinite(base_logits).all():
            raise FloatingPointError("base ensemble produced non-finite logits")
        base_best_index = np.argmax(base_logits, axis=1)
        if len(classes) == 1:
            base_unique = np.ones(len(keys), dtype=bool)
        else:
            ordered = np.sort(base_logits, axis=1)
            base_unique = (ordered[:, -1] - ordered[:, -2]) > 1e-12
        agreements = np.zeros(len(keys), dtype=np.int32)
        for coefficients, valid_member in zip(bootstrap, bootstrap_valid, strict=True):
            if not valid_member:
                continue
            bootstrap_logits = features @ coefficients.T
            if not np.isfinite(bootstrap_logits).all():
                raise FloatingPointError("bootstrap ensemble produced non-finite logits")
            best = np.argmax(bootstrap_logits, axis=1)
            if len(classes) == 1:
                bootstrap_unique = np.ones(len(keys), dtype=bool)
            else:
                ordered = np.sort(bootstrap_logits, axis=1)
                bootstrap_unique = (ordered[:, -1] - ordered[:, -2]) > 1e-12
            agreements += base_unique & bootstrap_unique & (best == base_best_index)
        part = keys.copy()
        part["best_candidate"] = classes[base_best_index]
        part["stability"] = agreements / len(bootstrap)
        part["bootstrap_replicates"] = len(bootstrap)
        part["valid_bootstrap_replicates"] = ensemble.valid_bootstrap_replicates
        pieces.append(part)
    if not pieces:
        return pd.DataFrame(
            columns=[
                "donor_id",
                "object_id",
                "best_candidate",
                "stability",
                "bootstrap_replicates",
                "valid_bootstrap_replicates",
            ]
        )
    out = pd.concat(pieces, ignore_index=True)
    if out.duplicated(["donor_id", "object_id"]).any():
        raise AssertionError("bundle stability produced duplicate cells")
    return out


__all__ = [
    "BUNDLE_METHOD_VERSION",
    "BUNDLE_SCHEMA_VERSION",
    "ClassErrorResult",
    "EVALUATION_BOOTSTRAP_UNIT",
    "EVALUATION_MULTIPLICITY_METHOD",
    "EVALUATION_REPORT_METHOD_VERSION",
    "FrozenModelEnsemble",
    "MINIMUM_VALID_BOOTSTRAP_FRACTION",
    "PRODUCTION_MAX_ABS_WEIGHT",
    "PROMOTION_DECISION",
    "PROMOTION_SCHEMA_VERSION",
    "MetricInterval",
    "PromotionGatePolicy",
    "PromotionReplayReceipt",
    "SCORE_SEMANTICS",
    "SplitEvaluationReport",
    "TEMPERATURE_MAX",
    "TEMPERATURE_MIN",
    "TypingModelBundle",
    "TypingModelPromotionManifest",
    "TargetEvaluationResult",
    "ensemble_stability_wide",
    "load_typing_model_promotion_manifest",
    "load_typing_model_bundle",
    "promotion_replay_code_sha256",
    "promotion_replay_environment_sha256",
]
