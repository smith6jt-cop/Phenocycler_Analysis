"""Build auditable held-out reports from immutable candidate assignments.

The public builder in this module accepts predictions and blinded review
artifacts, never precomputed metrics.  It validates their complete lineage,
joins cells by the donor/object key, applies the review sampling weights, and
uses two-stage donor and within-stratum survey bootstrap resampling for every
confidence interval.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TypeAlias

import numpy as np
import pandas as pd

from .artifacts import (
    RunManifest,
    dataframe_content_sha256,
    hash_identifiers,
    sha256_json,
)
from .candidate_evaluation import (
    CandidateEvaluationManifest,
    load_candidate_evaluation_manifest,
)
from .hierarchical_typing import (
    AMBIGUOUS,
    AUTHORITATIVE_ANCHOR,
    INFERRED,
    OTHER,
    UNAVAILABLE,
    TypingRegistry,
)
from .marker_registry import MarkerRegistry
from .refinement_contracts import (
    REVIEW_LABEL_COLUMNS,
    ReviewLabelProvenance,
    validate_review_labels,
    validate_review_sampling_against_assignments,
)
from .typing_model_bundle import (
    ClassErrorResult,
    EvaluatedDecisionRow,
    FrozenModelEnsemble,
    MetricInterval,
    PromotionGatePolicy,
    SplitEvaluationReport,
    TargetEvaluationResult,
    TypingModelBundle,
    evaluation_bootstrap_random_state,
)

_KEY_COLUMNS = ("donor_id", "object_id")
_ACCEPTED_STATUSES = frozenset({AUTHORITATIVE_ANCHOR, INFERRED})
_KNOWN_STATUSES = frozenset(
    {AUTHORITATIVE_ANCHOR, INFERRED, AMBIGUOUS, OTHER, UNAVAILABLE, "not_applicable"}
)


@dataclass(frozen=True)
class EvaluationTargetReview:
    """One prediction-blinded sample and its independently reviewed ledger."""

    sample_frame: pd.DataFrame
    reviewer_frame: pd.DataFrame
    review_labels: pd.DataFrame
    provenance: ReviewLabelProvenance

    def __post_init__(self) -> None:
        if not isinstance(self.sample_frame, pd.DataFrame):
            raise TypeError("target review sample_frame must be a pandas DataFrame")
        if not isinstance(self.reviewer_frame, pd.DataFrame):
            raise TypeError("target review reviewer_frame must be a pandas DataFrame")
        if not isinstance(self.review_labels, pd.DataFrame):
            raise TypeError("target review review_labels must be a pandas DataFrame")
        if not isinstance(self.provenance, ReviewLabelProvenance):
            raise TypeError("target review provenance must be ReviewLabelProvenance")

    @property
    def target_key(self) -> tuple[str, str | None]:
        parent = (
            None
            if self.provenance.level in {"broad", "specific"}
            else self.provenance.parent
        )
        return self.provenance.level, parent


@dataclass(frozen=True)
class _EvaluationTargetSpec:
    level: str
    parent: str | None
    classes: tuple[str, ...]


def _specific_target_spec(model_bundle: TypingModelBundle) -> _EvaluationTargetSpec:
    return _EvaluationTargetSpec(
        level="specific",
        parent=None,
        classes=tuple(model_bundle.fit_parameters["final_specific_classes"]),
    )


CandidateEvaluationSource: TypeAlias = CandidateEvaluationManifest | Path | str


def _positive_json_int(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < 1:
        raise ValueError(f"{name} must be positive")
    return result


def _random_seed(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError("random_state must be an integer")
    result = int(value)
    if result < 0:
        raise ValueError("random_state must be nonnegative")
    return result


def _confidence_level(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError("confidence_level must be numeric")
    result = float(value)
    if not np.isfinite(result) or result <= 0 or result >= 1:
        raise ValueError("confidence_level must lie in (0, 1)")
    return result


def _explicit_bool(series: pd.Series, *, name: str) -> pd.Series:
    if series.isna().any() or not series.isin([True, False, 0, 1]).all():
        raise ValueError(f"{name} must contain explicit booleans")
    return series.astype(bool)


def _load_and_validate_candidate(
    source: CandidateEvaluationSource,
    *,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    source_run_manifest: RunManifest,
    allow_locked_test: bool,
) -> CandidateEvaluationManifest:
    if isinstance(source, CandidateEvaluationManifest):
        source.validate_current(
            model_bundle=model_bundle,
            marker_registry=marker_registry,
            typing_registry=typing_registry,
            source_run_manifest=source_run_manifest,
            allow_locked_test=allow_locked_test,
        )
        return source
    return load_candidate_evaluation_manifest(
        source,
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=allow_locked_test,
    )


def _read_exact_assignments(
    manifest: CandidateEvaluationManifest,
    *,
    columns: Sequence[str],
) -> pd.DataFrame:
    requested_columns = tuple(dict.fromkeys((*_KEY_COLUMNS, *columns)))
    frames: list[pd.DataFrame] = []
    for donor_snapshot in manifest.output.donors:
        donor_frames = [
            pd.read_parquet(
                Path(fingerprint.path),
                columns=list(requested_columns),
            )
            for fingerprint in donor_snapshot.files
        ]
        frame = (
            pd.concat(donor_frames, ignore_index=True) if len(donor_frames) > 1 else donor_frames[0]
        )
        missing = sorted(set(_KEY_COLUMNS) - set(frame.columns))
        if missing:
            raise KeyError(
                f"candidate assignments for {donor_snapshot.donor_id} "
                f"are missing key columns: {missing}"
            )
        frame = frame.copy()
        for column in _KEY_COLUMNS:
            frame[column] = frame[column].astype("string").fillna("").str.strip()
        if (frame.loc[:, list(_KEY_COLUMNS)] == "").any().any():
            raise ValueError("candidate assignments contain blank donor/object keys")
        if set(frame["donor_id"]) != {donor_snapshot.donor_id}:
            raise ValueError(
                f"candidate assignment partition {donor_snapshot.donor_id} contains another donor"
            )
        if len(frame) != donor_snapshot.row_count:
            raise ValueError(
                f"candidate assignment partition {donor_snapshot.donor_id} "
                "row count differs from its snapshot"
            )
        if hash_identifiers(frame["object_id"]) != donor_snapshot.object_id_sha256:
            raise ValueError(
                f"candidate assignment partition {donor_snapshot.donor_id} "
                "object universe differs from its snapshot"
            )
        frames.append(frame)
    assignments = pd.concat(frames, ignore_index=True)
    if len(assignments) != manifest.output.total_rows:
        raise ValueError("candidate assignment row count differs from its snapshot")
    if assignments.duplicated(list(_KEY_COLUMNS)).any():
        raise ValueError("candidate assignments contain duplicate donor/object keys")
    if set(assignments["donor_id"]) != set(manifest.donors):
        raise ValueError("candidate assignments do not exactly cover manifest donors")
    return assignments


def _validate_assignment_lineage(
    assignments: pd.DataFrame,
    *,
    manifest: CandidateEvaluationManifest,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
) -> None:
    required = {
        "candidate_assignment_run_id": manifest.assignment_run_id,
        "typing_model_bundle_content_id": model_bundle.content_id,
        "typing_mode": "probabilistic_candidate",
        "marker_registry_fingerprint": marker_registry.fingerprint,
        "typing_rules_fingerprint": typing_registry.rules_fingerprint,
    }
    missing = sorted(set(required) - set(assignments.columns))
    if missing:
        raise KeyError(f"candidate assignments are missing lineage columns: {missing}")
    for column, expected in required.items():
        values = assignments[column].astype("string").fillna("")
        if not values.eq(expected).all():
            raise ValueError(
                f"candidate assignment {column} is not uniformly bound to {expected!r}"
            )


def _target_columns(level: str) -> dict[str, str]:
    if level == "broad":
        return {
            "prediction": "broad_label",
            "status": "broad_assignment_status",
            "reason": "broad_reason",
            "eligible": "broad_threshold_eligible",
        }
    if level == "subtype":
        return {
            "prediction": "subtype_label",
            "status": "subtype_assignment_status",
            "reason": "subtype_reason",
            "eligible": "subtype_threshold_eligible",
        }
    return {
        "prediction": "specific_type",
        "status": "specific_assignment_status",
        "reason": "",
        "eligible": "",
    }


def _join_target_frame(
    *,
    review: EvaluationTargetReview,
    ensemble: FrozenModelEnsemble | _EvaluationTargetSpec,
    assignments: pd.DataFrame,
    manifest: CandidateEvaluationManifest,
    model_bundle: TypingModelBundle,
    review_context_artifact_content_ids: frozenset[str],
    review_blinding_key: str | bytes,
    allow_locked_test: bool,
) -> pd.DataFrame:
    provenance = review.provenance
    if provenance.split_name != manifest.split_name:
        raise ValueError("target review and candidate assignment use different splits")
    if (
        provenance.review_context_artifact_content_id
        not in review_context_artifact_content_ids
    ):
        raise ValueError(
            f"target {review.target_key} review context is not an output artifact "
            "in the bound source run"
        )
    validate_review_sampling_against_assignments(
        review.sample_frame,
        assignments,
        splits=model_bundle.split_manifest,
        split_name=manifest.split_name,
        candidate_assignment_run_id=manifest.assignment_run_id,
        level=ensemble.level,
        parent=ensemble.parent,
        blinding_key=review_blinding_key,
    )
    labels = validate_review_labels(
        review.review_labels,
        provenance=provenance,
        splits=model_bundle.split_manifest,
        sample_frame=review.sample_frame,
        reviewer_frame=review.reviewer_frame,
        assignment_run_id=manifest.assignment_run_id,
        typing_model_bundle_content_id=model_bundle.content_id,
        candidate_evaluation_manifest_content_id=manifest.content_id,
        candidate_assignment_output_content_id=manifest.output_content_sha256,
        allow_locked_test=allow_locked_test,
    )
    if not provenance.complete_evidence_sufficient:
        raise ValueError(f"target {review.target_key} is not completely evidence-sufficient")

    sample = review.sample_frame.copy()
    if "_review_sample_record" in sample.columns:
        raise ValueError("review sample uses reserved column '_review_sample_record'")
    for column in _KEY_COLUMNS:
        sample[column] = sample[column].astype("string").fillna("").str.strip()
    sample["sampling_weight"] = pd.to_numeric(sample["sampling_weight"], errors="raise").astype(
        float
    )
    sample["_review_sample_record"] = sample.to_dict(orient="records")
    reviewer = review.reviewer_frame.copy()
    if "_reviewer_record" in reviewer.columns:
        raise ValueError("reviewer frame uses reserved column '_reviewer_record'")
    for column in _KEY_COLUMNS:
        reviewer[column] = reviewer[column].astype("string").fillna("").str.strip()
    reviewer["_reviewer_record"] = reviewer.to_dict(orient="records")
    sample = sample.merge(
        reviewer.loc[:, [*_KEY_COLUMNS, "_reviewer_record"]],
        on=list(_KEY_COLUMNS),
        how="inner",
        validate="1:1",
    )
    label_frame = labels.loc[:, list(REVIEW_LABEL_COLUMNS)].copy()
    overlap = sorted(
        (set(sample.columns) - set(_KEY_COLUMNS) - {"_review_sample_record"})
        & (set(REVIEW_LABEL_COLUMNS) - set(_KEY_COLUMNS))
    )
    if overlap:
        raise ValueError(f"review sample overlaps review-ledger columns: {overlap}")
    joined = sample.merge(
        label_frame,
        on=list(_KEY_COLUMNS),
        how="inner",
        validate="1:1",
    )
    columns = _target_columns(ensemble.level)
    required_assignment_columns = set(_KEY_COLUMNS) | {
        value for value in columns.values() if value
    }
    if ensemble.level == "subtype":
        required_assignment_columns.update(
            {
                "broad_label",
                "broad_assignment_status",
                "typing_subtype_model_source",
            }
        )
    elif ensemble.level == "specific":
        required_assignment_columns.update(
            {
                "broad_assignment_status",
                "broad_reason",
                "broad_threshold_eligible",
                "subtype_assignment_status",
                "subtype_reason",
                "subtype_threshold_eligible",
            }
        )
    missing = sorted(required_assignment_columns - set(assignments.columns))
    if missing:
        raise KeyError(
            f"candidate assignments for target {review.target_key} are missing columns: {missing}"
        )
    selected_assignment_columns = [
        *_KEY_COLUMNS,
        *sorted(required_assignment_columns - set(_KEY_COLUMNS)),
    ]
    joined = joined.merge(
        assignments.loc[:, selected_assignment_columns],
        on=list(_KEY_COLUMNS),
        how="left",
        validate="1:1",
        indicator=True,
    )
    missing_cells = joined["_merge"].ne("both")
    if missing_cells.any():
        examples = joined.loc[missing_cells, list(_KEY_COLUMNS)].head(5)
        raise ValueError(
            "candidate assignments are missing reviewed cells: "
            f"{examples.to_dict(orient='records')}"
        )
    joined = joined.drop(columns="_merge")

    classes = set(ensemble.classes)
    # ``Other`` is a reviewed, explicitly non-modeled negative.  Retaining it
    # is necessary to measure forced named calls as false positives.
    supported_references = classes | {"Other"}
    unsupported = sorted(set(joined["reference_label"]) - supported_references)
    if unsupported:
        raise ValueError(
            f"target {review.target_key} review contains unsupported labels: {unsupported}"
        )
    prediction = joined[columns["prediction"]].astype("string").fillna("").str.strip()
    status = joined[columns["status"]].astype("string").fillna("").str.strip()
    if ensemble.level == "specific":
        broad_status = (
            joined["broad_assignment_status"].astype("string").fillna("").str.strip()
        )
        subtype_status = (
            joined["subtype_assignment_status"].astype("string").fillna("").str.strip()
        )
        use_subtype = broad_status.isin(_ACCEPTED_STATUSES) & subtype_status.ne(
            "not_applicable"
        )
        broad_reason = joined["broad_reason"].astype("string").fillna("").str.strip()
        subtype_reason = (
            joined["subtype_reason"].astype("string").fillna("").str.strip()
        )
        reason = broad_reason.where(~use_subtype, subtype_reason)
        broad_eligible = _explicit_bool(
            joined["broad_threshold_eligible"], name="broad_threshold_eligible"
        )
        subtype_eligible = _explicit_bool(
            joined["subtype_threshold_eligible"], name="subtype_threshold_eligible"
        )
        eligible = broad_eligible | (use_subtype & subtype_eligible)
    else:
        reason = joined[columns["reason"]].astype("string").fillna("").str.strip()
        eligible = _explicit_bool(joined[columns["eligible"]], name=columns["eligible"])
    unknown_status = sorted(set(status) - _KNOWN_STATUSES)
    if unknown_status:
        raise ValueError(
            f"target {review.target_key} contains unknown assignment statuses: {unknown_status}"
        )
    if reason.eq("").any():
        raise ValueError(f"target {review.target_key} contains blank assignment reasons")

    if ensemble.level in {"broad", "specific"}:
        scope = pd.Series(True, index=joined.index, dtype=bool)
    else:
        broad_status = joined["broad_assignment_status"].astype("string").fillna("").str.strip()
        broad_label = joined["broad_label"].astype("string").fillna("").str.strip()
        subtype_source = (
            joined["typing_subtype_model_source"].astype("string").fillna("").str.strip()
        )
        scope = (
            broad_status.isin(_ACCEPTED_STATUSES)
            & broad_label.eq(str(ensemble.parent))
            & subtype_source.eq("probabilistic_bundle")
        )

    accepted = scope & status.isin(_ACCEPTED_STATUSES)
    invalid_accepted = accepted & ~prediction.isin(classes)
    if invalid_accepted.any():
        values = sorted(set(prediction[invalid_accepted]))
        raise ValueError(
            f"target {review.target_key} has accepted predictions outside its "
            f"modeled classes: {values}"
        )
    target_eligible = scope & eligible
    if ensemble.level != "specific":
        bad_eligible_status = target_eligible & ~status.isin({AMBIGUOUS, INFERRED})
        if bad_eligible_status.any():
            raise ValueError(
                f"target {review.target_key} threshold eligibility conflicts with assignment status"
            )
    if ensemble.level == "specific":
        inferred = accepted & (
            broad_status.eq(INFERRED) | (use_subtype & subtype_status.eq(INFERRED))
        )
    else:
        inferred = scope & status.eq(INFERRED)
    if (inferred & ~eligible).any():
        raise ValueError(f"target {review.target_key} has inferred calls outside the rescue lane")

    joined["_prediction"] = prediction
    joined["_status"] = status
    joined["_reason"] = reason
    joined["_target_scope"] = scope.astype(bool)
    joined["_called"] = accepted.astype(bool)
    joined["_rescue_eligible"] = target_eligible.astype(bool)
    joined["_rescue_called"] = inferred.astype(bool)
    return joined


def _ratio(numerator: float, denominator: float, *, zero_denominator: float) -> float:
    if denominator <= 0:
        return zero_denominator
    value = numerator / denominator
    return float(np.clip(value, 0.0, 1.0))


def _interval(
    estimate: float,
    bootstrap_values: np.ndarray,
    *,
    confidence_level: float,
) -> MetricInterval:
    if bootstrap_values.ndim != 1 or not len(bootstrap_values):
        raise ValueError("bootstrap metric vector must be non-empty and one-dimensional")
    if not np.isfinite(bootstrap_values).all():
        raise ValueError("bootstrap metric vector contains non-finite values")
    alpha = (1.0 - confidence_level) / 2.0
    lower, upper = np.quantile(
        bootstrap_values,
        [alpha, 1.0 - alpha],
        method="linear",
    )
    # Promotion's strict interval schema requires the observed estimate to be
    # contained.  Expanding only when a biased percentile interval excludes it
    # is conservative for both lower-bound and upper-bound gates.
    lower = min(float(lower), estimate)
    upper = max(float(upper), estimate)
    return MetricInterval(estimate=estimate, lower=lower, upper=upper)


def _target_seed(random_state: int, key: tuple[str, str | None]) -> int:
    digest = sha256_json({"random_state": random_state, "level": key[0], "parent": key[1]})
    return int(digest[:16], 16)


def _build_target_result(
    *,
    frame: pd.DataFrame,
    review: EvaluationTargetReview,
    ensemble: FrozenModelEnsemble,
    confidence_level: float,
    n_bootstrap: int,
    random_state: int,
) -> tuple[TargetEvaluationResult, str]:
    weight = frame["sampling_weight"].to_numpy(dtype=float)
    reference = frame["reference_label"].astype(str).to_numpy()
    prediction = frame["_prediction"].astype(str).to_numpy()
    called = frame["_called"].to_numpy(dtype=bool)
    eligible = frame["_rescue_eligible"].to_numpy(dtype=bool)
    rescue_called = frame["_rescue_called"].to_numpy(dtype=bool)
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
    component_values = [
        weight,
        weight * called,
        weight * wrong,
        weight * eligible,
        weight * rescue_called,
        weight * rescue_wrong,
    ]
    for label in ensemble.classes:
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
    donors = tuple(sorted(frame["donor_id"].astype(str).unique()))
    donor_index = {donor: index for index, donor in enumerate(donors)}
    donor_values = np.zeros((len(donors), values.shape[1]), dtype=float)
    for row, donor in enumerate(frame["donor_id"].astype(str)):
        donor_values[donor_index[donor]] += values[row]
    point = donor_values.sum(axis=0)
    position = {name: index for index, name in enumerate(component_names)}

    if point[position["called"]] <= 0:
        raise ValueError(f"target {review.target_key} has no accepted calls")
    if point[position["rescue_eligible"]] <= 0:
        raise ValueError(f"target {review.target_key} has no rescue-eligible cells")
    if point[position["rescue_called"]] <= 0:
        raise ValueError(f"target {review.target_key} has no rescue calls")
    for label in ensemble.classes:
        if point[position[f"{label}:positive"]] <= 0:
            raise ValueError(f"target {review.target_key} class {label!r} has no positive support")
        if point[position[f"{label}:negative"]] <= 0:
            raise ValueError(f"target {review.target_key} class {label!r} has no negative support")

    rng = np.random.default_rng(_target_seed(random_state, review.target_key))
    draws = rng.integers(0, len(donors), size=(n_bootstrap, len(donors)))
    multiplicities = np.zeros((n_bootstrap, len(donors)), dtype=float)
    for replicate, sampled in enumerate(draws):
        multiplicities[replicate] = np.bincount(sampled, minlength=len(donors))
    bootstrap = multiplicities @ donor_values

    metric_specs = {
        "coverage": ("called", "population", 0.0),
        "selective_error": ("wrong", "called", 1.0),
        "rescue_coverage": ("rescue_called", "rescue_eligible", 0.0),
        "rescue_selective_error": ("rescue_wrong", "rescue_called", 1.0),
    }
    intervals: dict[str, MetricInterval] = {}
    for name, (numerator, denominator, fallback) in metric_specs.items():
        numerator_index = position[numerator]
        denominator_index = position[denominator]
        estimate = _ratio(
            point[numerator_index],
            point[denominator_index],
            zero_denominator=fallback,
        )
        bootstrap_metric = np.asarray(
            [
                _ratio(numerator_value, denominator_value, zero_denominator=fallback)
                for numerator_value, denominator_value in zip(
                    bootstrap[:, numerator_index],
                    bootstrap[:, denominator_index],
                    strict=True,
                )
            ],
            dtype=float,
        )
        intervals[name] = _interval(
            estimate,
            bootstrap_metric,
            confidence_level=confidence_level,
        )

    class_errors: list[ClassErrorResult] = []
    for label in ensemble.classes:
        positive_index = position[f"{label}:positive"]
        negative_index = position[f"{label}:negative"]
        fp_index = position[f"{label}:fp"]
        fn_index = position[f"{label}:fn"]
        fp_estimate = _ratio(point[fp_index], point[negative_index], zero_denominator=1.0)
        fn_estimate = _ratio(point[fn_index], point[positive_index], zero_denominator=1.0)
        fp_bootstrap = np.asarray(
            [
                _ratio(numerator, denominator, zero_denominator=1.0)
                for numerator, denominator in zip(
                    bootstrap[:, fp_index],
                    bootstrap[:, negative_index],
                    strict=True,
                )
            ]
        )
        fn_bootstrap = np.asarray(
            [
                _ratio(numerator, denominator, zero_denominator=1.0)
                for numerator, denominator in zip(
                    bootstrap[:, fn_index],
                    bootstrap[:, positive_index],
                    strict=True,
                )
            ]
        )
        class_errors.append(
            ClassErrorResult(
                label=label,
                reference_positive_weight=float(point[positive_index]),
                reference_negative_weight=float(point[negative_index]),
                false_positive_weight=float(point[fp_index]),
                false_negative_weight=float(point[fn_index]),
                false_positive_rate=_interval(
                    fp_estimate,
                    fp_bootstrap,
                    confidence_level=confidence_level,
                ),
                false_negative_rate=_interval(
                    fn_estimate,
                    fn_bootstrap,
                    confidence_level=confidence_level,
                ),
            )
        )

    evaluated = frame.loc[
        :,
        [
            *_KEY_COLUMNS,
            "sampling_weight",
            "reference_label",
            "_prediction",
            "_status",
            "_reason",
            "_target_scope",
            "_called",
            "_rescue_eligible",
            "_rescue_called",
        ],
    ].copy()
    frame_fingerprint = dataframe_content_sha256(evaluated, sort_by=list(_KEY_COLUMNS))
    return (
        TargetEvaluationResult(
            level=ensemble.level,
            parent=ensemble.parent,
            classes=ensemble.classes,
            row_count=len(frame),
            donors=donors,
            population_weight=float(point[position["population"]]),
            called_weight=float(point[position["called"]]),
            wrong_call_weight=float(point[position["wrong"]]),
            rescue_eligible_weight=float(point[position["rescue_eligible"]]),
            rescue_called_weight=float(point[position["rescue_called"]]),
            rescue_wrong_call_weight=float(point[position["rescue_wrong"]]),
            coverage=intervals["coverage"],
            selective_error=intervals["selective_error"],
            rescue_coverage=intervals["rescue_coverage"],
            rescue_selective_error=intervals["rescue_selective_error"],
            class_errors=tuple(class_errors),
            review_provenance=review.provenance,
        ),
        frame_fingerprint,
    )


def _build_replayable_target_result(
    *,
    frame: pd.DataFrame,
    review: EvaluationTargetReview,
    ensemble: FrozenModelEnsemble | _EvaluationTargetSpec,
    confidence_level: float,
    n_bootstrap: int,
    random_state: int,
    multiplicity_target_count: int,
) -> TargetEvaluationResult:
    """Normalize scientific rows and delegate every metric to replay code."""

    rows = tuple(
        EvaluatedDecisionRow(
            donor_id=row["donor_id"],
            object_id=row["object_id"],
            sample_record=row["_review_sample_record"],
            reviewer_record=row["_reviewer_record"],
            reference_label=row["reference_label"],
            reviewer_1=row["reviewer_1"],
            reviewer_2=row["reviewer_2"],
            adjudicated=bool(row["adjudicated"]),
            root_cause=row["root_cause"],
            evidence_sufficient=bool(row["evidence_sufficient"]),
            prediction=row["_prediction"],
            assignment_status=row["_status"],
            assignment_reason=row["_reason"],
            target_scope=bool(row["_target_scope"]),
            called=bool(row["_called"]),
            rescue_eligible=bool(row["_rescue_eligible"]),
            rescue_called=bool(row["_rescue_called"]),
        )
        for row in frame.to_dict(orient="records")
    )
    return TargetEvaluationResult.from_evaluated_rows(
        level=ensemble.level,
        parent=ensemble.parent,
        classes=ensemble.classes,
        evaluated_rows=rows,
        confidence_level=confidence_level,
        n_bootstrap=n_bootstrap,
        random_state=random_state,
        multiplicity_target_count=multiplicity_target_count,
        review_provenance=review.provenance,
    )


def build_split_evaluation_report(
    *,
    candidate_evaluation: CandidateEvaluationSource,
    model_bundle: TypingModelBundle,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    source_run_manifest: RunManifest,
    target_reviews: Sequence[EvaluationTargetReview],
    promotion_policy: PromotionGatePolicy,
    review_blinding_key: str | bytes,
    random_state: int | None = None,
    allow_locked_test: bool = False,
) -> SplitEvaluationReport:
    """Derive one strict validation or explicitly unlocked test report.

    ``target_reviews`` must contain exactly one broad review, one all-cell
    end-to-end specific review, and one diagnostic review for every
    probabilistic subtype ensemble in ``model_bundle``.  All metrics,
    support weights, pass-relevant confidence bounds, and artifact identities
    are computed here; callers cannot attest any result values.
    """

    if not isinstance(promotion_policy, PromotionGatePolicy):
        raise TypeError("promotion_policy must be PromotionGatePolicy")
    if (
        model_bundle.split_manifest.promotion_policy_content_id
        != promotion_policy.content_id
    ):
        raise ValueError(
            "held-out evaluation policy was not preregistered in the donor split"
        )
    from .typing_model_bundle import _promotion_report_confidence_level

    confidence = _promotion_report_confidence_level(
        promotion_policy.confidence_level
    )
    bootstrap_count = promotion_policy.bootstrap_replicates
    if model_bundle.threshold_selection_provenance is None:
        raise ValueError(
            "held-out evaluation requires a calibration-frozen model bundle; "
            "preliminary fits are restricted to calibration threshold selection"
        )
    manifest = _load_and_validate_candidate(
        candidate_evaluation,
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=allow_locked_test,
    )
    if manifest.split_name not in {"validation", "locked_test"}:
        raise ValueError(
            "held-out evaluation reports require validation or locked_test; "
            "calibration assignments are reserved for threshold selection"
        )
    seed = evaluation_bootstrap_random_state(
        split_name=manifest.split_name,
        model_bundle_content_id=model_bundle.content_id,
        assignment_run_id=manifest.assignment_run_id,
    )
    if random_state is not None and _random_seed(random_state) != seed:
        raise ValueError(
            "random_state is fixed by immutable evaluation lineage; omit it or "
            f"supply {seed}"
        )
    assignment_columns = {
        "candidate_assignment_run_id",
        "typing_model_bundle_content_id",
        "typing_mode",
        "marker_registry_fingerprint",
        "typing_rules_fingerprint",
        *_target_columns("broad").values(),
        "specific_type",
        "specific_assignment_status",
        "subtype_assignment_status",
        "subtype_reason",
        "subtype_threshold_eligible",
    }
    if model_bundle.subtypes:
        assignment_columns.update(_target_columns("subtype").values())
        assignment_columns.update(
            {
                "broad_label",
                "broad_assignment_status",
                "typing_subtype_model_source",
            }
        )
    assignments = _read_exact_assignments(
        manifest,
        columns=tuple(sorted(assignment_columns)),
    )
    _validate_assignment_lineage(
        assignments,
        manifest=manifest,
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )

    if isinstance(target_reviews, (str, bytes)) or not isinstance(target_reviews, Sequence):
        raise TypeError("target_reviews must be a sequence")
    reviews = tuple(target_reviews)
    if any(not isinstance(value, EvaluationTargetReview) for value in reviews):
        raise TypeError("target_reviews must contain EvaluationTargetReview objects")
    observed_keys = [review.target_key for review in reviews]
    if len(set(observed_keys)) != len(observed_keys):
        raise ValueError("target_reviews contains duplicate model targets")
    observed = dict(zip(observed_keys, reviews, strict=True))
    ensembles = (
        model_bundle.broad,
        _specific_target_spec(model_bundle),
        *model_bundle.subtypes,
    )
    multiplicity_target_count = len(ensembles)
    expected = {(ensemble.level, ensemble.parent): ensemble for ensemble in ensembles}
    missing = sorted(set(expected) - set(observed), key=str)
    extra = sorted(set(observed) - set(expected), key=str)
    if missing or extra:
        raise ValueError(f"target review coverage is incomplete: missing={missing}, extra={extra}")

    target_results: list[TargetEvaluationResult] = []
    ingest_references = tuple(
        reference
        for reference in source_run_manifest.stages
        if reference.stage == "ingest"
    )
    if len(ingest_references) != 1:
        raise ValueError(
            "held-out review context requires exactly one source-run ingest output"
        )
    review_context_ids = frozenset(
        {ingest_references[0].output_content_sha256}
    )
    for ensemble in ensembles:
        key = ensemble.level, ensemble.parent
        review = observed[key]
        frame = _join_target_frame(
            review=review,
            ensemble=ensemble,
            assignments=assignments,
            manifest=manifest,
            model_bundle=model_bundle,
            review_context_artifact_content_ids=review_context_ids,
            review_blinding_key=review_blinding_key,
            allow_locked_test=allow_locked_test,
        )
        result = _build_replayable_target_result(
            frame=frame,
            review=review,
            ensemble=ensemble,
            confidence_level=confidence,
            n_bootstrap=bootstrap_count,
            random_state=seed,
            multiplicity_target_count=multiplicity_target_count,
        )
        target_results.append(result)

    # Recheck immutable bytes after reading to narrow the validation/read race.
    manifest.validate_current(
        model_bundle=model_bundle,
        marker_registry=marker_registry,
        typing_registry=typing_registry,
        source_run_manifest=source_run_manifest,
        allow_locked_test=allow_locked_test,
    )
    return SplitEvaluationReport.from_target_results(
        split_name=manifest.split_name,
        model_bundle_content_id=model_bundle.content_id,
        split_manifest_content_id=model_bundle.split_manifest.content_id,
        source_run_id=model_bundle.split_manifest.source_run_id,
        source_run_manifest=source_run_manifest.to_dict(),
        candidate_evaluation_manifest_content_id=manifest.content_id,
        candidate_evaluation_manifest=manifest.to_dict(),
        assignment_run_id=manifest.assignment_run_id,
        assignment_artifact_content_id=manifest.output_content_sha256,
        evaluation_donors=manifest.donors,
        promotion_policy=promotion_policy,
        targets=tuple(target_results),
    )


__all__ = [
    "CandidateEvaluationSource",
    "EvaluationTargetReview",
    "build_split_evaluation_report",
]
