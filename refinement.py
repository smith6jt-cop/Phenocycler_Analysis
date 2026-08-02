"""Detached review sampling and error-aware assignment diagnostics.

The helpers in this module support *audit* and *refinement* notebooks.  They do
not change a cell's authoritative assignment, infer biology from target class
frequencies, or turn a review-priority flag into a diagnosis.  In particular,
``candidate_fp`` and ``candidate_fn`` mean "worth reviewing for this failure
mode", not "known to be a false positive/negative".

Only NumPy and pandas are required so the same calculations can be reused by
headless validation jobs and interactive notebooks.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping, Sequence
from numbers import Integral
from typing import Any

import numpy as np
import pandas as pd

DEFAULT_SAMPLE_SALT = "phenocycler-refinement-v1"
DEFAULT_ABSTAIN_LABELS = ("Ambiguous", "Unavailable")
DEFAULT_NON_CALL_LABELS = (*DEFAULT_ABSTAIN_LABELS, "Other")
# Keep this detached helper dependency-light while requiring the exact score
# contract stamped by a promoted TypingModelBundle.  A generic claim such as
# ``posterior_probability`` is not sufficient provenance for calibration plots.
POSTERIOR_PROBABILITY_SEMANTICS = (
    "donor_class_balanced_temperature_scaled_softmax_probability"
)

_SAMPLE_AUDIT_COLUMNS = (
    "sampling_stratum",
    "stratum_population",
    "stratum_sampled",
    "sampling_weight",
)
_MISSING_PREDICTION = object()


def _canonical_scalar(value: Any) -> list[str]:
    """Return a type-aware, JSON-safe representation of one scalar.

    Type information prevents common identifiers such as integer ``1`` and
    string ``"1"`` from sharing a hash input.  Values used as identifiers are
    expected to have stable string representations.
    """

    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, pd.Timestamp):
        return ["pandas.Timestamp", value.isoformat()]
    if isinstance(value, bytes):
        return ["builtins.bytes", value.hex()]
    type_name = f"{type(value).__module__}.{type(value).__qualname__}"
    if isinstance(value, float):
        return [type_name, value.hex()]
    return [type_name, str(value)]


def _canonical_token(values: Sequence[Any]) -> str:
    return json.dumps(
        [_canonical_scalar(value) for value in values],
        ensure_ascii=True,
        separators=(",", ":"),
    )


def _stable_digest(values: Sequence[Any], *, salt: str) -> str:
    payload = _canonical_token(values).encode("utf-8")
    return hashlib.sha256(salt.encode("utf-8") + b"\0" + payload).hexdigest()


def _validate_positive_integer(value: Any, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or int(value) <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def deterministic_stratified_sample(
    frame: pd.DataFrame,
    strata_col: str,
    n_per_stratum: int | Mapping[Any, int],
    id_cols: Sequence[str] = ("donor_id", "object_id"),
    salt: str = DEFAULT_SAMPLE_SALT,
) -> pd.DataFrame:
    """Select a reproducible audit sample without depending on row order.

    Rows are ranked inside each stratum by a SHA-256 digest of the stratum and
    immutable identifier columns.  The returned table has a deterministic row
    order and includes inverse-probability weights.  If a stratum is smaller
    than its requested sample size, all rows are returned with weight one.

    ``n_per_stratum`` may be one positive integer for every stratum or a mapping
    from each observed stratum value to its requested positive integer.  Keys
    must be unique and identifiers/strata must be non-missing and non-blank.
    These constraints prevent an audit sample from silently changing when the
    input is reordered or re-indexed.
    """

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    if not isinstance(strata_col, str) or not strata_col:
        raise ValueError("strata_col must be a non-empty column name")
    id_cols = tuple(id_cols)
    if not id_cols or any(not isinstance(column, str) or not column for column in id_cols):
        raise ValueError("id_cols must contain at least one non-empty column name")
    if len(set(id_cols)) != len(id_cols):
        raise ValueError("id_cols contains duplicate column names")
    required = [strata_col, *id_cols]
    missing_columns = [column for column in required if column not in frame.columns]
    if missing_columns:
        raise KeyError(f"frame is missing required columns {missing_columns}")
    if frame[required].isna().any(axis=None):
        missing = [column for column in required if frame[column].isna().any()]
        raise ValueError(f"sampling keys/strata contain missing values in {missing}")
    blank = frame[required].apply(
        lambda column: column.map(
            lambda value: isinstance(value, str) and not value.strip()
        )
    )
    if blank.any(axis=None):
        columns = [column for column in required if blank[column].any()]
        raise ValueError(f"sampling keys/strata contain blank values in {columns}")
    if frame.duplicated(list(id_cols)).any():
        examples = frame.loc[frame.duplicated(list(id_cols), keep=False), list(id_cols)].head(3)
        raise ValueError(
            "sampling identifier columns must be unique; duplicate examples: "
            f"{examples.to_dict(orient='records')}"
        )
    if not isinstance(salt, str) or not salt:
        raise ValueError("salt must be a non-empty string")

    requested_by_stratum: Mapping[Any, int] | None
    if isinstance(n_per_stratum, Mapping):
        requested_by_stratum = n_per_stratum
    else:
        common_request = _validate_positive_integer(
            n_per_stratum, name="n_per_stratum"
        )
        requested_by_stratum = None

    if frame.empty:
        out = frame.copy()
        out["sampling_stratum"] = pd.Series(dtype="object")
        out["stratum_population"] = pd.Series(dtype="int64")
        out["stratum_sampled"] = pd.Series(dtype="int64")
        out["sampling_weight"] = pd.Series(dtype="float64")
        return out

    work = frame.copy()
    hidden = {
        "stratum": "__refinement_stratum_token__",
        "key": "__refinement_key_token__",
        "rank": "__refinement_rank_token__",
    }
    collisions = [name for name in hidden.values() if name in work.columns]
    if collisions:
        raise ValueError(f"frame uses reserved internal columns {collisions}")

    work[hidden["stratum"]] = [
        _canonical_token((value,)) for value in work[strata_col]
    ]
    work[hidden["key"]] = [
        _canonical_token(tuple(row))
        for row in work.loc[:, list(id_cols)].itertuples(index=False, name=None)
    ]
    work[hidden["rank"]] = [
        _stable_digest((stratum, *keys), salt=salt)
        for stratum, *keys in work.loc[:, [strata_col, *id_cols]].itertuples(
            index=False, name=None
        )
    ]

    sampled_groups: list[pd.DataFrame] = []
    for _, group in work.groupby(hidden["stratum"], sort=True, observed=True):
        stratum = group[strata_col].iloc[0]
        if requested_by_stratum is None:
            requested = common_request
        else:
            try:
                requested_value = requested_by_stratum[stratum]
            except (KeyError, TypeError) as exc:
                raise KeyError(
                    f"n_per_stratum has no sample size for stratum {stratum!r}"
                ) from exc
            requested = _validate_positive_integer(
                requested_value, name=f"n_per_stratum[{stratum!r}]"
            )
        population = len(group)
        sampled = min(requested, population)
        chosen = (
            group.sort_values(
                [hidden["rank"], hidden["key"]],
                kind="mergesort",
            )
            .head(sampled)
            .copy()
        )
        chosen["sampling_stratum"] = stratum
        chosen["stratum_population"] = population
        chosen["stratum_sampled"] = sampled
        chosen["sampling_weight"] = population / sampled
        sampled_groups.append(chosen)

    out = pd.concat(sampled_groups, ignore_index=True)
    out = out.drop(columns=list(hidden.values()))
    # Put audit fields at the end even if an input table already carried stale
    # versions of them; the values above are authoritative for this draw.
    base_columns = [column for column in frame.columns if column not in _SAMPLE_AUDIT_COLUMNS]
    return out.loc[:, [*base_columns, *_SAMPLE_AUDIT_COLUMNS]]


def _as_series(values: Sequence[Any], *, name: str) -> pd.Series:
    if isinstance(values, pd.DataFrame):
        raise TypeError(f"{name} must be one-dimensional")
    if isinstance(values, (str, bytes)) or np.isscalar(values):
        raise TypeError(f"{name} must be a one-dimensional sequence")
    array = np.asarray(values, dtype=object)
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return pd.Series(array, dtype="object")


def _validated_weights(weights: Sequence[float] | None, *, length: int) -> np.ndarray:
    if weights is None:
        return np.ones(length, dtype=float)
    if isinstance(weights, (str, bytes)) or np.isscalar(weights):
        raise TypeError("weights must be a one-dimensional sequence")
    try:
        result = np.asarray(weights, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("weights must be numeric") from exc
    if result.ndim != 1 or len(result) != length:
        raise ValueError(f"weights must be one-dimensional with length {length}")
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError("weights must be finite and non-negative")
    if result.sum() <= 0:
        raise ValueError("weights must have positive total weight")
    return result


def _validated_abstentions(abstain_labels: Sequence[Any]) -> tuple[Any, ...]:
    if isinstance(abstain_labels, (str, bytes)):
        labels = (abstain_labels,)
    else:
        labels = tuple(abstain_labels)
    if any(pd.isna(label) for label in labels):
        raise ValueError("abstain_labels cannot contain missing values")
    return labels


def _validated_threshold_lane(
    candidate_prediction: Sequence[Any] | None,
    threshold_eligible: Sequence[bool] | None,
    *,
    length: int,
    abstentions: Sequence[Any],
) -> tuple[pd.Series | None, np.ndarray | None]:
    """Validate an optional lane whose candidates are rethresholded."""

    if (candidate_prediction is None) != (threshold_eligible is None):
        raise ValueError(
            "candidate_prediction and threshold_eligible must be provided together"
        )
    if candidate_prediction is None:
        return None, None
    candidate = _as_series(candidate_prediction, name="candidate_prediction")
    eligible_series = _as_series(threshold_eligible, name="threshold_eligible")
    if len(candidate) != length or len(eligible_series) != length:
        raise ValueError(
            "candidate_prediction and threshold_eligible must match reference length"
        )
    if eligible_series.isna().any() or not eligible_series.isin(
        [True, False, 0, 1]
    ).all():
        raise ValueError("threshold_eligible must contain explicit booleans")
    eligible = eligible_series.astype(bool).to_numpy()
    missing_candidate = candidate.isna().to_numpy()
    abstaining_candidate = candidate.map(
        lambda value: False if pd.isna(value) else _is_in(value, abstentions)
    ).to_numpy(dtype=bool)
    if (eligible & (missing_candidate | abstaining_candidate)).any():
        raise ValueError(
            "every threshold-eligible row requires a non-abstaining candidate"
        )
    return candidate, eligible


def _is_in(value: Any, candidates: Sequence[Any]) -> bool:
    return any(value == candidate for candidate in candidates)


def _metric_labels(
    reference: pd.Series,
    prediction: pd.Series,
    labels: Sequence[Any] | None,
    abstentions: tuple[Any, ...],
) -> list[Any]:
    if labels is None:
        candidates = pd.unique(pd.concat([reference, prediction], ignore_index=True))
    elif isinstance(labels, (str, bytes)):
        candidates = np.asarray([labels], dtype=object)
    else:
        candidates = np.asarray(list(labels), dtype=object)

    result: list[Any] = []
    for label in candidates:
        if label is _MISSING_PREDICTION or pd.isna(label) or _is_in(label, abstentions):
            continue
        if not any(label == existing for existing in result):
            result.append(label)
    return sorted(result, key=lambda value: _canonical_token((value,)))


def _safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator / denominator) if denominator > 0 else np.nan


def _require_probability_semantics(score_semantics: str) -> None:
    """Refuse calibration diagnostics for scores that are not probabilities."""

    if score_semantics != POSTERIOR_PROBABILITY_SEMANTICS:
        raise ValueError(
            "probability calibration requires "
            f"score_semantics={POSTERIOR_PROBABILITY_SEMANTICS!r}; "
            f"received {score_semantics!r}"
        )


def _validated_probability_frame(
    probabilities: pd.DataFrame | Mapping[Any, Sequence[float]],
    *,
    length: int,
) -> pd.DataFrame:
    if isinstance(probabilities, pd.DataFrame):
        frame = probabilities.copy()
    elif isinstance(probabilities, Mapping):
        frame = pd.DataFrame(dict(probabilities))
    else:
        raise TypeError("probabilities must be a DataFrame or mapping of class scores")
    if len(frame) != length:
        raise ValueError(f"probabilities must have length {length}")
    if frame.columns.has_duplicates:
        raise ValueError("probabilities contains duplicate class columns")
    return frame.reset_index(drop=True)


def _calibration_labels(
    truth: pd.Series,
    probability_columns: Sequence[Any],
    labels: Sequence[Any] | None,
    abstentions: tuple[Any, ...],
) -> list[Any]:
    available = list(probability_columns)
    if labels is None:
        candidates = available
    elif isinstance(labels, (str, bytes)):
        candidates = [labels]
    else:
        candidates = list(labels)
    if not candidates:
        raise ValueError("at least one probability class is required")

    result: list[Any] = []
    for label in candidates:
        if pd.isna(label):
            raise ValueError("probability class labels cannot be missing")
        if _is_in(label, abstentions):
            raise ValueError("abstention labels cannot be probability classes")
        if not any(label == column for column in available):
            raise KeyError(f"probabilities is missing class column {label!r}")
        if not any(label == existing for existing in result):
            result.append(label)

    # Resolved references outside the model vocabulary (most importantly
    # ``Other``) remain valid one-vs-rest negatives for every modeled class.
    # Dropping them would hide forced-class false-positive probability mass and
    # make Brier/ECE estimates optimistically biased.
    return sorted(result, key=lambda value: _canonical_token((value,)))


def one_vs_rest_reliability_table(
    reference: Sequence[Any],
    probabilities: pd.DataFrame | Mapping[Any, Sequence[float]],
    labels: Sequence[Any] | None = None,
    weights: Sequence[float] | None = None,
    abstain_labels: Sequence[Any] = DEFAULT_ABSTAIN_LABELS,
    *,
    score_semantics: str,
    n_bins: int = 10,
) -> pd.DataFrame:
    """Return weighted one-vs-rest reliability bins and Brier diagnostics.

    Calibration is meaningful only for a declared posterior probability.  The
    explicit ``score_semantics`` argument prevents empirical tail evidence,
    margins, logits, or arbitrary ranking scores in ``[0, 1]`` from being
    silently treated as calibrated probabilities.

    Equal-width bins include zero in the first bin and one in the last.  Empty
    bins are retained with zero weight and undefined means, which gives plots a
    stable class-by-bin shape.  ``brier_score`` and
    ``expected_calibration_error`` are class-level weighted summaries repeated
    on each bin row for a self-contained audit table. Resolved reference labels
    such as ``Other`` need not have their own probability column: they are
    retained as negatives for every modeled class.
    """

    _require_probability_semantics(score_semantics)
    truth = _as_series(reference, name="reference")
    if truth.isna().any():
        raise ValueError("reference contains missing values")
    abstentions = _validated_abstentions(abstain_labels)
    if truth.map(lambda value: _is_in(value, abstentions)).any():
        raise ValueError(
            "reference contains abstention labels; adjudicate or exclude those rows first"
        )
    probability_frame = _validated_probability_frame(
        probabilities,
        length=len(truth),
    )
    try:
        probability_frame = probability_frame.apply(
            pd.to_numeric, errors="raise"
        ).astype(float)
    except (TypeError, ValueError) as exc:
        raise ValueError("all probability columns must be numeric") from exc
    probability_values = probability_frame.to_numpy(dtype=float)
    if (
        not np.isfinite(probability_values).all()
        or ((probability_values < 0) | (probability_values > 1)).any()
    ):
        raise ValueError("probabilities must be finite and lie in [0, 1]")
    if not np.allclose(probability_values.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError(
            "temperature-scaled softmax probability rows must sum to one"
        )
    class_labels = _calibration_labels(
        truth,
        probability_frame.columns,
        labels,
        abstentions,
    )
    sample_weight = _validated_weights(weights, length=len(truth))
    bin_count = _validate_positive_integer(n_bins, name="n_bins")
    total_weight = float(sample_weight.sum())

    rows: list[dict[str, Any]] = []
    for label in class_labels:
        score = probability_frame[label].to_numpy(dtype=float)
        outcome = truth.eq(label).to_numpy(dtype=float)
        brier = float(np.sum(sample_weight * np.square(score - outcome)) / total_weight)
        indices = np.minimum(
            np.floor(score * bin_count).astype(int),
            bin_count - 1,
        )
        label_rows: list[dict[str, Any]] = []
        ece = 0.0
        for bin_index in range(bin_count):
            selected = indices == bin_index
            selected_weight = float(sample_weight[selected].sum())
            selected_count = int(selected.sum())
            positive_weight = float(np.sum(sample_weight[selected] * outcome[selected]))
            if selected_weight > 0:
                mean_probability = float(
                    np.sum(sample_weight[selected] * score[selected]) / selected_weight
                )
                observed_fraction = positive_weight / selected_weight
                gap = observed_fraction - mean_probability
                ece += selected_weight / total_weight * abs(gap)
            else:
                mean_probability = np.nan
                observed_fraction = np.nan
                gap = np.nan
            label_rows.append(
                {
                    "label": label,
                    "bin_index": bin_index,
                    "bin_low": bin_index / bin_count,
                    "bin_high": (bin_index + 1) / bin_count,
                    "bin_count": selected_count,
                    "bin_weight": selected_weight,
                    "bin_positive_weight": positive_weight,
                    "mean_probability": mean_probability,
                    "observed_positive_fraction": observed_fraction,
                    "calibration_gap": gap,
                    "absolute_calibration_gap": abs(gap) if np.isfinite(gap) else np.nan,
                    "total_weight": total_weight,
                    "class_support_weight": float(
                        np.sum(sample_weight * outcome)
                    ),
                    "brier_score": brier,
                }
            )
        for row in label_rows:
            row["expected_calibration_error"] = float(ece)
        rows.extend(label_rows)

    columns = [
        "label",
        "bin_index",
        "bin_low",
        "bin_high",
        "bin_count",
        "bin_weight",
        "bin_positive_weight",
        "mean_probability",
        "observed_positive_fraction",
        "calibration_gap",
        "absolute_calibration_gap",
        "total_weight",
        "class_support_weight",
        "brier_score",
        "expected_calibration_error",
    ]
    return pd.DataFrame(rows, columns=columns)


def per_class_metrics(
    reference: Sequence[Any],
    prediction: Sequence[Any],
    labels: Sequence[Any] | None = None,
    weights: Sequence[float] | None = None,
    abstain_labels: Sequence[Any] = DEFAULT_ABSTAIN_LABELS,
) -> pd.DataFrame:
    """Compute weighted one-vs-rest error metrics for each biological class.

    An abstaining prediction is a false negative for its known reference class
    and is never promoted to an ``Ambiguous``/``Unavailable`` pseudo-class.
    Missing predictions are treated as abstentions.  Reference labels must be
    resolved, because an unresolved review cannot establish either a false
    positive or a false negative.

    ``TP``, ``FP``, ``FN``, ``TN``, ``support``, and ``calls`` are weighted
    totals (ordinary counts when ``weights`` is omitted).  Undefined ratios are
    returned as ``NaN`` rather than silently reported as zero.
    """

    truth = _as_series(reference, name="reference")
    called = _as_series(prediction, name="prediction")
    if len(truth) != len(called):
        raise ValueError("reference and prediction must have the same length")
    if truth.isna().any():
        raise ValueError("reference contains missing values")
    abstentions = _validated_abstentions(abstain_labels)
    unresolved_reference = truth.map(lambda value: _is_in(value, abstentions))
    if unresolved_reference.any():
        raise ValueError(
            "reference contains abstention labels; adjudicate or exclude those rows first"
        )
    sample_weight = _validated_weights(weights, length=len(truth))
    called.loc[called.isna()] = _MISSING_PREDICTION
    class_labels = _metric_labels(truth, called, labels, abstentions)

    rows: list[dict[str, Any]] = []
    for label in class_labels:
        truth_positive = truth.eq(label).to_numpy(dtype=bool)
        called_positive = called.eq(label).to_numpy(dtype=bool)
        tp = float(sample_weight[truth_positive & called_positive].sum())
        fp = float(sample_weight[~truth_positive & called_positive].sum())
        fn = float(sample_weight[truth_positive & ~called_positive].sum())
        tn = float(sample_weight[~truth_positive & ~called_positive].sum())
        precision = _safe_ratio(tp, tp + fp)
        recall = _safe_ratio(tp, tp + fn)
        specificity = _safe_ratio(tn, tn + fp)
        negative_predictive_value = _safe_ratio(tn, tn + fn)
        false_positive_rate = _safe_ratio(fp, fp + tn)
        false_negative_rate = _safe_ratio(fn, fn + tp)
        f1 = _safe_ratio(2 * tp, 2 * tp + fp + fn)
        rows.append(
            {
                "label": label,
                "TP": tp,
                "FP": fp,
                "FN": fn,
                "TN": tn,
                "precision": precision,
                "recall": recall,
                "specificity": specificity,
                "NPV": negative_predictive_value,
                "FPR": false_positive_rate,
                "FNR": false_negative_rate,
                "F1": f1,
                "support": tp + fn,
                "calls": tp + fp,
            }
        )
    columns = [
        "label",
        "TP",
        "FP",
        "FN",
        "TN",
        "precision",
        "recall",
        "specificity",
        "NPV",
        "FPR",
        "FNR",
        "F1",
        "support",
        "calls",
    ]
    return pd.DataFrame(rows, columns=columns)


def grouped_bootstrap_metrics(
    reference: Sequence[Any],
    prediction: Sequence[Any],
    groups: Sequence[Any],
    labels: Sequence[Any] | None = None,
    weights: Sequence[float] | None = None,
    abstain_labels: Sequence[Any] = DEFAULT_ABSTAIN_LABELS,
    *,
    replicates: int = 1_000,
    confidence: float = 0.95,
    seed: int = 0,
) -> pd.DataFrame:
    """Return donor/ROI-clustered bootstrap intervals for error metrics.

    Whole groups are sampled with replacement, preserving within-group cell
    dependence and any inverse-probability review weights.  The point estimate
    is calculated on the original sample.  Undefined replicate metrics are
    omitted from that metric's interval and reported through
    ``valid_replicates``; they are never silently replaced by zero.
    """

    truth = _as_series(reference, name="reference")
    called = _as_series(prediction, name="prediction")
    cluster = _as_series(groups, name="groups")
    if len(truth) != len(called) or len(truth) != len(cluster):
        raise ValueError("reference, prediction, and groups must have the same length")
    if cluster.isna().any():
        raise ValueError("groups contains missing values")
    n_replicates = _validate_positive_integer(replicates, name="replicates")
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between 0 and 1")
    if isinstance(seed, bool) or not isinstance(seed, Integral):
        raise ValueError("seed must be an integer")
    sample_weight = _validated_weights(weights, length=len(truth))
    cluster_tokens = cluster.map(lambda value: _canonical_token((value,)))
    group_codes, unique_groups = pd.factorize(cluster_tokens, sort=True)
    if len(unique_groups) < 2:
        raise ValueError("grouped bootstrap requires at least two groups")

    point = per_class_metrics(
        truth,
        called,
        labels=labels,
        weights=sample_weight,
        abstain_labels=abstain_labels,
    )
    metric_names = ("precision", "recall", "specificity", "NPV", "FPR", "FNR", "F1")
    class_labels = point["label"].tolist()
    values = {
        (label, metric): np.full(n_replicates, np.nan, dtype=float)
        for label in class_labels
        for metric in metric_names
    }
    rng = np.random.default_rng(int(seed))
    n_groups = len(unique_groups)
    for replicate in range(n_replicates):
        sampled_groups = rng.integers(0, n_groups, size=n_groups)
        multiplicity = np.bincount(sampled_groups, minlength=n_groups)
        replicate_weights = sample_weight * multiplicity[group_codes]
        replicate_metrics = per_class_metrics(
            truth,
            called,
            labels=class_labels,
            weights=replicate_weights,
            abstain_labels=abstain_labels,
        ).set_index("label")
        for label in class_labels:
            for metric in metric_names:
                values[(label, metric)][replicate] = replicate_metrics.loc[label, metric]

    alpha = (1.0 - float(confidence)) / 2.0
    point_index = point.set_index("label")
    rows = []
    for label in class_labels:
        for metric in metric_names:
            finite = values[(label, metric)]
            finite = finite[np.isfinite(finite)]
            low, high = (
                np.quantile(finite, [alpha, 1.0 - alpha])
                if len(finite)
                else (np.nan, np.nan)
            )
            rows.append(
                {
                    "label": label,
                    "metric": metric,
                    "estimate": float(point_index.loc[label, metric]),
                    "ci_low": float(low),
                    "ci_high": float(high),
                    "valid_replicates": int(len(finite)),
                    "replicates": n_replicates,
                    "groups": n_groups,
                    "confidence_level": float(confidence),
                }
            )
    return pd.DataFrame(rows)


def coverage_risk_curve(
    reference: Sequence[Any],
    prediction: Sequence[Any],
    confidence: Sequence[float],
    thresholds: Sequence[float] | None = None,
    weights: Sequence[float] | None = None,
    abstain_labels: Sequence[Any] = DEFAULT_ABSTAIN_LABELS,
    non_call_labels: Sequence[Any] = DEFAULT_NON_CALL_LABELS,
    *,
    candidate_prediction: Sequence[Any] | None = None,
    threshold_eligible: Sequence[bool] | None = None,
) -> pd.DataFrame:
    """Compute a weighted coverage-versus-error curve for selective calls.

    By default, a row is called only when its prediction is not an abstention
    and its finite confidence is at least the threshold.

    ``candidate_prediction`` plus ``threshold_eligible`` activates the
    production-aware mode: only eligible inference-lane rows are rethresholded;
    authoritative anchors, ``Other``, unavailable cells, and conflict
    abstentions retain their existing disposition. ``Other`` is a valid
    reference negative but is not counted as a named call. This is required when one
    frozen probabilistic gate is being varied inside a larger hierarchical
    decision contract. The selective error is undefined (``NaN``) when no
    cells are called.
    """

    truth = _as_series(reference, name="reference")
    called = _as_series(prediction, name="prediction")
    if len(truth) != len(called):
        raise ValueError("reference and prediction must have the same length")
    if truth.isna().any():
        raise ValueError("reference contains missing values")
    abstentions = _validated_abstentions(abstain_labels)
    non_calls = _validated_abstentions(non_call_labels)
    if truth.map(lambda value: _is_in(value, abstentions)).any():
        raise ValueError(
            "reference contains abstention labels; adjudicate or exclude those rows first"
        )
    candidate, eligible = _validated_threshold_lane(
        candidate_prediction,
        threshold_eligible,
        length=len(truth),
        abstentions=abstentions,
    )
    try:
        score = np.asarray(confidence, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("confidence must be numeric") from exc
    if score.ndim != 1 or len(score) != len(truth):
        raise ValueError(
            f"confidence must be one-dimensional with length {len(truth)}"
        )
    finite_score = np.isfinite(score)
    invalid_score = finite_score & ((score < 0) | (score > 1))
    if invalid_score.any() or np.isinf(score).any():
        raise ValueError("finite confidence values must lie in [0, 1]")
    sample_weight = _validated_weights(weights, length=len(truth))

    if thresholds is None:
        finite_values = score[
            finite_score if eligible is None else finite_score & eligible
        ]
        threshold_values = np.unique(
            np.concatenate(([0.0], finite_values, [1.0]))
        )
    else:
        if isinstance(thresholds, (str, bytes)) or np.isscalar(thresholds):
            threshold_values = np.asarray([thresholds], dtype=float)
        else:
            try:
                threshold_values = np.asarray(thresholds, dtype=float)
            except (TypeError, ValueError) as exc:
                raise ValueError("thresholds must be numeric") from exc
        if threshold_values.ndim != 1 or not np.isfinite(threshold_values).all():
            raise ValueError("thresholds must be a finite one-dimensional sequence")
        if ((threshold_values < 0) | (threshold_values > 1)).any():
            raise ValueError("thresholds must lie in [0, 1]")
        threshold_values = np.unique(threshold_values)

    called.loc[called.isna()] = _MISSING_PREDICTION
    total_weight = float(sample_weight.sum())

    rows = []
    for threshold in threshold_values:
        thresholded = called.copy()
        if eligible is None:
            score_pass = finite_score & (score >= threshold)
        else:
            score_pass = eligible & finite_score & (score >= threshold)
            thresholded.iloc[np.flatnonzero(eligible & ~score_pass)] = (
                _MISSING_PREDICTION
            )
            thresholded.iloc[np.flatnonzero(score_pass)] = candidate.iloc[
                np.flatnonzero(score_pass)
            ].to_numpy()
        is_non_call = thresholded.map(
            lambda value: value is _MISSING_PREDICTION
            or _is_in(value, non_calls)
        ).to_numpy(dtype=bool)
        selected = ~is_non_call & (
            score_pass if eligible is None else np.ones(len(truth), dtype=bool)
        )
        correct = thresholded.eq(truth).to_numpy(dtype=bool)
        selected_weight = float(sample_weight[selected].sum())
        correct_weight = float(sample_weight[selected & correct].sum())
        rescue_eligible_weight = (
            float(sample_weight[eligible].sum()) if eligible is not None else np.nan
        )
        rescue_called = (
            eligible & ~is_non_call
            if eligible is not None
            else np.zeros(len(truth), dtype=bool)
        )
        rescue_called_weight = (
            float(sample_weight[rescue_called].sum())
            if eligible is not None
            else np.nan
        )
        rescue_correct_weight = (
            float(sample_weight[rescue_called & correct].sum())
            if eligible is not None
            else np.nan
        )
        rows.append(
            {
                "threshold": float(threshold),
                "coverage": selected_weight / total_weight,
                "selective_error": (
                    1.0 - correct_weight / selected_weight
                    if selected_weight > 0
                    else np.nan
                ),
                "correct_weight": correct_weight,
                "called_weight": selected_weight,
                "total_weight": total_weight,
                "rescue_eligible_weight": rescue_eligible_weight,
                "rescue_called_weight": rescue_called_weight,
                "rescue_coverage": (
                    rescue_called_weight / rescue_eligible_weight
                    if eligible is not None and rescue_eligible_weight > 0
                    else np.nan
                ),
                "rescue_selective_error": (
                    1.0 - rescue_correct_weight / rescue_called_weight
                    if eligible is not None and rescue_called_weight > 0
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(
        rows,
        columns=[
            "threshold",
            "coverage",
            "selective_error",
            "correct_weight",
            "called_weight",
            "total_weight",
            "rescue_eligible_weight",
            "rescue_called_weight",
            "rescue_coverage",
            "rescue_selective_error",
        ],
    )


def grouped_bootstrap_coverage_risk(
    reference: Sequence[Any],
    prediction: Sequence[Any],
    score: Sequence[float],
    groups: Sequence[Any],
    thresholds: Sequence[float] | None = None,
    weights: Sequence[float] | None = None,
    abstain_labels: Sequence[Any] = DEFAULT_ABSTAIN_LABELS,
    non_call_labels: Sequence[Any] = DEFAULT_NON_CALL_LABELS,
    *,
    candidate_prediction: Sequence[Any] | None = None,
    threshold_eligible: Sequence[bool] | None = None,
    replicates: int = 1_000,
    confidence: float = 0.95,
    seed: int = 0,
) -> pd.DataFrame:
    """Add donor/ROI-clustered intervals to a coverage-versus-risk curve.

    Whole groups are sampled with replacement while every replicate is
    evaluated at the point curve's fixed thresholds.  This makes uncertainty
    bands comparable along the curve and preserves within-group cell
    dependence.  A replicate with no positive total weight is omitted from
    both intervals; selective-error intervals also omit replicates with no
    calls at that threshold.  The corresponding valid-replicate counts make
    those omissions explicit.

    Rows are canonically ordered before aggregation and group identifiers are
    canonically sorted before random draws, so a fixed seed gives the same
    result after an input-row permutation.
    """

    truth = _as_series(reference, name="reference")
    called = _as_series(prediction, name="prediction")
    cluster = _as_series(groups, name="groups")
    if len(truth) != len(called) or len(truth) != len(cluster):
        raise ValueError(
            "reference, prediction, score, and groups must have the same length"
        )
    try:
        score_values = np.asarray(score, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError("score must be numeric") from exc
    if score_values.ndim != 1 or len(score_values) != len(truth):
        raise ValueError(f"score must be one-dimensional with length {len(truth)}")
    if cluster.isna().any():
        raise ValueError("groups contains missing values")
    abstentions = _validated_abstentions(abstain_labels)
    candidate, eligible = _validated_threshold_lane(
        candidate_prediction,
        threshold_eligible,
        length=len(truth),
        abstentions=abstentions,
    )

    n_replicates = _validate_positive_integer(replicates, name="replicates")
    if not np.isfinite(confidence) or not 0 < confidence < 1:
        raise ValueError("confidence must lie strictly between 0 and 1")
    if isinstance(seed, bool) or not isinstance(seed, Integral):
        raise ValueError("seed must be an integer")
    sample_weight = _validated_weights(weights, length=len(truth))

    # Deterministic row ordering eliminates floating-point aggregation changes
    # caused solely by an input permutation.  Duplicate canonical rows are
    # interchangeable because they make identical contributions to the curve.
    cluster_tokens = cluster.map(lambda value: _canonical_token((value,)))
    row_tokens = np.asarray(
        [
            _canonical_token(
                (
                    group_token,
                    truth_value,
                    called_value,
                    score_value,
                    weight_value,
                    None if candidate is None else candidate_value,
                    None if eligible is None else eligible_value,
                )
            )
            for (
                group_token,
                truth_value,
                called_value,
                score_value,
                weight_value,
                candidate_value,
                eligible_value,
            ) in zip(
                cluster_tokens,
                truth,
                called,
                score_values,
                sample_weight,
                (
                    np.full(len(truth), None, dtype=object)
                    if candidate is None
                    else candidate
                ),
                (
                    np.full(len(truth), None, dtype=object)
                    if eligible is None
                    else eligible
                ),
                strict=True,
            )
        ],
        dtype=str,
    )
    order = np.argsort(row_tokens, kind="stable")
    truth = truth.iloc[order].reset_index(drop=True)
    called = called.iloc[order].reset_index(drop=True)
    score_values = score_values[order]
    sample_weight = sample_weight[order]
    cluster_tokens = cluster_tokens.iloc[order].reset_index(drop=True)
    if candidate is not None:
        candidate = candidate.iloc[order].reset_index(drop=True)
        eligible = eligible[order]

    group_codes, unique_groups = pd.factorize(cluster_tokens, sort=True)
    n_groups = len(unique_groups)
    if n_groups < 2:
        raise ValueError("grouped bootstrap requires at least two groups")

    point = coverage_risk_curve(
        truth,
        called,
        score_values,
        thresholds=thresholds,
        weights=sample_weight,
        abstain_labels=abstain_labels,
        non_call_labels=non_call_labels,
        candidate_prediction=candidate,
        threshold_eligible=eligible,
    )
    threshold_values = point["threshold"].to_numpy(dtype=float)
    n_thresholds = len(threshold_values)
    coverage_values = np.full((n_replicates, n_thresholds), np.nan, dtype=float)
    error_values = np.full((n_replicates, n_thresholds), np.nan, dtype=float)

    rng = np.random.default_rng(int(seed))
    for replicate in range(n_replicates):
        sampled_groups = rng.integers(0, n_groups, size=n_groups)
        multiplicity = np.bincount(sampled_groups, minlength=n_groups)
        replicate_weights = sample_weight * multiplicity[group_codes]
        if replicate_weights.sum() <= 0:
            continue
        replicate_curve = coverage_risk_curve(
            truth,
            called,
            score_values,
            thresholds=threshold_values,
            weights=replicate_weights,
            abstain_labels=abstain_labels,
            non_call_labels=non_call_labels,
            candidate_prediction=candidate,
            threshold_eligible=eligible,
        )
        coverage_values[replicate] = replicate_curve["coverage"].to_numpy(
            dtype=float
        )
        error_values[replicate] = replicate_curve["selective_error"].to_numpy(
            dtype=float
        )

    alpha = (1.0 - float(confidence)) / 2.0

    def interval(values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        low = np.full(n_thresholds, np.nan, dtype=float)
        high = np.full(n_thresholds, np.nan, dtype=float)
        valid = np.zeros(n_thresholds, dtype=int)
        for index in range(n_thresholds):
            finite = values[:, index]
            finite = finite[np.isfinite(finite)]
            valid[index] = len(finite)
            if len(finite):
                low[index], high[index] = np.quantile(
                    finite,
                    [alpha, 1.0 - alpha],
                )
        return low, high, valid

    coverage_low, coverage_high, coverage_valid = interval(coverage_values)
    error_low, error_high, error_valid = interval(error_values)
    result = point.copy()
    result["coverage_ci_low"] = coverage_low
    result["coverage_ci_high"] = coverage_high
    result["selective_error_ci_low"] = error_low
    result["selective_error_ci_high"] = error_high
    result["coverage_valid_replicates"] = coverage_valid
    result["selective_error_valid_replicates"] = error_valid
    result["replicates"] = n_replicates
    result["groups"] = n_groups
    result["confidence_level"] = float(confidence)
    return result


def _normalised_text(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series("", index=frame.index, dtype="object")
    return frame[column].fillna("").astype(str).str.strip().str.lower()


def _numeric_column(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series(np.nan, index=frame.index, dtype=float)
    return pd.to_numeric(frame[column], errors="coerce")


def _candidate_present(frame: pd.DataFrame, column: str | None) -> pd.Series:
    if column is None or column not in frame:
        return pd.Series(True, index=frame.index, dtype=bool)
    text = frame[column].fillna("").astype(str).str.strip().str.lower()
    return ~text.isin({"", "nan", "none"})


def challenge_flags(
    assignments: pd.DataFrame,
    low_margin: float = 0.25,
    high_probability: float = 0.8,
) -> pd.DataFrame:
    """Add deterministic strata for false-positive/negative challenge review.

    The flags use only standard assignment audit fields.  They deliberately
    over-sample plausible failure modes while retaining ``random_baseline`` on
    every row, which is necessary to estimate error rates instead of reviewing
    only suspicious cells.

    False-positive candidates are accepted calls with a low model margin or an
    accepted status paired with an explicitly contradictory reason.  False-
    negative candidates are high-probability abstentions, anchor conflicts,
    ``Other`` calls that still retain a ranked candidate, geometry-QC excluded
    objects, or parent-accepted cells with an unresolved subtype.  These are
    screening heuristics: image/marker review or an independent reference is
    required to establish a real error.
    """

    if not isinstance(assignments, pd.DataFrame):
        raise TypeError("assignments must be a pandas DataFrame")
    for value, name in (
        (low_margin, "low_margin"),
        (high_probability, "high_probability"),
    ):
        if not np.isfinite(value) or value < 0 or value > 1:
            raise ValueError(f"{name} must be finite and lie in [0, 1]")

    out = assignments.copy()
    accepted_statuses = {"anchor", "inferred", "accepted"}
    abstained_statuses = {"ambiguous", "unavailable"}
    contradictory_pattern = (
        r"conflict|failed|unavailable|ineligible|no_identity_evidence|"
        r"required_models_unavailable|parent_not_classified|"
        r"all_required_models_valid_and_negative|terminal_parent"
    )

    level_specs: list[tuple[str, str | None, str | None, str | None, str | None]] = []
    if "broad_assignment_status" in out:
        level_specs.append(
            (
                "broad_assignment_status",
                "broad_margin",
                "broad_best_probability",
                "broad_best_candidate",
                "broad_reason",
            )
        )
    if "subtype_assignment_status" in out:
        level_specs.append(
            (
                "subtype_assignment_status",
                "subtype_margin",
                "subtype_best_probability",
                "subtype_best_candidate",
                "subtype_reason",
            )
        )
    elif "specific_assignment_status" in out:
        level_specs.append(
            (
                "specific_assignment_status",
                next(
                    (name for name in ("specific_margin", "subtype_margin") if name in out),
                    None,
                ),
                next(
                    (
                        name
                        for name in ("specific_best_probability", "subtype_best_probability")
                        if name in out
                    ),
                    None,
                ),
                next(
                    (
                        name
                        for name in ("best_specific_type", "specific_best_candidate")
                        if name in out
                    ),
                    None,
                ),
                "specific_reason" if "specific_reason" in out else None,
            )
        )
    if "assignment_status" in out and not level_specs:
        generic_margin = next(
            (name for name in ("margin", "specific_margin", "broad_margin") if name in out),
            None,
        )
        generic_probability = next(
            (
                name
                for name in ("confidence", "best_probability", "broad_best_probability")
                if name in out
            ),
            None,
        )
        generic_candidate = next(
            (
                name
                for name in ("best_specific_type", "best_broad_type", "best_candidate")
                if name in out
            ),
            None,
        )
        generic_reason = next(
            (name for name in ("assignment_reason", "reason", "broad_reason") if name in out),
            None,
        )
        level_specs.append(
            (
                "assignment_status",
                generic_margin,
                generic_probability,
                generic_candidate,
                generic_reason,
            )
        )

    accepted_low_margin = pd.Series(False, index=out.index, dtype=bool)
    status_reason_contradiction = pd.Series(False, index=out.index, dtype=bool)
    abstained_high_probability = pd.Series(False, index=out.index, dtype=bool)
    other_with_candidate = pd.Series(False, index=out.index, dtype=bool)
    anchor_conflict = pd.Series(False, index=out.index, dtype=bool)

    for status_col, margin_col, probability_col, candidate_col, reason_col in level_specs:
        status = _normalised_text(out, status_col)
        accepted = status.isin(accepted_statuses)
        if margin_col is not None:
            margin = _numeric_column(out, margin_col)
            accepted_low_margin |= accepted & margin.notna() & (margin < low_margin)

        reason = (
            _normalised_text(out, reason_col)
            if reason_col is not None
            else pd.Series("", index=out.index, dtype="object")
        )
        contradictory = reason.str.contains(contradictory_pattern, regex=True, na=False)
        status_reason_contradiction |= accepted & contradictory
        anchor_conflict |= reason.str.contains(
            r"anchor.*conflict|conflicting.*anchor", regex=True, na=False
        )

        if probability_col is not None:
            probability = _numeric_column(out, probability_col)
            high_candidate = probability.notna() & (probability >= high_probability)
            abstained_high_probability |= status.isin(abstained_statuses) & high_candidate
            other_with_candidate |= (
                status.eq("other")
                & probability.notna()
                & _candidate_present(out, candidate_col)
            )

    for anchor_column in (
        "broad_anchor_classes",
        "subtype_anchor_classes",
        "anchor_classes",
    ):
        if anchor_column in out:
            anchor_text = out[anchor_column].fillna("").astype(str)
            anchor_conflict |= anchor_text.str.count(r"\|").ge(1)

    geometry_qc_ineligible = pd.Series(False, index=out.index, dtype=bool)
    if "qc_analysis_eligible" in out:
        geometry = out["qc_analysis_eligible"]
        geometry_qc_ineligible = geometry.map(
            lambda value: (
                False
                if pd.isna(value)
                else str(value).strip().lower() in {"false", "0", "no"}
            )
        ).astype(bool)

    subtype_unresolved = pd.Series(False, index=out.index, dtype=bool)
    if "specific_type" in out:
        subtype_unresolved = out["specific_type"].fillna("").astype(str).str.endswith(
            "-unclassified"
        )
    elif "cell_type" in out:
        subtype_unresolved = out["cell_type"].fillna("").astype(str).str.endswith(
            "-unclassified"
        )

    out["candidate_fp"] = (
        accepted_low_margin | status_reason_contradiction
    ).astype(bool)
    out["candidate_fn"] = (
        abstained_high_probability
        | anchor_conflict
        | other_with_candidate
        | geometry_qc_ineligible
        | subtype_unresolved
    ).astype(bool)

    ordered_strata = [
        ("random_baseline", pd.Series(True, index=out.index, dtype=bool)),
        ("accepted_low_margin", accepted_low_margin),
        ("status_reason_contradiction", status_reason_contradiction),
        ("abstained_high_probability", abstained_high_probability),
        ("anchor_conflict", anchor_conflict),
        ("other_with_candidate", other_with_candidate),
        ("geometry_qc_ineligible", geometry_qc_ineligible),
        ("subtype_unresolved", subtype_unresolved),
    ]
    stratum_masks = [mask.to_numpy(dtype=bool) for _, mask in ordered_strata]
    review_strata = [
        ";".join(
            name
            for (name, _), mask in zip(ordered_strata, stratum_masks, strict=True)
            if mask[position]
        )
        for position in range(len(out))
    ]
    out["review_strata"] = review_strata
    return out


__all__ = [
    "DEFAULT_ABSTAIN_LABELS",
    "DEFAULT_NON_CALL_LABELS",
    "DEFAULT_SAMPLE_SALT",
    "POSTERIOR_PROBABILITY_SEMANTICS",
    "challenge_flags",
    "coverage_risk_curve",
    "deterministic_stratified_sample",
    "grouped_bootstrap_coverage_risk",
    "grouped_bootstrap_metrics",
    "one_vs_rest_reliability_table",
    "per_class_metrics",
]
