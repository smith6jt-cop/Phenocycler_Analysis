"""Donor-local, frequency-independent marker calibration.

The estimator intentionally has no mixture model and no target-positive
population discovery.  Its only estimation population is an *explicit*
reference-positive mask supplied by the caller for the marker's validated,
mutually exclusive reference.  Target abundance elsewhere in the donor never
enters the model.

For each donor and marker the method:

1. intersects reference-positive cells with an explicit estimation-eligibility
   mask and valid coordinates/target measurements;
2. deterministically balances those controls over space and caps the sample at
   20,000 cells (10,000 minimum by default);
3. estimates background as median/Qn on ``log1p`` target intensity;
4. marks controls above robust z=6 as possible target contamination and
   invalidates the calibration when more than 5% are contaminated;
5. takes a fixed, registry-allocated empirical upper-tail order statistic;
6. bootstraps that threshold to define an indeterminate interval; and
7. applies the frozen donor-marker model to *all* measurable cells, including
   cells excluded from estimation.

No cohort values, target frequency, maximum, GMM, NNMF, Otsu threshold, or
data-dependent family allocation is used.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from enum import Enum
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from statsmodels.robust.scale import qn_scale

from .marker_registry import MarkerRegistry, MarkerSpec, load_registry


METHOD_VERSION = "donor-marker-calibration-v1"
MAX_CONTROL_CAP = 20_000
EVIDENCE_DELIMITER = "__"


class CalibrationStatus(str, Enum):
    """Disposition of one donor-marker estimator."""

    VALID = "valid"
    INSUFFICIENT_CONTROLS = "insufficient_controls"
    INSUFFICIENT_TAIL_RESOLUTION = "insufficient_tail_resolution"
    EXCESS_CONTROL_CONTAMINATION = "excess_control_contamination"
    INVALID_BACKGROUND = "invalid_background"


class EvidenceState(str, Enum):
    """Per-cell marker state.

    ``indeterminate`` means the measurement exists but the calibration is
    invalid or the value lies inside the bootstrap threshold interval.
    ``unavailable`` means the cell itself has no finite nonnegative target
    measurement.
    """

    POSITIVE = "positive"
    NEGATIVE = "negative"
    INDETERMINATE = "indeterminate"
    UNAVAILABLE = "unavailable"


@dataclass(frozen=True)
class CalibrationConfig:
    """Scientific and computational settings for the donor-local estimator."""

    control_cap: int = MAX_CONTROL_CAP
    min_controls: int = 10_000
    spatial_bins: int = 8
    contamination_z: float = 6.0
    max_contamination_fraction: float = 0.05
    bootstrap_replicates: int = 200
    bootstrap_confidence: float = 0.95
    random_seed: int = 0
    min_log_scale: float = 1e-9

    def __post_init__(self) -> None:
        if not (1 <= self.control_cap <= MAX_CONTROL_CAP):
            raise ValueError(
                f"control_cap must be in [1, {MAX_CONTROL_CAP}]"
            )
        if not (2 <= self.min_controls <= self.control_cap):
            raise ValueError("min_controls must be in [2, control_cap]")
        if self.spatial_bins < 1:
            raise ValueError("spatial_bins must be positive")
        if not np.isfinite(self.contamination_z) or self.contamination_z <= 0:
            raise ValueError("contamination_z must be finite and positive")
        if not (0.0 <= self.max_contamination_fraction < 1.0):
            raise ValueError(
                "max_contamination_fraction must be in [0, 1)"
            )
        if self.bootstrap_replicates < 20:
            raise ValueError("bootstrap_replicates must be at least 20")
        if not (0.0 < self.bootstrap_confidence < 1.0):
            raise ValueError("bootstrap_confidence must be between 0 and 1")
        if not np.isfinite(self.min_log_scale) or self.min_log_scale <= 0:
            raise ValueError("min_log_scale must be finite and positive")


@dataclass(frozen=True)
class DonorMarkerCalibration:
    """One fitted donor-marker model plus its empirical control distribution."""

    method_version: str
    registry_version: str
    registry_fingerprint: str
    donor_id: str
    marker: str
    reference: str
    kind: str
    uses: tuple[str, ...]
    broad_compartment: str | None
    families: tuple[str, ...]
    compartment: str
    spillover_policy: str
    status: CalibrationStatus
    status_reason: str
    marker_alpha: float
    n_cells: int
    n_estimation_eligible: int
    n_reference_positive: int
    n_control_candidates: int
    n_controls_selected: int
    n_controls_clean: int
    n_controls_contaminated: int
    contamination_fraction: float
    background_log1p_median: float
    background_log1p_qn: float
    threshold_log1p: float
    threshold_raw: float
    threshold_ci_low_log1p: float
    threshold_ci_high_log1p: float
    threshold_ci_low_raw: float
    threshold_ci_high_raw: float
    bootstrap_replicates: int
    control_cap: int
    min_controls: int
    spatial_bins: int
    control_indices: tuple[int, ...] = field(repr=False, compare=False)
    control_object_ids: tuple[str, ...] = field(repr=False, compare=False)
    clean_control_log1p: np.ndarray = field(repr=False, compare=False)

    @property
    def valid(self) -> bool:
        return self.status is CalibrationStatus.VALID

    def to_record(self) -> dict[str, object]:
        """Compact, one-row representation suitable for Parquet/CSV.

        The selected control values themselves remain in the in-memory model;
        their counts, robust summaries, and selection policy are recorded here.
        This keeps the long model table compact while per-cell empirical
        evidence is generated during the same calibrated application.
        """

        return {
            "method_version": self.method_version,
            "registry_version": self.registry_version,
            "registry_fingerprint": self.registry_fingerprint,
            "donor_id": self.donor_id,
            "marker": self.marker,
            "reference": self.reference,
            "kind": self.kind,
            "uses": json.dumps(self.uses),
            "broad_compartment": self.broad_compartment,
            "families": json.dumps(self.families),
            "compartment": self.compartment,
            "spillover_policy": self.spillover_policy,
            "status": self.status.value,
            "status_reason": self.status_reason,
            "marker_alpha": self.marker_alpha,
            "n_cells": self.n_cells,
            "n_estimation_eligible": self.n_estimation_eligible,
            "n_reference_positive": self.n_reference_positive,
            "n_control_candidates": self.n_control_candidates,
            "n_controls_selected": self.n_controls_selected,
            "n_controls_clean": self.n_controls_clean,
            "n_controls_contaminated": self.n_controls_contaminated,
            "contamination_fraction": self.contamination_fraction,
            "background_log1p_median": self.background_log1p_median,
            "background_log1p_qn": self.background_log1p_qn,
            "threshold_log1p": self.threshold_log1p,
            "threshold_raw": self.threshold_raw,
            "threshold_ci_low_log1p": self.threshold_ci_low_log1p,
            "threshold_ci_high_log1p": self.threshold_ci_high_log1p,
            "threshold_ci_low_raw": self.threshold_ci_low_raw,
            "threshold_ci_high_raw": self.threshold_ci_high_raw,
            "bootstrap_replicates": self.bootstrap_replicates,
            "control_cap": self.control_cap,
            "min_controls": self.min_controls,
            "spatial_bins": self.spatial_bins,
        }


@dataclass(frozen=True)
class DonorCalibrationResult:
    """All active marker models and long evidence for one donor."""

    donor_id: str
    models: tuple[DonorMarkerCalibration, ...]
    evidence_long: pd.DataFrame

    @property
    def model_long(self) -> pd.DataFrame:
        return models_to_long(self.models)

    def evidence_wide(
        self,
        *,
        value_columns: Sequence[str] = (
            "state",
            "log2_threshold_ratio",
            "empirical_tail_probability",
            "expression_probability",
            "empirical_tail_evidence",
        ),
    ) -> pd.DataFrame:
        return evidence_to_wide(self.evidence_long, value_columns=value_columns)


def _as_1d(values: Sequence[object] | np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(values)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return result


def _as_bool(values: Sequence[object] | np.ndarray, name: str) -> np.ndarray:
    raw = _as_1d(values, name)
    if pd.isna(raw).any():
        raise ValueError(f"{name} must not contain missing values")
    if raw.dtype.kind not in "b":
        allowed = np.isin(raw, [True, False, 0, 1])
        if not bool(np.all(allowed)):
            raise ValueError(f"{name} must contain only boolean values")
    return raw.astype(bool, copy=False)


def _validate_lengths(expected: int, **arrays: np.ndarray) -> None:
    bad = {name: len(values) for name, values in arrays.items() if len(values) != expected}
    if bad:
        raise ValueError(
            f"all cell arrays must have length {expected}; mismatches={bad}"
        )


def _stable_hash(values: np.ndarray, salt: str) -> np.ndarray:
    """Stable uint64 hash independent of Python's hash randomization."""

    strings = pd.Series(values.astype(str), dtype="string") + "|" + salt
    return pd.util.hash_pandas_object(strings, index=False).to_numpy(np.uint64)


def _spatial_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    if maximum <= minimum:
        return np.zeros(len(values), dtype=np.int64)
    scaled = (values - minimum) / (maximum - minimum)
    return np.minimum((scaled * n_bins).astype(np.int64), n_bins - 1)


def spatially_balanced_control_indices(
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    candidate_mask: Sequence[bool] | np.ndarray,
    object_ids: Sequence[object] | np.ndarray,
    *,
    cap: int,
    spatial_bins: int,
    salt: str,
) -> np.ndarray:
    """Select controls round-robin over occupied equal-width spatial tiles.

    Within each tile cells are ordered by a stable hash of object ID and the
    donor-marker salt.  The result is therefore deterministic and invariant to
    input row ordering.  Taking one cell per occupied tile before taking a
    second prevents a dense focus from monopolizing the capped control sample.
    """

    x_values = _as_1d(x, "x").astype(float, copy=False)
    y_values = _as_1d(y, "y").astype(float, copy=False)
    candidates = _as_bool(candidate_mask, "candidate_mask")
    ids = _as_1d(object_ids, "object_ids")
    _validate_lengths(
        len(candidates), x=x_values, y=y_values, object_ids=ids
    )
    if cap < 1:
        raise ValueError("cap must be positive")
    if spatial_bins < 1:
        raise ValueError("spatial_bins must be positive")
    if pd.isna(ids).any():
        raise ValueError("object_ids must not contain missing values")
    id_strings = ids.astype(str)
    if len(np.unique(id_strings)) != len(id_strings):
        raise ValueError("object_ids must be unique")

    candidate_indices = np.flatnonzero(
        candidates & np.isfinite(x_values) & np.isfinite(y_values)
    )
    if len(candidate_indices) <= cap:
        # Still return a canonical order so audits are row-order invariant.
        hashes = _stable_hash(id_strings[candidate_indices], salt)
        return candidate_indices[np.argsort(hashes, kind="stable")]

    cx = x_values[candidate_indices]
    cy = y_values[candidate_indices]
    tiles = _spatial_bins(cx, spatial_bins) + spatial_bins * _spatial_bins(
        cy, spatial_bins
    )
    hashes = _stable_hash(id_strings[candidate_indices], salt)

    # Establish stable within-tile ranks, then sort by rank first.  This is a
    # compact round-robin over tiles without a Python loop over cells.
    within_order = np.lexsort((hashes, tiles))
    sorted_tiles = tiles[within_order]
    rounds = np.empty(len(within_order), dtype=np.int64)
    starts = np.r_[0, np.flatnonzero(np.diff(sorted_tiles)) + 1]
    ends = np.r_[starts[1:], len(sorted_tiles)]
    for start, end in zip(starts, ends):
        rounds[within_order[start:end]] = np.arange(end - start)
    balanced_order = np.lexsort((hashes, tiles, rounds))
    return candidate_indices[balanced_order[:cap]]


def _robust_background(log_values: np.ndarray) -> tuple[float, float]:
    median = float(np.median(log_values))
    scale = float(qn_scale(log_values))
    return median, scale


def _empirical_order_index(n: int, alpha: float) -> int | None:
    """Conformal upper-tail order index, or ``None`` if n cannot resolve alpha."""

    if not (0.0 < alpha < 1.0):
        raise ValueError("marker_alpha must be between 0 and 1")
    # With p=(1 + number of controls >= x)/(n+1), no observation can have
    # p <= alpha when alpha < 1/(n+1).
    if alpha < 1.0 / (n + 1):
        return None
    index = int(np.ceil((1.0 - alpha) * (n + 1))) - 1
    return int(np.clip(index, 0, n - 1))


def empirical_upper_threshold(log_values: np.ndarray, alpha: float) -> float:
    """Fixed empirical upper-tail threshold using the finite-sample order rule."""

    values = np.asarray(log_values, dtype=float)
    if values.ndim != 1 or len(values) < 2:
        raise ValueError("log_values must contain at least two values")
    if not np.isfinite(values).all():
        raise ValueError("log_values must be finite")
    index = _empirical_order_index(len(values), alpha)
    if index is None:
        raise ValueError(
            f"{len(values)} controls cannot resolve marker alpha {alpha:g}"
        )
    return float(np.partition(values, index)[index])


def empirical_tail_probability(
    log_values: Sequence[float] | np.ndarray,
    sorted_control_log_values: Sequence[float] | np.ndarray,
) -> np.ndarray:
    """Finite-sample empirical upper-tail p-values.

    The +1 correction prevents zero p-values and matches the order statistic
    used by :func:`empirical_upper_threshold`.
    """

    values = np.asarray(log_values, dtype=float)
    controls = np.asarray(sorted_control_log_values, dtype=float)
    if controls.ndim != 1 or len(controls) < 2:
        raise ValueError("sorted controls must contain at least two values")
    if not np.isfinite(controls).all() or np.any(np.diff(controls) < 0):
        raise ValueError("sorted controls must be finite and nondecreasing")
    left = np.searchsorted(controls, values, side="left")
    return (1.0 + len(controls) - left) / (len(controls) + 1.0)


def _bootstrap_threshold_interval(
    clean_log_values: np.ndarray,
    *,
    alpha: float,
    replicates: int,
    confidence: float,
    seed: int,
) -> tuple[float, float]:
    n = len(clean_log_values)
    index = _empirical_order_index(n, alpha)
    if index is None:
        raise ValueError(f"{n} controls cannot resolve marker alpha {alpha:g}")
    rng = np.random.default_rng(seed)
    thresholds = np.empty(replicates, dtype=float)
    for iteration in range(replicates):
        sample_indices = rng.integers(0, n, size=n)
        sample = clean_log_values[sample_indices]
        thresholds[iteration] = np.partition(sample, index)[index]
    tail = (1.0 - confidence) / 2.0
    low, high = np.quantile(thresholds, [tail, 1.0 - tail])
    return float(low), float(high)


def _deterministic_seed(*parts: object, base_seed: int) -> int:
    text = "|".join(str(part) for part in (*parts, base_seed))
    digest = hashlib.blake2b(text.encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "little", signed=False)


def _resolve_spec_and_alpha(
    marker: str | MarkerSpec,
    registry: MarkerRegistry | None,
    marker_alpha: float | None,
) -> tuple[MarkerSpec, MarkerRegistry, float]:
    active_registry = load_registry() if registry is None else registry
    if isinstance(marker, str):
        spec = active_registry.marker(marker)
    else:
        spec = marker
        registered = active_registry.by_name.get(spec.name)
        if registered is not None and registered != spec:
            raise ValueError(
                f"marker spec for {spec.name} differs from registry "
                f"{active_registry.registry_version}"
            )
    if not spec.active or spec.reference is None or not spec.mutually_exclusive:
        raise ValueError(
            f"{spec.name} does not have an active mutually-exclusive calibration"
        )
    alpha = (
        active_registry.marker_alpha(spec)
        if marker_alpha is None
        else float(marker_alpha)
    )
    if not (0.0 < alpha < 1.0):
        raise ValueError("marker_alpha must be between 0 and 1")
    return spec, active_registry, alpha


def _empty_model(
    *,
    spec: MarkerSpec,
    registry: MarkerRegistry,
    donor_id: str,
    status: CalibrationStatus,
    reason: str,
    alpha: float,
    n_cells: int,
    n_estimation_eligible: int,
    n_reference_positive: int,
    n_candidates: int,
    config: CalibrationConfig,
    selected_indices: np.ndarray | None = None,
    selected_object_ids: np.ndarray | None = None,
    n_contaminated: int = 0,
    contamination_fraction: float = float("nan"),
    clean_log: np.ndarray | None = None,
    background_median: float = float("nan"),
    background_qn: float = float("nan"),
) -> DonorMarkerCalibration:
    indices = (
        np.empty(0, dtype=int)
        if selected_indices is None
        else np.asarray(selected_indices, dtype=int)
    )
    ids = (
        np.empty(0, dtype=str)
        if selected_object_ids is None
        else np.asarray(selected_object_ids).astype(str)
    )
    clean = (
        np.empty(0, dtype=float)
        if clean_log is None
        else np.sort(np.asarray(clean_log, dtype=float))
    )
    return DonorMarkerCalibration(
        method_version=registry.method_version,
        registry_version=registry.registry_version,
        registry_fingerprint=registry.fingerprint,
        donor_id=str(donor_id),
        marker=spec.name,
        reference=spec.reference or "",
        kind=spec.kind,
        uses=spec.uses,
        broad_compartment=spec.broad_compartment,
        families=spec.families,
        compartment=spec.compartment,
        spillover_policy=spec.spillover_policy,
        status=status,
        status_reason=reason,
        marker_alpha=alpha,
        n_cells=n_cells,
        n_estimation_eligible=n_estimation_eligible,
        n_reference_positive=n_reference_positive,
        n_control_candidates=n_candidates,
        n_controls_selected=len(indices),
        n_controls_clean=len(clean),
        n_controls_contaminated=n_contaminated,
        contamination_fraction=contamination_fraction,
        background_log1p_median=background_median,
        background_log1p_qn=background_qn,
        threshold_log1p=float("nan"),
        threshold_raw=float("nan"),
        threshold_ci_low_log1p=float("nan"),
        threshold_ci_high_log1p=float("nan"),
        threshold_ci_low_raw=float("nan"),
        threshold_ci_high_raw=float("nan"),
        bootstrap_replicates=config.bootstrap_replicates,
        control_cap=config.control_cap,
        min_controls=config.min_controls,
        spatial_bins=config.spatial_bins,
        control_indices=tuple(int(value) for value in indices),
        control_object_ids=tuple(str(value) for value in ids),
        clean_control_log1p=clean,
    )


def fit_marker_calibration(
    target_values: Sequence[float] | np.ndarray,
    reference_positive: Sequence[bool] | np.ndarray,
    estimation_eligible: Sequence[bool] | np.ndarray,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    object_ids: Sequence[object] | np.ndarray,
    *,
    donor_id: str,
    marker: str | MarkerSpec,
    registry: MarkerRegistry | None = None,
    marker_alpha: float | None = None,
    config: CalibrationConfig | None = None,
) -> DonorMarkerCalibration:
    """Fit one donor-marker model from explicit reference-positive controls.

    ``estimation_eligible`` and ``reference_positive`` are required arguments;
    there is no permissive default and no attempt to infer either mask from
    target intensity.
    """

    config = CalibrationConfig() if config is None else config
    spec, active_registry, alpha = _resolve_spec_and_alpha(
        marker, registry, marker_alpha
    )
    target = _as_1d(target_values, "target_values").astype(float, copy=False)
    reference_mask = _as_bool(reference_positive, "reference_positive")
    eligible = _as_bool(estimation_eligible, "estimation_eligible")
    x_values = _as_1d(x, "x").astype(float, copy=False)
    y_values = _as_1d(y, "y").astype(float, copy=False)
    ids = _as_1d(object_ids, "object_ids")
    n = len(target)
    _validate_lengths(
        n,
        reference_positive=reference_mask,
        estimation_eligible=eligible,
        x=x_values,
        y=y_values,
        object_ids=ids,
    )
    if pd.isna(ids).any():
        raise ValueError("object_ids must not contain missing values")
    id_strings = ids.astype(str)
    if len(np.unique(id_strings)) != n:
        raise ValueError("object_ids must be unique")

    measurable = np.isfinite(target) & (target >= 0)
    spatially_defined = np.isfinite(x_values) & np.isfinite(y_values)
    candidates = eligible & reference_mask & measurable & spatially_defined
    n_candidates = int(candidates.sum())
    common = dict(
        spec=spec,
        registry=active_registry,
        donor_id=str(donor_id),
        alpha=alpha,
        n_cells=n,
        n_estimation_eligible=int(eligible.sum()),
        n_reference_positive=int(reference_mask.sum()),
        n_candidates=n_candidates,
        config=config,
    )
    if n_candidates < config.min_controls:
        return _empty_model(
            **common,
            status=CalibrationStatus.INSUFFICIENT_CONTROLS,
            reason=(
                f"{n_candidates:,} eligible, measurable, spatially defined "
                f"reference-positive controls; need {config.min_controls:,}"
            ),
        )

    salt = (
        f"{active_registry.registry_version}|{donor_id}|{spec.name}|"
        f"{config.random_seed}"
    )
    selected = spatially_balanced_control_indices(
        x_values,
        y_values,
        candidates,
        id_strings,
        cap=config.control_cap,
        spatial_bins=config.spatial_bins,
        salt=salt,
    )
    selected_ids = id_strings[selected]
    control_log = np.log1p(target[selected])
    initial_median, initial_qn = _robust_background(control_log)
    if not np.isfinite(initial_median) or not np.isfinite(initial_qn):
        return _empty_model(
            **common,
            status=CalibrationStatus.INVALID_BACKGROUND,
            reason="median/Qn background estimate is not finite",
            selected_indices=selected,
            selected_object_ids=selected_ids,
        )

    scale_for_z = max(initial_qn, config.min_log_scale)
    contamination = (
        (control_log - initial_median) / scale_for_z
    ) > config.contamination_z
    n_contaminated = int(contamination.sum())
    contamination_fraction = n_contaminated / len(control_log)
    clean_log = control_log[~contamination]
    if contamination_fraction > config.max_contamination_fraction:
        return _empty_model(
            **common,
            status=CalibrationStatus.EXCESS_CONTROL_CONTAMINATION,
            reason=(
                f"{100 * contamination_fraction:.3f}% of selected controls "
                f"exceed robust z>{config.contamination_z:g}; allowed at most "
                f"{100 * config.max_contamination_fraction:.3f}%"
            ),
            selected_indices=selected,
            selected_object_ids=selected_ids,
            n_contaminated=n_contaminated,
            contamination_fraction=contamination_fraction,
            clean_log=clean_log,
            background_median=initial_median,
            background_qn=initial_qn,
        )
    if len(clean_log) < config.min_controls:
        return _empty_model(
            **common,
            status=CalibrationStatus.INSUFFICIENT_CONTROLS,
            reason=(
                f"{len(clean_log):,} clean controls remain after contamination "
                f"screen; need {config.min_controls:,}"
            ),
            selected_indices=selected,
            selected_object_ids=selected_ids,
            n_contaminated=n_contaminated,
            contamination_fraction=contamination_fraction,
            clean_log=clean_log,
            background_median=initial_median,
            background_qn=initial_qn,
        )

    background_median, background_qn = _robust_background(clean_log)
    if not np.isfinite(background_median) or not np.isfinite(background_qn):
        return _empty_model(
            **common,
            status=CalibrationStatus.INVALID_BACKGROUND,
            reason="clean-control median/Qn background estimate is not finite",
            selected_indices=selected,
            selected_object_ids=selected_ids,
            n_contaminated=n_contaminated,
            contamination_fraction=contamination_fraction,
            clean_log=clean_log,
            background_median=background_median,
            background_qn=background_qn,
        )

    order_index = _empirical_order_index(len(clean_log), alpha)
    if order_index is None:
        return _empty_model(
            **common,
            status=CalibrationStatus.INSUFFICIENT_TAIL_RESOLUTION,
            reason=(
                f"{len(clean_log):,} clean controls cannot resolve fixed "
                f"marker alpha {alpha:g}; need at least "
                f"{int(np.ceil(1.0 / alpha - 1)):,}"
            ),
            selected_indices=selected,
            selected_object_ids=selected_ids,
            n_contaminated=n_contaminated,
            contamination_fraction=contamination_fraction,
            clean_log=clean_log,
            background_median=background_median,
            background_qn=background_qn,
        )

    threshold_log = float(np.partition(clean_log, order_index)[order_index])
    seed = _deterministic_seed(
        active_registry.registry_version,
        donor_id,
        spec.name,
        base_seed=config.random_seed,
    )
    ci_low, ci_high = _bootstrap_threshold_interval(
        clean_log,
        alpha=alpha,
        replicates=config.bootstrap_replicates,
        confidence=config.bootstrap_confidence,
        seed=seed,
    )
    # The nominal estimate should always lie in the reported uncertainty
    # interval, even for a discrete/tied control distribution.
    ci_low = min(ci_low, threshold_log)
    ci_high = max(ci_high, threshold_log)
    threshold_raw = float(np.expm1(threshold_log))
    ci_low_raw = float(np.expm1(ci_low))
    ci_high_raw = float(np.expm1(ci_high))
    sorted_clean = np.sort(clean_log)

    return DonorMarkerCalibration(
        method_version=active_registry.method_version,
        registry_version=active_registry.registry_version,
        registry_fingerprint=active_registry.fingerprint,
        donor_id=str(donor_id),
        marker=spec.name,
        reference=spec.reference,
        kind=spec.kind,
        uses=spec.uses,
        broad_compartment=spec.broad_compartment,
        families=spec.families,
        compartment=spec.compartment,
        spillover_policy=spec.spillover_policy,
        status=CalibrationStatus.VALID,
        status_reason=(
            "fixed spatially balanced reference-positive controls passed "
            "contamination and empirical-tail checks"
        ),
        marker_alpha=alpha,
        n_cells=n,
        n_estimation_eligible=int(eligible.sum()),
        n_reference_positive=int(reference_mask.sum()),
        n_control_candidates=n_candidates,
        n_controls_selected=len(selected),
        n_controls_clean=len(sorted_clean),
        n_controls_contaminated=n_contaminated,
        contamination_fraction=contamination_fraction,
        background_log1p_median=background_median,
        background_log1p_qn=background_qn,
        threshold_log1p=threshold_log,
        threshold_raw=threshold_raw,
        threshold_ci_low_log1p=ci_low,
        threshold_ci_high_log1p=ci_high,
        threshold_ci_low_raw=ci_low_raw,
        threshold_ci_high_raw=ci_high_raw,
        bootstrap_replicates=config.bootstrap_replicates,
        control_cap=config.control_cap,
        min_controls=config.min_controls,
        spatial_bins=config.spatial_bins,
        control_indices=tuple(int(value) for value in selected),
        control_object_ids=tuple(str(value) for value in selected_ids),
        clean_control_log1p=sorted_clean,
    )


def apply_marker_calibration(
    target_values: Sequence[float] | np.ndarray,
    model: DonorMarkerCalibration,
    *,
    object_ids: Sequence[object] | np.ndarray | None = None,
) -> pd.DataFrame:
    """Apply a fitted calibration to every measurable target value.

    Estimation eligibility is intentionally absent from this signature.
    Ineligible-but-measurable cells receive the same calibrated evidence as
    every other cell.
    """

    target = _as_1d(target_values, "target_values").astype(float, copy=False)
    n = len(target)
    if object_ids is None:
        ids = np.arange(n).astype(str)
    else:
        ids = _as_1d(object_ids, "object_ids")
        _validate_lengths(n, object_ids=ids)
        if pd.isna(ids).any():
            raise ValueError("object_ids must not contain missing values")
        ids = ids.astype(str)

    measurable = np.isfinite(target) & (target >= 0)
    states = np.full(n, EvidenceState.UNAVAILABLE.value, dtype="<U13")
    states[measurable] = EvidenceState.INDETERMINATE.value
    ratio = np.full(n, np.nan, dtype=float)
    tail_probability = np.full(n, np.nan, dtype=float)
    expression_probability = np.full(n, np.nan, dtype=float)
    tail_evidence = np.full(n, np.nan, dtype=float)

    if model.valid:
        if len(model.clean_control_log1p) < 2:
            raise ValueError("valid model is missing its empirical control distribution")
        log_values = np.log1p(target[measurable])
        ratio[measurable] = (
            log_values - model.threshold_log1p
        ) / np.log(2.0)
        probabilities = empirical_tail_probability(
            log_values, model.clean_control_log1p
        )
        tail_probability[measurable] = probabilities
        expression_probability[measurable] = 1.0 - probabilities
        tail_evidence[measurable] = -np.log10(probabilities)

        stable_positive = (
            (log_values > model.threshold_ci_high_log1p)
            & (probabilities <= model.marker_alpha)
        )
        stable_negative = (
            (log_values <= model.threshold_ci_low_log1p)
            & (probabilities > model.marker_alpha)
        )
        measurable_indices = np.flatnonzero(measurable)
        states[measurable_indices[stable_positive]] = EvidenceState.POSITIVE.value
        states[measurable_indices[stable_negative]] = EvidenceState.NEGATIVE.value

    return pd.DataFrame(
        {
            "object_id": ids,
            "donor_id": model.donor_id,
            "marker": model.marker,
            "reference": model.reference,
            "measurement": target,
            "state": states,
            "log2_threshold_ratio": ratio,
            "empirical_tail_probability": tail_probability,
            "expression_probability": expression_probability,
            "empirical_tail_evidence": tail_evidence,
            "threshold_raw": model.threshold_raw,
            "threshold_ci_low_raw": model.threshold_ci_low_raw,
            "threshold_ci_high_raw": model.threshold_ci_high_raw,
            "calibration_status": model.status.value,
            "method_version": model.method_version,
            "registry_version": model.registry_version,
        }
    )


def calibrate_marker_frame(
    cells: pd.DataFrame,
    *,
    donor_id: str,
    marker: str | MarkerSpec,
    reference_positive: str | Sequence[bool] | np.ndarray,
    registry: MarkerRegistry | None = None,
    marker_alpha: float | None = None,
    config: CalibrationConfig | None = None,
    eligibility_column: str = "qc_estimation_eligible",
    x_column: str = "X_centroid",
    y_column: str = "Y_centroid",
    object_id_column: str = "object_id",
    value_column: str | None = None,
) -> tuple[DonorMarkerCalibration, pd.DataFrame]:
    """Fit and apply one registered marker from a donor cell frame."""

    active_registry = load_registry() if registry is None else registry
    spec = active_registry.marker(marker) if isinstance(marker, str) else marker
    target_column = spec.measurement_column if value_column is None else value_column
    required = {
        eligibility_column,
        x_column,
        y_column,
        object_id_column,
        target_column,
    }
    if isinstance(reference_positive, str):
        required.add(reference_positive)
    missing = sorted(required.difference(cells.columns))
    if missing:
        raise ValueError(f"cell frame is missing required columns: {missing}")
    reference_mask = (
        cells[reference_positive].to_numpy()
        if isinstance(reference_positive, str)
        else reference_positive
    )
    model = fit_marker_calibration(
        cells[target_column].to_numpy(),
        reference_mask,
        cells[eligibility_column].to_numpy(),
        cells[x_column].to_numpy(),
        cells[y_column].to_numpy(),
        cells[object_id_column].to_numpy(),
        donor_id=str(donor_id),
        marker=spec,
        registry=active_registry,
        marker_alpha=marker_alpha,
        config=config,
    )
    evidence = apply_marker_calibration(
        cells[target_column].to_numpy(),
        model,
        object_ids=cells[object_id_column].to_numpy(),
    )
    return model, evidence


def _reference_mask_for_marker(
    masks: Mapping[object, str | Sequence[bool] | np.ndarray],
    spec: MarkerSpec,
) -> str | Sequence[bool] | np.ndarray:
    for key in ((spec.name, spec.reference), spec.name, spec.reference):
        if key in masks:
            return masks[key]
    raise KeyError(
        f"no explicit reference-positive mask for "
        f"{spec.name} <- {spec.reference}; supply a mapping key of the pair, "
        "target, or reference"
    )


def calibrate_donor(
    cells: pd.DataFrame,
    *,
    donor_id: str,
    reference_positive: Mapping[
        object, str | Sequence[bool] | np.ndarray
    ],
    registry: MarkerRegistry | None = None,
    markers: Iterable[str] | None = None,
    value_columns: Mapping[str, str] | None = None,
    config: CalibrationConfig | None = None,
    eligibility_column: str = "qc_estimation_eligible",
    x_column: str = "X_centroid",
    y_column: str = "Y_centroid",
    object_id_column: str = "object_id",
) -> DonorCalibrationResult:
    """Calibrate active markers and return long model/evidence products.

    The reference mask mapping may be keyed by ``(target, reference)``, target,
    or reference (checked in that order).  This supports one shared validated
    reference mask without forcing it when a pair-specific mask is needed.
    """

    active_registry = load_registry() if registry is None else registry
    selected_specs = (
        active_registry.active_markers
        if markers is None
        else tuple(active_registry.marker(name) for name in markers)
    )
    models: list[DonorMarkerCalibration] = []
    evidence: list[pd.DataFrame] = []
    value_columns = {} if value_columns is None else dict(value_columns)
    for spec in selected_specs:
        if not spec.active:
            raise ValueError(f"{spec.name} does not have an active calibration")
        mask = _reference_mask_for_marker(reference_positive, spec)
        model, marker_evidence = calibrate_marker_frame(
            cells,
            donor_id=str(donor_id),
            marker=spec,
            reference_positive=mask,
            registry=active_registry,
            config=config,
            eligibility_column=eligibility_column,
            x_column=x_column,
            y_column=y_column,
            object_id_column=object_id_column,
            value_column=value_columns.get(spec.name),
        )
        models.append(model)
        evidence.append(marker_evidence)
    long = evidence_to_long(evidence)
    return DonorCalibrationResult(
        donor_id=str(donor_id),
        models=tuple(models),
        evidence_long=long,
    )


def models_to_long(
    models: Iterable[DonorMarkerCalibration],
) -> pd.DataFrame:
    """One compact row per donor-marker calibration."""

    records = [model.to_record() for model in models]
    if not records:
        return pd.DataFrame()
    return pd.DataFrame.from_records(records).sort_values(
        ["donor_id", "marker"], kind="stable", ignore_index=True
    )


def evidence_to_long(
    evidence: Iterable[pd.DataFrame] | pd.DataFrame,
) -> pd.DataFrame:
    """Validate/concatenate per-marker evidence into canonical long form."""

    frames = [evidence] if isinstance(evidence, pd.DataFrame) else list(evidence)
    if not frames:
        return pd.DataFrame()
    required = {
        "object_id",
        "donor_id",
        "marker",
        "state",
        "log2_threshold_ratio",
        "empirical_tail_probability",
        "expression_probability",
        "empirical_tail_evidence",
    }
    for index, frame in enumerate(frames):
        missing = sorted(required.difference(frame.columns))
        if missing:
            raise ValueError(
                f"evidence frame {index} is missing required columns: {missing}"
            )
    result = pd.concat(frames, ignore_index=True)
    keys = ["donor_id", "object_id", "marker"]
    if result.duplicated(keys).any():
        examples = result.loc[result.duplicated(keys, keep=False), keys].head()
        raise ValueError(
            "duplicate donor/object/marker evidence rows: "
            f"{examples.to_dict('records')}"
        )
    return result.sort_values(keys, kind="stable", ignore_index=True)


def evidence_to_wide(
    evidence_long: pd.DataFrame,
    *,
    value_columns: Sequence[str] = (
        "state",
        "log2_threshold_ratio",
        "empirical_tail_probability",
        "expression_probability",
        "empirical_tail_evidence",
    ),
    delimiter: str = EVIDENCE_DELIMITER,
) -> pd.DataFrame:
    """Pivot canonical long evidence to one row per donor/cell.

    Columns use ``<marker>__<metric>`` so marker names containing a single
    underscore remain unambiguous.
    """

    if not delimiter:
        raise ValueError("delimiter must be non-empty")
    required = {"donor_id", "object_id", "marker", *value_columns}
    missing = sorted(required.difference(evidence_long.columns))
    if missing:
        raise ValueError(f"long evidence is missing required columns: {missing}")
    keys = ["donor_id", "object_id", "marker"]
    if evidence_long.duplicated(keys).any():
        raise ValueError("long evidence has duplicate donor/object/marker rows")
    index_columns = ["donor_id", "object_id"]
    wide_parts: list[pd.DataFrame] = []
    for metric in value_columns:
        part = evidence_long.pivot(
            index=index_columns, columns="marker", values=metric
        )
        part = part.rename(
            columns={marker: f"{marker}{delimiter}{metric}" for marker in part.columns}
        )
        wide_parts.append(part)
    result = pd.concat(wide_parts, axis=1).reset_index()
    marker_order = sorted(evidence_long["marker"].astype(str).unique())
    ordered = index_columns + [
        f"{marker}{delimiter}{metric}"
        for marker in marker_order
        for metric in value_columns
    ]
    result = result.reindex(columns=ordered)
    result.columns.name = None
    return result


def evidence_from_wide(
    evidence_wide: pd.DataFrame,
    *,
    delimiter: str = EVIDENCE_DELIMITER,
) -> pd.DataFrame:
    """Convert a frame written by :func:`evidence_to_wide` back to long form."""

    id_columns = ["donor_id", "object_id"]
    missing_ids = sorted(set(id_columns).difference(evidence_wide.columns))
    if missing_ids:
        raise ValueError(f"wide evidence is missing ID columns: {missing_ids}")
    parsed: dict[str, dict[str, str]] = {}
    for column in evidence_wide.columns:
        if column in id_columns:
            continue
        if delimiter not in column:
            raise ValueError(
                f"wide evidence column does not contain {delimiter!r}: {column}"
            )
        marker, metric = column.rsplit(delimiter, 1)
        if not marker or not metric:
            raise ValueError(f"invalid wide evidence column: {column}")
        parsed.setdefault(marker, {})[metric] = column
    frames: list[pd.DataFrame] = []
    for marker, metrics in parsed.items():
        frame = evidence_wide[id_columns].copy()
        frame["marker"] = marker
        for metric, column in metrics.items():
            frame[metric] = evidence_wide[column].to_numpy()
        frames.append(frame)
    if not frames:
        return pd.DataFrame(columns=[*id_columns, "marker"])
    result = pd.concat(frames, ignore_index=True).sort_values(
        [*id_columns, "marker"], kind="stable", ignore_index=True
    )
    result.columns.name = None
    return result


__all__ = [
    "CalibrationConfig",
    "CalibrationStatus",
    "DonorCalibrationResult",
    "DonorMarkerCalibration",
    "EVIDENCE_DELIMITER",
    "EvidenceState",
    "MAX_CONTROL_CAP",
    "METHOD_VERSION",
    "apply_marker_calibration",
    "calibrate_donor",
    "calibrate_marker_frame",
    "empirical_tail_probability",
    "empirical_upper_threshold",
    "evidence_from_wide",
    "evidence_to_long",
    "evidence_to_wide",
    "fit_marker_calibration",
    "models_to_long",
    "spatially_balanced_control_indices",
]
