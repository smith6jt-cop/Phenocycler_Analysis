"""Explicit donor-local reference-positive seeds for marker calibration.

This module supplies the input that ``marker_calibration`` deliberately refuses
to infer.  For every distinct reference in the validated marker registry it
selects a fixed-size, high-rank seed set using only:

* that reference marker's registry-selected intensity;
* the explicit QC estimation-eligibility mask;
* cell coordinates; and
* stable object IDs.

Selection is deterministic and target-independent.  There is no target value,
target prevalence, phenotype annotation, mixture model, GMM, NNMF, or Otsu
threshold in this API.  Candidate controls must have a finite positive
reference measurement and finite coordinates.  Within each spatial tile,
reference intensity determines rank and object ID breaks ties; selection then
runs round-robin over tiles up to a versioned cap.  A donor/reference with
fewer than the calibration minimum fails closed to an all-false mask.

Externally supplied masks enter only through :func:`load_reference_masks`,
which validates the same canonical long schema, registry provenance, object
universe, fixed count, and rank invariants.  Annotation/class columns are
rejected.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .marker_calibration import MAX_CONTROL_CAP
from .marker_registry import MarkerRegistry, MarkerSpec, load_registry


MASK_SCHEMA_VERSION = 1
METHOD_VERSION = "reference-control-rank-v1"
EXTERNAL_METHOD_PATTERN = re.compile(r"^external-reference-controls-v[1-9][0-9]*$")

MASK_COLUMNS = (
    "mask_schema_version",
    "method_version",
    "registry_version",
    "registry_fingerprint",
    "source",
    "donor_id",
    "reference",
    "object_id",
    "reference_positive",
    "seed_rank",
)
FORBIDDEN_EXTERNAL_COLUMN_FRAGMENTS = (
    "annotation",
    "cell_type",
    "classification",
    "phenotype",
    "class_label",
)


class ReferenceControlStatus(str, Enum):
    """Disposition of one donor/reference seed set."""

    VALID = "valid"
    INSUFFICIENT_CONTROLS = "insufficient_controls"


@dataclass(frozen=True)
class ReferenceControlConfig:
    """Fixed-count seed-selection policy."""

    seed_cap: int = MAX_CONTROL_CAP
    min_controls: int = 10_000
    spatial_bins: int = 8
    random_seed: int = 0

    def __post_init__(self) -> None:
        if not (1 <= self.seed_cap <= MAX_CONTROL_CAP):
            raise ValueError(
                f"seed_cap must be in [1, {MAX_CONTROL_CAP}]"
            )
        if not (2 <= self.min_controls <= self.seed_cap):
            raise ValueError("min_controls must be in [2, seed_cap]")
        if self.spatial_bins < 1:
            raise ValueError("spatial_bins must be positive")


@dataclass(frozen=True)
class ReferenceControlResult:
    """Canonical long masks and one audit row per donor/reference."""

    donor_id: str
    masks_long: pd.DataFrame
    audit: pd.DataFrame

    def as_arrays(
        self, object_ids: Sequence[object] | np.ndarray
    ) -> dict[str, np.ndarray]:
        return reference_masks_as_arrays(
            self.masks_long,
            object_ids=object_ids,
            donor_id=self.donor_id,
        )


def registry_references(registry: MarkerRegistry | None = None) -> tuple[str, ...]:
    """Distinct active references in stable registry order."""

    active_registry = load_registry() if registry is None else registry
    seen: set[str] = set()
    references: list[str] = []
    for marker in active_registry.active_markers:
        reference = marker.reference
        if reference is not None and reference not in seen:
            references.append(reference)
            seen.add(reference)
    return tuple(references)


def reference_masks_to_wide(
    masks: pd.DataFrame,
    *,
    registry: MarkerRegistry | None = None,
) -> pd.DataFrame:
    """Compact canonical long masks to one row per donor/cell.

    The long schema remains the external interchange contract. Production
    persistence uses this lossless wide form so object IDs remain unique and a
    four-reference panel does not store each cell key four times.
    """

    active_registry = load_registry() if registry is None else registry
    _validate_mask_provenance(masks, registry=active_registry)
    keys = masks.loc[:, ["donor_id", "object_id"]].drop_duplicates(
        ignore_index=True
    )
    if keys.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("reference-mask keys are not unique")
    result = keys.copy()
    key_index = pd.MultiIndex.from_frame(keys)
    for reference in registry_references(active_registry):
        group = masks.loc[masks["reference"].astype(str) == reference].copy()
        if len(group) != len(keys):
            raise ValueError(
                f"{reference}: long mask does not cover every donor/object key"
            )
        group["donor_id"] = group["donor_id"].astype(str)
        group["object_id"] = group["object_id"].astype(str)
        group = group.set_index(["donor_id", "object_id"]).reindex(key_index)
        for field in ("reference_positive", "seed_rank", "source", "method_version"):
            result[f"{reference}__{field}"] = group[field].to_numpy()
    result["mask_schema_version"] = MASK_SCHEMA_VERSION
    result["registry_version"] = active_registry.registry_version
    result["registry_fingerprint"] = active_registry.fingerprint
    return result


def reference_masks_from_wide(
    masks_wide: pd.DataFrame,
    *,
    registry: MarkerRegistry | None = None,
) -> pd.DataFrame:
    """Expand the compact production schema to canonical long interchange."""

    active_registry = load_registry() if registry is None else registry
    required = {
        "donor_id",
        "object_id",
        "mask_schema_version",
        "registry_version",
        "registry_fingerprint",
    }
    missing = sorted(required.difference(masks_wide.columns))
    if missing:
        raise ValueError(f"wide reference masks are missing columns: {missing}")
    if masks_wide.duplicated(["donor_id", "object_id"]).any():
        raise ValueError("wide reference masks contain duplicate cells")
    frames = []
    for reference in registry_references(active_registry):
        columns = {
            field: f"{reference}__{field}"
            for field in (
                "reference_positive",
                "seed_rank",
                "source",
                "method_version",
            )
        }
        absent = sorted(set(columns.values()).difference(masks_wide.columns))
        if absent:
            raise ValueError(
                f"wide reference masks are missing {reference} columns: {absent}"
            )
        frames.append(
            pd.DataFrame(
                {
                    "mask_schema_version": masks_wide["mask_schema_version"],
                    "method_version": masks_wide[columns["method_version"]],
                    "registry_version": masks_wide["registry_version"],
                    "registry_fingerprint": masks_wide["registry_fingerprint"],
                    "source": masks_wide[columns["source"]],
                    "donor_id": masks_wide["donor_id"],
                    "reference": reference,
                    "object_id": masks_wide["object_id"],
                    "reference_positive": masks_wide[
                        columns["reference_positive"]
                    ],
                    "seed_rank": masks_wide[columns["seed_rank"]],
                }
            )
        )
    return pd.concat(frames, ignore_index=True).loc[:, MASK_COLUMNS]


def _as_1d(values: Sequence[object] | np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(values)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    return result


def _as_bool(values: Sequence[object] | np.ndarray, name: str) -> np.ndarray:
    raw = _as_1d(values, name)
    if pd.isna(raw).any():
        raise ValueError(f"{name} must not contain missing values")
    if raw.dtype.kind != "b":
        if not bool(np.all(np.isin(raw, [True, False, 0, 1]))):
            raise ValueError(f"{name} must contain only boolean values")
    return raw.astype(bool, copy=False)


def _validate_cell_vectors(
    reference_intensity: Sequence[float] | np.ndarray,
    estimation_eligible: Sequence[bool] | np.ndarray,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    object_ids: Sequence[object] | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    intensity = _as_1d(reference_intensity, "reference_intensity").astype(
        float, copy=False
    )
    eligible = _as_bool(estimation_eligible, "estimation_eligible")
    x_values = _as_1d(x, "x").astype(float, copy=False)
    y_values = _as_1d(y, "y").astype(float, copy=False)
    ids = _as_1d(object_ids, "object_ids")
    lengths = {
        "reference_intensity": len(intensity),
        "estimation_eligible": len(eligible),
        "x": len(x_values),
        "y": len(y_values),
        "object_ids": len(ids),
    }
    if len(set(lengths.values())) != 1:
        raise ValueError(f"all cell arrays must have equal length: {lengths}")
    if pd.isna(ids).any():
        raise ValueError("object_ids must not contain missing values")
    ids = ids.astype(str)
    if len(np.unique(ids)) != len(ids):
        raise ValueError("object_ids must be unique")
    return intensity, eligible, x_values, y_values, ids


def _stable_hash(values: np.ndarray, salt: str) -> np.ndarray:
    strings = pd.Series(values.astype(str), dtype="string") + "|" + salt
    return pd.util.hash_pandas_object(strings, index=False).to_numpy(np.uint64)


def _equal_width_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    minimum = float(np.min(values))
    maximum = float(np.max(values))
    if maximum <= minimum:
        return np.zeros(len(values), dtype=np.int64)
    scaled = (values - minimum) / (maximum - minimum)
    return np.minimum((scaled * n_bins).astype(np.int64), n_bins - 1)


def spatially_balanced_rank_indices(
    reference_intensity: Sequence[float] | np.ndarray,
    candidate_mask: Sequence[bool] | np.ndarray,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    object_ids: Sequence[object] | np.ndarray,
    *,
    cap: int,
    spatial_bins: int,
    salt: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return selected row indices and their spatial tile IDs.

    A single spatial bin reduces exactly to the top ``cap`` reference
    intensities.  With multiple bins, rank is first computed within each tile
    and equal rank rounds are interleaved across occupied tiles.  Stable object
    hashes resolve intensity ties.
    """

    intensity = _as_1d(reference_intensity, "reference_intensity").astype(
        float, copy=False
    )
    candidates = _as_bool(candidate_mask, "candidate_mask")
    x_values = _as_1d(x, "x").astype(float, copy=False)
    y_values = _as_1d(y, "y").astype(float, copy=False)
    ids = _as_1d(object_ids, "object_ids")
    if len({len(intensity), len(candidates), len(x_values), len(y_values), len(ids)}) != 1:
        raise ValueError("all control-selection arrays must have equal length")
    if not (1 <= cap <= MAX_CONTROL_CAP):
        raise ValueError(f"cap must be in [1, {MAX_CONTROL_CAP}]")
    if spatial_bins < 1:
        raise ValueError("spatial_bins must be positive")
    if pd.isna(ids).any() or len(np.unique(ids.astype(str))) != len(ids):
        raise ValueError("object_ids must be non-missing and unique")

    candidate_indices = np.flatnonzero(
        candidates
        & np.isfinite(intensity)
        & (intensity > 0)
        & np.isfinite(x_values)
        & np.isfinite(y_values)
    )
    if not len(candidate_indices):
        return np.empty(0, dtype=int), np.empty(0, dtype=int)
    candidate_intensity = intensity[candidate_indices]
    tile = _equal_width_bins(
        x_values[candidate_indices], spatial_bins
    ) + spatial_bins * _equal_width_bins(
        y_values[candidate_indices], spatial_bins
    )
    hashes = _stable_hash(ids.astype(str)[candidate_indices], salt)

    # Sort by tile, descending intensity, stable hash; convert position within
    # each tile to a round number.
    within_order = np.lexsort((hashes, -candidate_intensity, tile))
    ordered_tiles = tile[within_order]
    rounds = np.empty(len(within_order), dtype=np.int64)
    starts = np.r_[0, np.flatnonzero(np.diff(ordered_tiles)) + 1]
    ends = np.r_[starts[1:], len(ordered_tiles)]
    for start, end in zip(starts, ends):
        rounds[within_order[start:end]] = np.arange(end - start)

    # Rank round first, tile second, intensity third.  This takes the brightest
    # cell from every occupied tile before the second-brightest from any tile.
    order = np.lexsort((hashes, -candidate_intensity, tile, rounds))
    chosen = order[: min(cap, len(order))]
    return candidate_indices[chosen], tile[chosen]


def _selection_fingerprint(
    *,
    donor_id: str,
    reference: str,
    registry: MarkerRegistry,
    config: ReferenceControlConfig,
    selected_ids: np.ndarray,
) -> str:
    digest = hashlib.sha256()
    for part in (
        METHOD_VERSION,
        registry.registry_version,
        registry.fingerprint,
        donor_id,
        reference,
        config.seed_cap,
        config.min_controls,
        config.spatial_bins,
        config.random_seed,
    ):
        digest.update(str(part).encode("utf-8"))
        digest.update(b"\0")
    for object_id in selected_ids:
        digest.update(str(object_id).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def _long_mask(
    *,
    donor_id: str,
    reference: str,
    registry: MarkerRegistry,
    object_ids: np.ndarray,
    selected_indices: np.ndarray,
    source: str,
    method_version: str,
) -> pd.DataFrame:
    selected = np.zeros(len(object_ids), dtype=bool)
    selected[selected_indices] = True
    rank = pd.array(np.full(len(object_ids), pd.NA, dtype=object), dtype="Int64")
    if len(selected_indices):
        rank[selected_indices] = np.arange(1, len(selected_indices) + 1)
    return pd.DataFrame(
        {
            "mask_schema_version": MASK_SCHEMA_VERSION,
            "method_version": method_version,
            "registry_version": registry.registry_version,
            "registry_fingerprint": registry.fingerprint,
            "source": source,
            "donor_id": str(donor_id),
            "reference": reference,
            "object_id": object_ids.astype(str),
            "reference_positive": selected,
            "seed_rank": rank,
        }
    )


def build_reference_control(
    reference_intensity: Sequence[float] | np.ndarray,
    estimation_eligible: Sequence[bool] | np.ndarray,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    object_ids: Sequence[object] | np.ndarray,
    *,
    donor_id: str,
    reference: str | MarkerSpec,
    registry: MarkerRegistry | None = None,
    config: ReferenceControlConfig | None = None,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Build one canonical donor/reference mask and audit row."""

    active_registry = load_registry() if registry is None else registry
    config = ReferenceControlConfig() if config is None else config
    spec = (
        active_registry.marker(reference)
        if isinstance(reference, str)
        else reference
    )
    required_references = set(registry_references(active_registry))
    if spec.name not in required_references:
        raise ValueError(
            f"{spec.name} is not an active reference in "
            f"{active_registry.registry_version}"
        )
    intensity, eligible, x_values, y_values, ids = _validate_cell_vectors(
        reference_intensity, estimation_eligible, x, y, object_ids
    )
    candidate = (
        eligible
        & np.isfinite(intensity)
        & (intensity > 0)
        & np.isfinite(x_values)
        & np.isfinite(y_values)
    )
    salt = (
        f"{METHOD_VERSION}|{active_registry.registry_version}|{donor_id}|"
        f"{spec.name}|{config.random_seed}"
    )
    selected, selected_tiles = spatially_balanced_rank_indices(
        intensity,
        candidate,
        x_values,
        y_values,
        ids,
        cap=config.seed_cap,
        spatial_bins=config.spatial_bins,
        salt=salt,
    )
    n_candidates = int(candidate.sum())
    status = ReferenceControlStatus.VALID
    reason = (
        "fixed-count reference-intensity ranks selected deterministically "
        "across spatial tiles"
    )
    if len(selected) < config.min_controls:
        # Fail closed: a partial mask must never become a small, unstable
        # negative-control population downstream.
        status = ReferenceControlStatus.INSUFFICIENT_CONTROLS
        reason = (
            f"{len(selected):,} eligible finite-positive spatial controls; "
            f"need {config.min_controls:,}"
        )
        selected = np.empty(0, dtype=int)
        selected_tiles = np.empty(0, dtype=int)

    mask = _long_mask(
        donor_id=str(donor_id),
        reference=spec.name,
        registry=active_registry,
        object_ids=ids,
        selected_indices=selected,
        source="rank_seed",
        method_version=METHOD_VERSION,
    )
    candidate_values = intensity[candidate]
    selected_values = intensity[selected]
    fingerprint = _selection_fingerprint(
        donor_id=str(donor_id),
        reference=spec.name,
        registry=active_registry,
        config=config,
        selected_ids=ids[selected],
    )
    audit: dict[str, object] = {
        "mask_schema_version": MASK_SCHEMA_VERSION,
        "method_version": METHOD_VERSION,
        "registry_version": active_registry.registry_version,
        "registry_fingerprint": active_registry.fingerprint,
        "selection_fingerprint": fingerprint,
        "source": "rank_seed",
        "donor_id": str(donor_id),
        "reference": spec.name,
        "measurement_column": spec.measurement_column,
        "compartment": spec.compartment,
        "spillover_policy": spec.spillover_policy,
        "status": status.value,
        "status_reason": reason,
        "n_cells": len(ids),
        "n_qc_estimation_eligible": int(eligible.sum()),
        "n_control_candidates": n_candidates,
        "n_selected": len(selected),
        "seed_cap": config.seed_cap,
        "min_controls": config.min_controls,
        "spatial_bins": config.spatial_bins,
        "n_candidate_tiles": int(
            len(
                np.unique(
                    _equal_width_bins(x_values[candidate], config.spatial_bins)
                    + config.spatial_bins
                    * _equal_width_bins(y_values[candidate], config.spatial_bins)
                )
            )
            if n_candidates
            else 0
        ),
        "n_selected_tiles": int(len(np.unique(selected_tiles))),
        "candidate_intensity_min": (
            float(np.min(candidate_values)) if len(candidate_values) else np.nan
        ),
        "candidate_intensity_median": (
            float(np.median(candidate_values)) if len(candidate_values) else np.nan
        ),
        "candidate_intensity_max": (
            float(np.max(candidate_values)) if len(candidate_values) else np.nan
        ),
        "selected_intensity_min": (
            float(np.min(selected_values)) if len(selected_values) else np.nan
        ),
        "selected_intensity_median": (
            float(np.median(selected_values)) if len(selected_values) else np.nan
        ),
        "selected_intensity_max": (
            float(np.max(selected_values)) if len(selected_values) else np.nan
        ),
    }
    return mask, audit


def build_reference_controls(
    cells: pd.DataFrame,
    *,
    donor_id: str,
    registry: MarkerRegistry | None = None,
    references: Iterable[str] | None = None,
    value_columns: Mapping[str, str] | None = None,
    config: ReferenceControlConfig | None = None,
    eligibility_column: str = "qc_estimation_eligible",
    x_column: str = "X_centroid",
    y_column: str = "Y_centroid",
    object_id_column: str = "object_id",
) -> ReferenceControlResult:
    """Build all requested donor-local reference masks from a cell frame."""

    active_registry = load_registry() if registry is None else registry
    config = ReferenceControlConfig() if config is None else config
    requested = (
        registry_references(active_registry)
        if references is None
        else tuple(str(reference) for reference in references)
    )
    if len(requested) != len(set(requested)):
        raise ValueError("references must be unique")
    if not requested:
        raise ValueError("at least one reference is required")
    value_columns = {} if value_columns is None else dict(value_columns)

    common_required = {
        eligibility_column,
        x_column,
        y_column,
        object_id_column,
    }
    missing_common = sorted(common_required.difference(cells.columns))
    if missing_common:
        raise ValueError(
            f"cell frame is missing required columns: {missing_common}"
        )
    masks: list[pd.DataFrame] = []
    audits: list[dict[str, object]] = []
    for reference in requested:
        spec = active_registry.marker(reference)
        value_column = value_columns.get(reference, spec.measurement_column)
        if value_column not in cells:
            raise ValueError(
                f"cell frame is missing {reference} measurement column "
                f"{value_column!r}"
            )
        mask, audit = build_reference_control(
            cells[value_column].to_numpy(),
            cells[eligibility_column].to_numpy(),
            cells[x_column].to_numpy(),
            cells[y_column].to_numpy(),
            cells[object_id_column].to_numpy(),
            donor_id=str(donor_id),
            reference=spec,
            registry=active_registry,
            config=config,
        )
        masks.append(mask)
        audits.append(audit)
    return ReferenceControlResult(
        donor_id=str(donor_id),
        masks_long=pd.concat(masks, ignore_index=True),
        audit=pd.DataFrame.from_records(audits).sort_values(
            ["donor_id", "reference"], kind="stable", ignore_index=True
        ),
    )


def _read_mask_source(source: str | Path | pd.DataFrame) -> pd.DataFrame:
    if isinstance(source, pd.DataFrame):
        return source.copy()
    path = Path(source)
    suffix = path.suffix.lower()
    if suffix in {".parquet", ".pq"}:
        return pd.read_parquet(path)
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path)
    raise ValueError(
        f"unsupported reference-mask file type {suffix!r}; use CSV or Parquet"
    )


def _validate_mask_provenance(
    frame: pd.DataFrame,
    *,
    registry: MarkerRegistry,
) -> None:
    missing = sorted(set(MASK_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"reference masks are missing schema columns: {missing}")
    forbidden = sorted(
        column
        for column in frame.columns
        if any(
            fragment in column.lower()
            for fragment in FORBIDDEN_EXTERNAL_COLUMN_FRAGMENTS
        )
    )
    if forbidden:
        raise ValueError(
            "reference masks must not contain phenotype/annotation columns: "
            f"{forbidden}"
        )
    if frame.empty:
        raise ValueError("reference masks must not be empty")
    if set(pd.to_numeric(frame["mask_schema_version"], errors="coerce")) != {
        MASK_SCHEMA_VERSION
    }:
        raise ValueError(
            f"mask_schema_version must be {MASK_SCHEMA_VERSION}"
        )
    if set(frame["registry_version"].astype(str)) != {
        registry.registry_version
    }:
        raise ValueError("reference-mask registry_version does not match registry")
    if set(frame["registry_fingerprint"].astype(str)) != {
        registry.fingerprint
    }:
        raise ValueError(
            "reference-mask registry_fingerprint does not match registry"
        )
    for source_name, method_version in frame.groupby("source", sort=False)[
        "method_version"
    ]:
        versions = set(method_version.astype(str))
        if len(versions) != 1:
            raise ValueError(
                f"source {source_name!r} carries multiple method versions"
            )
        version = next(iter(versions))
        if source_name == "rank_seed":
            if version != METHOD_VERSION:
                raise ValueError(
                    f"rank_seed masks require method_version {METHOD_VERSION!r}"
                )
        elif source_name == "external":
            if EXTERNAL_METHOD_PATTERN.fullmatch(version) is None:
                raise ValueError(
                    "external masks require an explicit method version like "
                    "'external-reference-controls-v1'"
                )
        else:
            raise ValueError(
                f"unknown reference-mask source {source_name!r}"
            )


def validate_reference_masks(
    masks: pd.DataFrame,
    *,
    donor_id: str,
    object_ids: Sequence[object] | np.ndarray,
    registry: MarkerRegistry | None = None,
    references: Iterable[str] | None = None,
    config: ReferenceControlConfig | None = None,
) -> ReferenceControlResult:
    """Validate canonical masks, align their universe, and fail small sets closed."""

    active_registry = load_registry() if registry is None else registry
    config = ReferenceControlConfig() if config is None else config
    # Canonicalize the row index before group-wise positional updates; external
    # DataFrames commonly arrive with filtered/non-consecutive index labels.
    frame = masks.copy().reset_index(drop=True)
    _validate_mask_provenance(frame, registry=active_registry)
    ids = _as_1d(object_ids, "object_ids")
    if pd.isna(ids).any():
        raise ValueError("expected object_ids must not contain missing values")
    ids = ids.astype(str)
    if len(np.unique(ids)) != len(ids):
        raise ValueError("expected object_ids must be unique")
    frame["donor_id"] = frame["donor_id"].astype(str)
    frame["reference"] = frame["reference"].astype(str)
    frame["object_id"] = frame["object_id"].astype(str)
    if set(frame["donor_id"]) != {str(donor_id)}:
        raise ValueError(
            f"reference masks must contain donor {donor_id!r} only"
        )
    expected_references = (
        registry_references(active_registry)
        if references is None
        else tuple(str(value) for value in references)
    )
    if set(frame["reference"]) != set(expected_references):
        raise ValueError(
            "reference-mask reference set does not match expected set: "
            f"observed={sorted(frame['reference'].unique())}, "
            f"expected={sorted(expected_references)}"
        )
    allowed_references = set(registry_references(active_registry))
    unknown = sorted(set(expected_references).difference(allowed_references))
    if unknown:
        raise ValueError(f"markers are not active references: {unknown}")
    if frame.duplicated(["donor_id", "reference", "object_id"]).any():
        raise ValueError(
            "reference masks contain duplicate donor/reference/object rows"
        )
    expected_id_set = set(ids)
    for reference, group in frame.groupby("reference", sort=False):
        observed = set(group["object_id"])
        if observed != expected_id_set or len(group) != len(ids):
            raise ValueError(
                f"{reference}: mask object universe differs from donor cells "
                f"(observed={len(observed):,}, expected={len(ids):,})"
            )

    positive = _as_bool(
        frame["reference_positive"].to_numpy(), "reference_positive"
    )
    frame["reference_positive"] = positive
    rank_numeric = pd.to_numeric(frame["seed_rank"], errors="coerce")
    bad_false_rank = (~positive) & rank_numeric.notna().to_numpy()
    bad_true_rank = positive & rank_numeric.isna().to_numpy()
    if bad_false_rank.any() or bad_true_rank.any():
        raise ValueError(
            "seed_rank must be present exactly for reference-positive rows"
        )

    audits: list[dict[str, object]] = []
    for reference, group_index in frame.groupby(
        "reference", sort=False
    ).groups.items():
        loc = np.asarray(list(group_index), dtype=int)
        provenance = frame.loc[loc, ["source", "method_version"]].drop_duplicates()
        if len(provenance) != 1:
            raise ValueError(
                f"{reference}: every row must share one source/method_version"
            )
        group_positive = frame.loc[loc, "reference_positive"].to_numpy(bool)
        group_ranks = pd.to_numeric(
            frame.loc[loc, "seed_rank"], errors="coerce"
        )
        selected_n = int(group_positive.sum())
        if selected_n:
            selected_ranks = group_ranks[group_positive].to_numpy(dtype=float)
            if (
                (selected_ranks <= 0).any()
                or not np.allclose(selected_ranks, np.round(selected_ranks))
            ):
                raise ValueError(
                    f"{reference}: selected seed ranks must be positive integers"
                )
            observed_ranks = np.sort(selected_ranks.astype(int))
            expected_ranks = np.arange(1, selected_n + 1)
            if not np.array_equal(observed_ranks, expected_ranks):
                raise ValueError(
                    f"{reference}: selected seed ranks must be exactly "
                    f"1..{selected_n}"
                )
        if selected_n > config.seed_cap:
            raise ValueError(
                f"{reference}: {selected_n:,} controls exceed cap "
                f"{config.seed_cap:,}"
            )
        status = ReferenceControlStatus.VALID
        reason = "canonical versioned reference mask passed schema validation"
        if selected_n < config.min_controls:
            status = ReferenceControlStatus.INSUFFICIENT_CONTROLS
            reason = (
                f"{selected_n:,} supplied controls; need "
                f"{config.min_controls:,}; failed closed"
            )
            frame.loc[loc, "reference_positive"] = False
            frame.loc[loc, "seed_rank"] = pd.NA
            selected_n = 0
        spec = active_registry.marker(reference)
        source = str(frame.loc[loc[0], "source"])
        method_version = str(frame.loc[loc[0], "method_version"])
        audits.append(
            {
                "mask_schema_version": MASK_SCHEMA_VERSION,
                "method_version": method_version,
                "registry_version": active_registry.registry_version,
                "registry_fingerprint": active_registry.fingerprint,
                "selection_fingerprint": np.nan,
                "source": source,
                "donor_id": str(donor_id),
                "reference": reference,
                "measurement_column": spec.measurement_column,
                "compartment": spec.compartment,
                "spillover_policy": spec.spillover_policy,
                "status": status.value,
                "status_reason": reason,
                "n_cells": len(ids),
                "n_qc_estimation_eligible": np.nan,
                "n_control_candidates": np.nan,
                "n_selected": selected_n,
                "seed_cap": config.seed_cap,
                "min_controls": config.min_controls,
                "spatial_bins": np.nan,
                "n_candidate_tiles": np.nan,
                "n_selected_tiles": np.nan,
                "candidate_intensity_min": np.nan,
                "candidate_intensity_median": np.nan,
                "candidate_intensity_max": np.nan,
                "selected_intensity_min": np.nan,
                "selected_intensity_median": np.nan,
                "selected_intensity_max": np.nan,
            }
        )
    frame["seed_rank"] = pd.array(frame["seed_rank"], dtype="Int64")
    frame = frame.loc[:, MASK_COLUMNS].sort_values(
        ["donor_id", "reference", "object_id"],
        kind="stable",
        ignore_index=True,
    )
    return ReferenceControlResult(
        donor_id=str(donor_id),
        masks_long=frame,
        audit=pd.DataFrame.from_records(audits).sort_values(
            ["donor_id", "reference"], kind="stable", ignore_index=True
        ),
    )


def load_reference_masks(
    source: str | Path | pd.DataFrame,
    *,
    donor_id: str,
    object_ids: Sequence[object] | np.ndarray,
    registry: MarkerRegistry | None = None,
    references: Iterable[str] | None = None,
    config: ReferenceControlConfig | None = None,
) -> ReferenceControlResult:
    """Load CSV/Parquet/external frame only through canonical validation."""

    return validate_reference_masks(
        _read_mask_source(source),
        donor_id=str(donor_id),
        object_ids=object_ids,
        registry=registry,
        references=references,
        config=config,
    )


def reference_masks_as_arrays(
    masks_long: pd.DataFrame,
    *,
    object_ids: Sequence[object] | np.ndarray,
    donor_id: str,
) -> dict[str, np.ndarray]:
    """Align canonical long masks to a caller's object-ID order."""

    ids = _as_1d(object_ids, "object_ids")
    if pd.isna(ids).any():
        raise ValueError("object_ids must not contain missing values")
    ids = ids.astype(str)
    if len(np.unique(ids)) != len(ids):
        raise ValueError("object_ids must be unique")
    frame = masks_long.copy()
    required = {"donor_id", "reference", "object_id", "reference_positive"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"reference masks are missing columns: {missing}")
    frame["donor_id"] = frame["donor_id"].astype(str)
    frame["object_id"] = frame["object_id"].astype(str)
    if set(frame["donor_id"]) != {str(donor_id)}:
        raise ValueError(f"reference masks do not contain donor {donor_id!r} only")
    result: dict[str, np.ndarray] = {}
    expected = pd.Index(ids, name="object_id")
    for reference, group in frame.groupby("reference", sort=False):
        if group["object_id"].duplicated().any():
            raise ValueError(f"{reference}: duplicate object IDs")
        indexed = group.set_index("object_id")["reference_positive"]
        if set(indexed.index) != set(expected):
            raise ValueError(f"{reference}: object universe differs from caller")
        result[str(reference)] = indexed.reindex(expected).to_numpy(dtype=bool)
    return result


__all__ = [
    "EXTERNAL_METHOD_PATTERN",
    "MASK_COLUMNS",
    "MASK_SCHEMA_VERSION",
    "METHOD_VERSION",
    "ReferenceControlConfig",
    "ReferenceControlResult",
    "ReferenceControlStatus",
    "build_reference_control",
    "build_reference_controls",
    "load_reference_masks",
    "reference_masks_as_arrays",
    "reference_masks_from_wide",
    "reference_masks_to_wide",
    "registry_references",
    "spatially_balanced_rank_indices",
    "validate_reference_masks",
]
