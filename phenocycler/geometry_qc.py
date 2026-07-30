"""Pre-REDSEA geometry quality control from tabular geometry metrics.

This stage deliberately consumes no marker intensity and no REDSEA output.  A
geometry exporter/rasterizer can produce the small table described by
``GEOMETRY_METRIC_COLUMNS`` once, after which the canonical analysis universe is
fixed before any marker-specific parameter is estimated.

Three eligibility decisions are kept separate:

``analysis_eligible``
    The object is a defensible segmented cell and may receive marker/type calls.
``estimation_eligible``
    The object may additionally help estimate donor-marker background/models.
    A real cell with effectively no cytoplasm remains analysis-eligible but is
    not allowed to set compartment-dependent parameters.
``spillover_context_eligible``
    The cell polygon may contribute physical boundary context to REDSEA. By
    default a real but analysis-excluded polygon (for example a tiny fragment)
    may still contribute, while a missing/empty geometry or duplicate copy may
    not.

Results are frozen records, not a mutable input DataFrame with flags attached.
``to_frame()`` always returns a new convenience copy.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .artifacts import (
    ContractError,
    DuplicateObjectError,
    GeometryContractError,
    dataframe_content_sha256,
    sha256_json,
)


GEOMETRY_METRIC_COLUMNS: tuple[str, ...] = (
    "object_id",
    "cell_geometry_present",
    "nucleus_geometry_present",
    "centroid_x_um",
    "centroid_y_um",
    "cell_area_um2",
    "cell_area_px",
    "nucleus_area_px",
    "cell_solidity",
)

ANALYSIS_REASON_ORDER: tuple[str, ...] = (
    "missing_cell_geometry",
    "missing_nucleus_geometry",
    "invalid_centroid",
    "invalid_cell_area",
    "cell_area_below_min",
    "cell_area_above_max",
    "invalid_solidity",
    "low_solidity",
    "raster_empty",
    "raster_area_mismatch",
    "nucleus_raster_empty",
    "duplicate_detection",
)


class DuplicateGeometryError(GeometryContractError):
    """Near-identical tabular geometries violate the default fatal policy."""


@dataclass(frozen=True)
class GeometryQCConfig:
    """Frozen donor-local geometry policy."""

    pixel_size_um_x: float
    pixel_size_um_y: float
    min_cell_area_um2: float = 20.0
    max_cell_area_median_multiple: float = 4.0
    min_cell_solidity: float = 0.70
    min_raster_area_ratio: float = 0.50
    max_raster_area_ratio: float = 2.00
    no_cytoplasm_nc_ratio: float = 0.999
    require_nucleus_geometry: bool = True
    missing_geometry_policy: Literal["fatal", "exclude"] = "fatal"
    duplicate_policy: Literal["fatal", "deterministic_keep"] = "fatal"
    duplicate_radius_um: float = 1.0
    duplicate_area_tolerance: float = 0.15
    retain_analysis_excluded_for_spillover: bool = True

    def __post_init__(self) -> None:
        numeric_positive = {
            "pixel_size_um_x": self.pixel_size_um_x,
            "pixel_size_um_y": self.pixel_size_um_y,
            "min_cell_area_um2": self.min_cell_area_um2,
            "max_cell_area_median_multiple": self.max_cell_area_median_multiple,
            "duplicate_radius_um": self.duplicate_radius_um,
        }
        for name, value in numeric_positive.items():
            if not math.isfinite(value) or value <= 0:
                raise ContractError(f"{name} must be finite and > 0")
        if not 0 <= self.min_cell_solidity <= 1:
            raise ContractError("min_cell_solidity must be in [0, 1]")
        if not 0 < self.min_raster_area_ratio <= self.max_raster_area_ratio:
            raise ContractError(
                "raster-area bounds must satisfy 0 < min <= max"
            )
        if not 0 < self.no_cytoplasm_nc_ratio <= 1:
            raise ContractError("no_cytoplasm_nc_ratio must be in (0, 1]")
        if not 0 <= self.duplicate_area_tolerance < 1:
            raise ContractError("duplicate_area_tolerance must be in [0, 1)")
        if self.missing_geometry_policy not in {"fatal", "exclude"}:
            raise ContractError(
                "missing_geometry_policy must be 'fatal' or 'exclude'"
            )
        if self.duplicate_policy not in {"fatal", "deterministic_keep"}:
            raise ContractError(
                "duplicate_policy must be 'fatal' or 'deterministic_keep'"
            )

    @property
    def pixel_area_um2(self) -> float:
        return self.pixel_size_um_x * self.pixel_size_um_y

    @property
    def content_id(self) -> str:
        return sha256_json(asdict(self))


@dataclass(frozen=True)
class CellGeometryEligibility:
    """Immutable eligibility decision for one object."""

    donor_id: str
    object_id: str
    analysis_eligible: bool
    estimation_eligible: bool
    spillover_context_eligible: bool
    analysis_reasons: tuple[str, ...]
    estimation_reasons: tuple[str, ...]
    spillover_context_reasons: tuple[str, ...]
    duplicate_group: str | None
    duplicate_survivor: bool | None
    raster_area_ratio: float | None
    nucleus_cell_area_ratio: float | None

    @property
    def analysis_reason(self) -> str:
        return self.analysis_reasons[0] if self.analysis_reasons else ""

    @property
    def estimation_reason(self) -> str:
        return self.estimation_reasons[0] if self.estimation_reasons else ""

    @property
    def spillover_context_reason(self) -> str:
        return (
            self.spillover_context_reasons[0]
            if self.spillover_context_reasons
            else ""
        )


@dataclass(frozen=True)
class GeometryQCResult:
    """Immutable donor result plus exact input/config provenance."""

    donor_id: str
    records: tuple[CellGeometryEligibility, ...]
    donor_median_cell_area_um2: float
    donor_max_cell_area_um2: float
    input_content_sha256: str
    config_content_sha256: str
    content_id: str

    def to_frame(self) -> pd.DataFrame:
        """Return a new mutable convenience view; records remain immutable."""
        rows = []
        for record in self.records:
            row = asdict(record)
            row["analysis_reason"] = record.analysis_reason
            row["estimation_reason"] = record.estimation_reason
            row["spillover_context_reason"] = record.spillover_context_reason
            row["analysis_reasons"] = ";".join(record.analysis_reasons)
            row["estimation_reasons"] = ";".join(record.estimation_reasons)
            row["spillover_context_reasons"] = ";".join(
                record.spillover_context_reasons
            )
            rows.append(row)
        return pd.DataFrame(rows)

    def verify_identity(self) -> None:
        payload = {
            "donor_id": self.donor_id,
            "records": [asdict(record) for record in self.records],
            "donor_median_cell_area_um2": self.donor_median_cell_area_um2,
            "donor_max_cell_area_um2": self.donor_max_cell_area_um2,
            "input_content_sha256": self.input_content_sha256,
            "config_content_sha256": self.config_content_sha256,
        }
        if sha256_json(payload) != self.content_id:
            raise ContractError("geometry-QC result content_id is invalid")


def _as_bool(series: pd.Series, column: str) -> np.ndarray:
    if series.isna().any():
        raise GeometryContractError(f"{column} contains null")
    values = series.to_numpy()
    if values.dtype == bool:
        return values.astype(bool, copy=False)
    accepted = {
        True: True,
        False: False,
        1: True,
        0: False,
        "true": True,
        "false": False,
        "True": True,
        "False": False,
    }
    try:
        return np.asarray([accepted[value] for value in values], dtype=bool)
    except (KeyError, TypeError) as exc:
        raise GeometryContractError(
            f"{column} must contain explicit booleans"
        ) from exc


def _numeric(frame: pd.DataFrame, column: str) -> np.ndarray:
    return pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)


def _connected_duplicate_groups(
    object_ids: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    area: np.ndarray,
    eligible: np.ndarray,
    *,
    radius_um: float,
    area_tolerance: float,
) -> tuple[tuple[int, ...], ...]:
    candidate_idx = np.flatnonzero(
        eligible & np.isfinite(x) & np.isfinite(y) & np.isfinite(area) & (area > 0)
    )
    if len(candidate_idx) < 2:
        return ()
    tree = cKDTree(np.column_stack((x[candidate_idx], y[candidate_idx])))
    local_pairs = tree.query_pairs(radius_um, output_type="set")

    parent = {int(index): int(index) for index in candidate_idx}

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for left_local, right_local in local_pairs:
        left = int(candidate_idx[left_local])
        right = int(candidate_idx[right_local])
        relative_delta = abs(area[left] - area[right]) / max(
            area[left], area[right], 1e-12
        )
        if relative_delta <= area_tolerance:
            union(left, right)

    components: dict[int, list[int]] = {}
    for index in candidate_idx:
        components.setdefault(find(int(index)), []).append(int(index))
    groups = [
        tuple(sorted(indices, key=lambda index: str(object_ids[index])))
        for indices in components.values()
        if len(indices) > 1
    ]
    return tuple(
        sorted(groups, key=lambda group: tuple(str(object_ids[index]) for index in group))
    )


def _duplicate_assignments(
    groups: Sequence[Sequence[int]],
    object_ids: np.ndarray,
    cell_area_px: np.ndarray,
) -> tuple[dict[int, str], set[int], set[int]]:
    group_by_index: dict[int, str] = {}
    survivors: set[int] = set()
    losers: set[int] = set()
    for group in groups:
        identifiers = sorted(str(object_ids[index]) for index in group)
        group_id = sha256_json(identifiers)[:16]

        def survivor_key(index: int) -> tuple[float, str]:
            pixels = cell_area_px[index]
            # A duplicate that won the raster mask retains the most pixels.
            raster_rank = -pixels if np.isfinite(pixels) else math.inf
            return raster_rank, str(object_ids[index])

        survivor = min(group, key=survivor_key)
        survivors.add(survivor)
        for index in group:
            group_by_index[index] = group_id
            if index != survivor:
                losers.add(index)
    return group_by_index, survivors, losers


def _ordered_reasons(flags: Mapping[str, np.ndarray], index: int) -> tuple[str, ...]:
    return tuple(
        reason
        for reason in ANALYSIS_REASON_ORDER
        if reason in flags and bool(flags[reason][index])
    )


def evaluate_geometry_metrics(
    metrics: pd.DataFrame,
    *,
    donor_id: str,
    config: GeometryQCConfig,
) -> GeometryQCResult:
    """Evaluate one donor's tabular geometry metrics without REDSEA.

    Duplicate geometries and missing required geometry are fatal by default.
    The explicit ``deterministic_keep`` and ``exclude`` policies are provided for
    controlled recovery of known historical exports; neither is silent.
    """
    donor_id = str(donor_id).strip()
    if not donor_id:
        raise ContractError("donor_id is required")
    missing_columns = set(GEOMETRY_METRIC_COLUMNS).difference(metrics.columns)
    if missing_columns:
        raise GeometryContractError(
            f"geometry metrics are missing columns {sorted(missing_columns)}"
        )
    if len(metrics) == 0:
        raise GeometryContractError(f"donor {donor_id}: geometry metric table is empty")

    frame = metrics.loc[:, GEOMETRY_METRIC_COLUMNS].copy()
    if frame["object_id"].isna().any():
        raise GeometryContractError("object_id contains null")
    frame["object_id"] = frame["object_id"].astype(str)
    if (frame["object_id"] == "").any():
        raise GeometryContractError("object_id contains an empty value")
    if not frame["object_id"].is_unique:
        duplicates = frame.loc[
            frame["object_id"].duplicated(keep=False), "object_id"
        ].unique()
        raise DuplicateObjectError(
            f"duplicate object_id values: {sorted(map(str, duplicates))[:5]}"
        )
    # A stable order makes every decision and content hash independent of input row order.
    frame = frame.sort_values("object_id", kind="mergesort").reset_index(drop=True)

    object_ids = frame["object_id"].to_numpy(dtype=object)
    cell_present = _as_bool(frame["cell_geometry_present"], "cell_geometry_present")
    nucleus_present = _as_bool(
        frame["nucleus_geometry_present"], "nucleus_geometry_present"
    )
    x = _numeric(frame, "centroid_x_um")
    y = _numeric(frame, "centroid_y_um")
    cell_area_um2 = _numeric(frame, "cell_area_um2")
    cell_area_px = _numeric(frame, "cell_area_px")
    nucleus_area_px = _numeric(frame, "nucleus_area_px")
    solidity = _numeric(frame, "cell_solidity")

    missing_cell = ~cell_present
    missing_nucleus = config.require_nucleus_geometry & ~nucleus_present
    if config.missing_geometry_policy == "fatal":
        if missing_cell.any():
            examples = object_ids[np.flatnonzero(missing_cell)[:5]].tolist()
            raise GeometryContractError(
                f"donor {donor_id}: {int(missing_cell.sum())} objects lack cell geometry "
                f"(examples={examples})"
            )
        if missing_nucleus.any():
            examples = object_ids[np.flatnonzero(missing_nucleus)[:5]].tolist()
            raise GeometryContractError(
                f"donor {donor_id}: {int(missing_nucleus.sum())} objects lack nucleus geometry "
                f"(examples={examples})"
            )

    valid_centroid = np.isfinite(x) & np.isfinite(y)
    valid_reported_area = np.isfinite(cell_area_um2) & (cell_area_um2 > 0)
    valid_raster_area = np.isfinite(cell_area_px) & (cell_area_px > 0)
    valid_solidity = np.isfinite(solidity) & (solidity >= 0) & (solidity <= 1)
    valid_nucleus_raster = np.isfinite(nucleus_area_px) & (nucleus_area_px > 0)

    duplicate_candidates = cell_present & valid_reported_area & valid_raster_area
    groups = _connected_duplicate_groups(
        object_ids,
        x,
        y,
        cell_area_um2,
        duplicate_candidates,
        radius_um=config.duplicate_radius_um,
        area_tolerance=config.duplicate_area_tolerance,
    )
    if groups and config.duplicate_policy == "fatal":
        examples = [
            [str(object_ids[index]) for index in group]
            for group in groups[:5]
        ]
        raise DuplicateGeometryError(
            f"donor {donor_id}: {len(groups)} duplicate geometry groups "
            f"(examples={examples})"
        )
    duplicate_groups, duplicate_survivors, duplicate_losers = _duplicate_assignments(
        groups, object_ids, cell_area_px
    )

    provisional = (
        cell_present
        & (~missing_nucleus)
        & valid_reported_area
        & (cell_area_um2 >= config.min_cell_area_um2)
        & valid_solidity
        & (solidity >= config.min_cell_solidity)
    )
    if duplicate_losers:
        provisional[list(duplicate_losers)] = False
    if not provisional.any():
        raise GeometryContractError(
            f"donor {donor_id}: no geometry survives the fixed gates needed "
            "to estimate the donor-relative area ceiling"
        )
    median_area = float(np.median(cell_area_um2[provisional]))
    max_area = config.max_cell_area_median_multiple * median_area

    expected_area_um2 = cell_area_px * config.pixel_area_um2
    with np.errstate(invalid="ignore", divide="ignore"):
        raster_ratio = expected_area_um2 / cell_area_um2
        nucleus_cell_ratio = nucleus_area_px / cell_area_px

    flags: dict[str, np.ndarray] = {
        "missing_cell_geometry": missing_cell,
        "missing_nucleus_geometry": missing_nucleus,
        "invalid_centroid": ~valid_centroid,
        "invalid_cell_area": ~valid_reported_area,
        "cell_area_below_min": valid_reported_area
        & (cell_area_um2 < config.min_cell_area_um2),
        "cell_area_above_max": valid_reported_area & (cell_area_um2 > max_area),
        "invalid_solidity": ~valid_solidity,
        "low_solidity": valid_solidity & (solidity < config.min_cell_solidity),
        "raster_empty": ~valid_raster_area,
        "raster_area_mismatch": (
            valid_reported_area
            & valid_raster_area
            & (
                (raster_ratio < config.min_raster_area_ratio)
                | (raster_ratio > config.max_raster_area_ratio)
            )
        ),
        "nucleus_raster_empty": config.require_nucleus_geometry
        & nucleus_present
        & ~valid_nucleus_raster,
        "duplicate_detection": np.asarray(
            [index in duplicate_losers for index in range(len(frame))], dtype=bool
        ),
    }

    records = []
    for index, object_id in enumerate(object_ids):
        analysis_reasons = _ordered_reasons(flags, index)
        analysis_eligible = not analysis_reasons

        estimation_reasons = list(analysis_reasons)
        no_cytoplasm = (
            analysis_eligible
            and np.isfinite(nucleus_cell_ratio[index])
            and nucleus_cell_ratio[index] >= config.no_cytoplasm_nc_ratio
        )
        if no_cytoplasm:
            estimation_reasons.append("no_cytoplasm")
        estimation_eligible = analysis_eligible and not no_cytoplasm

        spillover_reasons = []
        if not cell_present[index]:
            spillover_reasons.append("missing_cell_geometry")
        if not valid_raster_area[index]:
            spillover_reasons.append("raster_empty")
        if index in duplicate_losers:
            spillover_reasons.append("duplicate_detection")
        if (
            not config.retain_analysis_excluded_for_spillover
            and not analysis_eligible
            and not spillover_reasons
        ):
            spillover_reasons.append("analysis_ineligible")
        spillover_context_eligible = not spillover_reasons

        in_duplicate_group = index in duplicate_groups
        records.append(
            CellGeometryEligibility(
                donor_id=donor_id,
                object_id=str(object_id),
                analysis_eligible=bool(analysis_eligible),
                estimation_eligible=bool(estimation_eligible),
                spillover_context_eligible=bool(spillover_context_eligible),
                analysis_reasons=tuple(analysis_reasons),
                estimation_reasons=tuple(estimation_reasons),
                spillover_context_reasons=tuple(spillover_reasons),
                duplicate_group=duplicate_groups.get(index),
                duplicate_survivor=(
                    bool(index in duplicate_survivors)
                    if in_duplicate_group
                    else None
                ),
                raster_area_ratio=(
                    float(raster_ratio[index])
                    if np.isfinite(raster_ratio[index])
                    else None
                ),
                nucleus_cell_area_ratio=(
                    float(nucleus_cell_ratio[index])
                    if np.isfinite(nucleus_cell_ratio[index])
                    else None
                ),
            )
        )

    input_hash = dataframe_content_sha256(frame, sort_by=("object_id",))
    payload = {
        "donor_id": donor_id,
        "records": [asdict(record) for record in records],
        "donor_median_cell_area_um2": median_area,
        "donor_max_cell_area_um2": max_area,
        "input_content_sha256": input_hash,
        "config_content_sha256": config.content_id,
    }
    result = GeometryQCResult(
        donor_id=donor_id,
        records=tuple(records),
        donor_median_cell_area_um2=median_area,
        donor_max_cell_area_um2=max_area,
        input_content_sha256=input_hash,
        config_content_sha256=config.content_id,
        content_id=sha256_json(payload),
    )
    result.verify_identity()
    return result


__all__ = [
    "ANALYSIS_REASON_ORDER",
    "CellGeometryEligibility",
    "DuplicateGeometryError",
    "GEOMETRY_METRIC_COLUMNS",
    "GeometryQCConfig",
    "GeometryQCResult",
    "evaluate_geometry_metrics",
]
