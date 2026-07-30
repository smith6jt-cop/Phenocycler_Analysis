"""Build the one authoritative per-cell expression matrix.

Earlier workflow versions exposed several competing versions of the same
marker: raw QuPath whole-cell means, multiple corrected compartment means,
normalised values, and binary gates.  This module makes the current selection
explicit and registry-driven.  Each marker appears exactly once:

* ``redsea_corrected`` markers come from the marker's declared REDSEA
  compartment;
* ``uncompensated`` markers come from the declared raw QuPath compartment; and
* markers not acquired in an image are present as null, never silently treated
  as zero.

The resulting wide table is intentionally economical.  One row is one cell and
one column is one selected marker intensity.  Calibration evidence is a
separate wide table so the pipeline does not materialise tens of millions of
duplicated long-form identity rows merely to pivot them back again.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .artifacts import ContractError, DuplicateObjectError
from .marker_registry import MarkerRegistry, MarkerSpec, load_registry


METHOD_VERSION = "selected-expression-v1"
KEY_COLUMNS = ("donor_id", "object_id")
CONTEXT_COLUMNS = (
    "image",
    "X_centroid",
    "Y_centroid",
    "cell_region",
    "islet_num",
)


@dataclass(frozen=True)
class ExpressionAvailability:
    """The declared source and availability of one donor-marker measurement."""

    donor_id: str
    marker: str
    available: bool
    source: str
    compartment: str
    spillover_policy: str
    calibration_status: str
    reason: str


def _indexed(
    frame: pd.DataFrame,
    *,
    name: str,
    expected_ids: pd.Index | None = None,
) -> pd.DataFrame:
    if "object_id" not in frame:
        raise ContractError(f"{name} is missing object_id")
    if frame["object_id"].isna().any():
        raise ContractError(f"{name}.object_id contains null")
    result = frame.copy()
    result["object_id"] = result["object_id"].astype(str)
    if result["object_id"].duplicated().any():
        examples = result.loc[
            result["object_id"].duplicated(keep=False), "object_id"
        ].unique()[:5]
        raise DuplicateObjectError(
            f"{name} contains duplicate object_id values: {examples.tolist()}"
        )
    result = result.set_index("object_id", drop=False)
    if expected_ids is not None:
        missing = expected_ids.difference(result.index, sort=False)
        extra = result.index.difference(expected_ids, sort=False)
        if len(missing) or len(extra):
            raise ContractError(
                f"{name} object universe differs: "
                f"missing={len(missing):,}, extra={len(extra):,}"
            )
        result = result.reindex(expected_ids)
    return result


def _qc_columns(qc: pd.DataFrame) -> tuple[str, str]:
    analysis = next(
        (
            column
            for column in ("analysis_eligible", "qc_analysis_eligible", "qc_pass")
            if column in qc
        ),
        None,
    )
    estimation = next(
        (
            column
            for column in ("estimation_eligible", "qc_estimation_eligible")
            if column in qc
        ),
        None,
    )
    if analysis is None or estimation is None:
        raise ContractError(
            "geometry QC must provide analysis_eligible and estimation_eligible"
        )
    return analysis, estimation


def _raw_candidates(spec: MarkerSpec) -> tuple[str, ...]:
    selected = f"{spec.compartment}__{spec.name}"
    candidates = (
        selected,
        f"Raw__{selected}",
    )
    if spec.compartment == "Cell":
        candidates += (spec.name,)
    return candidates


def _resolve_measurement(
    *,
    spec: MarkerSpec,
    acquired: bool,
    cells: pd.DataFrame,
    redsea: pd.DataFrame,
) -> tuple[np.ndarray, str, str]:
    """Return values, source column and availability reason."""

    n = len(cells)
    if not acquired:
        return (
            np.full(n, np.nan, dtype=float),
            "",
            "marker_not_acquired_in_panel",
        )

    if spec.spillover_policy in {"redsea_corrected", "redsea_passthrough"}:
        source = spec.measurement_column
        if source not in redsea:
            raise ContractError(
                f"acquired marker {spec.name!r} requires REDSEA column {source!r}"
            )
        values = pd.to_numeric(redsea[source], errors="coerce").to_numpy(dtype=float)
        return values, f"redsea:{source}", ""

    if spec.spillover_policy == "uncompensated":
        source = next((column for column in _raw_candidates(spec) if column in cells), None)
        if source is None:
            raise ContractError(
                f"acquired marker {spec.name!r} requires raw QuPath "
                f"{spec.compartment} mean; tried {list(_raw_candidates(spec))}"
            )
        values = pd.to_numeric(cells[source], errors="coerce").to_numpy(dtype=float)
        return values, f"qupath:{source}", ""

    raise ContractError(
        f"{spec.name}: unsupported spillover policy {spec.spillover_policy!r}"
    )


def build_selected_expression(
    cells: pd.DataFrame,
    redsea_compartments: pd.DataFrame,
    geometry_qc: pd.DataFrame,
    *,
    donor_id: str,
    acquired_markers: Iterable[str],
    registry: MarkerRegistry | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Assemble one donor's selected expression and availability tables.

    ``acquired_markers`` must come from the versioned QuPath image contract,
    rather than being inferred from whichever columns happen to survive a
    partial stage.  If a declared acquired measurement is missing, the build is
    fatal.  A legitimately unacquired marker is represented by an all-null
    column and an explicit availability reason.
    """

    donor = str(donor_id).strip()
    if not donor:
        raise ContractError("donor_id is required")
    active_registry = load_registry() if registry is None else registry
    acquired = {str(marker) for marker in acquired_markers}
    unknown = sorted(acquired.difference(active_registry.by_name))
    if unknown:
        raise ContractError(
            f"QuPath panel contains marker(s) absent from the registry: {unknown}"
        )

    cell_indexed = _indexed(cells, name="cells")
    ids = cell_indexed.index
    redsea_indexed = _indexed(
        redsea_compartments, name="redsea_compartments", expected_ids=ids
    )
    qc_indexed = _indexed(geometry_qc, name="geometry_qc", expected_ids=ids)
    analysis_column, estimation_column = _qc_columns(qc_indexed)

    out = pd.DataFrame(
        {
            "donor_id": donor,
            "object_id": ids.to_numpy(dtype=str),
            "qc_analysis_eligible": qc_indexed[analysis_column]
            .astype(bool)
            .to_numpy(),
            "qc_estimation_eligible": qc_indexed[estimation_column]
            .astype(bool)
            .to_numpy(),
        }
    )
    for column in CONTEXT_COLUMNS:
        if column in cell_indexed:
            out[column] = cell_indexed[column].to_numpy()

    availability: list[ExpressionAvailability] = []
    for spec in active_registry.markers:
        is_acquired = spec.name in acquired
        values, source, reason = _resolve_measurement(
            spec=spec,
            acquired=is_acquired,
            cells=cell_indexed,
            redsea=redsea_indexed,
        )
        out[spec.name] = values
        availability.append(
            ExpressionAvailability(
                donor_id=donor,
                marker=spec.name,
                available=is_acquired,
                source=source,
                compartment=spec.compartment,
                spillover_policy=spec.spillover_policy,
                calibration_status=spec.calibration_status,
                reason=reason,
            )
        )

    availability_frame = pd.DataFrame(
        [record.__dict__ for record in availability]
    )
    return out, availability_frame


def write_donor_partition(
    frame: pd.DataFrame,
    root: Path,
    *,
    donor_id: str,
    filename: str = "data.parquet",
) -> Path:
    """Write one immutable donor partition.

    Production stages construct new content-addressed run directories.  An
    existing file therefore indicates an orchestration error and is never
    overwritten here.
    """

    donor = str(donor_id)
    partition = Path(root) / f"donor_id={donor}"
    path = partition / filename
    if path.exists():
        raise FileExistsError(f"immutable donor artifact already exists: {path}")
    partition.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return path


def read_single_partition(root: Path, donor_id: str) -> pd.DataFrame:
    """Read the only Parquet file in one donor partition."""

    partition = Path(root) / f"donor_id={donor_id}"
    files = sorted(partition.glob("*.parquet"))
    if len(files) != 1:
        raise ContractError(
            f"expected exactly one parquet for donor {donor_id} in {root}; "
            f"found {len(files)}"
        )
    return pd.read_parquet(files[0])


__all__ = [
    "CONTEXT_COLUMNS",
    "ExpressionAvailability",
    "KEY_COLUMNS",
    "METHOD_VERSION",
    "build_selected_expression",
    "read_single_partition",
    "write_donor_partition",
]
