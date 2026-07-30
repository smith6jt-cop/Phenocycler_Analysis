"""Raster geometry stage that runs before REDSEA parameter estimation.

Only geometry is read here: the qptiff metadata supplies raster dimensions and
pixel calibration, while the QuPath GeoJSON supplies cell/nucleus polygons.
No marker channel and no REDSEA output is consulted.  This removes the former
cycle in which cell QC depended on a REDSEA raster even though REDSEA was
supposed to run only after QC.
"""

from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from .artifacts import ContractError, QuPathImageManifest
from .geometry_qc import GeometryQCConfig, GeometryQCResult, evaluate_geometry_metrics
from .redsea import qptiff_channels, rasterize_mask


METHOD_VERSION = "geometry-raster-qc-v1"
GEOMETRY_SOURCE_COLUMNS = (
    "object_id",
    "X_centroid",
    "Y_centroid",
    "cell_area",
    "cell_solidity",
)


def _counts_for_canonical_ids(
    label_image: np.ndarray,
    geometry_object_ids: Sequence[object],
    canonical_object_ids: Sequence[object],
) -> np.ndarray:
    labels = np.bincount(
        np.asarray(label_image).ravel(),
        minlength=len(geometry_object_ids) + 1,
    )
    positions = {str(object_id): index + 1 for index, object_id in enumerate(geometry_object_ids)}
    if len(positions) != len(geometry_object_ids):
        raise ContractError("geometry object IDs are not unique")
    missing = [
        str(object_id)
        for object_id in canonical_object_ids
        if str(object_id) not in positions
    ]
    if missing:
        raise ContractError(
            f"{len(missing):,} canonical cells lack geometry "
            f"(examples={missing[:5]})"
        )
    return np.asarray(
        [labels[positions[str(object_id)]] for object_id in canonical_object_ids],
        dtype=float,
    )


def _presence_for_canonical_ids(
    geometry_object_ids: Sequence[object],
    geometry_presence: Sequence[object],
    canonical_object_ids: Sequence[object],
) -> np.ndarray:
    if len(geometry_presence) != len(geometry_object_ids):
        raise ContractError(
            "nucleus geometry presence must align with geometry object IDs"
        )
    positions = {
        str(object_id): index
        for index, object_id in enumerate(geometry_object_ids)
    }
    if len(positions) != len(geometry_object_ids):
        raise ContractError("geometry object IDs are not unique")
    missing = [
        str(object_id)
        for object_id in canonical_object_ids
        if str(object_id) not in positions
    ]
    if missing:
        raise ContractError(
            f"{len(missing):,} canonical cells lack geometry "
            f"(examples={missing[:5]})"
        )
    return np.asarray(
        [
            bool(geometry_presence[positions[str(object_id)]])
            for object_id in canonical_object_ids
        ],
        dtype=bool,
    )


def geometry_metrics_from_masks(
    cells: pd.DataFrame,
    *,
    cell_mask: np.ndarray,
    cell_geometry_object_ids: Sequence[object],
    nucleus_mask: np.ndarray,
    nucleus_geometry_object_ids: Sequence[object],
    nucleus_geometry_presence: Sequence[object] | None = None,
) -> pd.DataFrame:
    """Create the exact input table consumed by :func:`evaluate_geometry_metrics`."""

    required = set(GEOMETRY_SOURCE_COLUMNS)
    missing = sorted(required.difference(cells.columns))
    if missing:
        raise ContractError(f"canonical cells are missing geometry columns: {missing}")
    if cells["object_id"].isna().any() or cells["object_id"].astype(str).duplicated().any():
        raise ContractError("canonical cell object_id values must be non-null and unique")
    ids = cells["object_id"].astype(str).to_numpy()
    if nucleus_geometry_presence is None:
        nucleus_geometry_presence = [True] * len(nucleus_geometry_object_ids)
    return pd.DataFrame(
        {
            "object_id": ids,
            "cell_geometry_present": True,
            "nucleus_geometry_present": _presence_for_canonical_ids(
                nucleus_geometry_object_ids,
                nucleus_geometry_presence,
                ids,
            ),
            "centroid_x_um": pd.to_numeric(cells["X_centroid"], errors="coerce"),
            "centroid_y_um": pd.to_numeric(cells["Y_centroid"], errors="coerce"),
            "cell_area_um2": pd.to_numeric(cells["cell_area"], errors="coerce"),
            "cell_area_px": _counts_for_canonical_ids(
                cell_mask, cell_geometry_object_ids, ids
            ),
            "nucleus_area_px": _counts_for_canonical_ids(
                nucleus_mask, nucleus_geometry_object_ids, ids
            ),
            "cell_solidity": pd.to_numeric(cells["cell_solidity"], errors="coerce"),
        }
    )


def _separate_nucleus_mask(
    nucleus_mask: np.ndarray,
    nucleus_ids: Sequence[object],
    cell_mask: np.ndarray,
    cell_ids: Sequence[object],
) -> np.ndarray:
    """Relabel a separate nucleus GeoJSON into the cell mask's label space."""

    cell_labels = {str(object_id): index + 1 for index, object_id in enumerate(cell_ids)}
    lookup = np.zeros(len(nucleus_ids) + 1, dtype=cell_mask.dtype)
    for index, object_id in enumerate(nucleus_ids, start=1):
        try:
            lookup[index] = cell_labels[str(object_id)]
        except KeyError as exc:
            raise ContractError(
                f"nucleus object {object_id!r} has no cell geometry"
            ) from exc
    relabelled = lookup[np.asarray(nucleus_mask)]
    # A nucleus polygon is usable only inside its own cell polygon.
    return np.where(relabelled == cell_mask, relabelled, 0)


def run_geometry_qc_for_image(
    cells: pd.DataFrame,
    image: QuPathImageManifest,
    *,
    config: GeometryQCConfig | None = None,
) -> GeometryQCResult:
    """Rasterize one contracted image and evaluate immutable QC eligibility."""

    image.validate_current(mode="fast")
    tf, _pages, _channels, shape, pixel_size = qptiff_channels(
        Path(image.qptiff.path), 0
    )
    tf.close()
    if not np.isclose(pixel_size, image.pixel_size_um_x, rtol=1e-6, atol=1e-9):
        raise ContractError(
            f"{image.image_id}: qptiff pixel size {pixel_size} differs from "
            f"manifest x calibration {image.pixel_size_um_x}"
        )
    if not np.isclose(
        image.pixel_size_um_x, image.pixel_size_um_y, rtol=1e-6, atol=1e-9
    ):
        raise ContractError(
            f"{image.image_id}: REDSEA currently requires isotropic pixels; "
            f"manifest has {image.pixel_size_um_x} x {image.pixel_size_um_y} um"
        )

    if image.combined_geometry_file:
        cell_mask, cell_ids, nucleus_mask, nucleus_presence = rasterize_mask(
            Path(image.cell_geojson.path),
            shape,
            1.0,
            with_nucleus=True,
            return_nucleus_presence=True,
        )
        nucleus_ids = cell_ids
    else:
        cell_mask, cell_ids = rasterize_mask(
            Path(image.cell_geojson.path), shape, 1.0
        )
        raw_nucleus_mask, nucleus_ids = rasterize_mask(
            Path(image.nucleus_geojson.path), shape, 1.0
        )
        nucleus_mask = _separate_nucleus_mask(
            raw_nucleus_mask, nucleus_ids, cell_mask, cell_ids
        )
        nucleus_presence = [True] * len(nucleus_ids)

    metrics = geometry_metrics_from_masks(
        cells,
        cell_mask=cell_mask,
        cell_geometry_object_ids=cell_ids,
        nucleus_mask=nucleus_mask,
        nucleus_geometry_object_ids=nucleus_ids,
        nucleus_geometry_presence=nucleus_presence,
    )
    qc_config = config or GeometryQCConfig(
        pixel_size_um_x=image.pixel_size_um_x,
        pixel_size_um_y=image.pixel_size_um_y,
    )
    return evaluate_geometry_metrics(
        metrics,
        donor_id=image.donor_id,
        config=qc_config,
    )


__all__ = [
    "GEOMETRY_SOURCE_COLUMNS",
    "METHOD_VERSION",
    "geometry_metrics_from_masks",
    "run_geometry_qc_for_image",
]
