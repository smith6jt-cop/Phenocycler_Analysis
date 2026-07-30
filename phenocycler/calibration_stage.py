"""Scalable donor calibration stage.

The estimator APIs operate on one marker at a time.  This runner keeps that
property and writes a wide per-cell result, avoiding the historical pattern of
materialising ``cells × markers`` long tables and immediately pivoting them
back to one row per cell.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .marker_calibration import (
    CalibrationConfig,
    DonorMarkerCalibration,
    apply_marker_calibration,
    fit_marker_calibration,
    models_to_long,
)
from .marker_registry import MarkerRegistry, load_registry
from .reference_controls import (
    ReferenceControlConfig,
    ReferenceControlResult,
    build_reference_controls,
    validate_reference_masks,
)


METHOD_VERSION = "calibration-stage-wide-v1"


@dataclass(frozen=True)
class CalibrationStageResult:
    donor_id: str
    reference_controls: ReferenceControlResult
    models: tuple[DonorMarkerCalibration, ...]
    model_table: pd.DataFrame
    evidence_wide: pd.DataFrame


def calibrate_expression(
    expression: pd.DataFrame,
    *,
    donor_id: str,
    registry: MarkerRegistry | None = None,
    calibration_config: CalibrationConfig | None = None,
    reference_config: ReferenceControlConfig | None = None,
    reference_controls: ReferenceControlResult | pd.DataFrame | None = None,
) -> CalibrationStageResult:
    """Select explicit reference seeds, fit each marker, and apply to all cells."""

    active_registry = load_registry() if registry is None else registry
    calibration_config = calibration_config or CalibrationConfig()
    reference_config = reference_config or ReferenceControlConfig(
        min_controls=calibration_config.min_controls,
        spatial_bins=calibration_config.spatial_bins,
    )
    required = {
        "object_id",
        "X_centroid",
        "Y_centroid",
        "qc_analysis_eligible",
        "qc_estimation_eligible",
    }
    missing = sorted(required.difference(expression.columns))
    if missing:
        raise ValueError(f"selected expression is missing columns: {missing}")
    missing_markers = sorted(
        marker.name
        for marker in active_registry.markers
        if marker.name not in expression
    )
    if missing_markers:
        raise ValueError(
            "selected expression must contain every registered marker "
            f"(null is allowed): {missing_markers}"
        )

    value_columns = {
        marker.name: marker.name for marker in active_registry.markers
    }
    if reference_controls is None:
        controls = build_reference_controls(
            expression,
            donor_id=str(donor_id),
            registry=active_registry,
            value_columns=value_columns,
            config=reference_config,
        )
    elif isinstance(reference_controls, ReferenceControlResult):
        controls = validate_reference_masks(
            reference_controls.masks_long,
            donor_id=str(donor_id),
            object_ids=expression["object_id"].to_numpy(),
            registry=active_registry,
            config=reference_config,
        )
    else:
        controls = validate_reference_masks(
            reference_controls,
            donor_id=str(donor_id),
            object_ids=expression["object_id"].to_numpy(),
            registry=active_registry,
            config=reference_config,
        )
    reference_arrays = controls.as_arrays(expression["object_id"].to_numpy())

    evidence_columns: dict[str, object] = {
        "donor_id": expression["donor_id"].to_numpy(),
        "object_id": expression["object_id"].to_numpy(),
        "qc_analysis_eligible": expression[
            "qc_analysis_eligible"
        ].to_numpy(dtype=bool),
    }
    models: list[DonorMarkerCalibration] = []
    for spec in active_registry.active_markers:
        model = fit_marker_calibration(
            expression[spec.name].to_numpy(),
            reference_arrays[spec.reference],
            expression["qc_estimation_eligible"].to_numpy(),
            expression["X_centroid"].to_numpy(),
            expression["Y_centroid"].to_numpy(),
            expression["object_id"].to_numpy(),
            donor_id=str(donor_id),
            marker=spec,
            registry=active_registry,
            config=calibration_config,
        )
        evidence = apply_marker_calibration(
            expression[spec.name].to_numpy(),
            model,
            object_ids=expression["object_id"].to_numpy(),
        )
        prefix = f"{spec.name}__"
        evidence_columns[prefix + "corrected_intensity"] = evidence[
            "measurement"
        ].to_numpy()
        evidence_columns[prefix + "state"] = pd.Categorical(
            evidence["state"],
            categories=("negative", "indeterminate", "positive", "unavailable"),
        )
        evidence_columns[prefix + "log2_threshold_ratio"] = evidence[
            "log2_threshold_ratio"
        ].to_numpy()
        evidence_columns[prefix + "empirical_tail_probability"] = evidence[
            "empirical_tail_probability"
        ].to_numpy()
        evidence_columns[prefix + "expression_probability"] = evidence[
            "expression_probability"
        ].to_numpy()
        evidence_columns[prefix + "empirical_tail_evidence"] = evidence[
            "empirical_tail_evidence"
        ].to_numpy()
        evidence_columns[prefix + "model_valid"] = (
            evidence["calibration_status"].eq("valid").to_numpy()
            & ~evidence["state"].isin(["indeterminate", "unavailable"]).to_numpy()
        )
        # Constant donor-marker model fields compress efficiently in Parquet
        # and keep the canonical evidence artifact self-describing even if a
        # convenience audit table is moved or regenerated.
        evidence_columns[prefix + "calibration_status"] = np.full(
            len(expression), model.status.value, dtype=object
        )
        evidence_columns[prefix + "reference"] = np.full(
            len(expression), model.reference, dtype=object
        )
        evidence_columns[prefix + "marker_alpha"] = np.full(
            len(expression), model.marker_alpha, dtype=float
        )
        evidence_columns[prefix + "threshold"] = np.full(
            len(expression), model.threshold_raw, dtype=float
        )
        evidence_columns[prefix + "threshold_ci_low"] = np.full(
            len(expression), model.threshold_ci_low_raw, dtype=float
        )
        evidence_columns[prefix + "threshold_ci_high"] = np.full(
            len(expression), model.threshold_ci_high_raw, dtype=float
        )
        evidence_columns[prefix + "n_controls"] = np.full(
            len(expression), model.n_controls_clean, dtype=np.int32
        )
        evidence_columns[prefix + "control_contamination_fraction"] = np.full(
            len(expression), model.contamination_fraction, dtype=float
        )
        models.append(model)

    evidence_wide = pd.DataFrame(evidence_columns)
    model_table = models_to_long(models)
    return CalibrationStageResult(
        donor_id=str(donor_id),
        reference_controls=controls,
        models=tuple(models),
        model_table=model_table,
        evidence_wide=evidence_wide,
    )


__all__ = [
    "CalibrationStageResult",
    "METHOD_VERSION",
    "calibrate_expression",
]
