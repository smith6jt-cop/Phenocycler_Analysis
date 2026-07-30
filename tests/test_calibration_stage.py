from __future__ import annotations

import numpy as np
import pandas as pd

from phenocycler.calibration_stage import calibrate_expression
from phenocycler.marker_calibration import CalibrationConfig
from phenocycler.marker_registry import load_registry
from phenocycler.reference_controls import ReferenceControlConfig


def test_stage_preserves_every_cell_and_emits_wide_marker_evidence():
    registry = load_registry()
    n = 40
    expression = pd.DataFrame(
        {
            "donor_id": "D",
            "object_id": [f"c{i}" for i in range(n)],
            "X_centroid": np.arange(n, dtype=float),
            "Y_centroid": np.arange(n, dtype=float) % 5,
            "qc_analysis_eligible": True,
            "qc_estimation_eligible": True,
        }
    )
    for marker in registry.markers:
        expression[marker.name] = np.nan
    # Four references need enough finite-positive seed values. Active targets
    # use a quiet baseline with a few clear positives.
    for reference in registry.reference_graph.values():
        expression[reference] = np.linspace(1, 10, n)
    for marker in registry.active_markers:
        expression[marker.name] = np.r_[np.ones(n - 3), [20.0, 30.0, 40.0]]

    result = calibrate_expression(
        expression,
        donor_id="D",
        registry=registry,
        calibration_config=CalibrationConfig(
            control_cap=20,
            min_controls=10,
            spatial_bins=2,
            bootstrap_replicates=20,
        ),
        reference_config=ReferenceControlConfig(
            seed_cap=20,
            min_controls=10,
            spatial_bins=2,
        ),
    )
    assert result.evidence_wide["object_id"].tolist() == expression["object_id"].tolist()
    assert result.evidence_wide["qc_analysis_eligible"].all()
    assert len(result.model_table) == len(registry.active_markers)
    for marker in registry.active_markers:
        assert f"{marker.name}__corrected_intensity" in result.evidence_wide
        assert f"{marker.name}__state" in result.evidence_wide
        assert f"{marker.name}__expression_probability" in result.evidence_wide
        assert f"{marker.name}__empirical_tail_evidence" in result.evidence_wide
        probability = result.evidence_wide[
            f"{marker.name}__expression_probability"
        ]
        assert probability.dropna().between(0.0, 1.0).all()
