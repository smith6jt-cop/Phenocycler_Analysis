from __future__ import annotations

import numpy as np
import pandas as pd

from phenocycler.marker_registry import load_registry
from phenocycler.typing_stage import type_calibrated_cells


def test_typing_stage_keeps_all_cells_and_provenance():
    registry = load_registry()
    frame = pd.DataFrame(
        {
            "donor_id": ["D", "D"],
            "object_id": ["immune", "uncertain"],
            "qc_analysis_eligible": [True, False],
        }
    )
    for marker in registry.active_markers:
        frame[f"{marker.name}__state"] = "negative"
        frame[f"{marker.name}__expression_probability"] = 0.0
        frame[f"{marker.name}__model_valid"] = True
    frame.loc[0, "CD3e__state"] = "positive"
    frame.loc[0, "CD3e__expression_probability"] = 0.999
    frame.loc[1, "CD68__state"] = "indeterminate"
    frame.loc[1, "CD68__expression_probability"] = 0.7
    frame.loc[1, "CD68__model_valid"] = False

    out = type_calibrated_cells(frame, marker_registry=registry).set_index("object_id")
    assert set(out.index) == {"immune", "uncertain"}
    assert out.loc["immune", "broad_type"] == "Immune"
    assert out.loc["immune", "assignment_status"] == "anchor"
    assert out.loc["uncertain", "best_broad_type"]
    assert out.loc["uncertain", "assignment_status"] == "unavailable"
    assert out.loc["uncertain", "broad_reason"] == "geometry_qc_ineligible"
    assert not out.loc["uncertain", "qc_analysis_eligible"]
    assert out["typing_rules_fingerprint"].str.len().eq(64).all()
