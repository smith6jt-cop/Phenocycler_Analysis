import numpy as np
import pandas as pd

from phenocycler.state_markers import annotate_proliferation, proliferation_threshold


def test_proliferation_is_an_orthogonal_state_annotation():
    cells = pd.DataFrame(
        {
            "object_id": [f"c{i}" for i in range(11)],
            "Ki67": [1.0] * 9 + [1000.0, 1000.0],
            "qc_analysis_eligible": [True] * 10 + [False],
            "qc_estimation_eligible": [True] * 10 + [False],
        }
    )
    out, models = annotate_proliferation(cells, donor_id="D", k=1.0, min_n=5)
    assert "cell_type" not in out
    assert out.loc[9, "Ki67__state"] == "positive"
    assert out.loc[10, "Ki67__state"] == "unavailable"
    assert not out.loc[10, "qc_analysis_eligible"]
    assert models.loc[0, "status"] == "valid"


def test_proliferation_missing_or_too_small_is_unavailable():
    assert np.isinf(proliferation_threshold([1, 2], min_n=3))
    out, models = annotate_proliferation(
        pd.DataFrame(
            {
                "object_id": ["a", "b"],
                "qc_analysis_eligible": True,
                "qc_estimation_eligible": True,
            }
        ),
        donor_id="D",
        markers=["PCNA"],
        min_n=1,
    )
    assert set(out["PCNA__state"]) == {"unavailable"}
    assert models.loc[0, "status"] == "marker_unavailable"
