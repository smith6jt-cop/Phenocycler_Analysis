from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.expression import build_selected_expression
from phenocycler.marker_registry import load_registry


def _frames():
    cells = pd.DataFrame(
        {
            "object_id": ["c2", "c1", "c3"],
            "image": ["img", "img", "img"],
            "X_centroid": [2.0, 1.0, 3.0],
            "Y_centroid": [12.0, 11.0, 13.0],
            "Nucleus__Ki67": [4.0, 5.0, 6.0],
        }
    )
    redsea = pd.DataFrame(
        {
            "object_id": ["c1", "c2", "c3"],
            "Cell__CD3e": [10.0, 20.0, 30.0],
            "Cell__E_cadherin": [1.0, 2.0, 3.0],
        }
    )
    qc = pd.DataFrame(
        {
            "object_id": ["c3", "c1", "c2"],
            "analysis_eligible": [True, False, True],
            "estimation_eligible": [True, False, False],
        }
    )
    return cells, redsea, qc


def test_selected_expression_has_one_registry_selected_value_per_marker():
    cells, redsea, qc = _frames()
    registry = load_registry()
    expression, availability = build_selected_expression(
        cells,
        redsea,
        qc,
        donor_id="D1",
        acquired_markers={"CD3e", "E_cadherin", "Ki67"},
        registry=registry,
    )

    # Cell order is authoritative and every joined source is aligned to it.
    assert expression["object_id"].tolist() == ["c2", "c1", "c3"]
    assert expression["CD3e"].tolist() == [20.0, 10.0, 30.0]
    assert expression["Ki67"].tolist() == [4.0, 5.0, 6.0]
    assert expression["qc_analysis_eligible"].tolist() == [True, False, True]

    # A panel-absent marker is explicit null, not zero or a dropped column.
    assert expression["SST"].isna().all()
    sst = availability.set_index("marker").loc["SST"]
    assert not sst["available"]
    assert sst["reason"] == "marker_not_acquired_in_panel"

    ki67 = availability.set_index("marker").loc["Ki67"]
    assert ki67["source"] == "qupath:Nucleus__Ki67"
    assert ki67["spillover_policy"] == "uncompensated"


def test_declared_acquired_marker_missing_from_its_source_is_fatal():
    cells, redsea, qc = _frames()
    with pytest.raises(ContractError, match="Cell__SST"):
        build_selected_expression(
            cells,
            redsea,
            qc,
            donor_id="D1",
            acquired_markers={"CD3e", "E_cadherin", "Ki67", "SST"},
        )


def test_expression_rejects_partial_object_universes():
    cells, redsea, qc = _frames()
    with pytest.raises(ContractError, match="object universe differs"):
        build_selected_expression(
            cells,
            redsea.iloc[:-1],
            qc,
            donor_id="D1",
            acquired_markers={"CD3e", "E_cadherin", "Ki67"},
        )
