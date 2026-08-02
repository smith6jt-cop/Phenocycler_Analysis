from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import ContractError
from phenocycler.geometry_stage import geometry_metrics_from_masks


def test_geometry_metrics_align_masks_to_canonical_object_order():
    cells = pd.DataFrame(
        {
            "object_id": ["b", "a"],
            "X_centroid": [2.0, 1.0],
            "Y_centroid": [4.0, 3.0],
            "cell_area": [2.0, 3.0],
            "cell_solidity": [0.9, 0.8],
        }
    )
    cell = np.array([[1, 1, 0], [2, 2, 2]], dtype=np.int32)
    nucleus = np.array([[1, 0, 0], [0, 2, 0]], dtype=np.int32)
    result = geometry_metrics_from_masks(
        cells,
        cell_mask=cell,
        cell_geometry_object_ids=["a", "b"],
        nucleus_mask=nucleus,
        nucleus_geometry_object_ids=["a", "b"],
    )
    assert result["object_id"].tolist() == ["b", "a"]
    assert result["cell_area_px"].tolist() == [3.0, 2.0]
    assert result["nucleus_area_px"].tolist() == [1.0, 1.0]


def test_geometry_metrics_preserve_embedded_nucleus_presence_by_object_id():
    cells = pd.DataFrame(
        {
            "object_id": ["b", "a"],
            "X_centroid": [2.0, 1.0],
            "Y_centroid": [4.0, 3.0],
            "cell_area": [2.0, 3.0],
            "cell_solidity": [0.9, 0.8],
        }
    )
    result = geometry_metrics_from_masks(
        cells,
        cell_mask=np.array([[1, 1, 0], [2, 2, 2]], dtype=np.int32),
        cell_geometry_object_ids=["a", "b"],
        nucleus_mask=np.array([[1, 0, 0], [0, 0, 0]], dtype=np.int32),
        nucleus_geometry_object_ids=["a", "b"],
        nucleus_geometry_presence=[True, False],
    )

    assert result["object_id"].tolist() == ["b", "a"]
    assert result["nucleus_geometry_present"].tolist() == [False, True]
    assert result["nucleus_area_px"].tolist() == [0.0, 1.0]


def test_geometry_metrics_reject_missing_canonical_geometry():
    cells = pd.DataFrame(
        {
            "object_id": ["missing"],
            "X_centroid": [1.0],
            "Y_centroid": [1.0],
            "cell_area": [2.0],
            "cell_solidity": [0.9],
        }
    )
    with pytest.raises(ContractError, match="lack geometry"):
        geometry_metrics_from_masks(
            cells,
            cell_mask=np.array([[1]], dtype=np.int32),
            cell_geometry_object_ids=["other"],
            nucleus_mask=np.array([[1]], dtype=np.int32),
            nucleus_geometry_object_ids=["other"],
        )
