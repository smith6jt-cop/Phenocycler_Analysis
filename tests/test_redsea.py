"""Unit tests for the REDSEA math (no image data required)."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from phenocycler import redsea


def test_compensate_subtract_only_analytic():
    """Two mutually-contacting cells, one channel: hand-checked subtract-only result.

    F = row-normalized contact = [[0,1],[1,0]]; F@edge swaps the edge terms.
    corrected_sum = data - F@edge = [100-5, 10-20] -> clip -> [95, 0]; /size(=10).
    """
    M = sp.csr_matrix(np.array([[0.0, 4.0], [4.0, 0.0]]))
    data = np.array([[100.0], [10.0]])
    edge = np.array([[20.0], [5.0]])
    sizes = np.array([10.0, 10.0])

    corrected, clamp_frac, n_iso = redsea.compensate(data, edge, sizes, M, comp_mode=0, alpha=1.0)

    assert np.allclose(corrected, [[9.5], [0.0]])
    assert np.allclose(clamp_frac, [0.5])       # one of two cells clamped to 0
    assert n_iso == 0


def test_compensate_reinforce_mode():
    """comp_mode=1 adds the cell's own boundary term back (subtract + reinforce)."""
    M = sp.csr_matrix(np.array([[0.0, 4.0], [4.0, 0.0]]))
    data = np.array([[100.0], [10.0]])
    edge = np.array([[20.0], [5.0]])
    sizes = np.array([10.0, 10.0])

    corrected, _, _ = redsea.compensate(data, edge, sizes, M, comp_mode=1, alpha=1.0)
    # corrected_sum = data + edge - F@edge = [100+20-5, 10+5-20] -> [115, 0] -> /10
    assert np.allclose(corrected, [[11.5], [0.0]])


def test_compensate_alpha_scales_subtraction():
    M = sp.csr_matrix(np.array([[0.0, 4.0], [4.0, 0.0]]))
    data = np.array([[100.0], [100.0]])
    edge = np.array([[40.0], [0.0]])
    sizes = np.array([10.0, 10.0])
    # F@edge = [edge[1], edge[0]] = [0, 40]; alpha=0.5 -> subtract [0, 20]
    corrected, _, _ = redsea.compensate(data, edge, sizes, M, comp_mode=0, alpha=0.5)
    assert np.allclose(corrected, [[10.0], [8.0]])


def test_compensate_isolated_cell_passthrough():
    """A cell with no contacts (row sum 0) is left unchanged and counted isolated."""
    M = sp.csr_matrix(np.array([[0.0, 4.0, 0.0],
                                [4.0, 0.0, 0.0],
                                [0.0, 0.0, 0.0]]))
    data = np.array([[100.0], [10.0], [77.0]])
    edge = np.array([[20.0], [5.0], [50.0]])
    sizes = np.array([10.0, 10.0, 10.0])
    corrected, _, n_iso = redsea.compensate(data, edge, sizes, M, comp_mode=0, alpha=1.0)
    assert n_iso == 1
    assert np.isclose(corrected[2, 0], 7.7)     # unchanged: 77/10


def test_contact_matrix_counts_8_adjacencies():
    """Two 2x2 label blocks side by side: 4 shared-boundary pixel adjacencies."""
    mask = np.array([[1, 1, 2, 2],
                     [1, 1, 2, 2]], dtype=np.int32)
    M = redsea.contact_matrix(mask, n=2)
    dense = M.toarray()
    assert dense[0, 1] == dense[1, 0]
    assert dense[0, 1] == 4      # 2 horizontal + 2 diagonal adjacencies
    assert dense[0, 0] == 0 and dense[1, 1] == 0


def test_contact_matrix_symmetric_and_empty():
    mask = np.array([[1, 0, 0, 2]], dtype=np.int32)   # not adjacent
    M = redsea.contact_matrix(mask, n=2)
    assert M.nnz == 0


def test_rasterize_mask_fills_polygon_with_uuid(tmp_path):
    """A single square polygon rasterizes to label 1 and keeps its UUID."""
    gj = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "id": "uuid-A",
            "geometry": {"type": "Polygon",
                         "coordinates": [[[1, 1], [1, 3], [3, 3], [3, 1], [1, 1]]]},
            "properties": {},
        }],
    }
    p = tmp_path / "cells__img.geojson"
    p.write_text(json.dumps(gj))
    mask, oids = redsea.rasterize_mask(p, shape=(5, 5), downsample=1.0)
    assert oids == ["uuid-A"]
    assert (mask == 1).any()            # polygon filled
    assert mask.max() == 1              # only one label
    # background preserved at the far corner
    assert mask[0, 0] == 0


def test_redsea_params_from_config_preserves_defaults():
    from phenocycler import load_config
    cfg = load_config()
    p = redsea.RedseaParams.from_config(cfg)
    assert p.comp_mode == 0            # subtract-only
    assert p.alpha == 1.0
    assert p.edge_radius == 0          # 1-px band
    assert p.gap_bridge == 1


def _redsea_io_config(tmp_path):
    return SimpleNamespace(
        cells_dir=tmp_path / "cells",
        cells_redsea_dir=tmp_path / "cells_redsea",
    )


def _write_raw_ids(cfg, donor, ids):
    path = cfg.cells_dir / f"donor_id={donor}"
    path.mkdir(parents=True)
    pd.DataFrame({"object_id": ids}).to_parquet(path / "data_0.parquet", index=False)


def test_write_corrected_persists_only_canonical_raw_cells(tmp_path):
    cfg = _redsea_io_config(tmp_path)
    _write_raw_ids(cfg, "1", ["cell-b", "cell-a"])

    redsea.write_corrected(
        cfg,
        "1",
        ["cell-a", "tiny-fragment", "cell-b"],
        ["CD3e"],
        np.array([[10.0], [20.0], [30.0]]),
        np.array([5.0, 1.0, 6.0]),
    )

    out = pd.read_parquet(cfg.cells_redsea_dir / "donor_id=1" / "data_0.parquet")
    assert out["object_id"].tolist() == ["cell-a", "cell-b"]
    assert out["CD3e"].tolist() == [10.0, 30.0]
    assert out["cell_area_px"].tolist() == [5.0, 6.0]


def test_write_corrected_fails_when_canonical_cell_is_missing(tmp_path):
    cfg = _redsea_io_config(tmp_path)
    _write_raw_ids(cfg, "1", ["cell-a", "missing-cell"])

    with pytest.raises(ValueError, match="canonical raw cells are missing"):
        redsea.write_corrected(
            cfg,
            "1",
            ["cell-a", "tiny-fragment"],
            ["CD3e"],
            np.array([[10.0], [20.0]]),
            np.array([5.0, 1.0]),
        )


def test_reconcile_existing_redsea_is_exact_and_schema_preserving(tmp_path):
    cfg = _redsea_io_config(tmp_path)
    _write_raw_ids(cfg, "1", ["cell-a", "cell-b"])
    path = cfg.cells_redsea_dir / "donor_id=1"
    path.mkdir(parents=True)
    before = pd.DataFrame(
        {
            "object_id": ["cell-a", "tiny-fragment", "cell-b"],
            "donor_id": ["1", "1", "1"],
            "cell_area_px": np.array([5.0, 1.0, 6.0], dtype=np.float64),
            "CD3e": np.array([10.0, 20.0, 30.0], dtype=np.float64),
        }
    )
    before.to_parquet(path / "data_0.parquet", index=False)

    redsea.reconcile_existing_donor("1", cfg)

    out = pd.read_parquet(path / "data_0.parquet")
    assert out["object_id"].tolist() == ["cell-a", "cell-b"]
    assert out["CD3e"].tolist() == [10.0, 30.0]
    assert out.dtypes.to_dict() == before.dtypes.to_dict()
