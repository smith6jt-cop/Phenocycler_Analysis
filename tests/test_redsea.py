"""Unit tests for the REDSEA math (no image data required)."""

from __future__ import annotations

import json

import numpy as np
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


def test_donor_image_rejects_a_partition_holding_two_sections(tmp_path):
    """One partition must be one section, and this is where that gets enforced.

    `donor_image` picks the GeoJSON used to mask every cell in a partition. If two images
    ever share a partition — the exact failure the section key exists to prevent — taking
    `.iloc[0]` would mask one section's cells with the other section's polygons: most would
    fall outside every polygon and be silently zeroed, and any that landed inside would get a
    neighbouring section's spillover correction. No stage downstream would report anything,
    so the guard belongs at the read.
    """
    import pandas as pd
    import pytest

    from phenocycler import load_config

    cfg = load_config()
    cfg.data_dir = tmp_path
    part = cfg.cells_dir / "donor_id=6539"
    part.mkdir(parents=True)

    one = pd.DataFrame({"image": ["6539_Scan1.er.qptiff - resolution #1"] * 3})
    one.to_parquet(part / "data_0.parquet", index=False)
    assert redsea.donor_image(cfg, "6539").startswith("6539_Scan1")

    mixed = pd.DataFrame({"image": ["6539_Scan1.er.qptiff - resolution #1",
                                    "6539pLN_Scan1.er.qptiff - resolution #1"]})
    mixed.to_parquet(part / "data_0.parquet", index=False)
    with pytest.raises(SystemExit, match="contains 2 images"):
        redsea.donor_image(cfg, "6539")
