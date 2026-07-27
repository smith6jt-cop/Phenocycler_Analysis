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


def test_donor_and_recipient_normalisation_differ_on_unequal_perimeters():
    """The published (donor) operator divides by the NEIGHBOUR's total contact, not the recipient's.

    Reference: redseapy builds ``comp*I - M[i,j]/B[j]`` and transposes it, then applies the result as
    ``transpose(dot(transpose(E), P))`` == ``P.T @ E`` == ``A @ E`` — so the transpose is undone and the
    net operator is the donor form. The two forms coincide only when every cell has the same total
    contact perimeter, which is why the symmetric two-cell fixtures above cannot tell them apart.
    """
    #   cell0 -10- cell1 -5- cell2   ->  B = [10, 15, 5]
    M = sp.csr_matrix(np.array([[0.0, 10.0, 0.0],
                                [10.0, 0.0, 5.0],
                                [0.0, 5.0, 0.0]]))
    edge = np.array([[100.0], [1000.0], [20.0]])
    data = np.full((3, 1), 1e6)          # large enough that nothing clips
    sizes = np.ones(3)

    donor, _, _ = redsea.compensate(data, edge, sizes, M, comp_mode=0, alpha=1.0, norm_form="donor")
    recip, _, _ = redsea.compensate(data, edge, sizes, M, comp_mode=0, alpha=1.0, norm_form="recipient")

    # donor: subtract sum_j M[i,j]/B[j] * edge[j]
    #   cell0: (10/15)*1000                      = 666.67
    #   cell1: (10/10)*100 + (5/5)*20            = 120
    #   cell2: (5/15)*1000                       = 333.33
    assert np.allclose((data - donor).ravel(), [2000.0 / 3, 120.0, 1000.0 / 3])
    # recipient: rows sum to 1 -> a contact-weighted AVERAGE of the neighbours' whole rim sums
    #   cell0: (10/10)*1000                      = 1000
    #   cell1: (10/15)*100 + (5/15)*20           = 73.33
    #   cell2: (5/5)*1000                        = 1000
    assert np.allclose((data - recip).ravel(), [1000.0, 220.0 / 3, 1000.0])

    # Mass conservation is the point: the donor form redistributes exactly each cell's own rim signal.
    assert np.isclose((data - donor).sum(), edge.sum())
    assert not np.isclose((data - recip).sum(), edge.sum())


def test_compensate_rejects_unknown_norm_form():
    M = sp.csr_matrix(np.array([[0.0, 4.0], [4.0, 0.0]]))
    with pytest.raises(ValueError, match="norm_form"):
        redsea.compensate(np.ones((2, 1)), np.ones((2, 1)), np.ones(2), M,
                          comp_mode=1, alpha=1.0, norm_form="rowsum")


def test_no_spillover_channels_pass_through_uncorrected():
    """Nuclear/ECM channels get alpha=0 AND comp_mode=0, so corrected == raw for those columns."""
    cols = ["E_cadherin", "DAPI", "CD31", "Collagen_IV", "Ki67"]
    params = redsea.RedseaParams(comp_mode=1, alpha=1.0, exclude_no_spillover=True)

    alpha, comp_mode, n_skipped = redsea.apply_no_spillover_mask(params, cols, params.alpha,
                                                                 params.comp_mode)
    assert n_skipped == 3                                    # DAPI, Collagen_IV, Ki67
    assert np.allclose(alpha, [1.0, 0.0, 1.0, 0.0, 0.0])
    assert np.allclose(comp_mode, [1, 0, 1, 0, 0])

    M = sp.csr_matrix(np.array([[0.0, 4.0], [4.0, 0.0]]))
    data = np.tile(np.array([[100.0], [10.0]]), (1, len(cols)))
    edge = np.tile(np.array([[20.0], [5.0]]), (1, len(cols)))
    sizes = np.array([10.0, 10.0])
    corrected, _, _ = redsea.compensate(data, edge, sizes, M, comp_mode, alpha)

    keep = np.array([c not in ("DAPI", "Collagen_IV", "Ki67") for c in cols])
    raw_mean = data / sizes[:, None]
    assert np.allclose(corrected[:, ~keep], raw_mean[:, ~keep])   # untouched
    assert not np.allclose(corrected[:, keep], raw_mean[:, keep])  # corrected


def test_no_spillover_mask_is_a_no_op_when_disabled():
    cols = ["E_cadherin", "DAPI"]
    params = redsea.RedseaParams(comp_mode=1, alpha=1.0, exclude_no_spillover=False)
    alpha, comp_mode, n_skipped = redsea.apply_no_spillover_mask(params, cols, 1.0, 1)
    assert (alpha, comp_mode, n_skipped) == (1.0, 1, 0)        # scalars preserved -> fast path kept


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
    assert p.comp_mode == 1            # subtract + reinforce (Bai et al. default, restored in v10)
    assert p.norm_form == "donor"      # published mass-conserving operator
    assert p.exclude_no_spillover is True
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
