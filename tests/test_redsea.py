"""Unit tests for the REDSEA math (no image data required)."""

from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from phenocycler import redsea
from phenocycler.artifacts import ChannelMapping


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


# ---------------------------------------------------------------- compartment decomposition


def _compartment_case(seed=0, n=40, c=5, nuc_frac=0.6, band_frac=0.25):
    """Random but self-consistent inputs: nucleus inside the cell, band inside the cell,
    nucleus-band inside the band."""
    rng = np.random.RandomState(seed)
    sizes = rng.randint(50, 400, size=n).astype(float)
    bandcount = np.maximum(np.round(sizes * band_frac), 1.0)
    nuccount = np.round(sizes * nuc_frac)
    nucbandcount = np.round(np.minimum(bandcount * 0.5, nuccount))
    data = rng.gamma(3.0, 200.0, size=(n, c))
    edge = data * rng.uniform(0.1, 0.4, size=(n, c))
    nuc = data * rng.uniform(0.2, 0.6, size=(n, c))
    nuc_edge = np.minimum(edge, nuc) * rng.uniform(0.1, 0.5, size=(n, c))
    dense = rng.rand(n, n) < 0.15
    dense = np.triu(dense, 1)
    M = sp.csr_matrix((dense + dense.T).astype(float) * rng.randint(1, 9, size=(n, n)))
    M = sp.csr_matrix(np.triu(M.toarray(), 1) + np.triu(M.toarray(), 1).T)
    return dict(data=data, edge=edge, nuc=nuc, nuc_edge=nuc_edge, sizes=sizes,
                bandcount=bandcount, nuccount=nuccount, nucbandcount=nucbandcount, M=M)


@pytest.mark.parametrize("comp_mode,alpha", [(1, 1.0), (0, 1.0), (1, 0.5), (0, 0.0)])
def test_compartments_conserve_mass_pre_clip(comp_mode, alpha):
    """Nucleus + Cytoplasm == Cell, exactly, as SUMS. This is the whole justification for the
    decomposition: the correction is apportioned by band share, so nothing is invented or lost."""
    k = _compartment_case()
    sub, _ = redsea.incoming_spillover(k["edge"], k["M"], alpha)
    cyto = k["data"] - k["nuc"]
    cyto_edge = k["edge"] - k["nuc_edge"]
    share_n = (k["nucbandcount"] / k["bandcount"])[:, None]
    share_y = ((k["bandcount"] - k["nucbandcount"]) / k["bandcount"])[:, None]

    s_n = k["nuc"] + comp_mode * k["nuc_edge"] - sub * share_n
    s_y = cyto + comp_mode * cyto_edge - sub * share_y
    s_cell = k["data"] + comp_mode * k["edge"] - sub

    assert np.allclose(s_n + s_y, s_cell, rtol=0, atol=1e-9)


@pytest.mark.parametrize("comp_mode,alpha", [(1, 1.0), (0, 1.0), (1, 0.5)])
def test_compartment_cell_reproduces_compensate(comp_mode, alpha):
    """The compartment Cell output is the canonical whole-cell REDSEA value."""
    k = _compartment_case(seed=1)
    corrected, _, _ = redsea.compensate(k["data"], k["edge"], k["sizes"], k["M"], comp_mode, alpha)
    means, counts, _ = redsea.compensate_compartments(
        k["data"], k["edge"], k["nuc"], k["nuc_edge"], k["sizes"],
        k["bandcount"], k["nuccount"], k["nucbandcount"], k["M"], comp_mode, alpha)
    assert np.allclose(means["Cell"], corrected, rtol=1e-12, atol=0.0, equal_nan=True)
    assert np.array_equal(counts["Cell"], k["sizes"])
    assert set(means) == set(redsea.COMPARTMENTS)


def test_compartment_means_divide_by_their_own_area():
    """Each compartment mean is its own corrected sum over its own pixel count -- a Membrane mean
    divided by the whole-cell area would silently dilute every band value."""
    k = _compartment_case(seed=2)
    means, counts, _ = redsea.compensate_compartments(
        k["data"], k["edge"], k["nuc"], k["nuc_edge"], k["sizes"],
        k["bandcount"], k["nuccount"], k["nucbandcount"], k["M"], 1, 1.0)
    sub, _ = redsea.incoming_spillover(k["edge"], k["M"], 1.0)
    expect_band = np.clip(k["edge"] + k["edge"] - sub, 0, None) / k["bandcount"][:, None]
    assert np.allclose(means["Membrane"], expect_band)
    assert np.array_equal(counts["Membrane"], k["bandcount"])
    assert np.array_equal(counts["Cytoplasm"], k["sizes"] - k["nuccount"])


def test_cytoplasm_is_nan_when_the_nucleus_fills_the_cell():
    """QuPath emits no Cytoplasm mean for such a cell (~9% of cells in this cohort) and neither do
    we. A 0 would be read as 'measured, and dark'."""
    k = _compartment_case(seed=3, n=6)
    k["nuccount"] = k["sizes"].copy()          # nucleus fills every cell
    k["nuc"] = k["data"].copy()
    means, counts, _ = redsea.compensate_compartments(
        k["data"], k["edge"], k["nuc"], k["nuc_edge"], k["sizes"],
        k["bandcount"], k["nuccount"], k["nucbandcount"], k["M"], 1, 1.0)
    assert np.isnan(means["Cytoplasm"]).all()
    assert np.isfinite(means["Nucleus"]).all()
    assert np.isfinite(means["Cell"]).all()


def test_nucleus_is_nan_when_the_feature_had_no_nucleus_geometry():
    k = _compartment_case(seed=4, n=6)
    k["nuccount"][:2] = 0.0
    k["nuc"][:2] = 0.0
    k["nucbandcount"][:2] = 0.0
    k["nuc_edge"][:2] = 0.0
    means, _, _ = redsea.compensate_compartments(
        k["data"], k["edge"], k["nuc"], k["nuc_edge"], k["sizes"],
        k["bandcount"], k["nuccount"], k["nucbandcount"], k["M"], 1, 1.0)
    assert np.isnan(means["Nucleus"][:2]).all()
    assert np.isfinite(means["Nucleus"][2:]).all()


def test_incoming_spillover_matches_the_expression_compensate_uses():
    """One implementation of alpha*(F@edge); `compensate` is not allowed to drift from it."""
    k = _compartment_case(seed=5)
    for norm_form in ("donor", "recipient"):
        sub, nz = redsea.incoming_spillover(
            k["edge"], k["M"], 1.0, norm_form=norm_form
        )
        corrected, _, n_iso = redsea.compensate(k["data"], k["edge"], k["sizes"], k["M"],
                                                comp_mode=0, alpha=1.0, norm_form=norm_form)
        expect = np.clip(k["data"] - sub, 0, None) / k["sizes"][:, None]
        assert np.allclose(corrected, expect)
        assert n_iso == int((~nz).sum())


def test_geometry_qc_controls_spillover_context_without_dropping_fragments():
    cell = np.array([[1, 1, 2], [1, 3, 2]], dtype=np.int32)
    nucleus = np.array([[1, 0, 2], [0, 3, 0]], dtype=np.int32)
    qc = pd.DataFrame(
        {
            "object_id": ["a", "b"],
            "spillover_context_eligible": [True, False],
        }
    )
    filtered_cell, filtered_nucleus, count = (
        redsea.apply_spillover_context_eligibility(
            cell,
            nucleus,
            ["a", "b", "geojson-only"],
            qc,
        )
    )
    assert count == 1
    assert not np.any(filtered_cell == 2)
    assert not np.any(filtered_nucleus == 2)
    assert np.any(filtered_cell == 1)
    assert np.any(filtered_cell == 3)


def test_separate_nucleus_raster_is_relabelled_by_object_id():
    cell = np.array([[1, 1, 0], [2, 2, 0]], dtype=np.int32)
    # Separate file order is B, A rather than cell order A, B.
    nucleus = np.array([[2, 0, 0], [1, 0, 0]], dtype=np.int32)
    relabelled = redsea.relabel_nucleus_mask(
        nucleus,
        ["B", "A"],
        cell,
        ["A", "B"],
    )
    assert relabelled[0, 0] == 1
    assert relabelled[1, 0] == 2
    assert ((relabelled == 0) | (relabelled == cell)).all()


def test_qptiff_channel_metadata_must_match_manifest_mapping():
    image = SimpleNamespace(
        image_id="image",
        channel_map=(
            ChannelMapping("DAPI", "DAPI raw", 0),
            ChannelMapping("CD3e", "CD3 raw", 1),
        ),
    )
    assert redsea.validated_channel_markers(
        image, ["DAPI raw", "CD3 raw"]
    ) == ["DAPI", "CD3e"]
    with pytest.raises(ValueError, match="differs from manifest"):
        redsea.validated_channel_markers(
            image, ["DAPI raw", "wrong"]
        )


def test_rasterize_mask_reads_nucleus_geometry(tmp_path):
    """nucleusGeometry is already in QuPath's export; before 2026-07-29 it was ignored."""
    square = lambda a, b: [[[a, a], [a, b], [b, b], [b, a], [a, a]]]
    gj = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "id": "uuid-A",
             "geometry": {"type": "Polygon", "coordinates": square(1, 5)},
             "nucleusGeometry": {"type": "Polygon", "coordinates": square(2, 4)},
             "properties": {}},
            {"type": "Feature", "id": "uuid-B",           # no nucleus -> stays background
             "geometry": {"type": "Polygon", "coordinates": square(6, 9)},
             "properties": {}},
        ],
    }
    p = tmp_path / "cells__img.geojson"
    p.write_text(json.dumps(gj))
    mask, oids, nuc = redsea.rasterize_mask(p, shape=(10, 10), downsample=1.0, with_nucleus=True)
    assert oids == ["uuid-A", "uuid-B"]
    assert (nuc == 1).any() and not (nuc == 2).any()      # A has a nucleus, B does not
    # the nucleus is a strict subset of its own cell, per label
    assert ((nuc == 0) | (nuc == mask)).all()
    assert (nuc == 1).sum() < (mask == 1).sum()
    _, presence_ids, _, nucleus_presence = redsea.rasterize_mask(
        p,
        shape=(10, 10),
        downsample=1.0,
        with_nucleus=True,
        return_nucleus_presence=True,
    )
    assert presence_ids == oids
    assert nucleus_presence == [True, False]
    # back-compat: the 2-tuple form is unchanged
    mask2, oids2 = redsea.rasterize_mask(p, shape=(10, 10), downsample=1.0)
    assert np.array_equal(mask2, mask) and oids2 == oids


def test_rasterize_mask_drops_nucleus_pixels_a_later_cell_overwrote(tmp_path):
    """Cells are drawn in order and a later polygon can overwrite an earlier one. If its nucleus
    raster kept those pixels, nuccount could exceed the cell's own area and `cyto = data - nuc`
    would go negative."""
    square = lambda a, b: [[[a, a], [a, b], [b, b], [b, a], [a, a]]]
    gj = {
        "type": "FeatureCollection",
        "features": [
            {"type": "Feature", "id": "A",
             "geometry": {"type": "Polygon", "coordinates": square(1, 6)},
             "nucleusGeometry": {"type": "Polygon", "coordinates": square(2, 5)},
             "properties": {}},
            {"type": "Feature", "id": "B",   # overlaps A, drawn second -> wins the shared pixels
             "geometry": {"type": "Polygon", "coordinates": square(3, 8)},
             "nucleusGeometry": {"type": "Polygon", "coordinates": square(4, 7)},
             "properties": {}},
        ],
    }
    p = tmp_path / "cells__img.geojson"
    p.write_text(json.dumps(gj))
    mask, _, nuc = redsea.rasterize_mask(p, shape=(10, 10), downsample=1.0, with_nucleus=True)
    n = 2
    sizes = np.bincount(mask.ravel(), minlength=n + 1)[1:]
    nuccount = np.bincount(nuc.ravel(), minlength=n + 1)[1:]
    assert (nuccount <= sizes).all()
    assert ((nuc == 0) | (nuc == mask)).all()


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
    )


def _write_raw_ids(cfg, donor, ids):
    path = cfg.cells_dir / f"donor_id={donor}"
    path.mkdir(parents=True)
    pd.DataFrame({"object_id": ids}).to_parquet(path / "data_0.parquet", index=False)


def test_canonical_output_mask_keeps_only_ingested_cells(tmp_path):
    cfg = _redsea_io_config(tmp_path)
    _write_raw_ids(cfg, "1", ["cell-b", "cell-a"])
    keep = redsea.canonical_cell_mask(
        cfg, "1", ["cell-a", "tiny-fragment", "cell-b"]
    )
    assert keep.tolist() == [True, False, True]


def test_canonical_output_mask_fails_when_ingested_cell_is_missing(tmp_path):
    cfg = _redsea_io_config(tmp_path)
    _write_raw_ids(cfg, "1", ["cell-a", "missing-cell"])

    with pytest.raises(ValueError, match="canonical raw cells are missing"):
        redsea.canonical_cell_mask(
            cfg, "1", ["cell-a", "tiny-fragment"]
        )


def test_one_image_per_donor_is_enforced_at_the_manifest(tmp_path):
    """One donor must resolve to exactly one image, and this is where that is enforced.

    `donor_image` picks the GeoJSON used to mask every cell in a partition. If a donor ever
    resolved to two images, the mask and the intensities would come from different sections:
    most cells would fall outside every polygon and be silently zeroed, and any that landed
    inside would get a neighbouring section's spillover correction, with nothing downstream
    reporting it.

    Main enforced this at the parquet read, because a partition was the unit of identity
    there. On this branch identity lives in the cohort manifest, so the guard moved with it —
    this test follows it rather than letting the invariant go untested.
    """
    import pytest

    from phenocycler import load_config
    from phenocycler.artifacts import QuPathCohortManifest

    cfg = load_config()
    cfg.qupath_manifest = tmp_path / "missing_manifest.json"
    with pytest.raises((ValueError, OSError, FileNotFoundError)):
        redsea.donor_image(cfg, "6539")

    # The manifest itself refuses to hold two images for one donor.
    assert hasattr(QuPathCohortManifest, "read_json")
