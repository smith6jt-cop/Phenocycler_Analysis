"""Cell polygons from both platforms — the loader `match_by_iou` was missing.

The two formats disagree about layout, id and — the trap — units. QuPath GeoJSON is in
full-resolution qptiff **pixels**; Xenium `cell_boundaries.parquet` is in **microns**. Getting
that wrong does not error, it puts one section's polygons a thousand microns from their own
cells, so the conversion is fitted from the data and validated rather than assumed.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.contract import coerce_cell_table, partition_path, write_cell_table

pytest.importorskip("geopandas")
pytest.importorskip("shapely")

from phenocycler.integration.boundaries import (  # noqa: E402
    BoundaryError, fit_pixel_scale, load_pheno_polygons, load_xenium_polygons,
    read_geojson_polygons, warp_polygons,
)
from phenocycler.integration.sameslide import match_by_iou, match_tables_by_iou  # noqa: E402
from phenocycler.integration.transform import from_rigid  # noqa: E402

UM_PER_PX = 0.325          # a realistic PhenoCycler pixel


def _square(cx, cy, half):
    return [[[cx - half, cy - half], [cx + half, cy - half],
             [cx + half, cy + half], [cx - half, cy + half], [cx - half, cy - half]]]


def _geojson(tmp_path, centres_px, half_px=15.0, ids=None, name="cells__img.geojson"):
    ids = ids or [f"c{i}" for i in range(len(centres_px))]
    feats = [{"type": "Feature", "id": i,
              "geometry": {"type": "Polygon", "coordinates": _square(cx, cy, half_px)},
              "properties": {}}
             for i, (cx, cy) in zip(ids, centres_px)]
    p = tmp_path / name
    p.write_text(json.dumps({"type": "FeatureCollection", "features": feats}))
    return p


def _xenium_boundaries(tmp_path, centres_um, half_um=5.0, ids=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    ids = ids or [f"c{i}" for i in range(len(centres_um))]
    rows = []
    for cid, (cx, cy) in zip(ids, centres_um):
        for vx, vy in _square(cx, cy, half_um)[0][:-1]:
            rows.append({"cell_id": cid, "vertex_x": vx, "vertex_y": vy})
    p = tmp_path / "cell_boundaries.parquet"
    pd.DataFrame(rows).to_parquet(p, index=False)
    return p


def _table(ids):
    """A minimal CellTable — the re-indexing helper only ever reads `df['cell_id']`."""
    from phenocycler.integration.contract import CellTable
    return CellTable(df=pd.DataFrame({"cell_id": list(ids)}), modality="xenium",
                     donor_id="6539", roi="panc")


def _cell_table(cfg, base, centres_um, ids, modality, image="img"):
    df = pd.DataFrame({
        "cell_id": ids,
        "x_um": [c[0] for c in centres_um],
        "y_um": [c[1] for c in centres_um],
        "lineage_native": ["Endocrine"] * len(ids),
        "lineage_common": ["Endocrine"] * len(ids),
        "image": [image] * len(ids),
    })
    out = coerce_cell_table(df, modality=modality, donor_id="6539", roi="panc")
    out["image"] = image
    write_cell_table(out, partition_path(base, "6539", "panc"))
    return out


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #

def test_geojson_polygons_keep_their_ids(tmp_path):
    path = _geojson(tmp_path, [(100.0, 100.0), (200.0, 300.0)], ids=["a", "b"])
    polys = read_geojson_polygons(path)
    assert list(polys.index) == ["a", "b"]
    assert polys.iloc[0].area == pytest.approx(30.0 * 30.0)


def test_geojson_holes_are_kept(tmp_path):
    """A cell with a hole has less area than its outline, and IoU is an area ratio."""
    outer = _square(100.0, 100.0, 20.0)[0]
    inner = _square(100.0, 100.0, 5.0)[0]
    feat = {"type": "Feature", "id": "a",
            "geometry": {"type": "Polygon", "coordinates": [outer, inner]}, "properties": {}}
    p = tmp_path / "cells__img.geojson"
    p.write_text(json.dumps({"type": "FeatureCollection", "features": [feat]}))

    polys = read_geojson_polygons(p)
    assert polys.iloc[0].area == pytest.approx(40.0 * 40.0 - 10.0 * 10.0)


def test_xenium_boundaries_group_vertices_per_cell(tmp_path):
    """Upstream is long format — one row per vertex — and already in microns."""
    path = _xenium_boundaries(tmp_path, [(10.0, 10.0), (50.0, 50.0)], half_um=4.0)
    polys = load_xenium_polygons(path)
    assert len(polys) == 2
    assert polys.iloc[0].area == pytest.approx(8.0 * 8.0)


def test_missing_xenium_boundaries_is_an_actionable_error(tmp_path):
    with pytest.raises(BoundaryError, match="cell_boundaries"):
        load_xenium_polygons(tmp_path)


# --------------------------------------------------------------------------- #
# The units trap
# --------------------------------------------------------------------------- #

def test_pixel_scale_is_fitted_from_the_cells_themselves(tmp_path):
    """The GeoJSON is in pixels and the contract in microns; the same ids join them, so the
    conversion is recoverable exactly rather than read from qptiff tags that may be absent."""
    rng = np.random.default_rng(0)
    px = rng.uniform(100, 4000, (40, 2))
    ids = [f"c{i}" for i in range(40)]
    polys = read_geojson_polygons(_geojson(tmp_path, px, ids=ids))
    cells = pd.DataFrame({"cell_id": ids, "x_um": px[:, 0] * UM_PER_PX,
                          "y_um": px[:, 1] * UM_PER_PX})

    sx, bx, sy, by = fit_pixel_scale(polys, cells)
    assert sx == pytest.approx(UM_PER_PX, rel=1e-6)
    assert sy == pytest.approx(UM_PER_PX, rel=1e-6)
    assert abs(bx) < 1e-6 and abs(by) < 1e-6


def test_an_implausible_scale_is_refused(tmp_path):
    """If the fit says 40 um/px, the two files are not the same section. Failing loudly beats
    an overlay that is confidently a thousand microns out."""
    ids = [f"c{i}" for i in range(20)]
    px = np.random.default_rng(1).uniform(100, 500, (20, 2))
    polys = read_geojson_polygons(_geojson(tmp_path, px, ids=ids))
    cells = pd.DataFrame({"cell_id": ids, "x_um": px[:, 0] * 40.0, "y_um": px[:, 1] * 40.0})
    with pytest.raises(BoundaryError, match="outside"):
        fit_pixel_scale(polys, cells)


def test_non_corresponding_cells_are_refused(tmp_path):
    """Right scale, wrong cells: the residual catches what the scale check cannot."""
    rng = np.random.default_rng(2)
    ids = [f"c{i}" for i in range(30)]
    px = rng.uniform(100, 4000, (30, 2))
    polys = read_geojson_polygons(_geojson(tmp_path, px, ids=ids))
    scrambled = rng.permutation(px)
    cells = pd.DataFrame({"cell_id": ids, "x_um": scrambled[:, 0] * UM_PER_PX,
                          "y_um": scrambled[:, 1] * UM_PER_PX})
    with pytest.raises(BoundaryError, match="residual|outside"):
        fit_pixel_scale(polys, cells)


def test_too_few_shared_ids_is_refused(tmp_path):
    polys = read_geojson_polygons(_geojson(tmp_path, [(1.0, 1.0)], ids=["a"]))
    cells = pd.DataFrame({"cell_id": ["z"], "x_um": [1.0], "y_um": [1.0]})
    with pytest.raises(BoundaryError, match="in both"):
        fit_pixel_scale(polys, cells)


def test_pheno_polygons_land_in_microns_on_their_own_cells(tmp_path):
    """End to end: load, fit, convert — the polygon must contain its cell's µm centroid."""
    cfg = load_config()
    cfg.data_dir = tmp_path
    cfg.geojson_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(3)
    px = rng.uniform(200, 3000, (25, 2))
    ids = [f"c{i}" for i in range(25)]
    um = [(x * UM_PER_PX, y * UM_PER_PX) for x, y in px]
    _geojson(cfg.geojson_dir, px, ids=ids, name="cells__img.geojson")
    _cell_table(cfg, cfg.cells_pheno_dir, um, ids, "phenocycler", image="img")

    polys = load_pheno_polygons(cfg, "6539", "panc")
    assert len(polys) == 25
    from shapely.geometry import Point
    for cid, (x, y) in zip(ids, um):
        assert polys.loc[cid].contains(Point(x, y)), f"{cid} polygon does not contain its cell"


# --------------------------------------------------------------------------- #
# Warping and matching
# --------------------------------------------------------------------------- #

def test_warping_moves_polygons_through_the_transform(tmp_path):
    polys = load_xenium_polygons(_xenium_boundaries(tmp_path, [(10.0, 10.0)], half_um=3.0))
    tf = from_rigid(0.0, (100.0, -50.0))
    moved = warp_polygons(polys, tf)
    c = moved.iloc[0].centroid
    assert (c.x, c.y) == pytest.approx((110.0, -40.0), abs=1e-6)
    assert moved.iloc[0].area == pytest.approx(polys.iloc[0].area)   # rigid preserves area


def test_iou_separates_the_same_cell_from_its_neighbour(tmp_path):
    """The whole reason to prefer polygons.

    Two segmentations of one cell overlap; a neighbour 12 µm away does not — even though in
    packed tissue that neighbour's centroid can be nearer than the partner outline's.
    """
    a = load_xenium_polygons(_xenium_boundaries(tmp_path / "a", [(0.0, 0.0)], half_um=6.0,
                                                ids=["same"]))
    b_dir = tmp_path / "b"
    b = load_xenium_polygons(_xenium_boundaries(b_dir, [(1.0, 0.0), (14.0, 0.0)], half_um=6.0,
                                                ids=["same", "neighbour"]))

    pairs = match_by_iou(a, b, min_iou=0.25)
    assert len(pairs) == 1
    assert b.index[int(pairs.iloc[0]["fixed_idx"])] == "same"
    assert pairs.iloc[0]["iou"] > 0.7


def test_iou_rejects_a_pure_neighbour(tmp_path):
    a = load_xenium_polygons(_xenium_boundaries(tmp_path / "a", [(0.0, 0.0)], half_um=5.0))
    b = load_xenium_polygons(_xenium_boundaries(tmp_path / "b", [(14.0, 0.0)], half_um=5.0))
    assert match_by_iou(a, b, min_iou=0.25).empty


def test_iou_pairs_are_reindexed_onto_the_cell_tables(tmp_path):
    """`match_by_iou` returns positions into the GeoSeries it was handed; every caller then
    indexes the *cell tables* with them. Those orders are not the same — the polygon loaders
    group and filter — so a missing re-index would silently pair the wrong cells.
    """
    ids = ["c0", "c1", "c2"]
    centres = [(0.0, 0.0), (40.0, 0.0), (80.0, 0.0)]
    mv = load_xenium_polygons(_xenium_boundaries(tmp_path / "m", centres, half_um=6.0,
                                                 ids=ids))
    fx = load_xenium_polygons(_xenium_boundaries(tmp_path / "f", centres, half_um=6.0,
                                                 ids=ids))
    # Cell tables in a different row order, and carrying a cell the polygons do not.
    moving = _table(["c2", "c1", "c0", "ghost"])
    fixed = _table(["c1", "c0", "c2"])

    pairs = match_tables_by_iou(mv, fx, moving, fixed, min_iou=0.25)
    assert len(pairs) == 3
    for _, r in pairs.iterrows():
        assert moving.df["cell_id"].iloc[int(r["moving_idx"])] == \
               fixed.df["cell_id"].iloc[int(r["fixed_idx"])]


def test_polygons_without_a_cell_row_are_dropped_not_raised(tmp_path):
    """In `auto` mode a KeyError here would be swallowed by the centroid fallback and read as
    'no polygons available', which is a different and much more alarming diagnosis."""
    mv = load_xenium_polygons(_xenium_boundaries(tmp_path / "m", [(0.0, 0.0), (40.0, 0.0)],
                                                 half_um=6.0, ids=["c0", "orphan"]))
    fx = load_xenium_polygons(_xenium_boundaries(tmp_path / "f", [(0.0, 0.0), (40.0, 0.0)],
                                                 half_um=6.0, ids=["c0", "orphan"]))
    pairs = match_tables_by_iou(mv, fx, _table(["c0"]), _table(["c0"]), min_iou=0.25)
    assert len(pairs) == 1
    assert int(pairs.iloc[0]["moving_idx"]) == 0
