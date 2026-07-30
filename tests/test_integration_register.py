"""Registration: can a known transform be recovered, and are the conventions right?

The headline tests apply a *known* rigid / affine / non-rigid warp to synthetic data and
assert that registration recovers it. Everything else in the spatial pipeline is downstream
of that being true, and several of the failure modes here are silent — a wrong displacement
sign doubles the misalignment instead of removing it, and a reflection the fitter cannot
represent comes back as a confident, wrong rotation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.rasterize import (GridSpec, grid_for, normalize_image,
                                               rasterize_points)
from phenocycler.integration.register import (dice, icp, nearest_neighbour_rmse,
                                              ransac_affine, refine_bspline, refine_ecc,
                                              register_image, register_point_set)
from phenocycler.integration.transform import (Transform, affine_from_pixels, affine_to_pixels,
                                               fit_affine_lstsq, fit_rigid_umeyama, from_rigid)


def _islets(n, seed, lo=300, hi=3700):
    r = np.random.default_rng(seed)
    xy = r.uniform(lo, hi, (n, 2))
    return pd.DataFrame({
        "structure_id": [f"s{i}" for i in range(n)],
        "x_um": xy[:, 0], "y_um": xy[:, 1],
        "n_cells": r.integers(25, 180, n).astype(float),
        "area_um2": r.uniform(3e3, 1.2e4, n),
        "eccentricity": r.uniform(0.1, 0.6, n),
    })


# --------------------------------------------------------------------------- #
# Transform algebra
# --------------------------------------------------------------------------- #

def test_umeyama_recovers_rotation_and_translation():
    rng = np.random.default_rng(0)
    truth = from_rigid(17.0, (250.0, -120.0), center_um=(2000, 2000))
    src = rng.uniform(0, 4000, (60, 2))
    est = fit_rigid_umeyama(src, truth.apply(src))
    assert est.rotation_deg == pytest.approx(17.0, abs=1e-6)
    assert est.scale == pytest.approx(1.0, abs=1e-9)
    assert np.abs(est.apply(src) - truth.apply(src)).max() < 1e-8


def test_umeyama_reflection_must_be_representable():
    """Serial sections are sometimes mounted face-down.

    Textbook Umeyama forces det=+1, which cannot represent that and silently returns a wrong
    proper rotation instead. `allow_reflection` is what makes flipped donors registrable.
    """
    rng = np.random.default_rng(1)
    truth = from_rigid(-73.0, (90.0, 410.0), mirror=True, center_um=(2000, 2000))
    src = rng.uniform(0, 4000, (40, 2))
    dst = truth.apply(src)

    forced = fit_rigid_umeyama(src, dst, allow_reflection=False)
    assert not forced.is_mirrored
    assert np.abs(forced.apply(src) - dst).max() > 100.0      # genuinely wrong

    ok = fit_rigid_umeyama(src, dst, allow_reflection=True)
    assert ok.is_mirrored
    assert np.abs(ok.apply(src) - dst).max() < 1e-6


def test_affine_lstsq_recovers_shear():
    rng = np.random.default_rng(2)
    A = np.array([[1.02, 0.06, 133.0], [-0.04, 0.98, -87.0]])
    src = rng.uniform(0, 4000, (50, 2))
    dst = src @ A[:, :2].T + A[:, 2]
    assert np.abs(fit_affine_lstsq(src, dst).apply(src) - dst).max() < 1e-8


def test_pixel_micron_affine_roundtrip():
    """The single conversion point between the two coordinate systems."""
    grid = GridSpec(-100.0, -50.0, 2.5, 800, 900)
    tf = from_rigid(23.0, (140.0, -75.0), center_um=(1000, 1000))
    back = affine_from_pixels(affine_to_pixels(tf, grid), grid)
    assert np.abs(back.affine - tf.affine).max() < 1e-9


def test_compose_and_invert():
    rng = np.random.default_rng(3)
    a = from_rigid(15.0, (100.0, 40.0))
    b = from_rigid(-40.0, (-70.0, 210.0), scale=1.05)
    pts = rng.uniform(0, 3000, (30, 2))
    assert np.abs(a.compose(b).apply(pts) - a.apply(b.apply(pts))).max() < 1e-8
    assert np.abs(a.invert().apply(a.apply(pts)) - pts).max() < 1e-8


def test_displacement_inverse_roundtrip():
    """`inverse_apply` must undo `apply` including the non-rigid part.

    `invert()` deliberately drops the displacement; using it for warping instead of
    `inverse_apply` makes the whole non-rigid stage a silent no-op.
    """
    grid = GridSpec(0.0, 0.0, 10.0, 400, 400)
    rows, cols = np.mgrid[0:400, 0:400]
    disp = np.stack([30 * np.sin(rows / 60.0), 20 * np.cos(cols / 70.0)]).astype(np.float32)
    tf = Transform(affine=from_rigid(9.0, (55.0, -20.0)).affine, displacement=disp, grid=grid)

    rng = np.random.default_rng(4)
    pts = rng.uniform(400, 3600, (300, 2))
    assert np.abs(tf.inverse_apply(tf.apply(pts)) - pts).max() < 1.0
    assert tf.max_displacement_um > 20.0


def test_transform_io_roundtrip(tmp_path):
    grid = GridSpec(0.0, 0.0, 10.0, 60, 60)
    disp = np.ones((2, 60, 60), dtype=np.float32) * 3.0
    tf = Transform(affine=from_rigid(11.0, (5.0, 6.0)).affine, displacement=disp, grid=grid,
                   meta={"track": "image"})
    tf.save(tmp_path)
    back = Transform.load(tmp_path)
    assert np.allclose(back.affine, tf.affine)
    assert np.allclose(back.displacement, tf.displacement)
    assert back.grid.pixel_um == grid.pixel_um
    assert back.meta["track"] == "image"


# --------------------------------------------------------------------------- #
# Track A — point set
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("angle,tx,ty,mirror,drop,jitter", [
    (12.0, 220.0, -140.0, False, 0.00, 5.0),
    (47.0, -310.0, 260.0, False, 0.25, 15.0),
    (-73.0, 90.0, 410.0, True, 0.15, 10.0),      # mounted flipped
    (155.0, -40.0, -330.0, True, 0.30, 20.0),    # flipped + a third of islets missing
])
def test_point_set_recovers_known_transform(angle, tx, ty, mirror, drop, jitter):
    """Islets appear and vanish between sections, so the fit must survive missing points."""
    rng = np.random.default_rng(11)
    fixed = _islets(40, 5)
    truth = from_rigid(angle, (tx, ty), mirror=mirror, center_um=(2000, 2000))

    mv_xy = truth.invert().apply(fixed.loc[:, ["x_um", "y_um"]].to_numpy())
    mv_xy = mv_xy + rng.normal(0, jitter, mv_xy.shape)
    moving = fixed.copy()
    moving["x_um"], moving["y_um"] = mv_xy[:, 0], mv_xy[:, 1]
    moving["n_cells"] *= rng.uniform(0.8, 1.25, len(moving))
    moving = moving[rng.random(len(moving)) >= drop].reset_index(drop=True)
    moving = pd.concat([moving, _islets(6, 99)], ignore_index=True)   # section-specific islets

    res = register_point_set(moving, fixed, inlier_tol_um=150.0)
    got = res.transform
    resid = np.linalg.norm(got.apply(mv_xy) - fixed.loc[:, ["x_um", "y_um"]].to_numpy(), axis=1)

    assert got.is_mirrored == mirror
    assert got.scale == pytest.approx(1.0, abs=0.05)
    assert np.median(resid) < 3 * max(jitter, 5.0)
    assert res.n_inliers >= 10


def test_point_set_declines_when_there_is_nothing_to_match():
    res = register_point_set(_islets(2, 1), _islets(2, 2))
    assert res.n_inliers == 0
    assert any("too few structures" in n for n in res.notes)


def test_ransac_rejects_implausible_scale():
    """Both instruments report calibrated microns, so a large scale change is a bad sample."""
    rng = np.random.default_rng(7)
    mv = rng.uniform(0, 1000, (20, 2))
    fx = mv * 3.0                                   # 3x scale: physically impossible here
    cands = [(i, i) for i in range(20)]
    tf, inliers = ransac_affine(mv, fx, cands, n_iter=500, rng=rng)
    assert tf is None or tf.scale <= 1.25


def test_icp_refines_a_rough_seed():
    rng = np.random.default_rng(9)
    fixed = _islets(35, 3).loc[:, ["x_um", "y_um"]].to_numpy()
    truth = from_rigid(8.0, (60.0, -35.0), center_um=(2000, 2000))
    moving = truth.invert().apply(fixed) + rng.normal(0, 3.0, fixed.shape)

    seed = from_rigid(5.0, (40.0, -20.0), center_um=(2000, 2000))
    before, _ = nearest_neighbour_rmse(seed.apply(moving), fixed, max_dist_um=400)
    after_tf = icp(moving, fixed, initial=seed, max_dist_um=400)
    after, _ = nearest_neighbour_rmse(after_tf.apply(moving), fixed, max_dist_um=400)
    assert after < before


# --------------------------------------------------------------------------- #
# Track B — image
# --------------------------------------------------------------------------- #

def _blob_image(seed=21, n_blobs=14, pixel_um=8.0):
    from scipy import ndimage

    r = np.random.default_rng(seed)
    pts = np.vstack([r.normal(c, 60, (250, 2)) for c in r.uniform(400, 3600, (n_blobs, 2))])
    grid = grid_for(pts, pixel_um, margin_um=400)
    img = ndimage.gaussian_filter(normalize_image(rasterize_points(pts, grid)), 2.0)
    return pts, grid, img


@pytest.mark.parametrize("angle,tx,ty", [(0.0, 150.0, -90.0), (23.0, -200.0, 120.0),
                                         (-41.0, 80.0, 260.0)])
def test_image_track_recovers_known_rigid_transform(angle, tx, ty):
    from scipy import ndimage

    pts, grid, fixed = _blob_image()
    truth = from_rigid(angle, (tx, ty), center_um=(2000, 2000))
    mv_pts = truth.invert().apply(pts)
    moving = ndimage.gaussian_filter(normalize_image(rasterize_points(mv_pts, grid)), 2.0)

    res = register_image(moving, fixed, grid, nonrigid=False, coarse_step_deg=3.0)
    got = res.transform
    err = abs(((got.rotation_deg - truth.rotation_deg + 180) % 360) - 180)
    resid = np.linalg.norm(got.apply(mv_pts) - pts, axis=1)

    assert err < 1.5
    assert np.median(resid) < 30.0
    assert "coarse" in res.stages


def test_bspline_improves_a_nonrigid_warp_and_respects_the_cap():
    """The displacement sign convention and the safety cap, together.

    SimpleITK returns a fixed->moving field while `Transform.displacement` is moving->fixed.
    Getting that backwards does not fail loudly — it doubles the misalignment — so the test
    asserts that alignment actually *improves*.
    """
    from scipy import ndimage

    _pts, grid, fixed = _blob_image(seed=31, n_blobs=16, pixel_um=10.0)
    H, W = grid.shape
    rows, cols = np.mgrid[0:H, 0:W]
    dx = 12.0 * np.sin(2 * np.pi * rows / H)
    dy = 12.0 * np.cos(2 * np.pi * cols / W)
    moving = ndimage.map_coordinates(fixed, [rows + dy, cols + dx], order=1, mode="constant")

    before = dice(moving > 0.05, fixed > 0.05)
    out = refine_bspline(moving, fixed, grid, max_disp_um=250.0, grid_spacing_um=400.0,
                         iterations=150)
    if out is None:
        pytest.skip("SimpleITK unavailable")
    field, field_grid = out
    tf = Transform(affine=np.array([[1.0, 0, 0], [0, 1.0, 0]]), displacement=field,
                   grid=field_grid)
    after = dice(tf.warp_mask(moving, field_grid) > 0.05, fixed > 0.05)
    assert after > before + 0.05, f"non-rigid stage did not improve alignment ({before}->{after})"

    capped = refine_bspline(moving, fixed, grid, max_disp_um=30.0, grid_spacing_um=400.0,
                            iterations=150)
    assert np.linalg.norm(capped[0], axis=0).max() <= 30.0 + 1e-4


def test_ecc_failure_is_not_fatal():
    """A donor whose ECC diverges must not abort a cohort run."""
    grid = GridSpec(0.0, 0.0, 10.0, 64, 64)
    blank = np.zeros(grid.shape, dtype=np.float32)
    tf, val = refine_ecc(blank, blank, grid, Transform(), motion="affine")
    assert np.isnan(val) or np.isfinite(val)


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def test_dice_bounds():
    a = np.zeros((10, 10), bool)
    a[2:8, 2:8] = True
    assert dice(a, a) == pytest.approx(1.0)
    assert dice(a, ~a) == pytest.approx(0.0)
    assert np.isnan(dice(np.zeros((4, 4), bool), np.zeros((4, 4), bool)))


def test_nearest_neighbour_rmse_uses_mutual_pairs():
    """One-way NN distance can be driven to zero by collapsing everything onto one point."""
    src = np.array([[0.0, 0.0], [1000.0, 1000.0]])
    dst = np.array([[0.0, 0.0]])
    rmse, n = nearest_neighbour_rmse(src, dst, max_dist_um=50.0)
    assert n == 1 and rmse == pytest.approx(0.0)
