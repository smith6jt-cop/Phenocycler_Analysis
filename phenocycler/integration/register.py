"""
S4 — registering two serial sections.

Two independent tracks, both run, best-by-QC wins and the choice is recorded:

**Track A — point set.** RANSAC over candidate islet correspondences, then ICP refinement.
Needs no images at all, so it works on any donor whose cell tables exist. Islets are natural
fiducials in pancreas: there are tens to hundreds of them, they are well separated, and they
carry size and shape features that constrain matching far better than position alone.

**Track B — image.** Coarse orientation search, then opencv ECC (Euclidean then affine),
then an optional SimpleITK B-spline. Runs on DAPI-to-DAPI where the morphology images are
mounted, and falls back to tissue/lineage-density rasters derived from the cell tables where
they are not.

Three things this design takes seriously that the prototype it replaces did not.

*Rotation and mirroring are not optional.* Serial sections are picked up and mounted by hand;
they arrive at arbitrary rotation and are sometimes flipped. A translation-only registration
(``phase_cross_correlation`` and nothing else) cannot represent either, so it will happily
return a confident-looking shift that is simply wrong. The coarse stage searches rotation
exhaustively and evaluates the mirrored hypothesis explicitly.

*Non-rigid deformation is real but bounded.* Tissue tears, folds and stretches during
sectioning and mounting. A B-spline absorbs that. But given enough freedom a B-spline will
also happily warp any two images into agreement, so displacement is capped at
``reg_max_disp_um`` and the cap is enforced, not merely configured.

*A registration that fails must say so.* Every result carries a QC block, and
``qc.py`` gates on it. A donor whose sections genuinely do not correspond — different blocks,
too large a gap, tissue destroyed — should drop out of structure-level analysis rather than
contribute a plausible-looking alignment to it.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, read_cell_table
from .manifest import load_manifest, paired_rows
from .rasterize import (GridSpec, grid_for, lineage_density, normalize_image,
                        rasterize_structures, tissue_mask)
from .structures import structures_path
from .transform import (Transform, affine_from_pixels, affine_to_pixels, fit_affine_lstsq,
                        fit_rigid_umeyama)

#: Lineages used as density channels when registering without morphology images. Endocrine
#: is effectively an islet map and Vascular a vessel map; both are structural and largely
#: disease-independent.
DENSITY_LINEAGES = ("Endocrine", "Vascular", "Exocrine")


@dataclass
class RegistrationResult:
    """A transform plus everything needed to judge whether to believe it."""

    transform: Transform
    track: str                       # 'point_set' | 'image'
    grid: Optional[GridSpec] = None
    n_inliers: int = 0
    n_candidates: int = 0
    rmse_um: float = float("nan")
    tissue_dice: float = float("nan")
    islet_dice: float = float("nan")
    mirrored: bool = False
    stages: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def inlier_fraction(self) -> float:
        return self.n_inliers / self.n_candidates if self.n_candidates else 0.0

    def score(self) -> float:
        """Single comparable number for picking between tracks.

        Weighted toward tissue Dice because it is the measure least able to be gamed: a
        transform can drive islet RMSE down by matching a handful of convenient pairs, but it
        cannot make two tissue outlines overlap unless it is genuinely close to right.
        """
        dice = 0.0 if np.isnan(self.tissue_dice) else self.tissue_dice
        idice = 0.0 if np.isnan(self.islet_dice) else self.islet_dice
        rmse_term = 0.0
        if np.isfinite(self.rmse_um):
            rmse_term = float(np.exp(-self.rmse_um / 150.0))
        return 0.5 * dice + 0.3 * idice + 0.2 * rmse_term

    def to_dict(self) -> dict:
        return {
            "track": self.track,
            "n_inliers": self.n_inliers,
            "n_candidates": self.n_candidates,
            "inlier_fraction": self.inlier_fraction,
            "rmse_um": self.rmse_um,
            "tissue_dice": self.tissue_dice,
            "islet_dice": self.islet_dice,
            "mirrored": self.mirrored,
            "score": self.score(),
            "stages": self.stages,
            "notes": self.notes,
            "transform": self.transform.to_dict(),
        }


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def dice(a: np.ndarray, b: np.ndarray) -> float:
    """Sorensen-Dice between two boolean masks."""
    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    denom = a.sum() + b.sum()
    if denom == 0:
        return float("nan")
    return float(2.0 * np.logical_and(a, b).sum() / denom)


def nearest_neighbour_rmse(src: np.ndarray, dst: np.ndarray,
                           max_dist_um: float = np.inf) -> tuple[float, int]:
    """RMSE of mutual-nearest-neighbour pairs, and how many there were.

    Mutual nearest rather than one-way: a one-way nearest-neighbour distance can be driven
    arbitrarily low by collapsing everything onto one dense cluster, which is exactly the
    degenerate solution an unconstrained optimiser will find.
    """
    if len(src) == 0 or len(dst) == 0:
        return float("nan"), 0
    from scipy.spatial import cKDTree

    ts, td = cKDTree(src), cKDTree(dst)
    d_sd, i_sd = td.query(src, k=1)
    _d_ds, i_ds = ts.query(dst, k=1)
    mutual = i_ds[i_sd] == np.arange(len(src))
    keep = mutual & (d_sd <= max_dist_um)
    if not keep.any():
        return float("nan"), 0
    return float(np.sqrt(np.mean(d_sd[keep] ** 2))), int(keep.sum())


# --------------------------------------------------------------------------- #
# Track A — point-set registration on structure centroids
# --------------------------------------------------------------------------- #

def _candidate_pairs(
    moving: pd.DataFrame,
    fixed: pd.DataFrame,
    *,
    area_ratio: float = 3.0,
    size_ratio: float = 3.0,
) -> list[tuple[int, int]]:
    """Plausible islet correspondences, filtered on size before geometry is considered.

    An islet's cell count and hull area are roughly conserved across a few microns of
    sectioning, so a 30-cell islet is not a candidate partner for a 300-cell one no matter
    how close they land. Filtering on that first shrinks the RANSAC search space by an order
    of magnitude and, more importantly, removes the correspondences most likely to produce a
    confidently wrong fit.
    """
    pairs = []
    m_n = moving["n_cells"].to_numpy(np.float64)
    f_n = fixed["n_cells"].to_numpy(np.float64)
    m_a = moving.get("area_um2", pd.Series(np.zeros(len(moving)))).to_numpy(np.float64)
    f_a = fixed.get("area_um2", pd.Series(np.zeros(len(fixed)))).to_numpy(np.float64)

    for i in range(len(moving)):
        for j in range(len(fixed)):
            rn = m_n[i] / max(f_n[j], 1e-9)
            if not (1 / size_ratio <= rn <= size_ratio):
                continue
            if m_a[i] > 0 and f_a[j] > 0:
                ra = m_a[i] / f_a[j]
                if not (1 / area_ratio <= ra <= area_ratio):
                    continue
            pairs.append((i, j))
    return pairs


def ransac_affine(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    candidates: Sequence[tuple[int, int]],
    *,
    n_iter: int = 4000,
    inlier_tol_um: float = 150.0,
    rng: Optional[np.random.Generator] = None,
    similarity_only: bool = True,
) -> tuple[Optional[Transform], np.ndarray]:
    """RANSAC a transform from noisy candidate correspondences.

    Fits a *similarity* (rotation, isotropic scale, translation, optional mirror) rather than
    a free affine by default. With 4 degrees of freedom instead of 6 it cannot shear a bad
    correspondence set into apparent agreement, and two serial sections of the same block
    genuinely differ by little more than a rigid motion at this stage — shear belongs in the
    later refinement, if anywhere.
    """
    rng = rng or np.random.default_rng(0)
    if len(candidates) < 3:
        return None, np.zeros(0, dtype=int)

    cand = np.asarray(candidates)
    src_all = moving_xy[cand[:, 0]]
    dst_all = fixed_xy[cand[:, 1]]

    best_tf, best_inliers = None, np.zeros(0, dtype=int)
    n_sample = 2 if similarity_only else 3

    for _ in range(n_iter):
        idx = rng.choice(len(cand), size=n_sample, replace=False)
        s, d = src_all[idx], dst_all[idx]
        if np.linalg.norm(s[0] - s[1]) < 1e-6 or np.linalg.norm(d[0] - d[1]) < 1e-6:
            continue
        try:
            tf = (fit_rigid_umeyama(s, d, allow_reflection=True) if similarity_only
                  else fit_affine_lstsq(s, d))
        except (ValueError, np.linalg.LinAlgError):
            continue

        # Reject physically implausible fits before scoring: both modalities report calibrated
        # microns, so a scale far from 1 means a bad sample, not a real observation.
        if not (0.8 <= tf.scale <= 1.25):
            continue

        resid = np.linalg.norm(tf.apply(src_all) - dst_all, axis=1)
        inliers = np.flatnonzero(resid <= inlier_tol_um)
        if len(inliers) > len(best_inliers):
            best_tf, best_inliers = tf, inliers

    if best_tf is None or len(best_inliers) < 3:
        return best_tf, best_inliers

    # Refit on all inliers; a 2-point sample fixes the model but is a poor estimate of it.
    try:
        refined = (fit_rigid_umeyama(src_all[best_inliers], dst_all[best_inliers],
                                     allow_reflection=True)
                   if similarity_only else
                   fit_affine_lstsq(src_all[best_inliers], dst_all[best_inliers]))
        resid = np.linalg.norm(refined.apply(src_all) - dst_all, axis=1)
        best_inliers = np.flatnonzero(resid <= inlier_tol_um)
        best_tf = refined
    except (ValueError, np.linalg.LinAlgError):
        pass
    return best_tf, best_inliers


def icp(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    initial: Optional[Transform] = None,
    *,
    max_iter: int = 60,
    max_dist_um: float = 250.0,
    tol_um: float = 0.05,
    similarity_only: bool = True,
) -> Transform:
    """Iterative closest point refinement over mutual nearest neighbours.

    Restricted to mutual pairs within ``max_dist_um``: plain closest-point ICP on two islet
    fields that differ in count (islets appear and vanish between sections) will drag
    unmatched points onto whatever is nearest and bias the fit.
    """
    from scipy.spatial import cKDTree

    tf = initial or Transform()
    tree_fixed = cKDTree(fixed_xy)
    prev = np.inf

    for _ in range(max_iter):
        cur = tf.apply(moving_xy)
        d, j = tree_fixed.query(cur, k=1)
        tree_cur = cKDTree(cur)
        _d2, i2 = tree_cur.query(fixed_xy, k=1)
        mutual = i2[j] == np.arange(len(cur))
        keep = mutual & (d <= max_dist_um)
        if keep.sum() < 3:
            break

        try:
            step = (fit_rigid_umeyama(cur[keep], fixed_xy[j[keep]]) if similarity_only
                    else fit_affine_lstsq(cur[keep], fixed_xy[j[keep]]))
            # NB: no reflection here. The seed already fixed handedness; letting each ICP
            # step flip would let it oscillate between the two solutions.
        except (ValueError, np.linalg.LinAlgError):
            break
        tf = step.compose(tf)

        rmse = float(np.sqrt(np.mean(d[keep] ** 2)))
        if abs(prev - rmse) < tol_um:
            break
        prev = rmse
    return tf


def register_point_set(
    moving_structures: pd.DataFrame,
    fixed_structures: pd.DataFrame,
    *,
    inlier_tol_um: float = 150.0,
    seed: int = 0,
) -> RegistrationResult:
    """Track A: align two islet fields."""
    res = RegistrationResult(transform=Transform(), track="point_set")
    if moving_structures is None or fixed_structures is None:
        res.notes.append("no structures available")
        return res
    if len(moving_structures) < 3 or len(fixed_structures) < 3:
        res.notes.append(
            f"too few structures to register ({len(moving_structures)} moving, "
            f"{len(fixed_structures)} fixed; need >= 3 each)")
        return res

    mv = moving_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    fx = fixed_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)

    cands = _candidate_pairs(moving_structures, fixed_structures)
    res.n_candidates = len(cands)
    if not cands:
        res.notes.append("no size-compatible islet pairs")
        return res

    rng = np.random.default_rng(seed)
    tf, inliers = ransac_affine(mv, fx, cands, inlier_tol_um=inlier_tol_um, rng=rng)
    if tf is None:
        res.notes.append("RANSAC found no consistent transform")
        return res
    res.stages.append("ransac")

    tf = icp(mv, fx, initial=tf, max_dist_um=inlier_tol_um * 1.5)
    res.stages.append("icp")

    rmse, n_pairs = nearest_neighbour_rmse(tf.apply(mv), fx, max_dist_um=inlier_tol_um * 2)
    res.transform = tf
    res.n_inliers = int(n_pairs)
    res.rmse_um = rmse
    res.mirrored = tf.is_mirrored
    return res


# --------------------------------------------------------------------------- #
# Track B — image / raster registration
# --------------------------------------------------------------------------- #

def coarse_orientation(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    angles_deg: Optional[Sequence[float]] = None,
    try_mirror: bool = True,
    coarse_pixel_um: float = 25.0,
) -> list[Transform]:
    """Candidate coarse transforms, ranked by correlation after a translation search.

    Exhaustive over rotation because there is no cheap way to recover a large rotation from
    a translation-only correlation, and serial sections routinely differ by tens of degrees.
    Each hypothesis gets its optimal translation from phase correlation, which is exact and
    fast, so the search is only over angle (and the mirror flag).

    Runs at ``coarse_pixel_um`` rather than the registration grid's own resolution. A real
    section is 10-20 mm across, so at the 2 um working resolution the raster is 5000-10000 px
    on a side and rotating it ~144 times (72 angles x 2 mirror hypotheses) is intractable.
    Orientation is a low-frequency property — a 25 um pixel resolves it perfectly well — and
    because the transform is expressed in **microns** it drops straight into the fine stage at
    full resolution with no conversion. That is precisely why ``Transform`` is defined in
    micron space rather than pixels.

    Returns several candidates rather than one: the best-correlating angle is not reliably
    the right one on low-contrast masks, so the fine stage refines the top few and picks by
    its own metric.
    """
    from scipy import ndimage
    from skimage.registration import phase_cross_correlation

    angles = list(angles_deg) if angles_deg is not None else list(np.arange(0, 360, 5.0))

    factor = max(1, int(round(coarse_pixel_um / grid.pixel_um)))
    if factor > 1:
        moving = ndimage.zoom(moving, 1.0 / factor, order=1)
        fixed = ndimage.zoom(fixed, 1.0 / factor, order=1)
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, grid.pixel_um * factor,
                        fixed.shape[0], fixed.shape[1])

    cx = grid.width / 2.0
    cy = grid.height / 2.0
    fixed_n = normalize_image(fixed)

    scored: list[tuple[float, Transform]] = []
    for mirror in ([False, True] if try_mirror else [False]):
        base = moving[:, ::-1] if mirror else moving
        for ang in angles:
            rot = ndimage.rotate(base, ang, reshape=False, order=1, mode="constant", cval=0.0)
            rot_n = normalize_image(rot)
            try:
                shift, _err, _phase = phase_cross_correlation(
                    fixed_n, rot_n, upsample_factor=4, normalization=None)
            except TypeError:      # older scikit-image signature
                shift, _err, _phase = phase_cross_correlation(fixed_n, rot_n, upsample_factor=4)
            shifted = ndimage.shift(rot_n, shift, order=1, mode="constant", cval=0.0)

            denom = np.linalg.norm(shifted) * np.linalg.norm(fixed_n)
            corr = float((shifted * fixed_n).sum() / denom) if denom > 0 else 0.0

            # Build the equivalent micron transform. scipy rotates about the array centre and
            # `shift` is (row, col); both conversions happen here so nothing downstream has to
            # know about either convention.
            th = np.radians(ang)
            c, s = np.cos(th), np.sin(th)
            lin = np.array([[c, s], [-s, c]])         # ndimage.rotate is clockwise-positive
            if mirror:
                lin = lin @ np.array([[-1.0, 0.0], [0.0, 1.0]])
            centre = np.array([cx, cy])
            t_px = centre - lin @ centre + np.array([shift[1], shift[0]])
            tf = affine_from_pixels(np.column_stack([lin, t_px]), grid)
            tf.meta.update({"coarse_angle_deg": ang, "coarse_mirror": mirror,
                            "coarse_corr": corr})
            scored.append((corr, tf))

    scored.sort(key=lambda t: -t[0])
    return [tf for _corr, tf in scored]


def refine_ecc(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    initial: Transform,
    *,
    motion: str = "affine",
    iterations: int = 300,
    eps: float = 1e-6,
    max_pixels: int = 4_000_000,
) -> tuple[Transform, float]:
    """Refine with opencv's ECC (enhanced correlation coefficient).

    ECC optimises correlation directly and is robust to the intensity differences that make
    cross-modal registration hard — a DAPI image and a protein channel do not share a linear
    intensity relationship, but they do share structure. Returns ``(transform, ecc_value)``;
    on failure to converge it returns the input transform and ``nan`` rather than raising, so
    one bad donor does not abort a cohort run.
    """
    import cv2
    from scipy import ndimage

    warp_mode = {"euclidean": cv2.MOTION_EUCLIDEAN, "affine": cv2.MOTION_AFFINE}[motion]
    f = normalize_image(fixed).astype(np.float32)
    m = normalize_image(moving).astype(np.float32)

    # Cap the working resolution. ECC is O(pixels) per iteration and runs hundreds of them,
    # so a full-resolution whole-section raster is not viable. Downsampling is safe because
    # the affine lives in micron space: `grid` is rescaled alongside the images, so the
    # recovered transform is identical whatever resolution it was estimated at.
    if f.size > max_pixels:
        shrink = float(np.sqrt(max_pixels / f.size))
        f = ndimage.zoom(f, shrink, order=1)
        m = ndimage.zoom(m, shrink, order=1)
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, grid.pixel_um / shrink,
                        f.shape[0], f.shape[1])

    # ECC estimates the fixed->moving warp; we hold moving->fixed, so seed with the inverse.
    warp = affine_to_pixels(initial.invert(), grid).astype(np.float32)
    if motion == "euclidean":
        # Force an exact rotation matrix; ECC rejects a seed that is not one.
        u, _s, vt = np.linalg.svd(warp[:, :2])
        warp[:, :2] = u @ vt

    criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, iterations, eps)
    try:
        value, warp = cv2.findTransformECC(f, m, warp, warp_mode, criteria, None, 5)
    except cv2.error:
        return initial, float("nan")

    fixed_to_moving = affine_from_pixels(np.asarray(warp, dtype=np.float64), grid)
    return fixed_to_moving.invert(), float(value)


def refine_bspline(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    max_disp_um: float = 200.0,
    grid_spacing_um: float = 400.0,
    iterations: int = 200,
    max_pixels: int = 2_000_000,
) -> Optional[np.ndarray]:
    """Non-rigid B-spline refinement; returns a ``(2, H, W)`` micron displacement field.

    Deliberately coarse: the control-point spacing is hundreds of microns, so the field can
    absorb the smooth stretch and shear of sectioning without being able to shrink-wrap one
    section's islets onto another's. The result is then hard-capped at ``max_disp_um``.

    That cap is the important part. A B-spline with enough freedom can align any two images,
    including two that have nothing to do with each other, and a "successful" registration of
    unrelated sections is far more damaging than an obvious failure.
    """
    try:
        import SimpleITK as sitk
    except ImportError:
        print("[register] SimpleITK not installed — skipping the non-rigid stage "
              "(pip install -e '.[integration]')", flush=True)
        return None

    from scipy import ndimage

    fa = normalize_image(fixed).astype(np.float32)
    ma = normalize_image(moving).astype(np.float32)
    if fa.size > max_pixels:
        shrink = float(np.sqrt(max_pixels / fa.size))
        fa = ndimage.zoom(fa, shrink, order=1)
        ma = ndimage.zoom(ma, shrink, order=1)
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, grid.pixel_um / shrink,
                        fa.shape[0], fa.shape[1])
    f = sitk.GetImageFromArray(fa)
    m = sitk.GetImageFromArray(ma)

    mesh = [max(1, int(round(s * grid.pixel_um / grid_spacing_um)))
            for s in (grid.width, grid.height)]
    try:
        tx = sitk.BSplineTransformInitializer(f, mesh)
        reg = sitk.ImageRegistrationMethod()
        reg.SetMetricAsMattesMutualInformation(numberOfHistogramBins=32)
        reg.SetMetricSamplingStrategy(reg.RANDOM)
        reg.SetMetricSamplingPercentage(0.2, seed=1)
        reg.SetInterpolator(sitk.sitkLinear)
        reg.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=iterations,
                                          convergenceMinimumValue=1e-6, convergenceWindowSize=10)
        reg.SetOptimizerScalesFromPhysicalShift()
        reg.SetInitialTransform(tx, inPlace=True)
        reg.SetShrinkFactorsPerLevel([4, 2, 1])
        reg.SetSmoothingSigmasPerLevel([2, 1, 0])
        reg.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()
        out_tx = reg.Execute(f, m)
    except Exception as exc:  # noqa: BLE001 - ITK raises a wide variety on bad input
        print(f"[register] B-spline stage failed ({exc}); keeping the affine result",
              flush=True)
        return None

    disp_filter = sitk.TransformToDisplacementFieldFilter()
    disp_filter.SetReferenceImage(f)
    disp_filter.SetOutputPixelType(sitk.sitkVectorFloat64)
    disp = sitk.GetArrayFromImage(disp_filter.Execute(out_tx))   # (H, W, 2), index units

    # Sign convention. SimpleITK's registration returns the transform used to *resample* the
    # moving image into the fixed frame, i.e. it maps fixed -> moving. `Transform.displacement`
    # is the opposite: moving -> fixed, added after the affine. So the field is negated here.
    # Getting this backwards does not fail loudly — it produces a transform that doubles the
    # misalignment instead of removing it — so it is asserted by tests/test_integration_register.
    # Vector components are (x, y); pixel size converts index units to microns.
    field = -np.stack([disp[..., 0], disp[..., 1]]).astype(np.float32) * grid.pixel_um
    bspline_grid = grid

    mag = np.linalg.norm(field, axis=0)
    over = mag > max_disp_um
    if over.any():
        scale = np.ones_like(mag)
        scale[over] = max_disp_um / mag[over]
        field = field * scale
        print(f"[register] capped {100*over.mean():.1f}% of the displacement field at "
              f"{max_disp_um:.0f} um", flush=True)
    return field, bspline_grid


def register_image(
    moving: np.ndarray,
    fixed: np.ndarray,
    grid: GridSpec,
    *,
    nonrigid: bool = True,
    max_disp_um: float = 200.0,
    n_coarse_candidates: int = 3,
    coarse_step_deg: float = 5.0,
) -> RegistrationResult:
    """Track B: coarse orientation -> ECC euclidean -> ECC affine -> optional B-spline."""
    res = RegistrationResult(transform=Transform(), track="image", grid=grid)

    candidates = coarse_orientation(
        moving, fixed, grid, angles_deg=np.arange(0, 360, coarse_step_deg))
    if not candidates:
        res.notes.append("coarse orientation search produced no candidates")
        return res
    res.stages.append("coarse")

    best, best_val = None, -np.inf
    for cand in candidates[:n_coarse_candidates]:
        tf, val = refine_ecc(moving, fixed, grid, cand, motion="euclidean")
        if not np.isfinite(val):
            continue
        if val > best_val:
            best, best_val = tf, val
    if best is None:
        # Every ECC seed diverged; keep the best correlating coarse candidate so the caller
        # still gets a usable starting point rather than the identity.
        res.transform = candidates[0]
        res.mirrored = candidates[0].is_mirrored
        res.notes.append("ECC did not converge from any coarse seed; using coarse only")
        return res
    res.stages.append("ecc_euclidean")

    tf, val = refine_ecc(moving, fixed, grid, best, motion="affine")
    if np.isfinite(val) and val >= best_val:
        best, best_val = tf, val
        res.stages.append("ecc_affine")

    res.transform = best
    res.mirrored = best.is_mirrored
    res.notes.append(f"ecc={best_val:.4f}")

    if nonrigid:
        warped = best.warp_mask(moving, grid)
        out = refine_bspline(warped, fixed, grid, max_disp_um=max_disp_um)
        if out is not None:
            field, field_grid = out
            res.transform = Transform(affine=best.affine, displacement=field,
                                      grid=field_grid, meta=dict(best.meta))
            res.stages.append("bspline")
    return res


# --------------------------------------------------------------------------- #
# Orchestration
# --------------------------------------------------------------------------- #

def _load_structures(cfg: PipelineConfig, modality: str, donor: str, roi: str,
                     kind: str = "islet") -> Optional[pd.DataFrame]:
    path = structures_path(cfg, modality, donor, roi, kind)
    if not path.exists():
        return None
    return pd.read_parquet(path)


def evaluate(
    result: RegistrationResult,
    moving_table: CellTable,
    fixed_table: CellTable,
    grid: GridSpec,
    moving_structures: Optional[pd.DataFrame] = None,
    fixed_structures: Optional[pd.DataFrame] = None,
    *,
    qc_pixel_um: float = 10.0,
) -> RegistrationResult:
    """Fill in the QC block for a candidate transform.

    Deliberately measured on quantities the transform was *not* optimised against wherever
    possible: the point-set track never sees the tissue mask, and the image track never sees
    islet centroids, so each track's headline metric is an out-of-sample check on the other.

    Evaluated on a deliberately coarse grid. Dice between two tissue outlines is a
    low-frequency measurement — at 10 um a whole-section mask is a few hundred pixels on a
    side, the metric is unchanged to three decimals, and the warp is two orders of magnitude
    cheaper than at the 2 um working resolution. Cost matters because this runs once per
    track per donor.
    """
    if qc_pixel_um > grid.pixel_um:
        factor = qc_pixel_um / grid.pixel_um
        grid = GridSpec(grid.origin_x_um, grid.origin_y_um, qc_pixel_um,
                        max(1, int(round(grid.height / factor))),
                        max(1, int(round(grid.width / factor))))

    mv_mask = tissue_mask(moving_table.xy(), grid)
    fx_mask = tissue_mask(fixed_table.xy(), grid)
    warped = result.transform.warp_mask(mv_mask, grid)
    result.tissue_dice = dice(warped, fx_mask)

    if moving_structures is not None and fixed_structures is not None \
            and len(moving_structures) and len(fixed_structures):
        mv = moving_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        fx = fixed_structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
        rmse, n = nearest_neighbour_rmse(result.transform.apply(mv), fx, max_dist_um=300.0)
        if result.track == "image" or not np.isfinite(result.rmse_um):
            result.rmse_um = rmse
            result.n_inliers = n
            result.n_candidates = result.n_candidates or min(len(mv), len(fx))

        # Islet Dice: dilate each centroid to a ~75 um disc so the metric measures "did the
        # islets land on each other" rather than "did two zero-area points coincide", which
        # would be 0 for any transform. `iterations` on a 3x3 element rather than one large
        # structuring element — same result, far cheaper.
        from scipy import ndimage

        r = max(1, int(round(75.0 / grid.pixel_um)))
        mv_islets = ndimage.binary_dilation(rasterize_structures(moving_structures, grid) > 0,
                                            iterations=r)
        fx_islets = ndimage.binary_dilation(rasterize_structures(fixed_structures, grid) > 0,
                                            iterations=r)
        result.islet_dice = dice(result.transform.warp_mask(mv_islets, grid), fx_islets)
    return result


def register_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    use_images: bool = True,
    nonrigid: Optional[bool] = None,
) -> Optional[RegistrationResult]:
    """Register one donor's two sections, running both tracks and keeping the better one."""
    nonrigid = cfg.reg_nonrigid if nonrigid is None else nonrigid

    fixed_is_pheno = cfg.fixed_modality == "phenocycler"
    fixed_root = cfg.cells_pheno_dir if fixed_is_pheno else cfg.cells_xen_dir
    moving_root = cfg.cells_xen_dir if fixed_is_pheno else cfg.cells_pheno_dir
    fixed_mod = "phenocycler" if fixed_is_pheno else "xenium"
    moving_mod = "xenium" if fixed_is_pheno else "phenocycler"

    try:
        fixed_table = read_cell_table(fixed_root, donor, roi, modality=fixed_mod)
        moving_table = read_cell_table(moving_root, donor, roi, modality=moving_mod)
    except FileNotFoundError as exc:
        print(f"[register] donor {donor}: {exc}", flush=True)
        return None

    grid = grid_for(fixed_table.xy(), cfg.reg_pixel_um, other=moving_table.xy())

    fx_struct = _load_structures(cfg, fixed_mod, donor, roi)
    mv_struct = _load_structures(cfg, moving_mod, donor, roi)

    results: list[RegistrationResult] = []

    # Track A — point set.
    ps = register_point_set(mv_struct, fx_struct,
                            inlier_tol_um=cfg.match_max_dist_um * 0.75)
    if ps.transform is not None and (ps.n_inliers or ps.stages):
        results.append(evaluate(ps, moving_table, fixed_table, grid, mv_struct, fx_struct))

    # Track B — image/raster.
    if use_images:
        mv_img, fx_img = _registration_rasters(cfg, moving_table, fixed_table, grid,
                                               moving_mod, fixed_mod, donor)
        img = register_image(mv_img, fx_img, grid, nonrigid=nonrigid,
                             max_disp_um=cfg.reg_max_disp_um)
        results.append(evaluate(img, moving_table, fixed_table, grid, mv_struct, fx_struct))

    if not results:
        print(f"[register] donor {donor}: no track produced a transform", flush=True)
        return None

    best = max(results, key=lambda r: r.score())
    best.grid = grid
    if len(results) > 1:
        other = [r for r in results if r is not best][0]
        best.notes.append(
            f"chosen over {other.track} (score {best.score():.3f} vs {other.score():.3f})")
    return best


def _registration_rasters(
    cfg: PipelineConfig,
    moving_table: CellTable,
    fixed_table: CellTable,
    grid: GridSpec,
    moving_mod: str,
    fixed_mod: str,
    donor: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Best available raster pair: DAPI-to-DAPI if the images are mounted, else densities.

    Falls back rather than failing. Most donors in a cohort will not have both morphology
    trees mounted at the same time, and a density-based registration is a genuinely useful
    result, not a placeholder.
    """
    from .rasterize import (find_qptiff, image_to_grid, read_qptiff_channel,
                            read_xenium_morphology)

    images: dict[str, Optional[np.ndarray]] = {"phenocycler": None, "xenium": None}

    qptiff = find_qptiff(cfg.images_dir, donor)
    if qptiff is not None:
        try:
            img, um = read_qptiff_channel(qptiff, channel=0, target_um=cfg.reg_pixel_um)
            images["phenocycler"] = image_to_grid(img, um, grid)
            print(f"[register] donor {donor}: PhenoCycler DAPI from {qptiff.name} "
                  f"({um:.2f} um/px)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[register] donor {donor}: could not read {qptiff.name}: {exc}", flush=True)

    man = load_manifest(cfg)
    row = man[(man["donor_id"] == donor) & (man["roi"] == moving_table.roi)]
    bundle = str(row.iloc[0]["xenium_bundle"]) if len(row) else ""
    if bundle and Path(bundle).exists():
        try:
            img, um = read_xenium_morphology(Path(bundle), channel=0,
                                             target_um=cfg.reg_pixel_um)
            images["xenium"] = image_to_grid(img, um, grid)
            print(f"[register] donor {donor}: Xenium DAPI from {Path(bundle).name} "
                  f"({um:.2f} um/px)", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"[register] donor {donor}: could not read Xenium morphology: {exc}",
                  flush=True)

    def fallback(table: CellTable) -> np.ndarray:
        stack = lineage_density(table, grid, DENSITY_LINEAGES)
        return normalize_image(stack.sum(axis=0))

    mv = images[moving_mod] if images[moving_mod] is not None else fallback(moving_table)
    fx = images[fixed_mod] if images[fixed_mod] is not None else fallback(fixed_table)
    if images[moving_mod] is None or images[fixed_mod] is None:
        print(f"[register] donor {donor}: registering on cell-density rasters "
              f"(DAPI unavailable for at least one modality)", flush=True)
    return normalize_image(mv), normalize_image(fx)


def registration_dir(cfg: PipelineConfig, donor: str, roi: str) -> Path:
    base = cfg.registration_dir / f"donor_id={donor}"
    if roi:
        base = base / f"roi={roi}"
    return base


def run_register(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: Optional[str] = None,
    tissue: Optional[str] = None,
    use_images: bool = True,
) -> pd.DataFrame:
    """Stage entry point: register every paired section and persist the transforms.

    One call covers the whole selection. Registration is per (donor, roi) — a donor's
    pancreas and each of its lymph nodes are separate sections with separate transforms.
    """
    import json

    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[register] no paired donors to register", flush=True)
        return pd.DataFrame()

    rows = []
    for _, r in man.iterrows():
        donor, roi = r["donor_id"], r["roi"]
        res = register_donor(cfg, donor, roi, use_images=use_images)
        if res is None:
            continue
        out = registration_dir(cfg, donor, roi)
        out.mkdir(parents=True, exist_ok=True)
        res.transform.save(out)
        (out / "qc.json").write_text(json.dumps(res.to_dict(), indent=2))

        rows.append({
            "donor_id": donor, "roi": roi, "tissue": r.get("tissue", ""),
            "track": res.track,
            "n_inliers": res.n_inliers, "n_candidates": res.n_candidates,
            "rmse_um": res.rmse_um, "tissue_dice": res.tissue_dice,
            "islet_dice": res.islet_dice, "mirrored": res.mirrored,
            "rotation_deg": res.transform.rotation_deg, "scale": res.transform.scale,
            "max_disp_um": res.transform.max_displacement_um,
            "score": res.score(), "stages": "|".join(res.stages),
        })
        print(f"[register] {donor}/{roi}: track={res.track} stages={'|'.join(res.stages)} "
              f"rot={res.transform.rotation_deg:+.1f}deg "
              f"tissue_dice={res.tissue_dice:.3f} islet_rmse={res.rmse_um:.1f}um",
              flush=True)

    return pd.DataFrame(rows)


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Register PhenoCycler <-> Xenium serial sections")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None,
                    help="one ROI (default: every ROI of --tissue, or all tissues)")
    ap.add_argument("--tissue", default=None, help="restrict to one tissue's ROIs")
    ap.add_argument("--no-images", action="store_true",
                    help="point-set track only (skip all image reads)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    summary = run_register(cfg, args.donors, roi=args.roi, tissue=args.tissue,
                           use_images=not args.no_images)
    if len(summary):
        print()
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
