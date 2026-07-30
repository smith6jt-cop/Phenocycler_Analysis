"""
S5 — the transform object, and applying it.

A ``Transform`` is an affine 2x3 in **microns** plus an optional non-rigid displacement
field. Keeping the affine in microns rather than pixels means the transform is meaningful
independently of whatever grid happened to be used to estimate it, which matters because the
coarse, fine and non-rigid stages may all run at different resolutions.

We define our own object rather than using ``spatialdata.transformations`` for two reasons.
First, the integration core has to stay runnable without spatialdata — a donor with no
Xenium zarr should still register from cell tables alone. Second, and more substantively, a
dense displacement field does not round-trip cleanly through an affine-oriented
transformation graph; serialising it as one of ``spatialdata``'s sequences would either lose
the non-rigid part or smuggle it in as an opaque blob. ``to_spatialdata()`` is provided so
the affine part *can* be handed to SpatialData consumers, with the loss made explicit.

Convention, stated once so it is never ambiguous downstream::

    moving (Xenium, by default) --T--> fixed (PhenoCycler, by default)

``[x', y', 1]^T = A @ [x, y, 1]^T`` then plus the displacement sampled at ``(x', y')``.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from .rasterize import GridSpec


@dataclass
class Transform:
    """A micron-space affine plus an optional displacement field.

    Attributes
    ----------
    affine
        2x3 matrix mapping moving microns to fixed microns.
    displacement
        ``(2, H, W)`` field in microns, defined on ``grid``, added *after* the affine.
        ``displacement[0]`` is dx, ``displacement[1]`` is dy.
    grid
        The grid the displacement field lives on. Required whenever displacement is set.
    """

    affine: np.ndarray = field(default_factory=lambda: np.array([[1.0, 0.0, 0.0],
                                                                 [0.0, 1.0, 0.0]]))
    displacement: Optional[np.ndarray] = None
    grid: Optional[GridSpec] = None
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.affine = np.asarray(self.affine, dtype=np.float64).reshape(2, 3)
        if self.displacement is not None:
            self.displacement = np.asarray(self.displacement, dtype=np.float32)
            if self.displacement.ndim != 3 or self.displacement.shape[0] != 2:
                raise ValueError("displacement must have shape (2, H, W)")
            if self.grid is None:
                raise ValueError("a displacement field requires the grid it is defined on")

    # -- properties ---------------------------------------------------------
    @property
    def is_identity(self) -> bool:
        return (self.displacement is None
                and np.allclose(self.affine, [[1, 0, 0], [0, 1, 0]]))

    @property
    def scale(self) -> float:
        """Isotropic scale factor implied by the linear part.

        Should be ~1.0 between two sections of the same block imaged on calibrated
        instruments — both modalities report microns. A scale far from 1 means something is
        wrong with the units, not that the tissue changed size.
        """
        return float(np.sqrt(abs(np.linalg.det(self.affine[:, :2]))))

    @property
    def rotation_deg(self) -> float:
        """Rotation implied by the linear part, in degrees."""
        return float(np.degrees(np.arctan2(self.affine[1, 0], self.affine[0, 0])))

    @property
    def translation_um(self) -> tuple[float, float]:
        return (float(self.affine[0, 2]), float(self.affine[1, 2]))

    @property
    def is_mirrored(self) -> bool:
        """Whether the linear part flips handedness.

        Not an error — serial sections genuinely get flipped when mounted, and a registration
        that cannot represent that will fail on those donors. Recording it makes the flip an
        observation rather than a silent sign error.
        """
        return bool(np.linalg.det(self.affine[:, :2]) < 0)

    @property
    def max_displacement_um(self) -> float:
        if self.displacement is None:
            return 0.0
        return float(np.nanmax(np.linalg.norm(self.displacement, axis=0)))

    # -- application --------------------------------------------------------
    def apply(self, xy_um: np.ndarray) -> np.ndarray:
        """Map (n, 2) moving-frame microns into the fixed frame."""
        xy = np.asarray(xy_um, dtype=np.float64)
        if xy.size == 0:
            return xy.reshape(-1, 2).copy()
        out = xy @ self.affine[:, :2].T + self.affine[:, 2]

        if self.displacement is not None and self.grid is not None:
            out = out + self._sample_displacement(out)
        return out

    def _sample_displacement(self, xy_um: np.ndarray) -> np.ndarray:
        """Bilinear sample of the displacement field at fixed-frame microns."""
        from scipy import ndimage

        px = self.grid.to_pixels(xy_um)
        coords = np.vstack([px[:, 1], px[:, 0]])          # (row, col) for map_coordinates
        dx = ndimage.map_coordinates(self.displacement[0], coords, order=1, mode="nearest")
        dy = ndimage.map_coordinates(self.displacement[1], coords, order=1, mode="nearest")
        return np.column_stack([dx, dy])

    def apply_to_frame(self, df, x_col: str = "x_um", y_col: str = "y_um",
                       out_cols: tuple[str, str] = ("x_um_reg", "y_um_reg")):
        """Write registered coordinates into a DataFrame (returns it, modified in place)."""
        xy = df.loc[:, [x_col, y_col]].to_numpy(np.float64)
        reg = self.apply(xy)
        df[out_cols[0]] = reg[:, 0].astype("float32")
        df[out_cols[1]] = reg[:, 1].astype("float32")
        return df

    def inverse_apply(self, xy_um: np.ndarray, *, n_iter: int = 5) -> np.ndarray:
        """Map (n, 2) fixed-frame microns back into the moving frame.

        The affine part inverts in closed form. The displacement field does not, so it is
        inverted by fixed-point iteration: :meth:`apply` computes ``u = A x`` then
        ``y = u + D(u)``, so recovering ``u`` from ``y`` means solving ``u = y - D(u)``,
        which converges quickly whenever the field is smooth relative to its own magnitude —
        which a coarse B-spline is by construction.

        This is what pull-based warping needs, and it is *not* the same thing as
        :meth:`invert`, which deliberately returns an affine-only transform.
        """
        y = np.asarray(xy_um, dtype=np.float64)
        if y.size == 0:
            return y.reshape(-1, 2).copy()

        u = y.copy()
        if self.displacement is not None and self.grid is not None:
            for _ in range(max(1, n_iter)):
                u = y - self._sample_displacement(u)

        lin_inv = np.linalg.inv(self.affine[:, :2])
        return (u - self.affine[:, 2]) @ lin_inv.T

    def warp_mask(self, mask: np.ndarray, grid: GridSpec) -> np.ndarray:
        """Warp a boolean/float raster from the moving frame into the fixed frame.

        Pull-based: for each fixed pixel, work out where it came from. Push-based warping
        leaves holes wherever the transform expands the image. Uses :meth:`inverse_apply`, so
        the displacement field is honoured — going through :meth:`invert` here would silently
        drop it and make the whole non-rigid stage a no-op.
        """
        from scipy import ndimage

        rows, cols = np.mgrid[0:grid.height, 0:grid.width]
        fixed_um = grid.to_microns(np.column_stack([cols.ravel(), rows.ravel()]))
        moving_um = self.inverse_apply(fixed_um)
        px = grid.to_pixels(moving_um)
        coords = np.vstack([px[:, 1], px[:, 0]])
        order = 0 if mask.dtype == bool else 1
        out = ndimage.map_coordinates(mask.astype(np.float32), coords, order=order,
                                      mode="constant", cval=0.0)
        out = out.reshape(grid.shape)
        return out > 0.5 if mask.dtype == bool else out.astype(np.float32)

    def invert(self) -> "Transform":
        """Affine-only inverse, as a :class:`Transform`.

        A dense displacement field has no closed-form inverse, so it is dropped here and the
        loss is recorded in ``meta``. Use :meth:`inverse_apply` when the non-rigid part
        matters (warping, mapping fixed-frame points back) — this method is for the cases
        that genuinely only need the affine, such as reporting the reverse alignment.
        """
        lin = self.affine[:, :2]
        inv_lin = np.linalg.inv(lin)
        inv_t = -inv_lin @ self.affine[:, 2]
        aff = np.column_stack([inv_lin, inv_t])
        meta = dict(self.meta)
        if self.displacement is not None:
            meta["inverse_is_affine_only"] = True
        return Transform(affine=aff, meta=meta)

    def compose(self, other: "Transform") -> "Transform":
        """``self ∘ other``: apply ``other`` first, then ``self``.

        Only the affine parts compose analytically; a displacement on either side is kept as
        this transform's, which is correct for the way the pipeline builds them up (coarse
        rigid, then affine, then a single non-rigid refinement at the end).
        """
        a, b = self.affine, other.affine
        lin = a[:, :2] @ b[:, :2]
        trans = a[:, :2] @ b[:, 2] + a[:, 2]
        return Transform(
            affine=np.column_stack([lin, trans]),
            displacement=self.displacement if self.displacement is not None
            else other.displacement,
            grid=self.grid if self.displacement is not None else other.grid,
            meta={**other.meta, **self.meta},
        )

    # -- serialisation ------------------------------------------------------
    def to_dict(self) -> dict:
        return {
            "affine": self.affine.tolist(),
            "grid": self.grid.as_dict() if self.grid is not None else None,
            "has_displacement": self.displacement is not None,
            "scale": self.scale,
            "rotation_deg": self.rotation_deg,
            "translation_um": list(self.translation_um),
            "mirrored": self.is_mirrored,
            "max_displacement_um": self.max_displacement_um,
            "meta": self.meta,
        }

    def save(self, directory: Path, *, stem: str = "transform") -> Path:
        """Write ``<stem>.json`` and, when present, ``<stem>_displacement.npz``."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        js = directory / f"{stem}.json"
        js.write_text(json.dumps(self.to_dict(), indent=2))
        if self.displacement is not None:
            np.savez_compressed(directory / f"{stem}_displacement.npz",
                                displacement=self.displacement)
        return js

    @classmethod
    def load(cls, directory: Path, *, stem: str = "transform") -> "Transform":
        directory = Path(directory)
        data = json.loads((directory / f"{stem}.json").read_text())
        grid = GridSpec.from_dict(data["grid"]) if data.get("grid") else None
        disp = None
        npz = directory / f"{stem}_displacement.npz"
        if data.get("has_displacement") and npz.exists():
            disp = np.load(npz)["displacement"]
        return cls(affine=np.array(data["affine"]), displacement=disp, grid=grid,
                   meta=data.get("meta", {}))

    def to_spatialdata(self):
        """The **affine part** as a ``spatialdata.transformations.Affine``.

        Provided so SpatialData/Xenium-Explorer consumers can use the alignment. Lossy by
        construction when a displacement field is present — that part cannot be expressed in
        SpatialData's transformation graph — so it warns rather than silently dropping it.
        """
        try:
            from spatialdata.transformations import Affine
        except ImportError as exc:  # pragma: no cover
            raise ImportError(
                "spatialdata is needed for to_spatialdata(); pip install -e '.[integration]'"
            ) from exc

        if self.displacement is not None:
            print("[transform] WARNING: exporting affine only — the non-rigid displacement "
                  "field cannot be represented in a SpatialData transformation and is being "
                  "dropped from this export", flush=True)
        mat = np.eye(3)
        mat[:2, :] = self.affine
        return Affine(mat, input_axes=("x", "y"), output_axes=("x", "y"))


def from_rigid(angle_deg: float, translation_um: tuple[float, float],
               scale: float = 1.0, mirror: bool = False,
               center_um: tuple[float, float] = (0.0, 0.0)) -> Transform:
    """Build a similarity transform (rotation, isotropic scale, optional mirror, translation).

    ``center_um`` rotates about a point rather than the origin, which is what a coarse
    alignment estimated from two centroids actually means.
    """
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    lin = scale * np.array([[c, -s], [s, c]])
    if mirror:
        lin = lin @ np.array([[-1.0, 0.0], [0.0, 1.0]])
    cx, cy = center_um
    t = np.array([translation_um[0], translation_um[1]]) + np.array([cx, cy]) - lin @ [cx, cy]
    return Transform(affine=np.column_stack([lin, t]))


def affine_from_pixels(affine_px: np.ndarray, grid: GridSpec) -> Transform:
    """Convert a 2x3 estimated on ``grid`` pixels into the micron-space equivalent.

    Registration back-ends (opencv ECC, phase correlation) work in pixels; the transform
    object is defined in microns. This is the single conversion point, so a mismatch cannot
    creep in from two places disagreeing about pixel size.
    """
    a = np.asarray(affine_px, dtype=np.float64).reshape(2, 3)
    p = grid.pixel_um
    o = np.array([grid.origin_x_um, grid.origin_y_um])

    lin = a[:, :2]
    t_px = a[:, 2]
    # microns -> pixels -> transform -> pixels -> microns
    t_um = p * t_px + o - lin @ o
    return Transform(affine=np.column_stack([lin, t_um]))


def affine_to_pixels(tf: Transform, grid: GridSpec) -> np.ndarray:
    """Inverse of :func:`affine_from_pixels` — a 2x3 in ``grid`` pixel space."""
    p = grid.pixel_um
    o = np.array([grid.origin_x_um, grid.origin_y_um])
    lin = tf.affine[:, :2]
    t_um = tf.affine[:, 2]
    t_px = (t_um - o + lin @ o) / p
    return np.column_stack([lin, t_px])


def fit_affine_lstsq(src_um: np.ndarray, dst_um: np.ndarray) -> Transform:
    """Least-squares affine through matched point pairs. Needs >= 3 non-collinear pairs."""
    src = np.asarray(src_um, dtype=np.float64)
    dst = np.asarray(dst_um, dtype=np.float64)
    if len(src) < 3:
        raise ValueError(f"need at least 3 point pairs to fit an affine, got {len(src)}")
    A = np.column_stack([src, np.ones(len(src))])
    sol, *_ = np.linalg.lstsq(A, dst, rcond=None)
    return Transform(affine=sol.T)


def fit_rigid_umeyama(src_um: np.ndarray, dst_um: np.ndarray,
                      *, allow_scale: bool = True,
                      allow_reflection: bool = False) -> Transform:
    """Closed-form similarity fit (Umeyama). Needs >= 2 pairs.

    Preferred over a free affine when the pairs are few or noisy: a similarity has 4 degrees
    of freedom against an affine's 6, so it cannot absorb a bad correspondence by shearing
    the tissue into agreement.

    ``allow_reflection`` matters here in a way it does not in most applications. Textbook
    Umeyama forces a *proper* rotation (det = +1) because a reflection is usually unphysical.
    Serial histology sections are the exception: a section can be picked up and mounted
    face-down, so the two modalities genuinely differ by a flip. With reflection forbidden the
    fit cannot represent that and silently returns a wrong proper rotation instead — so when
    enabled we compute both the proper and the reflected solution and keep whichever actually
    fits the data better.
    """
    src = np.asarray(src_um, dtype=np.float64)
    dst = np.asarray(dst_um, dtype=np.float64)
    if len(src) < 2:
        raise ValueError(f"need at least 2 point pairs, got {len(src)}")

    mu_s, mu_d = src.mean(axis=0), dst.mean(axis=0)
    sc, dc = src - mu_s, dst - mu_d
    cov = dc.T @ sc / len(src)
    U, D, Vt = np.linalg.svd(cov)
    var = (sc ** 2).sum() / len(src)

    def _build(sign: int) -> Transform:
        S = np.eye(2)
        if sign < 0:
            S[1, 1] = -1
        R = U @ S @ Vt
        scale = 1.0
        if allow_scale and var > 0:
            scale = float((D * np.diag(S)).sum() / var)
        t = mu_d - scale * R @ mu_s
        return Transform(affine=np.column_stack([scale * R, t]))

    proper_sign = 1 if np.linalg.det(U) * np.linalg.det(Vt) >= 0 else -1
    best = _build(proper_sign)
    if not allow_reflection:
        return best

    alt = _build(-proper_sign)
    err_best = float(np.sum((best.apply(src) - dst) ** 2))
    err_alt = float(np.sum((alt.apply(src) - dst) ** 2))
    return alt if err_alt < err_best else best
