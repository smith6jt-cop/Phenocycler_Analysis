"""
S3b — cell tables and morphology images onto a common micron grid.

Registration compares images, but a cell table is not an image. This module builds the
bridge: a shared, explicitly-defined micron grid that both modalities are resampled onto, so
"the same pixel" means "the same physical location" regardless of which instrument produced
the data.

Two kinds of raster, and the pipeline uses both:

**Derived masks** (from the cell table). Always available, no image I/O, no 35 GB TIFF. A
tissue mask plus per-lineage density maps are enough for a coarse rigid/affine alignment,
and for donors where the morphology images are not mounted they are the *only* option.

**Morphology channels** (from the instruments). Much richer, and what the fine and non-rigid
stages want. The DAPI channel is common to both — Xenium's multimodal segmentation kit
stains DAPI plus a boundary protein and an interior RNA channel, and PhenoCycler's qptiff
has DAPI as channel 0 — so DAPI-to-DAPI is the natural correspondence.

The Xenium reader deliberately goes to the **original bundle**, not the zarr. The Xenium
ingest sets ``morphology_focus=False, morphology_mip=False, aligned_images=False`` to skip
"the 35GB morphology TIFFs", so the zarr simply does not contain them. Reading the bundle's
``morphology_focus/`` directly is the only route, and doing it at a pyramid level rather
than full resolution keeps it tractable.

The ``GridSpec`` is the contract: an origin in microns, a pixel size in microns, and a shape.
Everything else in registration works in grid pixels and converts back through it, so there
is exactly one place where the micron/pixel boundary lives.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from .contract import CellTable


@dataclass(frozen=True)
class GridSpec:
    """A micron-anchored raster grid: ``pixel = (micron - origin) / pixel_um``.

    Frozen because a grid silently changing under a transform is the kind of bug that
    produces a plausible-looking, completely wrong overlay.
    """

    origin_x_um: float
    origin_y_um: float
    pixel_um: float
    height: int
    width: int

    @property
    def shape(self) -> tuple[int, int]:
        return (self.height, self.width)

    def to_pixels(self, xy_um: np.ndarray) -> np.ndarray:
        """(n, 2) microns -> (n, 2) float pixel coordinates, (col, row) order."""
        xy = np.asarray(xy_um, dtype=np.float64)
        return np.column_stack([
            (xy[:, 0] - self.origin_x_um) / self.pixel_um,
            (xy[:, 1] - self.origin_y_um) / self.pixel_um,
        ])

    def to_microns(self, xy_px: np.ndarray) -> np.ndarray:
        """(n, 2) pixel coordinates -> (n, 2) microns."""
        xy = np.asarray(xy_px, dtype=np.float64)
        return np.column_stack([
            xy[:, 0] * self.pixel_um + self.origin_x_um,
            xy[:, 1] * self.pixel_um + self.origin_y_um,
        ])

    def as_dict(self) -> dict:
        return {"origin_x_um": self.origin_x_um, "origin_y_um": self.origin_y_um,
                "pixel_um": self.pixel_um, "height": self.height, "width": self.width}

    @classmethod
    def from_dict(cls, d: dict) -> "GridSpec":
        return cls(float(d["origin_x_um"]), float(d["origin_y_um"]), float(d["pixel_um"]),
                   int(d["height"]), int(d["width"]))


def grid_for(
    xy_um: np.ndarray,
    pixel_um: float,
    *,
    margin_um: float = 100.0,
    other: Optional[np.ndarray] = None,
) -> GridSpec:
    """Build a grid covering ``xy_um`` (and optionally a second point set).

    Passing ``other`` makes both modalities share one grid, which is what registration needs:
    a rigid transform between two rasters is only meaningful if the rasters agree on what a
    pixel means and where the origin is.
    """
    pts = [np.asarray(xy_um, dtype=np.float64)]
    if other is not None and len(other):
        pts.append(np.asarray(other, dtype=np.float64))
    allpts = np.vstack(pts)

    x0 = float(np.nanmin(allpts[:, 0])) - margin_um
    y0 = float(np.nanmin(allpts[:, 1])) - margin_um
    x1 = float(np.nanmax(allpts[:, 0])) + margin_um
    y1 = float(np.nanmax(allpts[:, 1])) + margin_um
    width = max(1, int(np.ceil((x1 - x0) / pixel_um)))
    height = max(1, int(np.ceil((y1 - y0) / pixel_um)))
    return GridSpec(x0, y0, float(pixel_um), height, width)


def rasterize_points(
    xy_um: np.ndarray,
    grid: GridSpec,
    *,
    weights: Optional[np.ndarray] = None,
) -> np.ndarray:
    """Point density on the grid, as float32. A plain 2-D histogram, done with bincount.

    Not a synthetic-disc segmentation. The prototype this pipeline replaces stamped a
    radius-10 disc at every centroid and then registered *that* against protein intensity —
    two images with almost no shared structure. A density map at a sensible pixel size
    (a few microns) carries the tissue architecture that registration can actually use,
    without pretending to know cell shape.
    """
    xy = np.asarray(xy_um, dtype=np.float64)
    if len(xy) == 0:
        return np.zeros(grid.shape, dtype=np.float32)

    px = grid.to_pixels(xy)
    col = np.floor(px[:, 0]).astype(np.int64)
    row = np.floor(px[:, 1]).astype(np.int64)
    ok = (col >= 0) & (col < grid.width) & (row >= 0) & (row < grid.height)
    col, row = col[ok], row[ok]

    w = (np.ones(ok.sum(), dtype=np.float64) if weights is None
         else np.asarray(weights, dtype=np.float64)[ok])
    flat = np.bincount(row * grid.width + col, weights=w,
                       minlength=grid.height * grid.width)
    return flat.reshape(grid.shape).astype(np.float32)


def tissue_mask(
    xy_um: np.ndarray,
    grid: GridSpec,
    *,
    close_um: float = 60.0,
    min_area_um2: float = 50_000.0,
) -> np.ndarray:
    """Binary tissue footprint from cell positions.

    Density -> binarise -> morphological close (bridge the gaps between cells) -> fill holes
    -> drop small components. The result is the outline of the section, which is the most
    robust thing to align first: it survives disease state, marker panel and section gap.
    """
    from scipy import ndimage

    dens = rasterize_points(xy_um, grid)
    mask = dens > 0

    radius = max(1, int(round(close_um / grid.pixel_um)))
    struct = _disk(radius)
    mask = ndimage.binary_closing(mask, structure=struct)
    mask = ndimage.binary_fill_holes(mask)

    min_px = max(1, int(min_area_um2 / (grid.pixel_um ** 2)))
    labels, n = ndimage.label(mask)
    if n:
        sizes = np.bincount(labels.ravel())
        sizes[0] = 0
        keep = np.isin(labels, np.flatnonzero(sizes >= min_px))
        mask = keep
    return mask.astype(bool)


def _disk(radius: int) -> np.ndarray:
    r = int(max(1, radius))
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def lineage_density(
    table: CellTable,
    grid: GridSpec,
    lineages: Sequence[str],
    *,
    smooth_um: float = 25.0,
    registered: bool = False,
) -> np.ndarray:
    """(len(lineages), H, W) smoothed per-lineage density stack.

    Gives registration something structural to work with when no images are available: the
    endocrine channel is essentially an islet map, the vascular channel a vessel map. Because
    both modalities are mapped to the *common* vocabulary first, the channels mean the same
    thing on both sides.
    """
    from scipy import ndimage

    xy = table.xy(registered=registered)
    lc = table.df["lineage_common"].to_numpy()
    sigma = max(0.5, smooth_um / grid.pixel_um)

    out = np.zeros((len(lineages), grid.height, grid.width), dtype=np.float32)
    for i, lineage in enumerate(lineages):
        sel = lc == lineage
        if not sel.any():
            continue
        out[i] = ndimage.gaussian_filter(rasterize_points(xy[sel], grid), sigma)
    return out


def rasterize_structures(
    structures: pd.DataFrame,
    grid: GridSpec,
    *,
    weight_by: Optional[str] = "n_cells",
) -> np.ndarray:
    """Structure centroids as a weighted density map — a landmark image.

    Sparser and cleaner than the cell-level density: a few dozen islets rather than a million
    cells, which is exactly what a coarse alignment wants.
    """
    if structures is None or structures.empty:
        return np.zeros(grid.shape, dtype=np.float32)
    xy = structures.loc[:, ["x_um", "y_um"]].to_numpy(np.float64)
    w = (structures[weight_by].to_numpy(np.float64)
         if weight_by and weight_by in structures.columns else None)
    return rasterize_points(xy, grid, weights=w)


# --------------------------------------------------------------------------- #
# Morphology images
# --------------------------------------------------------------------------- #

#: Xenium multimodal segmentation channel order in morphology_focus/, per 10x. Channel 0 is
#: DAPI, which is the one that corresponds to PhenoCycler's nuclear channel.
XENIUM_MORPHOLOGY_CHANNELS = ("DAPI", "boundary", "interior_rna", "interior_protein")


def _pick_pyramid_level(levels_shapes: Sequence[tuple[int, int]], base_um_per_px: float,
                        target_um: float) -> int:
    """Choose the pyramid level closest to ``target_um`` per pixel, never finer.

    Reading level 0 of a whole-slide multi-channel image is what makes naive implementations
    fall over — a 40-channel PhenoCycler slide is tens of GB of uint16. Registration does not
    need or want that resolution; a couple of microns per pixel is plenty and is what
    ``reg_pixel_um`` asks for.
    """
    if not levels_shapes:
        return 0
    h0 = levels_shapes[0][0]
    best, best_err = 0, np.inf
    for i, (h, _w) in enumerate(levels_shapes):
        if h <= 0:
            continue
        um_per_px = base_um_per_px * (h0 / h)
        if um_per_px > target_um * 1.6:      # too coarse; would lose real structure
            continue
        err = abs(um_per_px - target_um)
        if err < best_err:
            best, best_err = i, err
    return best


def read_qptiff_channel(
    path: Path,
    channel: int = 0,
    *,
    target_um: float = 2.0,
    default_um_per_px: float = 0.5,
) -> tuple[np.ndarray, float]:
    """One channel of a PhenoCycler ``.qptiff`` at approximately ``target_um`` per pixel.

    Returns ``(image, um_per_px)``. Channel 0 is DAPI in Akoya's channel order.

    Pyramid-aware and single-channel by construction: it reads the level whose resolution is
    closest to the target, and only the requested channel, so memory stays bounded no matter
    how large the slide is.
    """
    import tifffile

    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        levels = getattr(series, "levels", [series])
        shapes = []
        for lv in levels:
            shp = tuple(lv.shape)
            shapes.append((shp[-2], shp[-1]))

        um_per_px = default_um_per_px
        try:
            page = tif.pages[0]
            tags = page.tags
            if "XResolution" in tags:
                num, den = tags["XResolution"].value
                if num:
                    px_per_unit = num / den
                    unit = tags["ResolutionUnit"].value if "ResolutionUnit" in tags else 3
                    # 3 == centimetre, 2 == inch (TIFF ResolutionUnit)
                    per_cm = px_per_unit if int(unit) == 3 else px_per_unit / 2.54
                    if per_cm > 0:
                        um_per_px = 10_000.0 / per_cm
        except Exception:  # noqa: BLE001 - resolution tags are frequently absent/odd
            pass

        level = _pick_pyramid_level(shapes, um_per_px, target_um)
        arr = levels[level].asarray()

    if arr.ndim == 3:
        n_ch = arr.shape[0]
        if channel >= n_ch:
            raise IndexError(f"{path} has {n_ch} channels; channel {channel} requested")
        arr = arr[channel]
    scale = shapes[0][0] / arr.shape[0]
    return np.asarray(arr, dtype=np.float32), um_per_px * scale


def read_xenium_morphology(
    bundle_dir: Path,
    channel: int = 0,
    *,
    target_um: float = 2.0,
) -> tuple[np.ndarray, float]:
    """One morphology channel from a Xenium bundle at approximately ``target_um`` per pixel.

    Reads ``<bundle>/morphology_focus/`` directly rather than going through the SpatialData
    zarr, because the ingest explicitly excludes morphology images from the zarr. Handles
    both the multi-file layout (``morphology_focus_0000.ome.tif`` ... one per channel, the
    multimodal-segmentation layout) and the older single multi-channel file.

    Xenium morphology is 0.2125 um/px at full resolution.
    """
    import tifffile

    bundle = Path(bundle_dir)
    focus = bundle / "morphology_focus"
    candidates: list[Path] = []
    if focus.is_dir():
        candidates = sorted(focus.glob("morphology_focus_*.ome.tif")) or \
            sorted(focus.glob("*.ome.tif")) or sorted(focus.glob("*.tif"))
    elif focus.with_suffix(".ome.tif").exists():
        candidates = [focus.with_suffix(".ome.tif")]
    if not candidates:
        for name in ("morphology.ome.tif", "morphology_mip.ome.tif", "morphology_focus.ome.tif"):
            p = bundle / name
            if p.exists():
                candidates = [p]
                break
    if not candidates:
        raise FileNotFoundError(
            f"no morphology image under {bundle}. The Xenium SpatialData ingest sets "
            f"morphology_focus=False, so the zarr does not carry it either — point at the "
            f"original output-XETG... bundle."
        )

    base_um = 0.2125
    multi_file = len(candidates) > 1
    path = candidates[channel] if (multi_file and channel < len(candidates)) else candidates[0]

    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        levels = getattr(series, "levels", [series])
        shapes = [(tuple(lv.shape)[-2], tuple(lv.shape)[-1]) for lv in levels]
        level = _pick_pyramid_level(shapes, base_um, target_um)
        arr = levels[level].asarray()

    if arr.ndim == 3 and not multi_file:
        n_ch = arr.shape[0]
        arr = arr[min(channel, n_ch - 1)]
    elif arr.ndim == 3:
        arr = arr[0]

    scale = shapes[0][0] / arr.shape[0]
    return np.asarray(arr, dtype=np.float32), base_um * scale


def image_to_grid(image: np.ndarray, um_per_px: float, grid: GridSpec,
                  *, origin_um: tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    """Resample a morphology image onto the common micron grid.

    ``origin_um`` is where the image's (0, 0) pixel sits in the modality's micron frame —
    zero for both Xenium (whose cell coordinates share the morphology origin) and a
    QuPath-exported qptiff, but explicit so a cropped image can still be placed correctly.
    """
    from scipy import ndimage

    zoom = um_per_px / grid.pixel_um
    if abs(zoom - 1.0) > 1e-6:
        image = ndimage.zoom(image, zoom, order=1)

    out = np.zeros(grid.shape, dtype=np.float32)
    off_col = int(round((origin_um[0] - grid.origin_x_um) / grid.pixel_um))
    off_row = int(round((origin_um[1] - grid.origin_y_um) / grid.pixel_um))

    src_r0, src_c0 = max(0, -off_row), max(0, -off_col)
    dst_r0, dst_c0 = max(0, off_row), max(0, off_col)
    h = min(image.shape[0] - src_r0, grid.height - dst_r0)
    w = min(image.shape[1] - src_c0, grid.width - dst_c0)
    if h > 0 and w > 0:
        out[dst_r0:dst_r0 + h, dst_c0:dst_c0 + w] = image[src_r0:src_r0 + h, src_c0:src_c0 + w]
    return out


def normalize_image(image: np.ndarray, *, low: float = 1.0, high: float = 99.5) -> np.ndarray:
    """Percentile-stretch to [0, 1].

    Percentiles rather than min/max because a single hot pixel (dust, a saturated edge)
    otherwise compresses the whole dynamic range and destroys the correlation registration
    depends on.
    """
    img = np.asarray(image, dtype=np.float32)
    finite = img[np.isfinite(img)]
    if finite.size == 0:
        return np.zeros_like(img)
    lo, hi = np.percentile(finite, [low, high])
    if hi <= lo:
        return np.zeros_like(img)
    return np.clip((img - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def find_qptiff(images_dir: Path, donor_id: str) -> Optional[Path]:
    """Locate a donor's qptiff by the leading-digits convention the pipeline already uses."""
    images_dir = Path(images_dir)
    if not images_dir.exists():
        return None
    hits = sorted(images_dir.glob(f"{donor_id}*.qptiff"))
    if hits:
        return hits[0]
    for p in sorted(images_dir.glob("*.qptiff")):
        m = re.match(r"(\d+)", p.name)
        if m and m.group(1) == str(donor_id):
            return p
    return None
