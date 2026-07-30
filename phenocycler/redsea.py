#!/usr/bin/env python3
"""Scalable pixel-level REDSEA lateral-spillover correction.

The implementation preserves the REDSEA operator while using sparse contacts
and channel-at-a-time qptiff reads:

    data[i,ch]     = whole-cell SUM of channel ch over cell i's pixels   (np.bincount)
    edge[i,ch]     = SUM over cell i's boundary band (1-px inner rim)     (np.bincount)
    contact[i,j]   = # 8-connected pixel adjacencies between cells i,j    (sparse COO)
    F              = donor_normalize(contact)        F[i,j]=contact[i,j]/colsum_j
    corrected_sum  = clip( data + comp_mode*edge - alpha*(F @ edge), 0 )
    corrected_mean = corrected_sum / cell_area_px

The mandatory output resolves that same operator into
Nucleus/Cytoplasm/Membrane/Cell means. Nucleus and Cytoplasm conserve the
whole-cell result exactly before non-negativity clipping. The content-addressed
pipeline is the production entry point; this module contains the numerical
implementation and a thin stage-only command.
"""

from __future__ import annotations

import argparse
import functools
import glob
import json
import os
import re
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import tifffile
from skimage.draw import polygon as draw_polygon
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, diamond

from . import marker_taxonomy
from .config import PipelineConfig, load_config
from .gpu import get_backend
from .parallel import map_donors

# qptiff channel name -> parquet column name (matches upstream SANITIZE)
SANITIZE = lambda s: re.sub(r"[ /\-]", "_", s)


def log(msg):
    print(msg, flush=True)


@dataclass
class RedseaParams:
    """The per-run REDSEA knobs (picklable; carried into worker processes)."""
    downsample: float = 1.0
    edge_radius: int = 0
    comp_mode: int = 1
    alpha: float = 1.0
    gap_bridge: int = 1
    # "donor" == the published Bai et al. operator (mass-conserving); "recipient" == the pre-v10 form.
    norm_form: str = "donor"
    # Skip compensation for markers with no lateral membrane spillover (marker_taxonomy.NO_SPILLOVER).
    exclude_no_spillover: bool = True
    keep_mask: bool = False
    save_intermediates: bool = False
    use_gpu: bool = False
    gpu_device: int = 0

    @classmethod
    def from_config(cls, cfg: PipelineConfig, **overrides) -> "RedseaParams":
        p = cls(
            downsample=cfg.redsea_downsample,
            edge_radius=cfg.redsea_edge_radius,
            comp_mode=cfg.redsea_comp_mode,
            alpha=cfg.redsea_alpha,
            gap_bridge=cfg.redsea_gap_bridge,
            norm_form=cfg.redsea_norm_form,
            exclude_no_spillover=cfg.redsea_exclude_no_spillover,
            use_gpu=cfg.use_gpu,
            gpu_device=cfg.gpu_device,
        )
        for k, v in overrides.items():
            setattr(p, k, v)
        return p


# --------------------------------------------------------------------------- discovery
def contracted_image(cfg: PipelineConfig, donor: str):
    """Return the donor's exact image record from the cohort manifest."""

    from .artifacts import QuPathCohortManifest

    cohort = QuPathCohortManifest.read_json(cfg.qupath_manifest)
    matches = [image for image in cohort.images if image.donor_id == str(donor)]
    if len(matches) != 1:
        raise ValueError(
            f"QuPath cohort contract has {len(matches)} images for donor {donor}"
        )
    return matches[0]


def donor_image(cfg: PipelineConfig, donor: str) -> str:
    """The exact QuPath ``Image`` value contracted for this donor."""

    return contracted_image(cfg, donor).image_id


def resolve_paths(cfg: PipelineConfig, donor: str):
    image = contracted_image(cfg, donor)
    image.validate_current(mode="fast")
    return (
        image,
        Path(image.cell_geojson.path),
        Path(image.nucleus_geojson.path),
        image.combined_geometry_file,
        Path(image.qptiff.path),
    )


def validated_channel_markers(image, channel_names) -> list[str]:
    """Validate qptiff metadata against the manifest and return canonical markers."""

    names = [str(name) for name in channel_names]
    declared = tuple(
        sorted(image.channel_map, key=lambda item: item.channel_index)
    )
    if len(declared) != len(names):
        raise ValueError(
            f"{image.image_id}: manifest declares {len(declared)} channels "
            f"but qptiff contains {len(names)}"
        )
    if tuple(channel.channel_index for channel in declared) != tuple(
        range(len(names))
    ):
        raise ValueError(
            f"{image.image_id}: channel indices do not exactly cover qptiff"
        )
    mismatches = [
        (
            channel.channel_index,
            channel.channel_name,
            names[channel.channel_index],
        )
        for channel in declared
        if channel.channel_name != names[channel.channel_index]
    ]
    if mismatches:
        raise ValueError(
            f"{image.image_id}: qptiff channel metadata differs from manifest "
            f"(examples={mismatches[:5]})"
        )
    return [channel.marker for channel in declared]


# --------------------------------------------------------------------------- qptiff
def qptiff_channels(qptiff: Path, level: int):
    """Return (tiff, level_pages, channel_names, (H,W), pixel_size_um)."""
    tf = tifffile.TiffFile(str(qptiff))
    series = tf.series[0]
    lv = series.levels[level]
    pages = lv.pages
    names = []
    for pg in pages:
        d = pg.description or ""
        m = re.search(r"<Biomarker>([^<]+)</Biomarker>", d) or re.search(r"<Name>([^<]+)</Name>", d)
        names.append(m.group(1).strip() if m else f"ch{len(names)}")
    H, W = lv.shape[-2], lv.shape[-1]
    px_um = None
    try:
        p0 = series.levels[0].pages[0]
        xr = p0.tags["XResolution"].value
        unit = p0.tags["ResolutionUnit"].value  # 3==cm, 2==inch
        ppu = xr[0] / xr[1]
        per_cm = ppu if int(unit) == 3 else ppu / 2.54
        px_um = 1e4 / per_cm * (2 ** level)
    except Exception:
        pass
    return tf, pages, names, (H, W), px_um


def read_channel(pages, ch, dtype=np.float32):
    """One full-res channel as a 2D array (float32 keeps bincount weights light)."""
    return np.asarray(pages[ch].asarray(), dtype=dtype)


# --------------------------------------------------------------------------- rasterize
def _draw_geometry(target, geom, label, shape, s):
    """Rasterize one (Multi)Polygon into `target` at `label`; carve interior rings. Returns #holes."""
    n_holes = 0
    polys = geom["coordinates"] if geom["type"] == "MultiPolygon" else [geom["coordinates"]]
    for poly in polys:                           # poly = [exterior, hole1, ...]
        ext = np.asarray(poly[0], dtype=np.float64)
        if ext.shape[0] < 3:
            continue
        cc = ext[:, 0] * s                        # x -> col
        rr = ext[:, 1] * s                        # y -> row
        pr, pc = draw_polygon(rr, cc, shape=shape)
        target[pr, pc] = label
        for hole in poly[1:]:                     # carve interior rings back to background
            h = np.asarray(hole, dtype=np.float64)
            if h.shape[0] < 3:
                continue
            n_holes += 1
            hr, hc = draw_polygon(h[:, 1] * s, h[:, 0] * s, shape=shape)
            target[hr, hc] = 0
    return n_holes


def rasterize_mask(geojson: Path, shape, downsample: float, with_nucleus: bool = False):
    """Rasterize per-cell polygons into an int32 instance mask. label = feature_index+1.
    Returns (mask, object_ids[list]) where object_ids[label-1] is the cell's UUID.

    ``with_nucleus`` additionally rasterizes each feature's ``nucleusGeometry`` -- which QuPath's
    GeoJSON export already carries alongside the cell polygon, and which this function ignored until
    2026-07-29 -- into a second label image with the SAME labels, and returns
    ``(mask, object_ids, nucleus_mask)``. The nucleus raster is intersected with the final cell raster
    (``nuc == mask`` or 0) so that N is a strict subset of Cell per label: the compartment
    decomposition in :func:`compensate_compartments` needs ``nuc <= data`` and
    ``nuccount <= sizes`` to hold exactly, and cell polygons drawn later can otherwise overwrite an
    earlier cell's pixels while its nucleus raster still claims them. Cells whose nucleusGeometry is
    absent (segmentation fragments; every CANONICAL cell has one) get an empty nucleus and therefore
    a NaN Nucleus mean, never a silent zero."""
    log(f"  [raster] loading {geojson.name} ...")
    feats = json.load(open(geojson))["features"]
    n = len(feats)
    mask = np.zeros(shape, dtype=np.int32)
    nuc_mask = np.zeros(shape, dtype=np.int32) if with_nucleus else None
    oids = [None] * n
    s = 1.0 / downsample
    n_holes = 0
    n_nuc = 0
    t0 = time.time()
    for i, ft in enumerate(feats):
        oids[i] = ft.get("id")
        label = i + 1
        n_holes += _draw_geometry(mask, ft["geometry"], label, shape, s)
        if with_nucleus:
            ng = ft.get("nucleusGeometry")
            if ng is not None:
                n_nuc += 1
                _draw_geometry(nuc_mask, ng, label, shape, s)
        if (i + 1) % 100000 == 0:
            log(f"  [raster] {i+1:,}/{n:,} ({(time.time()-t0):.0f}s)")
    log(f"  [raster] {n:,} cells, {n_holes} holes carved, {time.time()-t0:.0f}s")
    if not with_nucleus:
        return mask, oids
    # Keep only nucleus pixels the cell raster still agrees with (see docstring).
    np.putmask(nuc_mask, nuc_mask != mask, 0)
    log(f"  [raster] nuclei: {n_nuc:,}/{n:,} features carried a nucleusGeometry")
    return mask, oids, nuc_mask


def relabel_nucleus_mask(
    nucleus_mask: np.ndarray,
    nucleus_object_ids,
    cell_mask: np.ndarray,
    cell_object_ids,
) -> np.ndarray:
    """Relabel a separate nucleus raster into the cell raster label space."""

    cell_labels = {
        str(object_id): index
        for index, object_id in enumerate(cell_object_ids, start=1)
    }
    if len(cell_labels) != len(cell_object_ids):
        raise ValueError("cell geometry object IDs are not unique")
    if len(set(map(str, nucleus_object_ids))) != len(nucleus_object_ids):
        raise ValueError("nucleus geometry object IDs are not unique")
    lookup = np.zeros(len(nucleus_object_ids) + 1, dtype=cell_mask.dtype)
    for index, object_id in enumerate(nucleus_object_ids, start=1):
        try:
            lookup[index] = cell_labels[str(object_id)]
        except KeyError as exc:
            raise ValueError(
                f"nucleus object {object_id!r} has no cell geometry"
            ) from exc
    relabelled = lookup[np.asarray(nucleus_mask)]
    return np.where(relabelled == cell_mask, relabelled, 0)


# --------------------------------------------------------------------------- contact matrix
def contact_matrix(mask: np.ndarray, n: int):
    """8-connected cell-cell adjacency counts as a symmetric sparse CSR (NxN)."""
    H, W = mask.shape
    keys = []
    for dy, dx in [(0, 1), (1, 0), (1, 1), (1, -1)]:   # 4 vectors cover all 8 neighbours
        if dx >= 0:
            a = mask[: H - dy, : W - dx]
            b = mask[dy:, dx:]
        else:
            a = mask[: H - dy, -dx:]
            b = mask[dy:, : W + dx]
        v = (a > 0) & (b > 0) & (a != b)
        if not v.any():
            continue
        aa = a[v].astype(np.int64)
        bb = b[v].astype(np.int64)
        lo = np.minimum(aa, bb)
        hi = np.maximum(aa, bb)
        keys.append(lo * (np.int64(n) + 1) + hi)
    if not keys:
        return sp.csr_matrix((n, n), dtype=np.float64)
    keys = np.concatenate(keys)
    uk, cnt = np.unique(keys, return_counts=True)
    lo = (uk // (n + 1)).astype(np.int64)
    hi = (uk % (n + 1)).astype(np.int64)
    M = sp.coo_matrix((cnt.astype(np.float64), (lo - 1, hi - 1)), shape=(n, n)).tocsr()
    M = M + M.T                                    # symmetric
    return M


def contact_matrix_gapbridge(mask, n, grow):
    """Like contact_matrix, but first fills thin (<=grow px) background gaps between cells
    with the nearest cell label, so neighbours separated by a 1-2 px gap still count as
    contacting. The grown mask is used ONLY for the contact graph; sums/edge use the original."""
    if grow <= 0:
        return contact_matrix(mask, n)
    from scipy.ndimage import distance_transform_edt
    bg = mask == 0
    dist, idx = distance_transform_edt(bg, return_distances=True, return_indices=True)
    grown = mask[tuple(idx)]
    grown[bg & (dist > grow + 0.5)] = 0            # only bridge gaps within `grow` px
    M = contact_matrix(grown, n)
    del grown, idx, dist
    return M


# --------------------------------------------------------------------------- compensate
def apply_no_spillover_mask(params, cols, alpha, comp_mode):
    """Zero ``alpha`` and ``comp_mode`` for channels REDSEA must not touch.

    Returns ``(alpha, comp_mode, n_skipped)``. Both are promoted to length-C vectors only when at least
    one channel is excluded, so a run with ``exclude_no_spillover=False`` keeps the scalar fast path and
    stays byte-identical to the pre-v10 behaviour.

    Setting both terms to 0 makes ``corrected_sum = data + 0*edge - 0*(F @ edge) = data`` for that
    channel — a pure passthrough. The column is still written, so every donor
    retains the registry-defined schema.
    """
    if not getattr(params, "exclude_no_spillover", True):
        return alpha, comp_mode, 0
    keep = marker_taxonomy.spillover_corrected_mask(cols)
    n_skipped = int((~keep).sum())
    if n_skipped == 0:
        return alpha, comp_mode, 0
    alpha_vec = np.where(keep, np.broadcast_to(np.asarray(alpha, dtype=np.float64), (len(cols),)), 0.0)
    comp_vec = np.where(keep, np.broadcast_to(np.asarray(comp_mode), (len(cols),)), 0)
    log(f"  [no-spillover] {n_skipped} channel(s) pass through uncorrected: "
        f"{[c for c, k in zip(cols, keep) if not k]}")
    return alpha_vec, comp_vec, n_skipped


COMPARTMENTS: tuple[str, ...] = ("Nucleus", "Cytoplasm", "Membrane", "Cell")


def contact_operator(M, norm_form="donor"):
    """F from the contact matrix, plus the has-neighbours mask. See :func:`compensate` for the forms."""
    if norm_form not in ("donor", "recipient"):
        raise ValueError(f"norm_form must be 'donor' or 'recipient', got {norm_form!r}")
    rowsum = np.asarray(M.sum(1)).ravel()
    inv = np.zeros_like(rowsum)
    nz = rowsum > 0
    inv[nz] = 1.0 / rowsum[nz]
    F = (M @ sp.diags(inv)) if norm_form == "donor" else (sp.diags(inv) @ M)
    return F, nz


def incoming_spillover(edge, M, alpha, norm_form="donor"):
    """``alpha * (F @ edge)`` -- signal received from neighbours, [N,C].

    This is the whole of REDSEA's subtraction term, and it is received at cell i's BOUNDARY BAND:
    ``F @ edge`` redistributes each neighbour's rim sum across the neighbours it touches, and cells
    touch only along their rims. That locality is what makes the compartment decomposition in
    :func:`compensate_compartments` a derivation rather than a modelling choice.

    Returns (sub[N,C], has_neighbours[N])."""
    F, nz = contact_operator(M, norm_form)
    return np.asarray(alpha, dtype=np.float64) * (F @ edge), nz


def compensate_compartments(data, edge, nuc, nuc_edge, sizes, bandcount, nuccount, nucbandcount,
                            M, comp_mode, alpha, norm_form="donor"):
    """Compartment-resolved REDSEA: the same operator as :func:`compensate`, resolved onto
    Nucleus / Cytoplasm / Membrane / Cell instead of collapsed to the whole cell.

    Every term REDSEA adds or removes is a BOUNDARY-BAND quantity -- the ``+comp_mode*edge``
    reinforcement is the cell's own rim sum, and the ``-alpha*(F@edge)`` subtraction is spillover
    received at that rim (see :func:`incoming_spillover`). So the correction apportions across
    compartments by each compartment's share of the band, and nothing has to be assumed about how
    signal distributes inside the cell::

        S_N    = nuc   + comp_mode*nuc_edge  - S_in*(nucbandcount /bandcount)
        S_Y    = cyto  + comp_mode*cyto_edge - S_in*(cytobandcount/bandcount)
        S_B    = edge  + comp_mode*edge      - S_in            # membrane == the whole band
        S_Cell = data  + comp_mode*edge      - S_in

    with ``cyto = data - nuc`` and ``cyto_edge = edge - nuc_edge``. Nucleus and Cytoplasm partition
    the cell (QuPath's own convention -- verified on this cohort: ``Cell`` is reproduced from
    ``Nucleus``/``Cytoplasm`` and their areas at r=0.999993-0.999997, median relative error <=0.01%
    over 332,770 cells of donor 6539), so::

        S_N + S_Y == S_Cell        exactly, pre-clip

    i.e. the decomposition is mass-conserving and reproduces the whole-cell value already written to
    the canonical whole-cell output by construction. NOTE the Membrane band
    overlaps both Nucleus and Cytoplasm; it is not a fourth partition member,
    matching QuPath's overlapping membrane convention.

    Post-clip the identity holds except where the >=0 clamp binds on a compartment but not on the
    whole cell; the ``Cell`` entry is clipped on its own sum, so it stays byte-identical to
    :func:`compensate`.

    Undefined compartments come back NaN, never 0: Cytoplasm where the nucleus fills the cell
    (~9% of cells here, matching QuPath's own Cytoplasm NaN rate) and Nucleus where a fragment
    carried no nucleusGeometry.

    Returns (means: dict[str, ndarray[N,C]], counts: dict[str, ndarray[N]], n_isolated)."""
    sub, nz = incoming_spillover(edge, M, alpha, norm_form=norm_form)
    cyto = data - nuc
    cyto_edge = edge - nuc_edge
    cytocount = sizes - nuccount
    cytobandcount = bandcount - nucbandcount

    # Band shares. A cell with no band pixels has no contacts either, so its `sub` row is 0 and the
    # share it is multiplied by is irrelevant -- 0 keeps it finite.
    safe_band = np.where(bandcount > 0, bandcount, 1.0)[:, None]
    nuc_share = np.where(bandcount[:, None] > 0, nucbandcount[:, None] / safe_band, 0.0)
    cyto_share = np.where(bandcount[:, None] > 0, cytobandcount[:, None] / safe_band, 0.0)

    sums = {
        "Nucleus": nuc + comp_mode * nuc_edge - sub * nuc_share,
        "Cytoplasm": cyto + comp_mode * cyto_edge - sub * cyto_share,
        "Membrane": edge + comp_mode * edge - sub,
        "Cell": data + comp_mode * edge - sub,
    }
    counts = {"Nucleus": nuccount, "Cytoplasm": cytocount,
              "Membrane": bandcount, "Cell": sizes}
    means = {}
    for name in COMPARTMENTS:
        s = np.clip(sums[name], 0, None)
        c = counts[name]
        with np.errstate(invalid="ignore", divide="ignore"):
            m = s / np.where(c > 0, c, np.nan)[:, None]
        m[c <= 0, :] = np.nan
        means[name] = m
    return means, counts, int((~nz).sum())


def compensate(data, edge, sizes, M, comp_mode, alpha, backend=None, norm_form="donor"):
    """REDSEA compensation. corrected_sum = data + comp_mode*edge - alpha*(F @ edge), clamp>=0.
      comp_mode=1 : subtract + reinforce (redseapy default, and ours since 2026-07-25)
      comp_mode=0 : subtract only (remove incoming neighbour spillover; no reinforce)

    ``norm_form`` selects how the contact matrix is normalized into ``F``:
      "donor"     : F[i,j] = contact[i,j] / B[j]   — the share of NEIGHBOUR j's boundary that faces i.
                    This is the published operator (Bai et al. 2021). It is mass-conserving:
                    ``sum_i F[i,j] * edge[j] == edge[j]``, i.e. j's rim signal is redistributed among
                    its neighbours in proportion to shared perimeter and nothing is invented.
      "recipient" : F[i,j] = contact[i,j] / B[i]   — rows sum to 1, so cell i subtracts a contact-
                    weighted AVERAGE of its neighbours' whole rim sums regardless of how much perimeter
                    it actually shares with them. Not mass-conserving; over-subtracts. This was the
                    default through METHOD v9 and is retained only to reproduce those artifacts.

    Why the published form was not obvious: ``redseapy`` builds ``comp*I - M[i,j]/B[j]`` and then
    transposes it (``redsea.py:219``), which looks like it yields the recipient form — but it applies the
    result as ``transpose(dot(transpose(E), P))`` (``redsea.py:294``), i.e. ``P.T @ E = A @ E``, undoing
    the transpose. The net operator is the donor form.

    Since ``M`` is symmetric, ``B = M.sum(0) == M.sum(1)``; the two forms differ whenever neighbouring
    cells have unequal total contact perimeter, and coincide in a regular tiling.

    ``alpha``     : scalar OR length-C vector (per-channel subtraction strength), broadcast over the
                    columns of (F @ edge) [N×C]. A 0 entry disables subtraction for that channel.
    ``comp_mode`` : scalar OR length-C vector, broadcast over the columns of ``edge``. Setting both
                    ``alpha[c] = 0`` and ``comp_mode[c] = 0`` makes channel c a pure passthrough
                    (corrected == raw), which is how markers with no lateral spillover are excluded
                    (see ``marker_taxonomy.NO_SPILLOVER``).
    ``backend`` (optional): a :class:`phenocycler.gpu.Backend`; when on GPU the sparse ``F @ edge``
                matmul runs on the device.
    Returns (corrected_mean[N,C], clamp_frac[C], n_isolated)."""
    if norm_form not in ("donor", "recipient"):
        raise ValueError(f"norm_form must be 'donor' or 'recipient', got {norm_form!r}")

    if backend is not None and backend.on_gpu:  # pragma: no cover - GPU only
        xp = backend.xp
        d = xp.asarray(data)
        e = xp.asarray(edge)
        s = xp.asarray(sizes)
        Mg = backend.csr_from_scipy(M)
        rowsum = xp.asarray(np.asarray(M.sum(1)).ravel())
        inv = xp.zeros_like(rowsum)
        nz = rowsum > 0
        inv[nz] = 1.0 / rowsum[nz]
        F = Mg @ backend.diags(inv) if norm_form == "donor" else backend.diags(inv) @ Mg
        corrected_sum = d + comp_mode * e - alpha * (F @ e)
        xp.clip(corrected_sum, 0, None, out=corrected_sum)
        corrected = corrected_sum / s[:, None]
        corrected = backend.asnumpy(corrected)
        sizes_np = np.asarray(sizes)
        corrected[sizes_np == 0, :] = np.nan
        clamp_frac = backend.asnumpy((corrected_sum == 0).mean(axis=0))
        return corrected, clamp_frac, int(backend.asnumpy((~nz).sum()))

    rowsum = np.asarray(M.sum(1)).ravel()
    inv = np.zeros_like(rowsum)
    nz = rowsum > 0
    inv[nz] = 1.0 / rowsum[nz]
    # donor form: divide column j by neighbour j's total contact perimeter B[j] (published, conserving).
    # recipient form: divide row i by i's own total contact perimeter B[i] (rows sum to 1; legacy).
    F = (M @ sp.diags(inv)) if norm_form == "donor" else (sp.diags(inv) @ M)
    corrected_sum = data + comp_mode * edge - np.asarray(
        alpha, dtype=np.float64
    ) * (F @ edge)
    np.clip(corrected_sum, 0, None, out=corrected_sum)
    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = corrected_sum / sizes[:, None]
    corrected[sizes == 0, :] = np.nan
    clamp_frac = (corrected_sum == 0).mean(axis=0)
    return corrected, clamp_frac, int((~nz).sum())


# --------------------------------------------------------------------------- intermediates / io
def save_intermediates(cfg: PipelineConfig, donor, oids, cols, data, edge, sizes, M, bandcount=None,
                       nuc=None, nuc_edge=None, nuccount=None, nucbandcount=None):
    cfg.inter_dir.mkdir(parents=True, exist_ok=True)
    p = cfg.inter_dir / f"{donor}.npz"
    M = M.tocsr()
    # The compartment decomposition needs band counts to turn corrected sums
    # into membrane means.
    extra = {}
    for name, arr in (("bandcount", bandcount), ("nuccount", nuccount),
                      ("nucbandcount", nucbandcount)):
        if arr is not None:
            extra[name] = np.asarray(arr, dtype=np.float32)
    for name, arr in (("nuc", nuc), ("nuc_edge", nuc_edge)):
        if arr is not None:
            extra[name] = np.asarray(arr, dtype=np.float32)
    np.savez_compressed(p, oids=np.array(oids, dtype=object), cols=np.array(cols),
                        data=data.astype(np.float32), edge=edge.astype(np.float32),
                        sizes=sizes, m_data=M.data, m_indices=M.indices, m_indptr=M.indptr,
                        m_shape=np.array(M.shape), **extra)
    log(f"[{donor}] saved intermediates -> {p}"
        + (f" (+{', '.join(sorted(extra))})" if extra else ""))


def load_intermediates(cfg: PipelineConfig, donor):
    """Returns (oids, cols, data, edge, sizes, M, bandcount, nuc, nuc_edge, nuccount, nucbandcount).

    Missing compartment arrays are reported by the caller as an incompatible
    cache; the current workflow never reconstructs or guesses them."""
    z = np.load(cfg.inter_dir / f"{donor}.npz", allow_pickle=True)
    M = sp.csr_matrix((z["m_data"], z["m_indices"], z["m_indptr"]), shape=tuple(z["m_shape"]))
    get = lambda k: z[k].astype(np.float64) if k in z.files else None
    return (list(z["oids"]), list(z["cols"]), z["data"].astype(np.float64),
            z["edge"].astype(np.float64), z["sizes"].astype(np.float64), M,
            get("bandcount"), get("nuc"), get("nuc_edge"), get("nuccount"), get("nucbandcount"))

def canonical_cell_mask(cfg: PipelineConfig, donor: str, oids) -> np.ndarray:
    """Return the GeoJSON rows that belong to the canonical raw-cell table.

    REDSEA must retain every GeoJSON object while building the mask and contact graph,
    including tiny segmentation fragments excluded from ``data/cells``. Persisting
    those fragments is unsafe, though: they have no canonical morphology or spatial
    metadata and create a larger downstream cell universe. Every canonical raw cell
    must occur exactly once in the GeoJSON; GeoJSON-only rows are deliberately omitted
    only when writing the corrected parquet.
    """
    files = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    if len(files) != 1:
        raise ValueError(
            f"expected one canonical raw parquet for donor {donor}, found {len(files)}"
        )

    canonical_values = pd.read_parquet(files[0], columns=["object_id"])["object_id"]
    geojson_values = pd.Series(oids, dtype="string")
    if canonical_values.isna().any() or geojson_values.isna().any():
        raise ValueError(f"donor {donor}: null object_id in canonical raw or GeoJSON cells")
    canonical = pd.Index(canonical_values.astype(str), name="object_id")
    geojson_ids = pd.Index(geojson_values.astype(str), name="object_id")
    if not canonical.is_unique:
        raise ValueError(f"donor {donor}: duplicate object_id in canonical raw cells")
    if not geojson_ids.is_unique:
        raise ValueError(f"donor {donor}: duplicate object_id in GeoJSON cells")

    missing = canonical.difference(geojson_ids, sort=False)
    if len(missing):
        preview = ", ".join(map(str, missing[:5]))
        raise ValueError(
            f"donor {donor}: {len(missing):,} canonical raw cells are missing from "
            f"the REDSEA GeoJSON ({preview})"
        )

    keep = geojson_ids.isin(canonical)
    n_extra = int((~keep).sum())
    if n_extra:
        log(
            f"[{donor}] output key reconciliation: retaining {int(keep.sum()):,} "
            f"canonical cells; omitting {n_extra:,} GeoJSON-only fragments"
        )
    return np.asarray(keep, dtype=bool)


def apply_spillover_context_eligibility(
    cell_mask: np.ndarray,
    nucleus_mask: np.ndarray,
    geometry_object_ids,
    geometry_qc: pd.DataFrame,
) -> tuple[np.ndarray, np.ndarray, int]:
    """Remove explicitly ineligible polygons from REDSEA physical context.

    GeoJSON-only fragments absent from canonical cells remain available as
    boundary context. For canonical cells, the immutable geometry-QC decision
    is authoritative and is applied before boundary/contact construction.
    """

    required = {"object_id", "spillover_context_eligible"}
    missing = sorted(required.difference(geometry_qc.columns))
    if missing:
        raise ValueError(f"geometry QC is missing columns: {missing}")
    qc = geometry_qc.loc[
        :, ["object_id", "spillover_context_eligible"]
    ].copy()
    if qc["object_id"].isna().any():
        raise ValueError("geometry QC object_id contains null")
    qc["object_id"] = qc["object_id"].astype(str)
    if qc["object_id"].duplicated().any():
        raise ValueError("geometry QC object_id contains duplicates")
    flags = qc["spillover_context_eligible"]
    if flags.isna().any() or not flags.isin([True, False, 0, 1]).all():
        raise ValueError(
            "spillover_context_eligible must contain explicit booleans"
        )

    label_by_id = {
        str(object_id): index
        for index, object_id in enumerate(geometry_object_ids, start=1)
    }
    if len(label_by_id) != len(geometry_object_ids):
        raise ValueError("GeoJSON geometry object IDs are not unique")
    missing_geometry = sorted(set(qc["object_id"]).difference(label_by_id))
    if missing_geometry:
        raise ValueError(
            f"{len(missing_geometry):,} geometry-QC cells are absent from "
            f"the REDSEA raster (examples={missing_geometry[:5]})"
        )
    excluded_labels = np.asarray(
        [
            label_by_id[object_id]
            for object_id, eligible in zip(qc["object_id"], flags.astype(bool))
            if not eligible
        ],
        dtype=np.int64,
    )
    if len(excluded_labels) == 0:
        return cell_mask, nucleus_mask, 0

    keep_label = np.ones(len(geometry_object_ids) + 1, dtype=bool)
    keep_label[excluded_labels] = False
    filtered_cell = np.where(keep_label[cell_mask], cell_mask, 0)
    filtered_nucleus = np.where(
        keep_label[nucleus_mask], nucleus_mask, 0
    )
    return filtered_cell, filtered_nucleus, len(excluded_labels)


def write_corrected_compartments(cfg: PipelineConfig, donor, oids, cols, means, counts):
    """Write the compartment-resolved REDSEA means for one donor.

    Schema: ``object_id, donor_id, {cell,band,nucleus,cytoplasm}_area_px`` then
    ``{Nucleus,Cytoplasm,Membrane,Cell}__<marker>`` for every channel.  This is
    the sole REDSEA dataset; a second whole-cell copy is intentionally not
    written because ``Cell__<marker>`` already is that value."""
    keep = canonical_cell_mask(cfg, str(donor), oids)
    object_ids = np.asarray(oids, dtype=object)[keep]
    from .marker_registry import load_registry

    registry = load_registry(cfg.marker_registry)
    unknown = sorted(set(cols).difference(registry.by_name))
    if unknown:
        raise ValueError(
            f"donor {donor}: qptiff channels absent from marker registry: {unknown}"
        )
    channel_index = {marker: index for index, marker in enumerate(cols)}

    output_columns = {
        "object_id": object_ids,
        "donor_id": np.full(len(object_ids), str(donor), dtype=object),
    }
    for name, key in (("cell", "Cell"), ("band", "Membrane"),
                      ("nucleus", "Nucleus"), ("cytoplasm", "Cytoplasm")):
        output_columns[f"{name}_area_px"] = np.asarray(counts[key])[keep]
    for compartment in COMPARTMENTS:
        block = np.asarray(means[compartment])[keep]
        for spec in registry.markers:
            index = channel_index.get(spec.name)
            output_columns[f"{compartment}__{spec.name}"] = (
                block[:, index]
                if index is not None
                else np.full(len(object_ids), np.nan)
            )
    out = pd.DataFrame(output_columns)

    out_dir = cfg.redsea_dir / f"donor_id={donor}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_dir / "data_0.parquet", index=False)
    log(f"[{donor}] wrote {out_dir/'data_0.parquet'} "
        f"({len(out):,} rows x {len(out.columns)} cols)")


# --------------------------------------------------------------------------- per donor
def process_donor(donor: str, cfg: PipelineConfig, params: RedseaParams):
    t_all = time.time()
    backend = get_backend(params.use_gpu, params.gpu_device)
    (
        image,
        cell_geojson,
        nucleus_geojson,
        combined_geometry,
        qptiff,
    ) = resolve_paths(cfg, donor)
    level = int(round(np.log2(params.downsample)))
    if 2 ** level != int(params.downsample):
        raise SystemExit(f"[err] --downsample must be a power of 2 (got {params.downsample})")
    log(f"[{donor}] image={image.image_id}")

    tf, pages, ch_names, (H, W), px_um = qptiff_channels(qptiff, level)
    cols = validated_channel_markers(image, ch_names)
    log(f"[{donor}] level={level} dims={H}x{W} channels={len(cols)} px_um={px_um} "
        f"gpu={backend.on_gpu}")

    # 1) Rasterize the mandatory cell and nucleus geometries.
    if combined_geometry:
        mask, oids, nuc_mask = rasterize_mask(
            cell_geojson, (H, W), params.downsample, with_nucleus=True
        )
    else:
        mask, oids = rasterize_mask(
            cell_geojson, (H, W), params.downsample
        )
        raw_nucleus_mask, nucleus_ids = rasterize_mask(
            nucleus_geojson, (H, W), params.downsample
        )
        nuc_mask = relabel_nucleus_mask(
            raw_nucleus_mask,
            nucleus_ids,
            mask,
            oids,
        )
    from .expression import read_single_partition

    mask, nuc_mask, n_context_excluded = (
        apply_spillover_context_eligibility(
            mask,
            nuc_mask,
            oids,
            read_single_partition(cfg.geometry_qc_dir, str(donor)),
        )
    )
    if n_context_excluded:
        log(
            f"[{donor}] excluded {n_context_excluded:,} geometry-QC labels "
            "from spillover context"
        )
    n = len(oids)
    sizes = np.bincount(mask.ravel(), minlength=n + 1)[1:].astype(np.float64)  # cell area in px
    n_empty = int((sizes == 0).sum())
    log(f"[{donor}] rasterized {n:,} cells | {n_empty} with 0 px | total fg px={int(sizes.sum()):,}")

    # 2) boundary band (within edge_radius px of a different label), 8-connected
    t0 = time.time()
    inner = find_boundaries(mask, mode="inner", connectivity=2, background=0)
    band = binary_dilation(inner, diamond(params.edge_radius)) & (mask > 0)
    band_idx = np.flatnonzero(band.ravel())
    mask_flat = mask.ravel()
    mask_band = mask_flat[band_idx]
    # The compartment decomposition divides the corrected band sum by this
    # count to obtain a membrane mean.
    bandcount = np.bincount(
        mask_band, minlength=n + 1
    )[1:].astype(np.float64)
    log(f"[{donor}] boundary band: {band_idx.size:,} px "
        f"({100*band_idx.size/mask.size:.1f}% of image), {time.time()-t0:.0f}s")
    del inner, band

    # 2b) nucleus pixel bookkeeping. nuc_mask is already intersected with `mask`, so a band pixel is
    # nuclear iff its nucleus label is nonzero; cytoplasm counts follow by subtraction.
    nuc_flat = nuc_mask.ravel()
    nuc_mask_band = nuc_flat[band_idx]
    nuccount = np.bincount(
        nuc_flat, minlength=n + 1
    )[1:].astype(np.float64)
    nucbandcount = np.bincount(
        nuc_mask_band, minlength=n + 1
    )[1:].astype(np.float64)
    n_nonuc = int(((sizes > 0) & (nuccount <= 0)).sum())
    n_full = int(((sizes > 0) & (nuccount >= sizes)).sum())
    log(
        f"[{donor}] nucleus: {int(nuccount.sum()):,} px | "
        f"{n_nonuc:,} cells with no nucleus | {n_full:,} "
        f"({100*n_full/max(int((sizes>0).sum()),1):.2f}%) whose nucleus "
        "fills the cell (Cytoplasm -> NaN)"
    )
    del nuc_mask

    # 3) per-channel whole-cell sums + boundary sums (streamed one channel at a time), plus the
    #    matching nucleus and nucleus-band sums when compartments are wanted. The channel READ
    #    dominates the loop, so the two extra bincounts cost far less than a second qptiff pass.
    nch = len(cols)
    data = np.zeros((n, nch), dtype=np.float64)
    edge = np.zeros((n, nch), dtype=np.float64)
    nuc = np.zeros((n, nch), dtype=np.float64)
    nuc_edge = np.zeros((n, nch), dtype=np.float64)
    t0 = time.time()
    if backend.on_gpu:  # pragma: no cover - GPU only
        mask_flat_g = backend.asarray(mask_flat)
        mask_band_g = backend.asarray(mask_band)
        band_idx_g = backend.asarray(band_idx)
        nuc_flat_g = backend.asarray(nuc_flat)
        nuc_mask_band_g = backend.asarray(nuc_mask_band)
    for k in range(nch):
        ch = read_channel(pages, k).ravel()
        if backend.on_gpu:  # pragma: no cover - GPU only
            chg = backend.asarray(ch)
            data[:, k] = backend.asnumpy(backend.bincount(mask_flat_g, weights=chg, minlength=n + 1))[1:]
            edge[:, k] = backend.asnumpy(
                backend.bincount(mask_band_g, weights=chg[band_idx_g], minlength=n + 1))[1:]
            nuc[:, k] = backend.asnumpy(
                backend.bincount(nuc_flat_g, weights=chg, minlength=n + 1))[1:]
            nuc_edge[:, k] = backend.asnumpy(
                backend.bincount(
                    nuc_mask_band_g,
                    weights=chg[band_idx_g],
                    minlength=n + 1,
                )
            )[1:]
        else:
            data[:, k] = np.bincount(mask_flat, weights=ch, minlength=n + 1)[1:]
            edge[:, k] = np.bincount(mask_band, weights=ch[band_idx], minlength=n + 1)[1:]
            nuc[:, k] = np.bincount(
                nuc_flat, weights=ch, minlength=n + 1
            )[1:]
            nuc_edge[:, k] = np.bincount(
                nuc_mask_band, weights=ch[band_idx], minlength=n + 1
            )[1:]
        del ch
        if (k + 1) % 10 == 0 or k == nch - 1:
            log(f"[{donor}] channels {k+1}/{nch} ({time.time()-t0:.0f}s)")
    tf.close()

    # 4) contact matrix (gap-bridged: cells here are not gapless)
    t0 = time.time()
    M = contact_matrix_gapbridge(mask, n, params.gap_bridge)
    log(f"[{donor}] contact matrix (gap_bridge={params.gap_bridge}): {M.nnz:,} contacts, "
        f"{time.time()-t0:.0f}s")
    if params.save_intermediates:
        save_intermediates(cfg, donor, oids, cols, data, edge, sizes, M, bandcount=bandcount,
                           nuc=nuc, nuc_edge=nuc_edge, nuccount=nuccount,
                           nucbandcount=nucbandcount)

    # 5) Apply the fixed scalar correction, with registry passthrough channels.
    alpha = params.alpha
    comp_mode = params.comp_mode
    alpha, comp_mode, n_skipped = apply_no_spillover_mask(params, cols, alpha, comp_mode)
    corrected, clamp_frac, n_isolated = compensate(
        data, edge, sizes, M, comp_mode, alpha,
        backend=backend, norm_form=params.norm_form)
    log(f"[{donor}] compensate(mode={params.comp_mode}, norm={params.norm_form}, "
        f"alpha={'vector' if np.ndim(alpha) else alpha}, "
        f"uncorrected_channels={n_skipped}): {n_isolated} isolated cells")

    # 5b) the same operator resolved onto Nucleus/Cytoplasm/Membrane/Cell. Its `Cell` output must
    #     reproduce `compensate` above -- that equality is what makes the two datasets one result.
    means, counts, _ = compensate_compartments(
        data,
        edge,
        nuc,
        nuc_edge,
        sizes,
        bandcount,
        nuccount,
        nucbandcount,
        M,
        comp_mode,
        alpha,
        norm_form=params.norm_form,
    )
    both = np.isfinite(means["Cell"]) & np.isfinite(corrected)
    if not np.array_equal(
        np.isfinite(means["Cell"]), np.isfinite(corrected)
    ) or not np.allclose(
        means["Cell"][both], corrected[both], rtol=1e-9, atol=0.0
    ):
        delta = float(
            np.max(
                np.abs(means["Cell"][both] - corrected[both]),
                initial=0.0,
            )
        )
        raise ValueError(
            f"donor {donor}: compartment Cell does not reproduce the "
            f"whole-cell REDSEA value (max abs delta={delta})"
        )
    log(f"[{donor}] compartments: Cell reproduces whole-cell correction")

    # Alignment is an input contract, not a best-effort diagnostic. A key/channel mismatch here
    # invalidates every corrected value, so fail the donor before any output is written.
    uncorr = np.where(sizes[:, None] > 0, data / np.where(sizes[:, None] == 0, 1, sizes[:, None]), np.nan)
    alignment_check(cfg, donor, oids, cols, uncorr, sizes, px_um)

    # 6) write the sole corrected dataset (one row per object_id).
    write_corrected_compartments(cfg, donor, oids, cols, means, counts)

    if params.keep_mask:
        cfg.mask_dir.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(cfg.mask_dir / f"{donor}.tif", mask, compression="zlib")
        log(f"[{donor}] persisted mask -> {cfg.mask_dir/f'{donor}.tif'}")
    del mask

    log(f"[{donor}] DONE in {time.time()-t_all:.0f}s | top-clamped: " +
        ", ".join(f"{cols[i]}={clamp_frac[i]*100:.1f}%" for i in np.argsort(clamp_frac)[-5:][::-1]))
    return donor


def alignment_check(cfg, donor, oids, cols, uncorr_mean, sizes, px_um):
    """Validate rerasterized means and object keys against the QuPath measurements.

    REDSEA used to print this comparison and continue even when it could not be performed. Production
    now requires complete object-key coverage plus a strong positive relationship for every available
    probe channel. The broad ratio bound catches channel/unit mismatches while allowing the small
    raster-boundary differences expected between QuPath and this rasterizer.
    """
    f = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))[0]
    probe = [c for c in ["DAPI", "INS", "Pan_Cytokeratin", "CD3e", "Vimentin", "CD31"] if c in cols]
    if len(probe) < 2:
        raise ValueError(f"donor {donor}: fewer than two alignment probe channels are available")
    ref = pd.read_parquet(f, columns=["object_id", "cell_area"] + probe)
    mine = pd.DataFrame({"object_id": oids, "_area_px": sizes})
    for c in probe:
        mine[c + "_mine"] = uncorr_mean[:, cols.index(c)]
    j = ref.merge(mine, on="object_id", how="inner")
    if len(j) != len(ref):
        raise ValueError(
            f"donor {donor}: REDSEA/QuPath key alignment matched {len(j):,}/{len(ref):,} cells"
        )
    log(f"[{donor}] alignment vs parquet ({len(j):,}/{len(ref):,} matched on object_id):")
    failures = []
    for c in probe:
        a = j[c].to_numpy(np.float64)
        b = j[c + "_mine"].to_numpy(np.float64)
        ok = np.isfinite(a) & np.isfinite(b)
        r = np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 10 else np.nan
        ratio = np.nanmedian(b[ok] / np.where(a[ok] == 0, np.nan, a[ok]))
        log(f"    {c:16s} pearson_r={r:.3f}  median(mine/parquet)={ratio:.3f}")
        if ok.sum() <= 10 or not np.isfinite(r) or r < 0.80:
            failures.append(f"{c}: r={r:.3f}")
        if not np.isfinite(ratio) or not 0.25 <= ratio <= 4.0:
            failures.append(f"{c}: median ratio={ratio:.3f}")
    if failures:
        raise ValueError(
            f"donor {donor}: REDSEA alignment contract failed ({'; '.join(failures)})"
        )


def recompensate_from_intermediates(donor: str, cfg: PipelineConfig, params: RedseaParams):
    """Reload a complete current cache and apply the configured operator."""
    backend = get_backend(params.use_gpu, params.gpu_device)
    (oids, cols, data, edge, sizes, M, bandcount,
     nuc, nuc_edge, nuccount, nucbandcount) = load_intermediates(cfg, donor)
    mask = None
    if params.gap_bridge > 0 and (cfg.mask_dir / f"{donor}.tif").exists():
        mask = tifffile.imread(cfg.mask_dir / f"{donor}.tif")
        M = contact_matrix_gapbridge(mask, len(oids), params.gap_bridge)
        log(f"[{donor}] recomputed gap-bridged contacts: {M.nnz:,}")
    del mask
    alpha = params.alpha
    comp_mode = params.comp_mode
    alpha, comp_mode, n_skipped = apply_no_spillover_mask(params, cols, alpha, comp_mode)
    corrected, clamp_frac, n_iso = compensate(
        data, edge, sizes, M, comp_mode, alpha,
        backend=backend, norm_form=params.norm_form)
    log(f"[{donor}] recompensate(mode={params.comp_mode}, norm={params.norm_form}, "
        f"alpha={'vector' if np.ndim(alpha) else alpha}, "
        f"uncorrected_channels={n_skipped}): {n_iso} isolated")
    cached = (nuc, nuc_edge, nuccount, nucbandcount, bandcount)
    if any(array is None for array in cached):
        missing = [
            name
            for name, array in zip(
                ("nuc", "nuc_edge", "nuccount", "nucbandcount", "bandcount"),
                cached,
            )
            if array is None
        ]
        raise ValueError(
            f"donor {donor}: intermediates predate canonical compartment "
            f"output (missing {missing}); rerun the full REDSEA pass"
        )
    means, counts, _ = compensate_compartments(
        data,
        edge,
        nuc,
        nuc_edge,
        sizes,
        bandcount,
        nuccount,
        nucbandcount,
        M,
        comp_mode,
        alpha,
        norm_form=params.norm_form,
    )
    write_corrected_compartments(cfg, donor, oids, cols, means, counts)
    return donor


# --------------------------------------------------------------------------- driver
def run_redsea(cfg: PipelineConfig, donors: list[str], params: RedseaParams,
               *, from_intermediates: bool = False, n_jobs: int = 1) -> list:
    """Run REDSEA over a list of donors, optionally across processes."""
    from .cohort import ensure_eligible_donors

    donors = list(ensure_eligible_donors(donors, context="REDSEA analysis"))
    worker = recompensate_from_intermediates if from_intermediates else process_donor
    fn = functools.partial(worker, cfg=cfg, params=params)
    results = map_donors(fn, donors, n_jobs=n_jobs, ordered=True, on_error="raise")
    if results != donors:
        raise RuntimeError(
            f"REDSEA returned an incomplete donor set: expected {donors}, got {results}"
        )
    return results


def main(argv=None):
    """Run only the manifest-controlled REDSEA stage.

    Scientific parameters come from the content-addressed pipeline config; the
    numerical module no longer exposes an independent experimental CLI.
    """

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--jobs", type=int)
    args = parser.parse_args(argv)
    command = ["run"]
    if args.config is not None:
        command.extend(["--config", str(args.config)])
    if args.jobs is not None:
        command.extend(["--jobs", str(args.jobs)])
    command.extend(["--only", "redsea", "--no-export"])
    from .pipeline import main as pipeline_main

    return pipeline_main(command)


if __name__ == "__main__":
    raise SystemExit(main())
