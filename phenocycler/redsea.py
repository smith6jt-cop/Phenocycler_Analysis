#!/usr/bin/env python3
"""
Step 2 — scalable pixel-level REDSEA (lateral spillover correction).

Faithful port of ``scripts/senior/redsea_full.py`` (Islet-Explorer-Senior).  It
reimplements the REDSEA math (Bai et al. 2021, *Front. Immunol.*) scalably — the
reference ``redseapy`` does not scale (dense cellNum² matrix, whole-image reads,
pure-python pixel loops):

    data[i,ch]     = whole-cell SUM of channel ch over cell i's pixels   (np.bincount)
    edge[i,ch]     = SUM over cell i's boundary band (1-px inner rim)     (np.bincount)
    contact[i,j]   = # 8-connected pixel adjacencies between cells i,j    (sparse COO)
    F              = donor_normalize(contact)        F[i,j]=contact[i,j]/colsum_j
    corrected_sum  = clip( data + comp_mode*edge - alpha*(F @ edge), 0 )
    corrected_mean = corrected_sum / cell_area_px

Defaults (METHOD v10, 2026-07-25): ``comp_mode=1`` (subtract + reinforce, the Bai et al. default),
``norm_form="donor"`` (the published mass-conserving operator), ``alpha=1.0``, ``edge_radius=0``
(the 1-px cell rim), ``gap_bridge=1``, and nuclear/ECM channels passed through uncorrected
(``marker_taxonomy.NO_SPILLOVER``).  The pre-v10 defaults — ``comp_mode=0``, ``norm_form="recipient"``,
every channel corrected — are still reachable by config/CLI to reproduce v9 artifacts; see
``compensate`` for why they over-subtract.  Mask source is per-cell GeoJSON exported from QuPath;
the qptiff is read one channel at a time.

Additions vs. upstream:
  * config-driven paths (no hardcoded ``/home/smith6jt/...``);
  * per-donor multiprocessing (``--jobs`` / ``n_jobs``) — every donor is
    independent, so this scales near-linearly and is the primary speed-up;
  * an optional CuPy backend (``--gpu`` / ``use_gpu``) for the per-channel
    ``bincount`` accumulation and the sparse ``F @ edge`` matmul.  GPU is meant
    for ``n_jobs=1`` (one donor on the device at a time); combine GPU with a
    SLURM array for cohort-scale parallelism instead of a local process pool.

    python -m phenocycler.redsea --donor 6539            # smoke test (smallest)
    python -m phenocycler.redsea --all --jobs 4          # 4 donors in parallel
    python -m phenocycler.redsea --all --gpu             # CuPy accumulation/matmul
    python -m phenocycler.redsea --all --reconcile-existing
                                                            # filter existing outputs to raw-cell IDs
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
from typing import Optional

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
    # -- optional per-channel alpha + receptivity ("brightness") gate: the abandoned "keratin
    #    lever", kept ADDITIVE and OFF by default.  When these stay at their defaults the
    #    compensation is bit-identical to scalar-alpha REDSEA. --
    alpha_per_channel: Optional[str] = None   # 'Name=val,Name2=val2'; unlisted channels use `alpha`
    gate_brightness: bool = False             # enable the receptivity gate
    gate_channels: Optional[str] = None       # comma-separated aggressor channels (default keratins+Vim)
    gate_kappa: float = 1.0                   # gate aggressiveness

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
def donor_image(cfg: PipelineConfig, donor: str) -> str:
    """The QuPath image string for this donor (from the cells parquet)."""
    f = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    if not f:
        raise SystemExit(f"[err] no {cfg.cells_dir}/donor_id={donor}")
    return str(pd.read_parquet(f[0], columns=["image"]).image.iloc[0])


def resolve_paths(cfg: PipelineConfig, donor: str):
    img = donor_image(cfg, donor)
    tag = re.sub(r"[^A-Za-z0-9._-]", "_", img)          # matches the Groovy sanitization
    geojson = cfg.geojson_dir / f"cells__{tag}.geojson"
    qptiff = cfg.images_dir / img.split(" - ")[0].strip()  # '...qptiff - resolution #1' -> '...qptiff'
    if not geojson.exists():
        raise SystemExit(f"[err] missing geojson: {geojson}\n      run export_cells_geojson.groovy for {img}")
    if not qptiff.exists():
        raise SystemExit(f"[err] missing qptiff: {qptiff}")
    return img, geojson, qptiff


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
def rasterize_mask(geojson: Path, shape, downsample: float):
    """Rasterize per-cell polygons into an int32 instance mask. label = feature_index+1.
    Returns (mask, object_ids[list]) where object_ids[label-1] is the cell's UUID."""
    log(f"  [raster] loading {geojson.name} ...")
    feats = json.load(open(geojson))["features"]
    n = len(feats)
    mask = np.zeros(shape, dtype=np.int32)
    oids = [None] * n
    s = 1.0 / downsample
    n_holes = 0
    t0 = time.time()
    for i, ft in enumerate(feats):
        oids[i] = ft.get("id")
        label = i + 1
        g = ft["geometry"]
        gtype = g["type"]
        polys = g["coordinates"] if gtype == "MultiPolygon" else [g["coordinates"]]
        for poly in polys:                       # poly = [exterior, hole1, ...]
            ext = np.asarray(poly[0], dtype=np.float64)
            if ext.shape[0] < 3:
                continue
            cc = ext[:, 0] * s                    # x -> col
            rr = ext[:, 1] * s                    # y -> row
            pr, pc = draw_polygon(rr, cc, shape=shape)
            mask[pr, pc] = label
            for hole in poly[1:]:                 # carve interior rings back to background
                h = np.asarray(hole, dtype=np.float64)
                if h.shape[0] < 3:
                    continue
                n_holes += 1
                hr, hc = draw_polygon(h[:, 1] * s, h[:, 0] * s, shape=shape)
                mask[hr, hc] = 0
        if (i + 1) % 100000 == 0:
            log(f"  [raster] {i+1:,}/{n:,} ({(time.time()-t0):.0f}s)")
    log(f"  [raster] {n:,} cells, {n_holes} holes carved, {time.time()-t0:.0f}s")
    return mask, oids


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
DEFAULT_GATE_CHANNELS = ["Pan_Cytokeratin", "Ker8_18", "Keratin_5", "Vimentin"]


def parse_alpha_per_channel(spec, cols, default):
    """Build a length-C per-channel alpha vector from 'Name=val,Name2=val2' (unlisted -> `default`).
    Warns on any name not in `cols`. Returns None if spec is falsy (caller uses the scalar `default`)."""
    if not spec:
        return None
    a = np.full(len(cols), float(default), dtype=np.float64)
    ci = {c: k for k, c in enumerate(cols)}
    for tok in str(spec).split(","):
        tok = tok.strip()
        if not tok:
            continue
        name, _, val = tok.partition("=")
        name = name.strip()
        if name not in ci:
            log(f"  [alpha] WARNING: channel '{name}' not in the panel — ignored")
            continue
        a[ci[name]] = float(val)
    return a


def apply_no_spillover_mask(params, cols, alpha, comp_mode):
    """Zero ``alpha`` and ``comp_mode`` for channels REDSEA must not touch.

    Returns ``(alpha, comp_mode, n_skipped)``. Both are promoted to length-C vectors only when at least
    one channel is excluded, so a run with ``exclude_no_spillover=False`` keeps the scalar fast path and
    stays byte-identical to the pre-v10 behaviour.

    Setting both terms to 0 makes ``corrected_sum = data + 0*edge - 0*(F @ edge) = data`` for that
    channel — a pure passthrough. The column is still WRITTEN, so the parquet schema is unchanged
    (``load_validation_sample`` reads every ``QUPATH_MARKERS`` entry and ``reconcile_existing_donor``
    compares Arrow schemas for equality).
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


def build_gate(params, cols):
    """Config for the receptivity (directional) gate, or None when disabled. `channels` is a bool[C]
    mask of the aggressor channels to gate (default the keratins + Vimentin); other channels always
    subtract fully (w=1).  `params` is a :class:`RedseaParams` (or any object exposing the gate_*
    attributes)."""
    if not getattr(params, "gate_brightness", False):
        return None
    want = (params.gate_channels.split(",") if getattr(params, "gate_channels", None)
            else DEFAULT_GATE_CHANNELS)
    want = {c.strip() for c in want if c.strip()}
    chan = np.array([c in want for c in cols], dtype=bool)
    missing = want - set(cols)
    if missing:
        log(f"  [gate] WARNING: gate channels not in panel — ignored: {sorted(missing)}")
    log(f"  [gate] receptivity gate ON: kappa={params.gate_kappa} "
        f"channels={[c for c in cols if c in want]}")
    return {"channels": chan, "kappa": float(params.gate_kappa)}


def receptivity_weight(data, edge, sizes, bandcount, gate):
    """Per-cell/per-channel weight w[i,c] in [0,1] scaling the spillover subtraction, from each cell's
    OWN spillover signature — no global reference (self-calibrating per cell).

    A spillover VICTIM has its signal concentrated at the shared boundary: the 1-px band is much brighter
    than the interior (neighbour keratin bleeding in). A real SOURCE is uniform: band ~= interior. So the
    fractional excess of band-over-interior density is the receptivity:
        band_density  = edge / bandcount                        (mean channel in the boundary band)
        int_density   = (data - edge) / (sizes - bandcount)     (mean channel in the interior)
        excess        = clip((band_density - int_density) / band_density, 0, 1)     # 0 uniform .. 1 rim-only
        w             = clip(kappa * excess, 0, 1)
    Victims (band >> interior) -> w->1 -> subtract the neighbour spillover; real acinar/vessel/fibroblast
    cells (uniform bright) -> w->0 -> PROTECTED. `kappa` scales aggressiveness (higher cleans more). Cells
    with no band (isolated) already have F-row 0, so their subtraction is 0 regardless. Non-gated channels
    get w=1 (unchanged full subtraction)."""
    if bandcount is None:
        raise ValueError("receptivity gate needs per-cell bandcount — regenerate intermediates with the "
                         "current code (--save-intermediates) or run with --keep-mask so it recomputes")
    bc = np.maximum(bandcount, 1.0)[:, None]
    ia = np.maximum(sizes - bandcount, 1.0)[:, None]
    band_density = edge / bc
    int_density = np.clip((data - edge) / ia, 0, None)
    with np.errstate(invalid="ignore", divide="ignore"):
        excess = np.where(band_density > 0, (band_density - int_density) / band_density, 0.0)
    excess = np.clip(excess, 0.0, 1.0)
    w = np.clip(gate["kappa"] * excess, 0.0, 1.0)
    return np.where(gate["channels"][None, :], w, 1.0)


def compensate(data, edge, sizes, M, comp_mode, alpha, bandcount=None, gate=None, backend=None,
               norm_form="donor"):
    """REDSEA compensation. corrected_sum = data + comp_mode*edge - (alpha [* w]) * (F @ edge), clamp>=0.
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
    ``gate``      : None (w=1) or a receptivity-gate config from :func:`build_gate` (needs ``bandcount``)
                    that protects bright-interior cells channel-wise.
    ``backend`` (optional): a :class:`phenocycler.gpu.Backend`; when on GPU the sparse ``F @ edge``
                matmul runs on the device. The GPU fast path is taken ONLY for a scalar-alpha, ungated
                call; gated / per-channel-alpha calls fall through to the CPU branch (correct, since
                data/edge/sizes/M are already NumPy there).

    A scalar-``alpha``, ungated, ``norm_form="recipient"`` call on CPU is BYTE-IDENTICAL to the original
    single-alpha REDSEA — it executes the exact same expression.
    Returns (corrected_mean[N,C], clamp_frac[C], n_isolated)."""
    if norm_form not in ("donor", "recipient"):
        raise ValueError(f"norm_form must be 'donor' or 'recipient', got {norm_form!r}")

    if (backend is not None and backend.on_gpu
            and gate is None and np.ndim(alpha) == 0):  # pragma: no cover - GPU only
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
    if gate is None and np.ndim(alpha) == 0:
        # DEFAULT path — identical expression to the original single-alpha REDSEA (bit-for-bit).
        corrected_sum = data + comp_mode * edge - alpha * (F @ edge)
    else:
        # opt-in "keratin lever": per-channel alpha vector and/or the receptivity gate.
        sub = np.asarray(alpha, dtype=np.float64) * (F @ edge)   # scalar or (C,) broadcasts over cols
        if gate is not None:
            sub = sub * receptivity_weight(data, edge, sizes, bandcount, gate)
        corrected_sum = data + comp_mode * edge - sub
    np.clip(corrected_sum, 0, None, out=corrected_sum)
    with np.errstate(invalid="ignore", divide="ignore"):
        corrected = corrected_sum / sizes[:, None]
    corrected[sizes == 0, :] = np.nan
    clamp_frac = (corrected_sum == 0).mean(axis=0)
    return corrected, clamp_frac, int((~nz).sum())


# --------------------------------------------------------------------------- intermediates / io
def save_intermediates(cfg: PipelineConfig, donor, oids, cols, data, edge, sizes, M, bandcount=None):
    cfg.inter_dir.mkdir(parents=True, exist_ok=True)
    p = cfg.inter_dir / f"{donor}.npz"
    M = M.tocsr()
    # bandcount is persisted ONLY when the receptivity gate is on (it is None otherwise), so a
    # default run writes exactly the same .npz keys as before.
    extra = {} if bandcount is None else {"bandcount": np.asarray(bandcount, dtype=np.float32)}
    np.savez_compressed(p, oids=np.array(oids, dtype=object), cols=np.array(cols),
                        data=data.astype(np.float32), edge=edge.astype(np.float32),
                        sizes=sizes, m_data=M.data, m_indices=M.indices, m_indptr=M.indptr,
                        m_shape=np.array(M.shape), **extra)
    log(f"[{donor}] saved intermediates -> {p}" + ("" if bandcount is None else " (+bandcount)"))


def load_intermediates(cfg: PipelineConfig, donor):
    """Returns (oids, cols, data, edge, sizes, M, bandcount). `bandcount` is None for .npz files
    written without the gate (the default), in which case a gated re-run recomputes it from the mask."""
    z = np.load(cfg.inter_dir / f"{donor}.npz", allow_pickle=True)
    M = sp.csr_matrix((z["m_data"], z["m_indices"], z["m_indptr"]), shape=tuple(z["m_shape"]))
    bandcount = z["bandcount"].astype(np.float64) if "bandcount" in z.files else None
    return (list(z["oids"]), list(z["cols"]), z["data"].astype(np.float64),
            z["edge"].astype(np.float64), z["sizes"].astype(np.float64), M, bandcount)


def bandcount_from_mask(mask, n, edge_radius):
    """Recompute per-cell boundary-band pixel counts from a persisted mask (for old intermediates that
    predate bandcount). Mirrors the band construction in process_donor."""
    inner = find_boundaries(mask, mode="inner", connectivity=2, background=0)
    band = binary_dilation(inner, diamond(edge_radius)) & (mask > 0)
    return np.bincount(mask.ravel()[np.flatnonzero(band.ravel())], minlength=n + 1)[1:].astype(np.float64)


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


def write_corrected(cfg: PipelineConfig, donor, oids, cols, corrected, sizes):
    keep = canonical_cell_mask(cfg, str(donor), oids)
    object_ids = np.asarray(oids, dtype=object)[keep]
    corrected = np.asarray(corrected)[keep]
    sizes = np.asarray(sizes)[keep]

    out = pd.DataFrame({"object_id": object_ids, "donor_id": donor, "cell_area_px": sizes})
    for k, c in enumerate(cols):
        out[c] = corrected[:, k]
    out_dir = cfg.cells_redsea_dir / f"donor_id={donor}"
    out_dir.mkdir(parents=True, exist_ok=True)
    out.to_parquet(out_dir / "data_0.parquet", index=False)
    log(f"[{donor}] wrote {out_dir/'data_0.parquet'} ({len(out):,} rows)")


def reconcile_existing_donor(donor: str, cfg: PipelineConfig) -> str:
    """Atomically restrict an existing REDSEA parquet to canonical raw object IDs."""
    import duckdb
    import pyarrow.parquet as pq

    raw_files = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    redsea_files = sorted(
        glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet"))
    )
    if len(raw_files) != 1 or len(redsea_files) != 1:
        raise ValueError(
            f"donor {donor}: expected one raw and one REDSEA parquet; found "
            f"{len(raw_files)} raw and {len(redsea_files)} REDSEA"
        )

    raw_path = Path(raw_files[0])
    redsea_path = Path(redsea_files[0])
    temp_path = redsea_path.with_name(f".{redsea_path.stem}.reconciled.parquet")
    temp_path.unlink(missing_ok=True)
    original_schema = pq.ParquetFile(redsea_path).schema_arrow

    def sql_path(path: Path) -> str:
        return str(path).replace("'", "''")

    con = duckdb.connect()
    try:
        raw_q = sql_path(raw_path)
        redsea_q = sql_path(redsea_path)
        temp_q = sql_path(temp_path)
        stats = con.execute(
            f"""
            WITH raw AS (
              SELECT CAST(object_id AS VARCHAR) AS object_id
              FROM read_parquet('{raw_q}', hive_partitioning=false)
            ), redsea AS (
              SELECT CAST(object_id AS VARCHAR) AS object_id
              FROM read_parquet('{redsea_q}', hive_partitioning=false)
            )
            SELECT
              (SELECT count(*) FROM raw),
              (SELECT count(DISTINCT object_id) FROM raw),
              (SELECT count(*) FROM raw WHERE object_id IS NULL),
              (SELECT count(*) FROM redsea),
              (SELECT count(DISTINCT object_id) FROM redsea),
              (SELECT count(*) FROM redsea WHERE object_id IS NULL),
              (SELECT count(*) FROM raw ANTI JOIN redsea USING (object_id))
            """
        ).fetchone()
        raw_n, raw_unique, raw_null, redsea_n, redsea_unique, redsea_null, raw_missing = stats
        if raw_null or redsea_null or raw_unique != raw_n or redsea_unique != redsea_n:
            raise ValueError(
                f"donor {donor}: object_id must be non-null and unique "
                f"(raw {raw_unique:,}/{raw_n:,}, REDSEA {redsea_unique:,}/{redsea_n:,})"
            )
        if raw_missing:
            raise ValueError(
                f"donor {donor}: {raw_missing:,} canonical raw cells are missing from REDSEA"
            )
        if raw_n == redsea_n:
            log(f"[{donor}] existing REDSEA keys already match {raw_n:,} canonical cells")
            return donor

        con.execute("SET preserve_insertion_order=true")
        con.execute(
            f"""
            COPY (
              SELECT redsea.*
              FROM read_parquet('{redsea_q}', hive_partitioning=false) AS redsea
              SEMI JOIN (
                SELECT CAST(object_id AS VARCHAR) AS object_id
                FROM read_parquet('{raw_q}', hive_partitioning=false)
              ) AS raw
              ON CAST(redsea.object_id AS VARCHAR) = raw.object_id
            ) TO '{temp_q}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )
        reconciled_n, reconciled_unique = con.execute(
            f"""
            SELECT count(*), count(DISTINCT CAST(object_id AS VARCHAR))
            FROM read_parquet('{temp_q}', hive_partitioning=false)
            """
        ).fetchone()
        if reconciled_n != raw_n or reconciled_unique != raw_n:
            raise ValueError(
                f"donor {donor}: reconciled REDSEA has {reconciled_n:,} rows / "
                f"{reconciled_unique:,} IDs; expected {raw_n:,}"
            )
        if not pq.ParquetFile(temp_path).schema_arrow.equals(original_schema):
            raise ValueError(f"donor {donor}: reconciliation changed the REDSEA schema")
        os.replace(temp_path, redsea_path)
        log(
            f"[{donor}] reconciled existing REDSEA: {redsea_n:,} -> {raw_n:,} "
            f"canonical rows"
        )
        return donor
    finally:
        con.close()
        temp_path.unlink(missing_ok=True)


def reconcile_existing(cfg: PipelineConfig, donors: list[str], *, n_jobs: int = 1) -> list:
    """Reconcile retained REDSEA outputs in parallel, one atomic file replacement per donor."""
    from .cohort import ensure_eligible_donors

    donors = list(
        ensure_eligible_donors(donors, context="REDSEA reconciliation")
    )
    fn = functools.partial(reconcile_existing_donor, cfg=cfg)
    return map_donors(fn, donors, n_jobs=n_jobs, ordered=True, on_error="raise")


# --------------------------------------------------------------------------- per donor
def process_donor(donor: str, cfg: PipelineConfig, params: RedseaParams):
    t_all = time.time()
    backend = get_backend(params.use_gpu, params.gpu_device)
    img, geojson, qptiff = resolve_paths(cfg, donor)
    level = int(round(np.log2(params.downsample)))
    if 2 ** level != int(params.downsample):
        raise SystemExit(f"[err] --downsample must be a power of 2 (got {params.downsample})")
    log(f"[{donor}] image={img}")

    tf, pages, ch_names, (H, W), px_um = qptiff_channels(qptiff, level)
    cols = [SANITIZE(c) for c in ch_names]
    gate = build_gate(params, cols)          # None unless --gate-brightness (default: off)
    log(f"[{donor}] level={level} dims={H}x{W} channels={len(cols)} px_um={px_um} "
        f"gpu={backend.on_gpu}")

    # 1) rasterize instance mask
    mask, oids = rasterize_mask(geojson, (H, W), params.downsample)
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
    # per-cell band-pixel count — only the receptivity gate needs it, so compute it only then
    bandcount = None
    if gate is not None:
        bandcount = np.bincount(mask_band, minlength=n + 1)[1:].astype(np.float64)
    log(f"[{donor}] boundary band: {band_idx.size:,} px "
        f"({100*band_idx.size/mask.size:.1f}% of image), {time.time()-t0:.0f}s")
    del inner, band

    # 3) per-channel whole-cell sums + boundary sums (streamed one channel at a time)
    nch = len(cols)
    data = np.zeros((n, nch), dtype=np.float64)
    edge = np.zeros((n, nch), dtype=np.float64)
    t0 = time.time()
    if backend.on_gpu:  # pragma: no cover - GPU only
        mask_flat_g = backend.asarray(mask_flat)
        mask_band_g = backend.asarray(mask_band)
        band_idx_g = backend.asarray(band_idx)
    for k in range(nch):
        ch = read_channel(pages, k).ravel()
        if backend.on_gpu:  # pragma: no cover - GPU only
            chg = backend.asarray(ch)
            data[:, k] = backend.asnumpy(backend.bincount(mask_flat_g, weights=chg, minlength=n + 1))[1:]
            edge[:, k] = backend.asnumpy(
                backend.bincount(mask_band_g, weights=chg[band_idx_g], minlength=n + 1))[1:]
        else:
            data[:, k] = np.bincount(mask_flat, weights=ch, minlength=n + 1)[1:]
            edge[:, k] = np.bincount(mask_band, weights=ch[band_idx], minlength=n + 1)[1:]
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
        save_intermediates(cfg, donor, oids, cols, data, edge, sizes, M, bandcount=bandcount)

    # 5) compensate (scalar alpha + no gate by default; optional per-channel vector + brightness gate)
    alpha = parse_alpha_per_channel(params.alpha_per_channel, cols, params.alpha)
    alpha = params.alpha if alpha is None else alpha
    comp_mode = params.comp_mode
    alpha, comp_mode, n_skipped = apply_no_spillover_mask(params, cols, alpha, comp_mode)
    corrected, clamp_frac, n_isolated = compensate(
        data, edge, sizes, M, comp_mode, alpha,
        bandcount=bandcount, gate=gate, backend=backend, norm_form=params.norm_form)
    log(f"[{donor}] compensate(mode={params.comp_mode}, norm={params.norm_form}, "
        f"alpha={'vector' if np.ndim(alpha) else alpha}, gate={'on' if gate else 'off'}, "
        f"uncorrected_channels={n_skipped}): {n_isolated} isolated cells")

    # alignment sanity vs the existing parquet (correlation of UNcorrected mean)
    uncorr = np.where(sizes[:, None] > 0, data / np.where(sizes[:, None] == 0, 1, sizes[:, None]), np.nan)
    try:
        alignment_check(cfg, donor, oids, cols, uncorr, sizes, px_um)
    except Exception as e:  # non-fatal sanity check
        log(f"[{donor}] alignment check skipped: {e}")

    # 6) write corrected parquet (one row per object_id)
    write_corrected(cfg, donor, oids, cols, corrected, sizes)

    if params.keep_mask:
        cfg.mask_dir.mkdir(parents=True, exist_ok=True)
        tifffile.imwrite(cfg.mask_dir / f"{donor}.tif", mask, compression="zlib")
        log(f"[{donor}] persisted mask -> {cfg.mask_dir/f'{donor}.tif'}")
    del mask

    log(f"[{donor}] DONE in {time.time()-t_all:.0f}s | top-clamped: " +
        ", ".join(f"{cols[i]}={clamp_frac[i]*100:.1f}%" for i in np.argsort(clamp_frac)[-5:][::-1]))
    return donor


def alignment_check(cfg, donor, oids, cols, uncorr_mean, sizes, px_um):
    """Correlate my UNcorrected full-res mean against the existing QuPath parquet means."""
    f = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))[0]
    probe = [c for c in ["DAPI", "INS", "Pan_Cytokeratin", "CD3e", "Vimentin", "CD31"] if c in cols]
    ref = pd.read_parquet(f, columns=["object_id", "cell_area"] + probe)
    mine = pd.DataFrame({"object_id": oids, "_area_px": sizes})
    for c in probe:
        mine[c + "_mine"] = uncorr_mean[:, cols.index(c)]
    j = ref.merge(mine, on="object_id", how="inner")
    log(f"[{donor}] alignment vs parquet ({len(j):,}/{len(ref):,} matched on object_id):")
    for c in probe:
        a = j[c].to_numpy(np.float64)
        b = j[c + "_mine"].to_numpy(np.float64)
        ok = np.isfinite(a) & np.isfinite(b)
        r = np.corrcoef(a[ok], b[ok])[0, 1] if ok.sum() > 10 else np.nan
        ratio = np.nanmedian(b[ok] / np.where(a[ok] == 0, np.nan, a[ok]))
        log(f"    {c:16s} pearson_r={r:.3f}  median(mine/parquet)={ratio:.3f}")


def recompensate_from_intermediates(donor: str, cfg: PipelineConfig, params: RedseaParams):
    """Reload cached intermediates and just (re)compensate + write (fast alpha sweeps)."""
    backend = get_backend(params.use_gpu, params.gpu_device)
    oids, cols, data, edge, sizes, M, bandcount = load_intermediates(cfg, donor)
    mask = None
    if params.gap_bridge > 0 and (cfg.mask_dir / f"{donor}.tif").exists():
        mask = tifffile.imread(cfg.mask_dir / f"{donor}.tif")
        M = contact_matrix_gapbridge(mask, len(oids), params.gap_bridge)
        log(f"[{donor}] recomputed gap-bridged contacts: {M.nnz:,}")
    alpha = parse_alpha_per_channel(params.alpha_per_channel, cols, params.alpha)
    alpha = params.alpha if alpha is None else alpha
    gate = build_gate(params, cols)
    if gate is not None and bandcount is None:   # old .npz predates bandcount -> recompute from mask
        if mask is None and (cfg.mask_dir / f"{donor}.tif").exists():
            mask = tifffile.imread(cfg.mask_dir / f"{donor}.tif")
        if mask is None:
            raise SystemExit("gate needs bandcount but intermediates predate it and no mask exists — "
                             "regenerate with --save-intermediates --keep-mask")
        bandcount = bandcount_from_mask(mask, len(oids), params.edge_radius)
        log(f"[{donor}] recomputed bandcount from mask")
    del mask
    comp_mode = params.comp_mode
    alpha, comp_mode, n_skipped = apply_no_spillover_mask(params, cols, alpha, comp_mode)
    corrected, clamp_frac, n_iso = compensate(
        data, edge, sizes, M, comp_mode, alpha,
        bandcount=bandcount, gate=gate, backend=backend, norm_form=params.norm_form)
    log(f"[{donor}] recompensate(mode={params.comp_mode}, norm={params.norm_form}, "
        f"alpha={'vector' if np.ndim(alpha) else alpha}, gate={'on' if gate else 'off'}, "
        f"uncorrected_channels={n_skipped}): {n_iso} isolated")
    write_corrected(cfg, donor, oids, cols, corrected, sizes)
    return donor


# --------------------------------------------------------------------------- driver
def run_redsea(cfg: PipelineConfig, donors: list[str], params: RedseaParams,
               *, from_intermediates: bool = False, n_jobs: int = 1) -> list:
    """Run REDSEA over a list of donors, optionally across processes."""
    from .cohort import ensure_eligible_donors

    donors = list(ensure_eligible_donors(donors, context="REDSEA analysis"))
    worker = recompensate_from_intermediates if from_intermediates else process_donor
    fn = functools.partial(_safe_worker, worker=worker, cfg=cfg, params=params)
    return map_donors(fn, donors, n_jobs=n_jobs, ordered=True, on_error="log")


def _safe_worker(donor, worker, cfg, params):
    """Top-level wrapper so SystemExit (a skipped donor) doesn't kill the pool."""
    try:
        return worker(donor, cfg, params)
    except SystemExit as e:
        log(f"[{donor}] SKIP: {e}")
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--donor", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--jobs", type=int, default=None, help="per-donor process pool size")
    ap.add_argument("--gpu", action="store_true", help="use the CuPy backend (n_jobs=1 recommended)")
    ap.add_argument("--downsample", type=float, default=None, help="power of 2; 1=full res")
    ap.add_argument("--edge-radius", type=int, default=None,
                    help="boundary band dilation radius (px); 0 = the 1-px cell rim")
    ap.add_argument("--comp-mode", type=int, default=None, choices=[0, 1],
                    help="1=subtract+reinforce (default, Bai et al.); 0=subtract only (pre-v10)")
    ap.add_argument("--norm-form", default=None, choices=["donor", "recipient"],
                    help="contact normalisation: donor=contact[i,j]/B[j] (published, mass-conserving, "
                         "default); recipient=contact[i,j]/B[i] (pre-v10, rows sum to 1)")
    ap.add_argument("--no-exclude-no-spillover", action="store_true",
                    help="correct EVERY channel including nuclear/ECM (pre-v10 behaviour); by default "
                         "marker_taxonomy.NO_SPILLOVER channels pass through uncorrected")
    ap.add_argument("--alpha", type=float, default=None, help="spillover subtraction strength")
    ap.add_argument("--alpha-per-channel", default=None,
                    help="per-channel subtraction strengths, e.g. "
                         "'Pan_Cytokeratin=2,Ker8_18=2,Keratin_5=1.5,Vimentin=1.5'. Channels not listed "
                         "use --alpha; default None -> the scalar --alpha for every channel (bit-identical).")
    ap.add_argument("--gate-brightness", action="store_true",
                    help="enable the receptivity (brightness) gate: protect cells bright in their OWN "
                         "interior (real acinar keratin / vessel-fibroblast Vimentin) and fully clean dim "
                         "victims. Off by default (w=1 -> bit-identical output).")
    ap.add_argument("--gate-channels", default=None,
                    help="comma-separated channels to brightness-gate (default: the keratins + Vimentin).")
    ap.add_argument("--gate-kappa", type=float, default=1.0,
                    help="receptivity-gate aggressiveness; higher -> clean more cells (w saturates faster).")
    ap.add_argument("--gap-bridge", type=int, default=None,
                    help="bridge background gaps <= N px in the contact graph (0 disables)")
    ap.add_argument("--keep-mask", action="store_true", help="persist the int32 instance mask")
    ap.add_argument("--save-intermediates", action="store_true",
                    help="dump data/edge/sizes/contact for fast --from-intermediates re-runs")
    ap.add_argument("--from-intermediates", action="store_true",
                    help="skip raster+channels; reload saved intermediates and just (re)compensate")
    ap.add_argument("--reconcile-existing", action="store_true",
                    help="atomically filter retained REDSEA parquets to canonical data/cells object IDs; "
                         "does not recompute REDSEA")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    if a.gpu:
        cfg.use_gpu = True

    over = {}
    for k, src in [("downsample", a.downsample), ("edge_radius", a.edge_radius),
                   ("comp_mode", a.comp_mode), ("alpha", a.alpha), ("gap_bridge", a.gap_bridge),
                   ("norm_form", a.norm_form)]:
        if src is not None:
            over[k] = src
    if a.no_exclude_no_spillover:
        over["exclude_no_spillover"] = False
    params = RedseaParams.from_config(cfg, keep_mask=a.keep_mask,
                                      save_intermediates=a.save_intermediates,
                                      alpha_per_channel=a.alpha_per_channel,
                                      gate_brightness=a.gate_brightness,
                                      gate_channels=a.gate_channels,
                                      gate_kappa=a.gate_kappa, **over)

    if a.all:
        donors = cfg.discover_donors()
    elif a.donor:
        from .cohort import ensure_eligible_donors

        donors = list(
            ensure_eligible_donors([a.donor], context="REDSEA analysis")
        )
    else:
        raise SystemExit("specify --donor <id> or --all")

    # GPU + a local process pool would contend on one device during REDSEA
    # computation. Existing-file reconciliation is CPU/I/O-only and remains parallel.
    n_jobs = cfg.n_jobs
    if not a.reconcile_existing and params.use_gpu and n_jobs != 1:
        log("[warn] --gpu with --jobs>1 would contend on one device; forcing jobs=1 "
            "(use a SLURM array for GPU cohort parallelism)")
        n_jobs = 1

    if a.reconcile_existing:
        if a.from_intermediates:
            raise SystemExit("--reconcile-existing cannot be combined with --from-intermediates")
        reconcile_existing(cfg, donors, n_jobs=n_jobs)
    else:
        run_redsea(cfg, donors, params, from_intermediates=a.from_intermediates, n_jobs=n_jobs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
