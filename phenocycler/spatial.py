#!/usr/bin/env python3
"""
REDSEA before/after **spatial** diagnostic — draw the cell segmentation masks over the marker channel
and fill each cell by its raw (before) vs REDSEA (after) mean and their Δ, for the most- vs least-changed
exemplars. Shows the actual lateral spillover in tissue: a dim cell next to a bright neighbour loses the
false signal that bled across their shared boundary (MOST-changed), while an isolated / uniform-neighbourhood
cell is left essentially untouched (Δ≈0, LEAST-changed).

Reads ONLY the ROI's tiles from the OME-TIFF (whole-channel reads are fatally slow off a spinning disk) and
ONLY the ROI's strips from the striped int32 label mask. Falls back to the non-spatial raw-vs-corrected
scatter/hist when a donor has no OME image or instance mask.

"before" = the raw QuPath mean (``cfg.cells_dir``); "after" = the REDSEA mean (``cfg.cells_redsea_dir``),
joined on ``object_id`` (never coordinates); corrected ≤ raw (subtract-only) and the two align at r ≈ 0.998.

Thin entry point for the notebooks::

    from phenocycler.spatial import plot_redsea_before_after
    plot_redsea_before_after(cfg, "6374", display_marker="auto", row_mode="both")

Ported from ``restore_threshold_tracing.ipynb`` Step A (globals threaded through as arguments).
"""

from __future__ import annotations

import glob
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import matplotlib as mpl
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries


# log10(MFI+1), NaN/neg-safe (REDSEA can go slightly < 0) — used for readable intensity plots.
def L(a):
    return np.log10(np.clip(np.asarray(a, float), 0, None) + 1)


_RC = {"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 12, "figure.facecolor": "white"}


# --------------------------------------------------------------------------- #
# Image I/O — tiled OME reads + striped mask reads (importable, no globals)
# --------------------------------------------------------------------------- #
def ome_px_um(tf):
    """µm/px from the OME-TIFF XResolution tag (ResolutionUnit=3 -> cm)."""
    xr = tf.series[0].levels[0].pages[0].tags["XResolution"].value
    return 1e4 / (xr[0] / xr[1])


def read_mask_band(page, y0, y1):
    """Rows [y0,y1) of the striped int32 label mask -> (band, band_y_origin). Reads only the ROI strips."""
    rps = page.rowsperstrip; W = page.shape[-1]; fh = page.parent.filehandle
    s_lo, s_hi = y0 // rps, (y1 - 1) // rps
    band = np.zeros(((s_hi - s_lo + 1) * rps, W), np.int32)
    for si in range(s_lo, s_hi + 1):
        fh.seek(page.dataoffsets[si]); seg = fh.read(page.databytecounts[si])
        a = np.asarray(page.decode(seg, si)[0]).reshape(-1, W)
        band[(si - s_lo) * rps:(si - s_lo) * rps + a.shape[0], :a.shape[1]] = a
    return band, s_lo * rps


def read_ome_region(tf, level, ch, y0, y1, x0, x1):
    """Decode ONLY the tiles overlapping [y0:y1, x0:x1] of one channel at ``level`` (level-local px).
    Avoids whole-channel reads (fatally slow off the /data HDD). Validated exact vs asarray()+crop."""
    lv = tf.series[0].levels[level]; kf = lv.keyframe; pg = lv.pages[ch]
    th, tw = int(kf.tilelength), int(kf.tilewidth)
    Hl, Wl = int(kf.shape[-2]), int(kf.shape[-1])
    y0, x0 = max(0, y0), max(0, x0); y1, x1 = min(Hl, y1), min(Wl, x1)
    tpr = (Wl + tw - 1) // tw
    tr0, tr1 = y0 // th, (y1 - 1) // th
    tc0, tc1 = x0 // tw, (x1 - 1) // tw
    out = np.zeros(((tr1 - tr0 + 1) * th, (tc1 - tc0 + 1) * tw), kf.dtype)
    fh = pg.parent.filehandle
    for tr in range(tr0, tr1 + 1):
        for tc in range(tc0, tc1 + 1):
            idx = tr * tpr + tc
            cnt = int(pg.databytecounts[idx])
            if cnt == 0:
                continue
            fh.seek(int(pg.dataoffsets[idx])); seg = fh.read(cnt)
            a = np.asarray(kf.decode(seg, idx)[0]).reshape(th, tw)
            out[(tr - tr0) * th:(tr - tr0 + 1) * th, (tc - tc0) * tw:(tc - tc0 + 1) * tw] = a
    return out[y0 - tr0 * th:y1 - tr0 * th, x0 - tc0 * tw:x1 - tc0 * tw]


def pick_level(tf, forced):
    if forced is not None:
        return forced
    return min(1, len(tf.series[0].levels) - 1)       # tiled reads are cheap -> just a display resolution


def build_roi(tf, mpg, px_um, cx_um, cy_um, ch, level, half_um):
    """Assemble a square ROI (gray marker image + int32 label mask) around a centroid (µm)."""
    sub = 2 ** level
    H, W = mpg.shape[-2], mpg.shape[-1]
    hw = int(round(half_um / px_um))
    y0, y1 = max(0, int(cy_um / px_um) - hw), min(H, int(cy_um / px_um) + hw)
    x0, x1 = max(0, int(cx_um / px_um) - hw), min(W, int(cx_um / px_um) + hw)
    gray = read_ome_region(tf, level, ch, y0 // sub, y1 // sub, x0 // sub, x1 // sub)
    band, by0 = read_mask_band(mpg, y0, y1)
    m = band[y0 - by0:y1 - by0, x0:x1]
    if sub > 1:
        m = m[::sub, ::sub]
    hh = min(gray.shape[0], m.shape[0]); ww = min(gray.shape[1], m.shape[1])
    return gray[:hh, :ww], m[:hh, :ww], sub


def fill_by_label(roi_mask, col_values):
    """Fill each cell in the ROI label mask by its per-cell value (label i -> col_values[i-1]; 0 -> NaN)."""
    lut = np.concatenate([[np.nan], col_values.astype(np.float32)])
    return lut[roi_mask]


# --------------------------------------------------------------------------- #
# Per-donor change table + exemplar selection
# --------------------------------------------------------------------------- #
def load_change(cfg, donor):
    """RAW (cells/, 'before') and COR (cells_redsea/, 'after') aligned so row i == mask label i+1.
    Returns (markers, object_ids, RAW, COR, D=after-before, meta[centroids/region/islet/area])."""
    fcor = glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet"))[0]
    markers = [c for c in pq.ParquetFile(fcor).schema_arrow.names
               if c not in ("object_id", "donor_id", "cell_area_px")]
    t = time.time()
    cortab = pq.ParquetFile(fcor).read(columns=["object_id"] + markers).to_pandas()
    oids = cortab["object_id"].astype(str).to_numpy()
    COR = cortab[markers].to_numpy(np.float32); del cortab
    fcell = glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet"))[0]
    cell_cols = list(pq.ParquetFile(fcell).schema_arrow.names)
    meta_want = [c for c in ("X_centroid", "Y_centroid", "cell_region", "islet_num", "cell_area") if c in cell_cols]
    raw = pq.ParquetFile(fcell).read(columns=["object_id"] + markers + meta_want).to_pandas()
    raw["object_id"] = raw["object_id"].astype(str)
    raw = raw.set_index("object_id").reindex(oids)                     # align to COR / mask-label order
    RAW = raw[markers].to_numpy(np.float32)
    meta = raw[meta_want].reset_index(drop=True)
    print(f"  loaded {len(oids):,} cells x {len(markers)} markers in {time.time() - t:.1f}s")
    return markers, oids, RAW, COR, COR - RAW, meta                    # D = after - before (<= 0; subtract-only)


def select_exemplars(D, meta, n_most=3, n_least=2, *, RAW=None, markers=None, dapi_marker="DAPI", dapi_pct=90):
    sad = np.abs(D).sum(1)                                             # NaN if any marker is NaN
    valid = np.isfinite(sad)
    region = meta["cell_region"].to_numpy() if "cell_region" in meta else np.array(["core"] * len(sad))
    area = meta["cell_area"].to_numpy(float) if "cell_area" in meta else np.full(len(sad), np.nan)
    a_lo, a_hi = (np.nanpercentile(area, 25), np.nanpercentile(area, 95)) if np.isfinite(area).any() else (-np.inf, np.inf)
    # Exclude highly-DAPI+ objects: a nucleus-dominated segmentation (the outline hugs the nucleus, not the
    # whole cell) reads as a DAPI-bright cell and is a poor whole-cell exemplar. Drop cells whose RAW DAPI is
    # above the `dapi_pct` percentile (set dapi_pct >= 100 to disable).
    if RAW is not None and markers is not None and dapi_marker in markers and 0 < dapi_pct < 100:
        dapi = RAW[:, markers.index(dapi_marker)].astype(float)
        not_dapi = np.isfinite(dapi) & (dapi <= np.nanpercentile(dapi, dapi_pct))
    else:
        not_dapi = np.ones(len(sad), bool)
    dense = np.isin(region, ["core", "peri"]) & valid & not_dapi & (area >= a_lo) & (area <= a_hi)
    cand = np.where(dense)[0]
    most = cand[np.argsort(sad[cand])[::-1][:n_most]]                  # biggest absolute cleanup
    core_pos = np.where(np.isin(region, ["core"]) & valid & not_dapi & (sad > 0) & (area >= a_lo) & (area <= a_hi))[0]
    thr = np.percentile(sad[core_pos], 7)                             # low pct (NOT the min: isolated cells are Δ≡0)
    near = core_pos[np.abs(sad[core_pos] - thr) < (thr * 0.5 + 1e-6)]
    least = near[np.argsort(sad[near])[:n_least]] if len(near) >= n_least else core_pos[np.argsort(sad[core_pos])[:n_least]]
    return list(most), list(least), sad


# --------------------------------------------------------------------------- #
# Panel drawing
# --------------------------------------------------------------------------- #
def _imshow_gray(ax, gray, mask=None):
    # Anchor black at 0 (no signal -> black) and white at the 99th-percentile of the IN-MASK (tissue) pixels.
    # This keeps low-expression markers dark and bright markers bright — instead of an independent per-panel
    # [1,99] stretch that amplifies background/low-signal noise to full contrast — without the over-saturation
    # of scaling raw pixels to the (narrower) per-cell-mean range.
    fg = gray[mask > 0] if (mask is not None and np.any(mask > 0)) else np.asarray(gray).ravel()
    hi = float(np.percentile(fg, 99)) if fg.size else 1.0
    ax.imshow(gray, cmap="gray", vmin=0, vmax=(hi if hi > 0 else 1.0), interpolation="nearest")


def _outlines(ax, roi_mask, exemplar_label=None):
    b = find_boundaries(roi_mask, mode="inner", connectivity=2, background=0)
    ax.imshow(np.ma.masked_where(~b, b), cmap=mpl.colors.ListedColormap(["#26c6da"]),
              vmin=0, vmax=1, interpolation="nearest", alpha=0.9)
    if exemplar_label is not None:
        be = find_boundaries(roi_mask == exemplar_label, mode="inner", connectivity=2, background=0)
        ax.imshow(np.ma.masked_where(~be, be), cmap=mpl.colors.ListedColormap(["#ffd400"]),
                  vmin=0, vmax=1, interpolation="nearest")


def _scalebar(ax, roi_mask, px_um, sub, um=20):
    bar = um / (px_um * sub); x = roi_mask.shape[1] * 0.06; y = roi_mask.shape[0] * 0.92
    ax.plot([x, x + bar], [y, y], "-", color="w", lw=3, solid_capstyle="butt")
    ax.text(x + bar / 2, y - roi_mask.shape[0] * 0.03, f"{um}µm", color="w", ha="center", va="bottom", fontsize=9)


def _resolve_display(markers, meanabs, display_marker, n_marker_rows):
    order = np.argsort(meanabs)[::-1]; top = [markers[i] for i in order]; inert = markers[int(np.argmin(meanabs))]
    if isinstance(display_marker, str) and display_marker != "auto":
        return {"named": display_marker}
    if isinstance(display_marker, (list, tuple)):
        names = []
        for x in display_marker:
            names += (top[:n_marker_rows - 1] + [inert]) if x == "auto" else [x]
        return {"list": list(dict.fromkeys(names))}
    return {"auto": (top, inert)}


def _panels(fig, axes, r, gray, roi_mask, before, after, px_um, sub, lab=None, scalebar=True):
    delta = after - before
    vals = np.concatenate([before[np.isfinite(before)], after[np.isfinite(after)]])
    vmin, vmax = (np.percentile(vals, [1, 99]) if vals.size else (0, 1))
    if vmax <= vmin:
        vmax = vmin + 1
    M = (np.nanpercentile(np.abs(delta), 98) if np.isfinite(delta).any() else 1.0) or 1.0
    vir = plt.cm.viridis.copy(); vir.set_bad("0.12"); rdb = plt.cm.RdBu_r.copy(); rdb.set_bad("0.12")
    a0, a1, a2, a3 = axes[r]
    _imshow_gray(a0, gray, roi_mask); _outlines(a0, roi_mask, lab)
    if scalebar:
        _scalebar(a0, roi_mask, px_um, sub)
    im_v = a1.imshow(np.ma.masked_invalid(before), cmap=vir, vmin=vmin, vmax=vmax, interpolation="nearest")
    a2.imshow(np.ma.masked_invalid(after), cmap=vir, vmin=vmin, vmax=vmax, interpolation="nearest")
    im_d = a3.imshow(np.ma.masked_invalid(delta), cmap=rdb, vmin=-M, vmax=M, interpolation="nearest")
    for a in (a1, a2, a3):
        _outlines(a, roi_mask, lab)
    for a in axes[r]:
        a.set_xticks([]); a.set_yticks([])
    fig.colorbar(im_v, ax=[a1, a2], fraction=0.03, pad=0.01, label=("marker mean" if r == 0 else None))
    fig.colorbar(im_d, ax=[a3], fraction=0.05, pad=0.01, label=("Δ" if r == 0 else None))
    if r == 0:
        a0.set_title("marker + segmentation"); a1.set_title("BEFORE (raw)")
        a2.set_title("AFTER (REDSEA)"); a3.set_title("Δ (after − before)")


def _render_cell(tf, mpg, px_um, markers, RAW, COR, D, meta, most, least, sad, meanabs, *,
                 donor, display_marker, n_most, n_least, adaptive_zoom, zoom_k, half_um_floor,
                 half_um, cell_level, level, n_marker_rows, fig_dir, save_png):
    rows = [("MOST", i) for i in most] + [("LEAST", i) for i in least]
    lvl = pick_level(tf, cell_level if cell_level is not None else level)
    disp = _resolve_display(markers, meanabs, display_marker, n_marker_rows)

    def marker_for(i):
        if "named" in disp:
            return disp["named"]
        if "list" in disp:
            idxs = [markers.index(m) for m in disp["list"] if m in markers]
            return markers[idxs[int(np.argmax(np.abs(D[i, idxs])))]]
        return markers[int(np.nanargmax(np.abs(D[i])))]

    PAN = 2.9
    fig, axes = plt.subplots(len(rows), 4, figsize=(4 * PAN + 2.4, len(rows) * PAN + 0.8), squeeze=False)
    for r, (kind, i) in enumerate(rows):
        ch = markers.index(marker_for(i))
        cx, cy = float(meta.X_centroid.iloc[i]), float(meta.Y_centroid.iloc[i])
        half = (max(half_um_floor, zoom_k * np.sqrt(float(meta.cell_area.iloc[i]) / np.pi))
                if (adaptive_zoom and "cell_area" in meta and np.isfinite(meta.cell_area.iloc[i])) else half_um)
        gray, roi_mask, sub = build_roi(tf, mpg, px_um, cx, cy, ch, lvl, half)
        before = fill_by_label(roi_mask, RAW[:, ch]); after = fill_by_label(roi_mask, COR[:, ch])
        _panels(fig, axes, r, gray, roi_mask, before, after, px_um, sub, lab=i + 1)
        axes[r][0].set_ylabel(f"{kind}\n{markers[ch]}\nraw {RAW[i,ch]:.0f}→{COR[i,ch]:.0f}\nΣ|Δ|={sad[i]:.0f}",
                              rotation=0, ha="right", va="center", labelpad=12, fontsize=10)
    fig.suptitle(f"{donor} — REDSEA spillover correction: most- vs least-changed cells "
                 f"(tight zoom per cell; display marker = per-cell largest |Δ|)", fontsize=15)
    if save_png:
        p = fig_dir / f"{donor}_redsea_before_after_cells.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); print("  saved", p)
    plt.show()


def _render_marker(tf, mpg, px_um, markers, RAW, COR, meta, most, meanabs, *,
                   donor, display_marker, marker_half_um, level, n_marker_rows, fig_dir, save_png):
    disp = _resolve_display(markers, meanabs, display_marker, n_marker_rows)
    if "named" in disp:
        mks = [disp["named"]]
    elif "list" in disp:
        mks = disp["list"]
    else:
        top, inert = disp["auto"]; mks = top[:n_marker_rows - 1] + [inert]
    i0 = most[0]; cx, cy = float(meta.X_centroid.iloc[i0]), float(meta.Y_centroid.iloc[i0])
    lvl = pick_level(tf, 0 if level is None else level); PAN = 3.0   # full-res: 2x-downsampled levels wash out / brighten
    fig, axes = plt.subplots(len(mks), 4, figsize=(4 * PAN + 2.4, len(mks) * PAN + 0.8), squeeze=False)
    for r, mk in enumerate(mks):
        ch = markers.index(mk)
        gray, roi_mask, sub = build_roi(tf, mpg, px_um, cx, cy, ch, lvl, marker_half_um)
        before = fill_by_label(roi_mask, RAW[:, ch]); after = fill_by_label(roi_mask, COR[:, ch])
        _panels(fig, axes, r, gray, roi_mask, before, after, px_um, sub, lab=None, scalebar=(r == 0))
        axes[r][0].set_ylabel(f"{mk}\nmean|Δ|={meanabs[ch]:.1f}", rotation=0, ha="right", va="center", labelpad=12, fontsize=10)
    try:
        isl = int(meta.islet_num.iloc[i0])
    except (ValueError, TypeError):
        isl = "?"
    fig.suptitle(f"{donor} — REDSEA correction across markers (one field around islet {isl})", fontsize=15)
    if save_png:
        p = fig_dir / f"{donor}_redsea_before_after_markers.png"; fig.savefig(p, dpi=150, bbox_inches="tight"); print("  saved", p)
    plt.show()


def plot_redsea_scatter_hist(cfg, donor, marker=None):
    """REDSEA before/after **correlation** diagnostic — always available (no OME/mask needed): the
    raw-vs-REDSEA scatter (log-log, with the y=x line) + the distribution-shift histogram for one marker,
    plus the Pearson r and mean shift. ``marker=None`` picks ``Pan_Cytokeratin`` (else the first marker).
    Corrected ≤ raw (subtract-only); heavily-corrected markers decorrelate more as more spillover is removed."""
    donor = str(donor)
    cor_files = glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet"))
    if not cor_files:
        print(f"[scatter skipped] no REDSEA output for donor {donor} under {cfg.cells_redsea_dir}. "
              f"Run the REDSEA cell first.")
        return
    if marker is None:
        cols = pq.ParquetFile(cor_files[0]).schema_arrow.names
        marker = "Pan_Cytokeratin" if "Pan_Cytokeratin" in cols else next(
            m for m in cols if m not in ("object_id", "donor_id", "cell_area_px"))
    _raw = pq.ParquetFile(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet"))[0]).read(
        columns=["object_id", marker]).to_pandas().rename(columns={marker: "raw"})
    _cor = pq.ParquetFile(cor_files[0]).read(
        columns=["object_id", marker]).to_pandas().rename(columns={marker: "cor"})
    _raw.object_id = _raw.object_id.astype(str); _cor.object_id = _cor.object_id.astype(str)
    rc = _raw.merge(_cor, on="object_id").dropna()
    rr = np.corrcoef(rc.raw, rc.cor)[0, 1]
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    pi = np.random.default_rng(0).choice(len(rc), min(30000, len(rc)), replace=False)
    ax[0].scatter(L(rc.raw.values[pi]), L(rc.cor.values[pi]), s=2, c="0.5")
    lim = [0, float(L(rc.raw).max())]; ax[0].plot(lim, lim, "r--", lw=1, label="y = x (no change)")
    ax[0].set_xlabel(f"log10(raw {marker}+1)"); ax[0].set_ylabel(f"log10(REDSEA {marker}+1)"); ax[0].legend(fontsize=8)
    ax[0].set_title(f"raw vs REDSEA   r = {rr:.4f}")
    ax[1].hist(L(rc.raw), bins=120, color="0.7", density=True, label="raw")
    ax[1].hist(L(rc.cor), bins=120, histtype="step", color="C0", lw=1.5, density=True, label="REDSEA-corrected")
    ax[1].set_xlabel(f"log10({marker}+1)"); ax[1].set_ylabel("density"); ax[1].legend(fontsize=8); ax[1].set_title("distribution shift")
    fig.suptitle(f"{donor} — REDSEA spillover correction ({marker})"); fig.tight_layout(); plt.show()
    print(f"{marker}: raw {rc.raw.mean():.1f} -> REDSEA {rc.cor.mean():.1f}   r={rr:.4f}   "
          f"{100*(rc.cor <= rc.raw + 1e-6).mean():.1f}% of cells corrected <= raw")


def plot_redsea_before_after(cfg, donor, *, display_marker="auto", row_mode="both",
                             n_most=3, n_least=2, adaptive_zoom=True, zoom_k=3.5, half_um_floor=16.0,
                             half_um=45, cell_level=0, marker_half_um=130, n_marker_rows=4, level=None,
                             ome_dir=Path("/data/islet_ome_tiff"), fig_dir=None, save_png=True,
                             fallback_marker="Pan_Cytokeratin", dapi_pct=90):
    """REDSEA before/after spatial diagnostic for one donor. Renders the mask-overlay panels (cell and/or
    marker layout) if the OME image + instance mask exist, else the raw-vs-corrected scatter/hist fallback.

    ``ome_dir`` — folder of per-donor pyramidal OME-TIFFs (default the viewer store ``/data/islet_ome_tiff``).
    The instance mask is read from ``cfg.mask_dir/<donor>.tif``. ``display_marker`` "auto" | a name | a list.
    ``row_mode`` "both" | "cell" | "marker"."""
    donor = str(donor)
    if not glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet")):
        print(f"[diagnostic skipped] no REDSEA output yet for donor {donor} under {cfg.cells_redsea_dir}.\n"
              f"    This is a before/AFTER comparison — run the REDSEA cell below first, then re-run this cell.")
        return
    ome = Path(ome_dir) / f"{donor}.ome.tiff"
    mask = cfg.mask_dir / f"{donor}.tif"
    fig_dir = Path(fig_dir) if fig_dir else (cfg.data_dir.parent / "figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    if not (ome.exists() and mask.exists()):
        markers = pq.ParquetFile(glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet"))[0]).schema_arrow.names
        mk = fallback_marker if fallback_marker in markers else next(
            (m for m in markers if m not in ("object_id", "donor_id", "cell_area_px")), None)
        print(f"[spatial skipped] no OME image / instance mask for {donor}\n"
              f"    expected: {ome}\n              {mask}\n"
              f"    (regenerate the mask with the QuPath cell export + `python -m phenocycler.redsea "
              f"--donor {donor} --keep-mask`). Falling back to the raw-vs-corrected scatter/hist.")
        with plt.rc_context(_RC):
            plot_redsea_scatter_hist(cfg, donor, mk)
        return

    import tifffile
    tf = tifffile.TiffFile(str(ome)); mpg = tifffile.TiffFile(str(mask)).series[0].pages[0]
    px_um = ome_px_um(tf)
    print(f"{donor}: px_um={px_um:.4f}  OME levels={len(tf.series[0].levels)}")
    markers, oids, RAW, COR, D, meta = load_change(cfg, donor)
    meanabs = np.nansum(np.abs(D), 0) / np.isfinite(D).sum(0)
    most, least, sad = select_exemplars(D, meta, n_most, n_least, RAW=RAW, markers=markers, dapi_pct=dapi_pct)
    print("  MOST :", [(markers[int(np.nanargmax(np.abs(D[i])))], round(float(sad[i]))) for i in most])
    print("  LEAST:", [round(float(sad[i]), 1) for i in least])
    with plt.rc_context(_RC):
        if row_mode in ("cell", "both"):
            _render_cell(tf, mpg, px_um, markers, RAW, COR, D, meta, most, least, sad, meanabs,
                         donor=donor, display_marker=display_marker, n_most=n_most, n_least=n_least,
                         adaptive_zoom=adaptive_zoom, zoom_k=zoom_k, half_um_floor=half_um_floor,
                         half_um=half_um, cell_level=cell_level, level=level, n_marker_rows=n_marker_rows,
                         fig_dir=fig_dir, save_png=save_png)
        if row_mode in ("marker", "both"):
            _render_marker(tf, mpg, px_um, markers, RAW, COR, meta, most, meanabs,
                           donor=donor, display_marker=display_marker, marker_half_um=marker_half_um,
                           level=level, n_marker_rows=n_marker_rows, fig_dir=fig_dir, save_png=save_png)
