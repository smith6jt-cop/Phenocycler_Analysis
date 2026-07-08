#!/usr/bin/env python3
"""
Step 3 — RESTORE per-image intensity normalization of mutually-exclusive markers.

Faithful port of ``scripts/senior/restore_normalize.py`` (Islet-Explorer-Senior).
Applies RESTORE (Chang Lab / OHSU; vendored under ``external/RESTORE``) to the
REDSEA-corrected single-cell data.  **Pan_Cytokeratin** (the abundant ~91% acinar
majority) is the universal negative control; the lone ``CD3e <- CD163`` exception
keeps a clean exclusive immune reference (see ``DEFAULT_MARKER_PAIRS`` in
``config.py``).  Three 2-cluster models (KMeans, GMM, SSC) split background vs.
signal per image×pair; the threshold is ``mean + 3σ`` of the target intensity in
the target-negative cluster.  SSC is the default gating model.

Operates on RAW / REDSEA-corrected QuPath mean intensities (RESTORE's idx_select
uses an absolute >50 floor and clusters on raw 2D intensities — do NOT log the
input).  SSC needs ``spams`` (``pip install spams-bin``).

Pipeline (config-driven paths; the REDSEA-corrected branch):
    <cells_redsea>/donor_id=*/data_0.parquet   (REDSEA-corrected MFI)
      -> per-image full-image idx_select, subsample-after cap for clustering
      -> Normalization(save_figs=False).run()  -> <restore_dir>/threshs.pkl
      -> flatten                                -> <restore_thresholds_csv>
      -> robust cohort-median guard on the chosen model's thresholds
      -> apply per-image thresholds to ALL cells (parallel over donors):
           <m>_pos   = raw >= thr
           <m>_norm  = raw / thr
           <m>_log2r = log2((raw + 1) / thr)
         -> <restore_gated_dir>/donor_id=<id>/data_0.parquet
      -> QC plots + reference-exclusivity audit -> <restore_dir>/qc/

    python -m phenocycler.restore                       # full run
    python -m phenocycler.restore --limit-scenes 1 --skip-apply   # dry run
    python -m phenocycler.restore --model GMM
"""

from __future__ import annotations

import argparse
import functools
import glob
import pickle
import re
import sys
import types
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import PipelineConfig, load_config, DEFAULT_MARKER_PAIRS, EXTRA_MARKER_PAIRS
from .parallel import map_donors

MODEL_COLORS = {"KMeans": "magenta", "GMM": "blue", "SSC": "green"}


# --------------------------------------------------------------------------- #
# Vendored RESTORE import (headless)
# --------------------------------------------------------------------------- #
def _install_headless_stubs():
    """RESTORE/RESTORE.py imports holoviews, selenium, and tqdm.notebook at module
    load and builds hv.* objects during threshold computation even when
    save_figs=False. Inject lightweight stubs so the vendored source runs headless
    and *unmodified* — we never render or save figures."""
    import tqdm

    class _HVObj:  # supports .opts(), .cols(), and the `*` overlay operator
        def opts(self, *a, **k):
            return self

        def cols(self, *a, **k):
            return self

        def __mul__(self, other):
            return self

        def __rmul__(self, other):
            return self

    def _factory(*a, **k):
        return _HVObj()

    hv = types.ModuleType("holoviews")
    for name in ("Curve", "Scatter", "VLine", "Overlay", "Layout"):
        setattr(hv, name, _factory)
    hv.save = lambda *a, **k: None
    hv.extension = lambda *a, **k: None
    sys.modules.setdefault("holoviews", hv)

    sel = types.ModuleType("selenium")
    swd = types.ModuleType("selenium.webdriver")
    sel.webdriver = swd
    sys.modules.setdefault("selenium", sel)
    sys.modules.setdefault("selenium.webdriver", swd)

    tn = types.ModuleType("tqdm.notebook")
    tn.tqdm = tqdm.tqdm
    sys.modules.setdefault("tqdm.notebook", tn)


def patch_restore_sigma(Normalization, sigma_map, default=3.0):
    """Override the target marker's sigma_weight (threshold = neg_mean + sigma·neg_σ)
    per target marker, leaving the vendored source untouched."""
    import functools as _ft
    if not sigma_map:
        return
    for fname in ("get_GMM_thresh", "get_clustering_thresh"):
        orig = getattr(Normalization, fname)

        def make(orig):
            @_ft.wraps(orig)
            def wrapper(self, data, marker_pair, model, sigma_weight, color, batch):
                sw = sigma_map.get(marker_pair[0], default)
                return orig(self, data, marker_pair, model, sw, color, batch)
            return wrapper
        setattr(Normalization, fname, make(orig))
    print(f"[sigma] per-marker sigma_weight override: {sigma_map} (default {default})")


def patch_restore_cluster_cap(Normalization, cap, seed=0):
    """Subsample-AFTER-idx_select: feed FULL per-image data to idx_select, then cap the
    resulting cloud to `cap` cells right before clustering so SSC/GMM/KMeans stay tractable."""
    import functools as _ft
    if not cap:
        return
    rng = np.random.default_rng(seed)
    for fname in ("get_GMM_thresh", "get_clustering_thresh"):
        orig = getattr(Normalization, fname)

        def make(orig):
            @_ft.wraps(orig)
            def wrapper(self, data, marker_pair, model, sigma_weight, color, batch):
                if data is not None and len(data) > cap:
                    data = data[rng.choice(len(data), size=cap, replace=False)]
                return orig(self, data, marker_pair, model, sigma_weight, color, batch)
            return wrapper
        setattr(Normalization, fname, make(orig))
    print(f"[cluster-cap] idx_select cloud capped to {cap:,} cells before clustering")


def patch_restore_no_kde(Normalization):
    """Neutralize RESTORE's figure-only O(n^2) gaussian_kde (a no-op with save_figs=False)."""
    R = sys.modules[Normalization.__module__]
    R.gaussian_kde = lambda dataT: (lambda q: np.zeros(np.asarray(q).shape[-1]))
    print("[kde-stub] figure-only gaussian_kde neutralized")


def import_restore(vendor: Path):
    _install_headless_stubs()
    if str(vendor) not in sys.path:
        sys.path.insert(0, str(vendor))
    try:
        from RESTORE import Normalization, nested_dict  # noqa
    except ModuleNotFoundError as e:
        if "spams" in str(e):
            raise SystemExit("[fatal] SSC model needs `spams`. Install:  pip install spams-bin")
        raise SystemExit(
            f"[fatal] cannot import RESTORE from {vendor} ({e}). "
            "Clone it:  git clone https://github.com/smith6jt-cop/RESTORE.git external/RESTORE")
    return Normalization, nested_dict


# --------------------------------------------------------------------------- #
# Data
# --------------------------------------------------------------------------- #
def donor_files(cells_dir, limit=None, donors=None):
    files = sorted(glob.glob(str(Path(cells_dir) / "donor_id=*" / "*.parquet")))
    pairs = [(re.search(r"donor_id=([^/]+)", f).group(1), f) for f in files]
    if donors:
        keep = {str(d) for d in donors}
        pairs = [p for p in pairs if p[0] in keep]
    return pairs[:limit] if limit else pairs


def build_restore_input(cells_dir, markers, subsample, seed, limit_scenes=None, donors=None):
    """Per-image (per-donor) subsample of raw MFI for threshold estimation.
    Returns a DataFrame with columns: scene (donor id), batch (0), + one per marker."""
    rng = np.random.default_rng(seed)
    frames = []
    for donor, f in donor_files(cells_dir, limit_scenes, donors):
        df = pd.read_parquet(f, columns=markers)
        n_raw = len(df)
        # REDSEA leaves a few all-marker-NaN edge cells; RESTORE's np.quantile is not NaN-aware.
        allnan = df.isna().all(axis=1)
        partial = df.isna().any(axis=1) & ~allnan
        if partial.any():
            print(f"[input] WARN {donor}: {int(partial.sum())} PARTIAL-NaN cells "
                  "(some markers valid) — dropped so the quantile stays finite")
        df = df[~df.isna().any(axis=1)]
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        df = df.astype("float64").reset_index(drop=True)
        df["scene"] = donor
        df["batch"] = 0
        frames.append(df)
        print(f"[input] image {donor}: {n_raw:,} cells ({n_raw - n:,} NaN dropped) -> {len(df):,} sampled")
    data = pd.concat(frames, ignore_index=True)
    print(f"[input] total RESTORE input: {len(data):,} cells / {data['scene'].nunique()} images")
    return data


# --------------------------------------------------------------------------- #
# Reference exclusivity QC
# --------------------------------------------------------------------------- #
def _idx_select(T, R, ratio_x=0.75, ratio_y=0.5):
    """Replicate RESTORE's idx_select prefilter (>50 floor + quantile box)."""
    return (((T > 50) & (R > np.quantile(R, ratio_y))) |
            ((R > 50) & (T > np.quantile(T, ratio_x))))


def reference_exclusivity_qc(cfg, cells_dir, pairs, subsample, seed, out_dir, donors=None):
    """Data-driven audit of each [target, reference] pair (svd ratio, pearson, abundance)."""
    rng = np.random.default_rng(seed)
    markers = sorted({m for pair in pairs for m in pair})
    frames = []
    for donor, f in donor_files(cells_dir, None, donors):
        df = pd.read_parquet(f, columns=["object_id"] + markers)
        df = df.dropna(subset=markers)
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        df = df.reset_index(drop=True)
        df["object_id"] = df["object_id"].astype(str)
        try:  # dist_islet lives in canonical data/cells (cells_redsea lacks it)
            raw = pd.read_parquet(cfg.cells_dir / f"donor_id={donor}" / "data_0.parquet",
                                  columns=["object_id", "dist_islet"])
            raw["object_id"] = raw["object_id"].astype(str)
            df = df.merge(raw, on="object_id", how="left")
        except Exception:
            df["dist_islet"] = np.nan
        frames.append(df)
    pool = pd.concat(frames, ignore_index=True)

    rows = []
    for target, ref in pairs:
        T = pool[target].to_numpy(np.float64); R = pool[ref].to_numpy(np.float64)
        sel = _idx_select(T, R)
        lt = np.log1p(T[sel]); lr = np.log1p(R[sel])
        M = np.vstack([(lt - lt.mean()) / (lt.std() or 1.0),
                       (lr - lr.mean()) / (lr.std() or 1.0)]).T
        s = np.linalg.svd(M, compute_uv=False)
        r = float(s[1] / s[0]) if s[0] > 0 else np.nan
        pear = float(np.corrcoef(np.log1p(T), np.log1p(R))[0, 1])
        refpos = R >= np.quantile(R, 0.75)
        bg = T[refpos]
        di = pool["dist_islet"].to_numpy(np.float64)[refpos]
        rows.append(dict(
            target=target, reference=ref,
            svd_ratio=round(r, 3), pearson_log=round(pear, 3),
            ref_frac_gt50=round(float((R > 50).mean()), 3),
            ref_frac_gt200=round(float((R > 200).mean()), 3),
            bg_median=round(float(np.median(bg)), 1),
            bg_p95=round(float(np.percentile(bg, 95)), 1),
            thr_mean3sd=round(float(bg.mean() + 3 * bg.std()), 1),
            refpos_dist_islet_med=round(float(np.nanmedian(di)), 1) if np.isfinite(di).any() else np.nan,
            flag_not_exclusive=bool(np.isfinite(r) and r < 0.4),
            flag_pos_corr=bool(pear > 0.3),
        ))
    qc = pd.DataFrame(rows)
    out_csv = Path(out_dir) / "reference_selection_qc.csv"
    qc.to_csv(out_csv, index=False)
    print(f"[ref-qc] wrote {out_csv}")
    print(qc.to_string(index=False))
    return qc


# --------------------------------------------------------------------------- #
# Thresholds
# --------------------------------------------------------------------------- #
def flatten_threshs(threshs):
    """threshs[batch][scene][target_marker][model] = value -> tidy DataFrame."""
    batches = [b for b in threshs if b != "global"] or list(threshs)
    rows = []
    for batch in batches:
        for scene in threshs[batch]:
            for marker in threshs[batch][scene]:
                for model, val in threshs[batch][scene][marker].items():
                    rows.append({"image": scene, "marker": marker, "model": model,
                                 "threshold": float(val)})
    return pd.DataFrame(rows)


def build_chosen_lut(thr_df, model, robust=True, factor=3.0):
    """Lookup (image, marker) -> threshold for the chosen model, with a robust
    cohort-median guard: any per-image threshold that is an outlier (> factor× or
    < 1/factor× the cohort median for that marker, cohort median over non-outliers)
    falls back to the cohort median.  Returns (lut, imputed)."""
    sub = thr_df[thr_df.model == model]
    lut = {(r.image, r.marker): float(r.threshold) for r in sub.itertuples()}
    imputed = {}
    if not robust:
        return lut, imputed
    for marker in sorted(sub.marker.unique()):
        vals = sub[sub.marker == marker].set_index("image").threshold.astype(float)
        med0 = vals.median()
        ok = vals[(vals <= factor * med0) & (vals >= med0 / factor)]
        med = ok.median() if len(ok) else med0
        for img, v in vals.items():
            if not np.isfinite(v) or v <= 0 or v > factor * med or v < med / factor:
                lut[(img, marker)] = float(med)
                imputed[(img, marker)] = (float(v), float(med))
    return lut, imputed


# --------------------------------------------------------------------------- #
# Apply to all cells (parallel over donors)
# --------------------------------------------------------------------------- #
def _apply_one_donor(donor, cells_dir, markers, lut, imputed, gated_dir):
    """Gate + normalize one donor's cells; write its partition; return frac rows."""
    files = sorted(glob.glob(str(Path(cells_dir) / f"donor_id={donor}" / "*.parquet")))
    if not files:
        return []
    df = pd.read_parquet(files[0], columns=["object_id"] + markers)
    out = pd.DataFrame({"object_id": df["object_id"].astype("string"), "donor_id": donor})
    frac_rows = []
    for m in markers:
        thr = lut.get((donor, m))
        if thr is None or not np.isfinite(thr) or thr <= 0:
            print(f"[apply] WARN image {donor} {m}: missing/invalid threshold ({thr}); skipped")
            continue
        raw = df[m].to_numpy(dtype="float64")
        pos = raw >= thr
        out[f"{m}_pos"] = pos
        out[f"{m}_norm"] = (raw / thr).astype("float32")
        out[f"{m}_log2r"] = np.log2((raw + 1.0) / thr).astype("float32")
        frac_rows.append({"image": donor, "marker": m, "threshold": float(thr),
                          "imputed": (donor, m) in imputed,
                          "n": int(len(raw)), "n_pos": int(pos.sum()),
                          "frac_pos": float(pos.mean())})
    dst = Path(gated_dir) / f"donor_id={donor}"
    dst.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dst / "data_0.parquet", index=False)
    print(f"[apply] image {donor}: wrote {len(out):,} rows -> {dst}/data_0.parquet", flush=True)
    return frac_rows


def apply_thresholds(cells_dir, markers, lut, imputed, gated_dir,
                     limit_scenes=None, donors=None, n_jobs=1):
    """Stream every donor parquet, gate + normalize, write partitions (parallel over donors)."""
    donor_ids = [d for d, _ in donor_files(cells_dir, limit_scenes, donors)]
    fn = functools.partial(_apply_one_donor, cells_dir=str(cells_dir), markers=markers,
                           lut=lut, imputed=imputed, gated_dir=str(gated_dir))
    results = map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True, on_error="log")
    rows = [r for res in results if res for r in res]
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# QC plots
# --------------------------------------------------------------------------- #
def qc_histograms(data, thr_df, markers, qc_dir):
    for scene in sorted(data["scene"].unique()):
        sd = data[data["scene"] == scene]
        ncol = 2
        nrow = int(np.ceil(len(markers) / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(6 * ncol, 4 * nrow))
        axes = np.atleast_1d(axes).ravel()
        for i, m in enumerate(markers):
            ax = axes[i]
            x = sd[m].to_numpy()
            hi = np.quantile(x, 0.995) if len(x) else 1.0
            ax.hist(np.clip(x, 0, hi), bins=80, color="0.75")
            ax.set_yscale("log")
            for model, col in MODEL_COLORS.items():
                r = thr_df[(thr_df.image == scene) & (thr_df.marker == m) & (thr_df.model == model)]
                if len(r):
                    t = float(r.threshold.iloc[0])
                    ax.axvline(t, color=col, lw=1.5, label=f"{model}={t:.0f}")
            ax.set_title(m); ax.set_xlabel("raw MFI"); ax.legend(fontsize=8)
        for j in range(len(markers), len(axes)):
            axes[j].axis("off")
        fig.suptitle(f"RESTORE thresholds — image {scene}")
        fig.tight_layout()
        fig.savefig(Path(qc_dir) / f"{scene}_hist.png", dpi=120)
        plt.close(fig)


def _annotated_heatmap(piv, title, cbar_label, path, fmt="{:.0f}", cmap="viridis"):
    fig, ax = plt.subplots(figsize=(1.0 * piv.shape[1] + 3, 0.6 * piv.shape[0] + 2))
    im = ax.imshow(piv.values.astype(float), aspect="auto", cmap=cmap)
    ax.set_xticks(range(piv.shape[1])); ax.set_xticklabels(piv.columns, rotation=45, ha="right")
    ax.set_yticks(range(piv.shape[0])); ax.set_yticklabels(piv.index)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if np.isfinite(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center", color="w", fontsize=7)
    fig.colorbar(im, ax=ax, label=cbar_label)
    ax.set_title(title); fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def qc_threshold_heatmap(thr_df, markers, model, qc_dir):
    piv = (thr_df[thr_df.model == model]
           .pivot(index="marker", columns="image", values="threshold").reindex(markers))
    _annotated_heatmap(piv, f"RESTORE per-image thresholds ({model})",
                       f"{model} threshold (raw MFI)", Path(qc_dir) / "threshold_heatmap.png")


def qc_positive_fractions(frac_df, markers, qc_dir):
    piv = (frac_df.pivot(index="marker", columns="image", values="frac_pos").reindex(markers) * 100.0)
    _annotated_heatmap(piv, "RESTORE positive fraction per image (chosen model)",
                       "% cells positive", Path(qc_dir) / "positive_fractions.png",
                       fmt="{:.1f}", cmap="magma")


# --------------------------------------------------------------------------- #
# Driver + Main
# --------------------------------------------------------------------------- #
def run_restore(cfg: PipelineConfig, *, cells_dir: Optional[Path] = None, model=None,
                subsample=None, marker_pairs=None, marker_sigma=None, limit_scenes=None,
                donors=None, skip_apply=False, robust=None, robust_factor=None,
                reuse_threshs=False, ref_qc=True, seed=None, n_jobs=None,
                out_dir: Optional[Path] = None, gated_dir: Optional[Path] = None,
                thresh_csv: Optional[Path] = None):
    """Run RESTORE fit + robust LUT + per-donor apply.  Reads the REDSEA-corrected
    cells by default; writes gated parquet + thresholds csv.

    ``out_dir`` / ``gated_dir`` / ``thresh_csv`` default (``None``) to the canonical
    ``cfg.restore_dir`` / ``cfg.restore_gated_dir`` / ``cfg.restore_thresholds_csv``, so the
    normal 10-marker path is byte-unchanged; :func:`run_restore_extra` overrides them to run
    the second (extra-marker) pass into the ``*_extra`` dirs without touching the validated
    10-marker gates."""
    cells_dir = Path(cells_dir) if cells_dir is not None else cfg.cells_redsea_dir
    model = model or cfg.restore_model
    subsample = cfg.restore_subsample if subsample is None else subsample
    robust = cfg.restore_robust if robust is None else robust
    robust_factor = cfg.restore_robust_factor if robust_factor is None else robust_factor
    seed = cfg.restore_seed if seed is None else seed
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs

    pairs = marker_pairs or [list(p) for p in DEFAULT_MARKER_PAIRS]
    markers = sorted({m for pair in pairs for m in pair})
    print(f"[config] pairs (target<-reference): {pairs}")
    print(f"[config] markers={markers}  chosen_model={model}  subsample={subsample}  cells_dir={cells_dir}")

    out_dir = cfg.restore_dir if out_dir is None else Path(out_dir)
    gated_dir = cfg.restore_gated_dir if gated_dir is None else Path(gated_dir)
    thresh_csv = cfg.restore_thresholds_csv if thresh_csv is None else Path(thresh_csv)
    qc_dir = out_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    Normalization, _ = import_restore(cfg.restore_vendor)
    if marker_sigma:
        patch_restore_sigma(Normalization, marker_sigma)
    patch_restore_cluster_cap(Normalization, subsample, seed)
    patch_restore_no_kde(Normalization)

    print("[step1] building RESTORE input (full images; idx_select runs on all cells) ...")
    data = build_restore_input(cells_dir, markers, None, seed, limit_scenes, donors)

    pkl = out_dir / "threshs.pkl"
    if reuse_threshs and pkl.exists():
        print(f"[step2] --reuse-threshs: loading existing {pkl}")
    else:
        print("[step2] running RESTORE (KMeans + GMM + SSC, save_figs=False) ...")
        norm = Normalization(data, pairs, save_dir=str(out_dir), save_figs=False)
        norm.run()
        print(f"[step2] wrote {pkl}")

    with open(pkl, "rb") as fh:
        threshs = pickle.load(fh)
    thr_df = flatten_threshs(threshs)
    thr_df["chosen"] = thr_df["model"] == model
    thr_df.sort_values(["marker", "image", "model"]).to_csv(thresh_csv, index=False)
    print(f"[step3] wrote {thresh_csv}  ({len(thr_df)} rows)")

    print("[step4] QC plots ...")
    qc_histograms(data, thr_df, markers, qc_dir)
    qc_threshold_heatmap(thr_df, markers, model, qc_dir)
    if ref_qc:
        print("[step4] reference-exclusivity audit ...")
        reference_exclusivity_qc(cfg, cells_dir, pairs, subsample, seed, out_dir, donors)
    print(f"[step4] wrote QC -> {qc_dir}")

    if skip_apply:
        print("[done] --skip-apply: thresholds + QC only")
        return thr_df

    chosen_lut, imputed = build_chosen_lut(thr_df, model, robust, robust_factor)
    if imputed:
        print(f"[robust] imputed {len(imputed)} degenerate {model} thresholds "
              f"(>{robust_factor}x cohort median) with the cohort median:")
        for (img, mk), (orig, med) in sorted(imputed.items()):
            print(f"  image {img} {mk}: {orig:.0f} -> {med:.0f}")
    elif robust:
        print(f"[robust] no outlier {model} thresholds to impute")

    print(f"[step5] applying thresholds to ALL cells (n_jobs={n_jobs}) ...")
    frac_df = apply_thresholds(cells_dir, markers, chosen_lut, imputed, gated_dir,
                               limit_scenes, donors, n_jobs=n_jobs)
    frac_df.sort_values(["marker", "image"]).to_csv(out_dir / "positive_fractions.csv", index=False)
    qc_positive_fractions(frac_df, markers, qc_dir)
    print(f"[step5] wrote {gated_dir} (partitioned) + {out_dir}/positive_fractions.csv")
    print("[done]")
    return thr_df


def run_restore_extra(cfg: PipelineConfig, *, donors=None):
    """Second RESTORE pass for the extra markers (CD99/B3TUBB/MPO <- Pan_Cytokeratin).

    Runs the SAME RESTORE machinery as :func:`run_restore` on ``EXTRA_MARKER_PAIRS`` with
    ``robust=False`` and ``ref_qc=False``, reading the REDSEA-corrected cells
    (``cfg.cells_redsea_dir``) and writing to the ``*_extra`` dirs
    (``cfg.restore_redsea_extra_dir`` / ``cfg.restore_gated_extra_dir`` /
    ``cfg.restore_thresholds_extra_csv``).  Keeping this a separate pass leaves the validated
    10-marker gates in ``restore_gated_redsea`` byte-identical.  Equivalent to the Senior
    command::

        restore_normalize.py \\
          --marker-pairs 'CD99:Pan_Cytokeratin,B3TUBB:Pan_Cytokeratin,MPO:Pan_Cytokeratin' \\
          --cells-dir data/cells_redsea --out-dir data/restore_redsea_extra \\
          --gated-dir data/restore_gated_redsea_extra \\
          --thresh-csv data/restore_thresholds_extra.csv --no-robust --no-ref-qc

    The gated output ({m}_pos / {m}_norm / {m}_log2r for CD99/B3TUBB/MPO) is merged into the
    broad-lineage call by ``object_id`` at assignment time."""
    pairs = [list(p) for p in EXTRA_MARKER_PAIRS]
    return run_restore(
        cfg,
        cells_dir=cfg.cells_redsea_dir,
        marker_pairs=pairs,
        donors=donors,
        robust=False,
        ref_qc=False,
        out_dir=cfg.restore_redsea_extra_dir,
        gated_dir=cfg.restore_gated_extra_dir,
        thresh_csv=cfg.restore_thresholds_extra_csv,
    )


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--cells-dir", type=Path, default=None,
                    help="input cells (default: the REDSEA-corrected cells_redsea)")
    ap.add_argument("--jobs", type=int, default=None, help="per-donor apply pool size")
    ap.add_argument("--subsample", type=int, default=None)
    ap.add_argument("--model", default=None, choices=["GMM", "KMeans", "SSC"])
    ap.add_argument("--marker-pairs", default=None, help="override pairs, e.g. 'CD3e:CD163,CD99:CD163'")
    ap.add_argument("--marker-sigma", default=None, help="per-target sigma, e.g. 'Pan_Cytokeratin:1.0'")
    ap.add_argument("--limit-scenes", type=int, default=None)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--skip-apply", action="store_true")
    ap.add_argument("--robust", action=argparse.BooleanOptionalAction, default=None)
    ap.add_argument("--robust-factor", type=float, default=None)
    ap.add_argument("--reuse-threshs", action="store_true")
    ap.add_argument("--ref-qc", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--extra", action="store_true",
                    help="run the SECOND RESTORE pass for the extra markers "
                         "(CD99/B3TUBB/MPO <- Pan_Cytokeratin) into the *_extra dirs; forces "
                         "robust off + ref-qc off so the 10-marker gates stay byte-identical "
                         "(other pair/model/sigma flags are ignored)")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    if a.extra:
        run_restore_extra(cfg, donors=a.donors)
        return 0
    pairs = ([p.split(":") for p in a.marker_pairs.split(",")] if a.marker_pairs else None)
    sigma = ({k: float(v) for k, v in (p.split(":") for p in a.marker_sigma.split(","))}
             if a.marker_sigma else None)
    run_restore(cfg, cells_dir=a.cells_dir, model=a.model, subsample=a.subsample,
                marker_pairs=pairs, marker_sigma=sigma, limit_scenes=a.limit_scenes,
                donors=a.donors, skip_apply=a.skip_apply, robust=a.robust,
                robust_factor=a.robust_factor, reuse_threshs=a.reuse_threshs,
                ref_qc=a.ref_qc, seed=a.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
