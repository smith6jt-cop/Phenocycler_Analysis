#!/usr/bin/env python3
"""
Step 3 — RESTORE per-image intensity normalization of mutually-exclusive markers.

RETIRED 2026-07-23 — the paired driver described below (``run_restore`` /
``run_restore_functional``) is DELETED along with the curated pair webs it consumed
(``MARKER_PAIRS`` / ``FUNCTIONAL_PAIRS``, see the tombstone in ``config.py``); ten of
those pairs referenced CD20, which is too sparse in pancreas to define a negative
control.  This module now serves only ``run_restore_proliferation`` plus the reusable
machinery below (``idx_select``, ``neg_stat``, ``fit_clusters``, the threshold LUT and
apply helpers, the QC plots), which the manuscript-faithful method reuses.  For Gate 2
use ``python -m phenocycler.restore_validation evaluate-locked``; see
``docs/restore_faithful_rebuild_plan.md``.  Everything from here to the Pipeline block
describes the RETIRED paired path and is kept for context only.

Applies RESTORE (Chang Lab / OHSU; vendored under ``external/RESTORE``) to the
REDSEA-corrected single-cell data.  Marker thresholds come from a curated web of
DIRECTIONAL mutually-exclusive pairs (``MARKER_PAIRS`` in ``config.py``): each pair
``[target, counterpart]`` thresholds the target using the counterpart's positive
cells as its negative population.  Three 2-cluster models (KMeans, GMM, SSC) split
background vs. signal per image×pair; the threshold is ``mean + 3σ`` of the target
intensity in the target-negative cluster.  SSC is the default gating model.

Operates on RAW / REDSEA-corrected QuPath mean intensities (RESTORE's idx_select
and clustering run on raw 2D intensities — do NOT log the input).  The co-positive
idx_select prefilter is replaced by the corrected mutually-exclusive single-positive
selection (see ``patch_restore_idx_floor``).  SSC needs ``spams`` (``pip install spams-bin``).

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
if not matplotlib.get_backend().startswith("module://"):  # keep a notebook's inline/widget backend; force Agg only headless
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import PipelineConfig, load_config, ENDOCRINE_MARKERS
from .marker_taxonomy import PROLIFERATION
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


# --------------------------------------------------------------------------- #
# Shared selection + threshold statistics — the SINGLE source of the idx-floor
# selection math, used both by the production ``patch_restore_idx_floor`` below
# and by the interactive threshold-tracing diagnostics (``phenocycler.diagnostics``).
# --------------------------------------------------------------------------- #
def idx_select(T, R, neg_q=0.5, ratio_x=0.75, ratio_y=0.5, nonzero_q=0.0):
    """Mutually-exclusive single-positive selection for a directional ``[target, reference]`` pair.

    Returns ``(mask, c0, c1)`` where ``mask`` keeps the two axis-hugging arms —
        REF+ / TGT<=c0   the target-NEGATIVE population that thresholds the target
        TGT+ / REF<=c1   signal
    and ``c0``/``c1`` are the per-axis near-zero ceilings = the ``neg_q`` quantile of T / R.
    ``neg_q <= 0`` reverts to RESTORE's vendored CO-POSITIVE ``>50`` prefilter (comparison escape hatch).
    ``nonzero_q > 0`` floors ``c0``/``c1`` at that quantile of the marker's NONZERO cells — a guard for
    REDSEA's zero-clip (if >=50% of a sparse marker is clamped to 0, ``quantile(T, neg_q)`` collapses to
    0 and the negative arm keeps only exact-zero cells). ``nonzero_q == 0`` is the exact current math.
    This is the exact selection the RESTORE idx-floor patch applies — kept in one place."""
    T = np.asarray(T, float); R = np.asarray(R, float)
    if neg_q <= 0:                                             # vendored co-positive >50 escape hatch
        f0 = f1 = 50.0
        mask = ((T > f0) & (R > np.quantile(R, ratio_y))) | ((R > f1) & (T > np.quantile(T, ratio_x)))
        return mask, f0, f1
    c0 = float(np.quantile(T, neg_q)); c1 = float(np.quantile(R, neg_q))   # OFF-marker near-zero ceilings
    if nonzero_q > 0:                       # guard REDSEA zero-clip collapse (>=50% zeros -> c0/c1 == 0)
        Tnz = T[T > 0]; Rnz = R[R > 0]
        if Tnz.size: c0 = max(c0, float(np.quantile(Tnz, nonzero_q)))
        if Rnz.size: c1 = max(c1, float(np.quantile(Rnz, nonzero_q)))
    mask = (((R > np.quantile(R, ratio_y)) & (T <= c0)) |     # REF+ / TGT-  = negative pop
            ((T > np.quantile(T, ratio_x)) & (R <= c1)))       # TGT+ / REF-  = signal
    return mask, c0, c1


def neg_stat(v, stat="mean3sd", k=3.0, pctile=99.865):
    """Threshold from the target-negative population's target-axis values (raw MFI).

    ``mean3sd`` is the production statistic (``mean + k·std``, ddof=0, matching vendored RESTORE);
    ``mad``/``logmad`` are sigma-consistent robust variants (median + k·1.4826·MAD, the latter in
    log space); ``pctile`` is coverage-matched to k·sigma (default 99.865 = 100·Φ(3))."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if v.size < 2:
        return np.nan
    if stat == "mean3sd":
        return float(v.mean() + k * v.std())
    if stat == "mad":
        m = np.median(v); return float(m + k * 1.4826 * np.median(np.abs(v - m)))
    if stat == "pctile":
        return float(np.percentile(v, pctile))
    if stat == "logmad":
        lv = np.log1p(np.clip(v, 0, None)); m = np.median(lv)   # clip: REDSEA can go slightly < 0
        return float(np.expm1(m + k * 1.4826 * np.median(np.abs(lv - m))))
    raise ValueError(f"unknown neg_stat: {stat!r}")


def fit_clusters(cloud, model="SSC", subsample=15000, seed=0, ssc_vendor=None):
    """Fit the RESTORE 2-cluster model on an idx_select cloud (columns = target, reference) and return
    ``(cloud, labels, neg)`` where ``neg`` is the cluster index with the largest std on the REFERENCE
    axis (reference-positive => target-negative — exactly the vendored negative-cluster selection).
    ``model`` in {"KMeans","GMM","SSC"}; SSC (default) needs the vendored ``ssc`` package (``spams``)."""
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture
    cloud = np.asarray(cloud, float)
    rng = np.random.default_rng(seed)
    if len(cloud) > subsample:
        cloud = cloud[rng.choice(len(cloud), subsample, replace=False)]
    if model == "KMeans":
        labels = KMeans(n_clusters=2, random_state=seed, n_init=10).fit_predict(cloud)
    elif model == "GMM":
        labels = GaussianMixture(n_components=2, n_init=5, random_state=seed).fit_predict(cloud)
    elif model == "SSC":
        vendor = Path(ssc_vendor) if ssc_vendor else (Path(__file__).resolve().parents[1]
                                                      / "external" / "RESTORE" / "python_code")
        if str(vendor) not in sys.path:
            sys.path.insert(0, str(vendor))
        try:
            from ssc.cluster.selfrepresentation import SparseSubspaceClusteringOMP  # NB stochastic + slow
        except ModuleNotFoundError as e:
            raise SystemExit(f"[fatal] SSC needs the vendored ssc pkg + spams ({e}); pip install spams-bin")
        m = SparseSubspaceClusteringOMP(n_clusters=2); m.fit(cloud); labels = m.labels_
    else:
        raise ValueError(f"unknown model: {model!r}")
    clusters = [cloud[labels == 0], cloud[labels == 1]]
    neg = int(np.argmax([c[:, 1].std() if len(c) else -1 for c in clusters]))
    return cloud, labels, neg


def r_otsu(x):
    """Pure Otsu valley on a 1-D array (``skimage.filters.threshold_otsu``); NaN if degenerate (<2 finite
    values or a flat array). Compute it on ``log10(MFI+1)`` for the honest background<->signal boundary.
    Lives here (not diagnostics.py) so the pipeline and ``phenocycler.diagnostics`` share one Otsu; the
    latter re-exports it. Lazy skimage import so ``restore`` loads without skimage installed."""
    from skimage.filters import threshold_otsu
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    return np.nan if x.size < 2 or x.min() == x.max() else float(threshold_otsu(x))


def bimodal_thresh(v, *, min_sep=2.0, min_frac_hi=0.02, min_n=200):
    """Standalone bimodal positivity cutoff for a single marker (NO mutually-exclusive counterpart).

    Fit a 2-component Gaussian on ``log10(v+1)`` of the NONZERO cells (drops the REDSEA zero-clip spike) and
    cut at the between-component crossover -> raw-MFI threshold. Returns ``+inf`` (ABSENT / unimodal) when the
    two components are not separated (``sep = |mu_hi - mu_lo| / pooled_sd < min_sep``), the positive component
    is negligible (``w_hi < min_frac_hi``), or too few cells. Used for LINEAGE-AGNOSTIC markers that have no
    RESTORE pair — proliferation (Ki67/PCNA) — and is the shared core of :func:`partner_low_thresh`."""
    from sklearn.mixture import GaussianMixture
    v = np.asarray(v, float); v = v[np.isfinite(v)]; v = v[v > 0]
    if v.size < min_n:
        return float("inf")
    vlog = np.log10(v + 1.0).reshape(-1, 1)
    g = GaussianMixture(2, n_init=3, random_state=0).fit(vlog)
    mu = g.means_.ravel(); sd = np.sqrt(g.covariances_.ravel()); w = g.weights_
    lo, hi = np.argsort(mu)
    sep = float((mu[hi] - mu[lo]) / np.sqrt(0.5 * (sd[lo] ** 2 + sd[hi] ** 2)))
    if sep < min_sep or w[hi] < min_frac_hi:                    # unimodal -> ABSENT (no manufactured cut)
        return float("inf")
    xs = np.linspace(mu[lo], mu[hi], 512)
    d_lo = w[lo] * np.exp(-0.5 * ((xs - mu[lo]) / sd[lo]) ** 2) / sd[lo]
    d_hi = w[hi] * np.exp(-0.5 * ((xs - mu[hi]) / sd[hi]) ** 2) / sd[hi]
    return float(10 ** xs[np.argmin(np.abs(d_lo - d_hi))] - 1.0)  # between-component crossover -> raw MFI


def logmean_ksigma(v, k=3.0, *, min_n=200):
    """Robust per-image upper-tail cutoff for a GRADED marker with no clean bimodality (proliferation
    Ki67/PCNA): ``mean + k*std`` of ``log10(nonzero+1)``, converted back to raw MFI. Per-donor adaptive
    (no fixed-fraction assumption); isolates the bright proliferating tail. ``+inf`` if too few cells."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]; v = v[v > 0]
    if v.size < min_n:
        return float("inf")
    lv = np.log10(v + 1.0)
    return float(10 ** (lv.mean() + k * lv.std()) - 1.0)


def partner_low_thresh(T, R, partner_low_q=0.5, *, islet=None, min_sep=3.0, min_frac_hi=0.03, min_n=200):
    """Threshold the TARGET among cells LOW for the mutually-exclusive REFERENCE (ideally within islets).

    For a directional ``[target, reference]`` pair (T=target, R=reference), among reference-LOW cells the
    target separates into background vs signal — the double-positive / co-expressing confound (RESTORE's
    problematic corner) is removed. This is the OPPOSITE conditioning to ``idx_select``'s negative arm
    (reference-HIGH). Because the non-islet tissue background dominates whole-slide and is broad, pass
    ``islet`` (a boolean mask of core/peri cells): within islets the reference-low target is cleanly bimodal.
    Delegates to :func:`bimodal_thresh` on the reference-low(-islet) subpopulation; returns ``+inf`` (ABSENT)
    when it is unimodal (e.g. a beta-loss donor's INS among GCG-low cells)."""
    T = np.asarray(T, float); R = np.asarray(R, float)
    sub = R <= np.quantile(R, partner_low_q)
    if islet is not None:
        sub = sub & np.asarray(islet, bool)
    return bimodal_thresh(T[sub], min_sep=min_sep, min_frac_hi=min_frac_hi, min_n=min_n)


def patch_restore_idx_floor(Normalization, neg_q=0.5, neg_stat_name="mean3sd",
                            nonzero_q=0.0, presence_min_sep=0.0, subsample=15000, seed=0,
                            partner_low_q=0.0):
    """Replace RESTORE's CO-POSITIVE idx_select prefilter with the correct MUTUALLY-EXCLUSIVE
    single-positive selection for the directional ``[target, reference]`` pair.

    The vendored idx_select selects cells HIGH IN BOTH markers — each of its two arms AND-requires
    the target AND the reference to be elevated, i.e. the co-positive top-right corner. For a
    mutually-exclusive pair those are the doublet/spillover/segmentation-error contamination — the
    one population that must NEVER define a negative control. Clustering that corner and taking its
    'negative' sub-cluster leaves the target-negative mean sitting out among the co-positive cells,
    so ``mean+3σ`` lands near the target p99 (on 6374 ND: INS → ~10,144 MFI → 0.6% of islet-core
    cells called INS⁺, i.e. an ND islet with no β-cells — false).

    Corrected: select the two axis-hugging single-positive arms —
        REF+ / TGT≈0  → the NEGATIVE population that thresholds the target
        TGT+ / REF≈0  → signal
    where 'near-zero' for the OFF marker is a per-(image, marker) quantile ``neg_q`` (DYNAMIC knob;
    default median). Downstream 2-cluster + max-reference-std still isolates the reference arm as the
    negative population. The vendored threshold core (and the sigma + cluster-cap patches) is
    otherwise unchanged; discarded figure objects are hv stubs (save_figs=False)."""
    R = sys.modules[Normalization.__module__]
    ratio_x, ratio_y = 0.75, 0.5

    @functools.wraps(Normalization.get_marker_pair_thresh)
    def wrapper(self, data, scene, marker_pair, batch):
        tmp = data[marker_pair].to_numpy()                       # col 0 = target, col 1 = reference
        sel, _c0, _c1 = idx_select(tmp[:, 0], tmp[:, 1], neg_q, ratio_x, ratio_y, nonzero_q)  # single-positive arms
        marker_pair_data = tmp[sel]
        xlabel = marker_pair[0]
        # Item 2 (opt-in, default off): per-image signal-presence guard. If the idx_select cloud has no
        # separated positive population, call the marker ABSENT for this image (thresh=+inf -> 0 positives)
        # rather than letting the 2-cluster fit emit a finite mu+3sigma on pure noise (e.g. beta-loss INS).
        if presence_min_sep > 0:
            cl, lab, neg = fit_clusters(marker_pair_data, "KMeans", subsample, seed)
            negv = cl[lab == neg][:, 0]; posv = cl[lab != neg][:, 0]
            nsd = float(negv.std()) if negv.size else 0.0
            sep = (float(posv.mean()) - float(negv.mean())) / nsd if (nsd > 0 and posv.size) else 0.0
            if sep < presence_min_sep:
                for nm in ("KMeans", "GMM", "SSC"):
                    self.threshs[batch][scene][xlabel][nm] = float("inf")
                print(f"[presence] {scene} {xlabel}: sep={sep:.2f} < {presence_min_sep} -> absent (0 positives)")
                return R.hv.Scatter(marker_pair_data), []
        # Partner-low (opt-in, default off): threshold the TARGET from the reference-LOW subpopulation
        # (where the target is cleanly bimodal, removing the double-positive confound), not the idx_select
        # arms. Deterministic (Otsu) -> one cutoff stored for all 3 model keys. +inf -> absent (beta-loss).
        # Partner-low is only valid where the reference-low target is cleanly bimodal WITHIN ISLETS, i.e.
        # for ENDOCRINE targets (the user's INS/GCG case); non-endocrine markers live in exocrine/immune
        # tissue, not islets, so they fall through to the normal negative-cluster threshold below.
        if partner_low_q > 0 and xlabel in ENDOCRINE_MARKERS:
            islet = (data["cell_region"].isin(("core", "peri")).to_numpy()
                     if "cell_region" in getattr(data, "columns", ()) else None)
            thr = partner_low_thresh(tmp[:, 0], tmp[:, 1], partner_low_q, islet=islet)
            for nm in ("KMeans", "GMM", "SSC"):
                self.threshs[batch][scene][xlabel][nm] = thr
            reg = "islet" if islet is not None else "whole-tissue(NO cell_region)"
            print(f"[partner-low] {scene} {xlabel}: ref-low q{partner_low_q} ({reg}) -> target thr={thr:.1f}")
            return R.hv.Scatter(marker_pair_data), []
        models = (("KMeans", "magenta", R.cluster.KMeans(n_clusters=2)),
                  ("GMM", "blue", R.mixture.GaussianMixture(n_components=2, n_init=10)),
                  ("SSC", "green", R.sr.SparseSubspaceClusteringOMP(n_clusters=2)))
        for name, color, model in models:
            if neg_stat_name == "mean3sd":   # vendored mu+3sigma (production); carries the sigma + cluster-cap patches
                fn = self.get_GMM_thresh if name == "GMM" else self.get_clustering_thresh
                thresh, _ = fn(marker_pair_data, marker_pair, model, 3, color, batch)
            else:                            # Item 1 (opt-in): robust negative-cluster statistic on the target axis
                cl, lab, neg = fit_clusters(marker_pair_data, name, subsample, seed)
                thresh = neg_stat(cl[lab == neg][:, 0], stat=neg_stat_name)
            self.threshs[batch][scene][xlabel][name] = thresh
        return R.hv.Scatter(marker_pair_data), []

    Normalization.get_marker_pair_thresh = wrapper
    extra = [f"neg_stat={neg_stat_name}"] if neg_stat_name != "mean3sd" else []
    if nonzero_q > 0: extra.append(f"nonzero_q={nonzero_q}")
    if presence_min_sep > 0: extra.append(f"presence_min_sep={presence_min_sep}")
    print(f"[idx-select] mutually-exclusive single-positive arms (REF+/TGT- neg | TGT+/REF- signal); "
          f"OFF-marker near-zero ceiling = per-(image,marker) quantile {neg_q} (was co-positive corner)"
          + (f"; {', '.join(extra)}" if extra else ""))


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
    from .cohort import ensure_eligible_donors, filter_eligible_donors

    files = sorted(glob.glob(str(Path(cells_dir) / "donor_id=*" / "*.parquet")))
    pairs = [(re.search(r"donor_id=([^/]+)", f).group(1), f) for f in files]
    if donors:
        keep = set(
            ensure_eligible_donors(donors, context="RESTORE analysis")
        )
        pairs = [p for p in pairs if p[0] in keep]
    else:
        eligible = set(filter_eligible_donors(donor for donor, _ in pairs))
        pairs = [pair for pair in pairs if pair[0] in eligible]
    return pairs[:limit] if limit else pairs


def build_restore_input(cells_dir, markers, subsample, seed, limit_scenes=None, donors=None, region_dir=None):
    """Per-image (per-donor) subsample of raw MFI for threshold estimation.
    Returns a DataFrame with columns: scene (donor id), batch (0), + one per marker. When ``region_dir`` is
    given, also attaches ``cell_region`` (joined on ``object_id`` from ``region_dir``, e.g. ``cfg.cells_dir``)
    so partner-low thresholding can restrict to islet (core/peri) cells."""
    rng = np.random.default_rng(seed)
    frames = []
    read_cols = (["object_id"] + markers) if region_dir is not None else markers
    for donor, f in donor_files(cells_dir, limit_scenes, donors):
        df = pd.read_parquet(f, columns=read_cols)
        if region_dir is not None:
            rf = glob.glob(str(Path(region_dir) / f"donor_id={donor}" / "*.parquet"))
            reg = (pd.read_parquet(rf[0], columns=["object_id", "cell_region"]) if rf
                   else pd.DataFrame({"object_id": [], "cell_region": []}))
            df["object_id"] = df["object_id"].astype(str); reg["object_id"] = reg["object_id"].astype(str)
            df = df.merge(reg, on="object_id", how="left").drop(columns="object_id")
        n_raw = len(df)
        # REDSEA leaves a few all-marker-NaN edge cells; RESTORE's np.quantile is not NaN-aware.
        allnan = df[markers].isna().all(axis=1)
        partial = df[markers].isna().any(axis=1) & ~allnan
        if partial.any():
            print(f"[input] WARN {donor}: {int(partial.sum())} PARTIAL-NaN cells "
                  "(some markers valid) — dropped so the quantile stays finite")
        df = df[~df[markers].isna().any(axis=1)]
        n = len(df)
        if subsample and n > subsample:
            df = df.iloc[rng.choice(n, size=subsample, replace=False)]
        region_col = df["cell_region"].reset_index(drop=True) if region_dir is not None else None
        df = df[markers].astype("float64").reset_index(drop=True)   # numeric markers only for the float cast
        df["scene"] = donor
        df["batch"] = 0
        if region_col is not None:
            df["cell_region"] = region_col
        frames.append(df)
        print(f"[input] image {donor}: {n_raw:,} cells ({n_raw - n:,} NaN dropped) -> {len(df):,} sampled")
    data = pd.concat(frames, ignore_index=True)
    print(f"[input] total RESTORE input: {len(data):,} cells / {data['scene'].nunique()} images")
    return data


# --------------------------------------------------------------------------- #
# Reference exclusivity QC
# --------------------------------------------------------------------------- #
def _idx_select(T, R, ratio_x=0.75, ratio_y=0.5):
    """AUDIT-ONLY co-expression probe (the vendored co-positive selection): its svd_ratio measures
    collinearity IN the co-expressing region (low ratio => the two markers co-vary => poor mutual
    exclusivity). This is intentionally NOT the corrected single-positive threshold selection used by
    ``patch_restore_idx_floor`` — here we WANT to look at the co-positive corner to score the pair."""
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
# RETIRED 2026-07-23 — `run_restore` and `run_restore_functional` are deleted with the curated pair webs
# they consumed (MARKER_PAIRS / FUNCTIONAL_PAIRS, see config.py). Ten of those pairs referenced CD20,
# which is too sparse in pancreas to define a negative-control population, and there is no mechanical
# substitute. The manuscript-faithful replacement lives in `restore_validation.py`; Gate 4 of
# docs/restore_faithful_rebuild_plan.md wires it into production. Everything above this banner --
# idx_select, neg_stat, fit_clusters, the threshold LUT and apply helpers, the QC plots -- is reusable
# machinery and stays.


def run_restore_proliferation(cfg: PipelineConfig, *, markers=None, donors=None, n_jobs=None, k=None):
    """STANDALONE binarization of lineage-agnostic PROLIFERATION markers (Ki67/PCNA). These have no
    mutually-exclusive counterpart AND are GRADED (not bimodal — a zero spike + one broad mode + a right
    tail, no valley), so RESTORE pairing / a bimodal cut do not apply. Instead take a robust per-image
    upper-tail cutoff :func:`logmean_ksigma` (``mean + k*std`` of ``log10(nonzero+1)``, per-donor adaptive,
    no fixed-fraction assumption), then apply ``raw >= thr`` to ALL cells and write ``{m}_pos/_norm/_log2r``
    to ``cfg.restore_gated_proliferation_dir`` + a thresholds/fractions csv. STATE annotations (a
    proliferating cell of whatever lineage), NOT a compartment call. ``k`` defaults to ``cfg.restore_proliferation_k``."""
    import pyarrow.parquet as _pq
    markers = list(markers or PROLIFERATION)
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    k = cfg.restore_proliferation_k if k is None else k
    cells_dir = cfg.cells_redsea_dir
    files = list(donor_files(cells_dir, None, donors))
    if not files:
        print("[proliferation] no donor parquet found"); return {}
    present = set(_pq.ParquetFile(files[0][1]).schema_arrow.names)
    markers = [m for m in markers if m in present]
    if not markers:
        print(f"[proliferation] none of {list(PROLIFERATION)} present in {cells_dir} — skipped"); return {}
    print(f"[proliferation] standalone log-mean+{k}sigma binarization (graded markers): {markers}")
    lut, rows = {}, []
    for donor, f in files:
        df = pd.read_parquet(f, columns=markers)
        for m in markers:
            thr = logmean_ksigma(df[m].to_numpy(dtype="float64"), k)
            lut[(donor, m)] = thr
            rows.append({"image": donor, "marker": m, "threshold": float(thr)})
            print(f"[proliferation] {donor} {m}: thr={thr:.1f}" + ("" if np.isfinite(thr) else "  (too few cells)"))
    gated_dir = cfg.restore_gated_proliferation_dir
    gated_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(cfg.data_dir / "restore_thresholds_proliferation.csv", index=False)
    frac_df = apply_thresholds(cells_dir, markers, lut, {}, gated_dir, None, donors, n_jobs=n_jobs)
    if len(frac_df):
        frac_df.sort_values(["marker", "image"]).to_csv(gated_dir / "proliferation_fractions.csv", index=False)
    print(f"[proliferation] wrote {gated_dir} (partitioned) + restore_thresholds_proliferation.csv")
    return lut


def main(argv=None):
    """CLI for the surviving standalone pass. The paired RESTORE driver was retired with its curated
    pair web (see the banner above); the manuscript-faithful replacement is
    ``python -m phenocycler.restore_validation evaluate-locked``."""
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=None, help="per-donor apply pool size")
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--proliferation", action="store_true",
                    help="STANDALONE binarization of proliferation markers (Ki67/PCNA) — no "
                         "mutually-exclusive pair; writes {m}_pos to restore_gated_proliferation")
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    if not a.proliferation:
        ap.error(
            "the paired RESTORE driver was retired with MARKER_PAIRS/FUNCTIONAL_PAIRS; "
            "use `python -m phenocycler.restore_validation evaluate-locked` for pair validation, "
            "or pass --proliferation for the standalone Ki67/PCNA pass"
        )
    run_restore_proliferation(cfg, donors=a.donors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
