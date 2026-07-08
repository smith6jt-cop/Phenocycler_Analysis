#!/usr/bin/env python3
"""
REDSEA reassessment — the diagnostic / acceptance yardstick (read-only, per-donor, islet-context-stratified).

Faithful port of ``scripts/senior/redsea_reassess_diag.py`` (Islet-Explorer-Senior), made config-driven: the
hardcoded absolute paths are replaced by a :class:`~phenocycler.config.PipelineConfig` and the module gains an
in-process ``run_reassess_diag(cfg)`` entry point. It stays READ-ONLY — it reads the pipeline outputs and
writes only its own diagnostics (CSV/PNG) under ``cfg.redsea_reassess_dir``; it modifies no pipeline data.
ALL analysis / metric logic (thresholds, strata, objective / guardrail / real-vs-false definitions) is
unchanged from the upstream script.

The Step-1 identity heatmap shows Endocrine cells lit up for the epithelial keratin block
(Pan_Cytokeratin/Ker8_18/Keratin_5) and for Vimentin — residual lateral spillover sitting on cells that are
*correctly* typed Endocrine. That contamination is SPATIAL: lone/scattered endocrine (the largest endocrine
stratum) sit inside Pan_CK+ acinar tissue and carry several-fold the keratin of large-islet-core endocrine.

This module quantifies that, per donor, stratified by islet context, and is the ACCEPTANCE YARDSTICK for any
change to REDSEA (`redsea_full.py`). It reports three things:

  OBJECTIVE (what we want to reduce):
    endocrine-typed mean CORRECTED keratin(3)+Vimentin, per stratum, as an absolute and as a RATIO to that
    donor's Epithelial reference; plus the impossible-co-expression rate = fraction of CONFIDENT endocrine
    cells that are Pan_Cytokeratin_pos (a confident β/α/δ cell should not be keratin-positive).

  GUARDRAILS (what must NOT regress under any lever):
    confident-acinar keratin retention (Epithelial, far-from-islet, Pan_CK+ cells keep their keratin);
    confident-vessel Vimentin retention (CD31+/SMA+ cells keep their Vimentin — CD31+SMA+Vim+ is a real
    vessel wall, not spillover); Fibroblast-compartment size; T1D hallmarks (β-loss INS⁺ + T-gain CD3e⁺
    across ND→T1D).

  REAL-vs-FALSE endocrine:
    split endocrine-typed cells into CONFIDENT (donor-percentile-high hormone `_norm`) vs MARGINAL (near the
    threshold) and report keratin load + morphology (solidity, max-diameter, N:C ratio) per stratum, so a
    segmentation doublet / spillover false-positive is never "fixed" with α.

Frequency-free by construction (repo working rule): "confident endocrine/acinar" is a per-donor PERCENTILE,
never a cohort constant, and every metric is reported per donor (6476 pancreatitis keeps its real extreme
CD3e — validated on its own data, not normalized toward the cohort).

The metric functions (`endocrine_objective`, `guardrails`, `realvsfalse`, `assign_stratum`) are importable so
`redsea_full_validate.py` applies the SAME definitions to α-swept, re-RESTORE'd data — one source of truth.

Inputs (read-only, all present, 22 donors):
  data/cells_redsea/donor_id=*        corrected MFI (Pan_Cytokeratin, Ker8_18, Keratin_5, Vimentin, ...)
  data/restore_gated_redsea/donor_id=* {m}_pos / {m}_norm for the 10 RESTORE markers
  data/cells/donor_id=*               cell_region, islet_num, dist_islet + morphology
  data/phenotype/broad/donor_id=*     broad_lineage (Epithelial reference, Fibroblast size)
Outputs: data/redsea_reassess/{endocrine_by_stratum,guardrails,realvsfalse}.csv + redsea_reassess_baseline.png

Run in `scvi-env`:  python -m phenocycler.reassess_diag                    # all donors under restore_gated_redsea
                    python -m phenocycler.reassess_diag --donors 6539 6521 6436
"""
from __future__ import annotations
import argparse
import glob
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config import PipelineConfig, load_config

# --- the aggressor channels we track on endocrine cells (corrected MFI column names) ---
KERATINS = ["Pan_Cytokeratin", "Ker8_18", "Keratin_5"]
SPILL = KERATINS + ["Vimentin"]
HORMONES = ["INS", "GCG", "SST"]

# --- stratification knobs ---
CONF_PCT = 75          # donor percentile of hormone strength defining a "confident" endocrine cell
LARGE_CORE = 60        # islet size (core-cell count) thresholds
SMALL_CORE = 15
FAR_DIST_PCT = 90      # donor percentile of dist_islet defining "far from islet" (confident acinar background)
STRATA = ["core_large", "core_med", "core_small", "peri", "lone_tissue"]
STATUS_ORDER = ["ND", "AAB", "T1D"]

# --- false-endocrine knobs (hormone _norm is threshold-relative: 1.0 == exactly at the RESTORE threshold) ---
CD99_BRIGHT = 3.0      # bright CD99 (>= 3x threshold) = PP/eps/EC endocrine we lack specific markers for
#                        NB: must equal marker_taxonomy.CD99_BRIGHT (the single source of truth for the gate).
WEAK_NORM = 2.0        # hormone _norm below this = near-threshold, noise-floor -> LIKELY FALSE
STRONG_NORM = 3.0      # hormone _norm at/above this = a real, separated positive -> real endocrine


def status_map(cfg: PipelineConfig) -> dict:
    m = pd.read_excel(cfg.donor_metadata)
    return dict(zip(m[cfg.metadata_donor_col].astype(str), m[cfg.metadata_status_col].astype(str)))


def _read(path_dir: Path, donor: str, cols):
    f = sorted(glob.glob(str(path_dir / f"donor_id={donor}" / "*.parquet")))
    if not f:
        return None
    df = pd.read_parquet(f[0], columns=cols)
    df["object_id"] = df["object_id"].astype(str)
    return df


def assign_stratum(cell_region, islet_size) -> np.ndarray:
    """core_large/med/small (by core-cell count of the cell's islet) · peri · lone_tissue."""
    reg = np.asarray(cell_region)
    sz = np.asarray(islet_size)
    out = np.full(len(reg), "lone_tissue", dtype="<U12")
    core = reg == "core"
    out[core & (sz > LARGE_CORE)] = "core_large"
    out[core & (sz > SMALL_CORE) & (sz <= LARGE_CORE)] = "core_med"
    out[core & (sz <= SMALL_CORE)] = "core_small"
    out[reg == "peri"] = "peri"
    return out


def load_donor(donor: str, cfg: PipelineConfig, corrected: pd.DataFrame | None = None,
               gated: pd.DataFrame | None = None) -> pd.DataFrame | None:
    """Merge corrected MFI + RESTORE _pos/_norm + cell context/morphology + broad lineage on object_id and
    add `islet_size`/`stratum`. `corrected`/`gated` override the on-disk α=1 versions (the validator passes
    α-swept / re-RESTORE'd frames keyed by object_id) — pass the SAME column names as the on-disk parquet."""
    red = corrected if corrected is not None else _read(
        cfg.cells_redsea_dir, donor, ["object_id"] + SPILL + ["INS", "GCG", "SST", "CD3e", "CD31", "SMA"])
    if red is None:
        return None
    red = red.copy(); red["object_id"] = red["object_id"].astype(str)
    gcols = ["object_id"] + [f"{m}_pos" for m in HORMONES + ["Pan_Cytokeratin", "CD3e", "CD31", "SMA", "Vimentin"]] \
        + [f"{m}_norm" for m in HORMONES]
    gat = gated if gated is not None else _read(cfg.restore_gated_dir, donor, gcols)
    if gat is None:
        return None
    gat = gat.copy(); gat["object_id"] = gat["object_id"].astype(str)
    ctx = _read(cfg.cells_dir, donor, ["object_id", "cell_region", "islet_num", "dist_islet",
                                       "cell_solidity", "cell_maxdiam", "nc_area_ratio"])
    bro = _read(cfg.broad_dir, donor, ["object_id", "broad_lineage", "epi_default"])
    if ctx is None or bro is None:
        return None
    df = red.merge(gat, on="object_id").merge(ctx, on="object_id").merge(bro, on="object_id", how="left")

    # optional bright-CD99 endocrine (PP/eps/EC) from the extra gating; NaN when not yet produced
    cd99 = _read(cfg.restore_gated_extra_dir, donor, ["object_id", "CD99_norm"])
    df = df.merge(cd99, on="object_id", how="left") if cd99 is not None else df.assign(CD99_norm=np.nan)

    core_sz = df.loc[df["cell_region"] == "core"].groupby("islet_num").size()
    df["islet_size"] = df["islet_num"].map(core_sz).fillna(0).astype(int)
    df["stratum"] = assign_stratum(df["cell_region"].to_numpy(), df["islet_size"].to_numpy())
    df["dist_islet"] = pd.to_numeric(df["dist_islet"], errors="coerce")
    return df


# --------------------------------------------------------------------------- shared derived masks
def _endocrine_typed(df) -> np.ndarray:
    """INS/GCG/SST _pos OR bright CD99 (>= CD99_BRIGHT x threshold) = PP/eps/EC endocrine."""
    horm = df["INS_pos"].to_numpy(bool) | df["GCG_pos"].to_numpy(bool) | df["SST_pos"].to_numpy(bool)
    c = df["CD99_norm"].to_numpy(float) if "CD99_norm" in df.columns else np.full(len(df), np.nan)
    cd99_bright = np.isfinite(c) & (c >= CD99_BRIGHT)
    return horm | cd99_bright


def _hormone_strength(df) -> np.ndarray:
    h = df[[f"{m}_norm" for m in HORMONES]].to_numpy(float)
    with warnings.catch_warnings():           # all-NaN rows -> NaN (not an error)
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmax(h, axis=1)


def _conf_threshold(hstr, endo) -> float:
    v = hstr[endo]
    v = v[np.isfinite(v)]
    return float(np.percentile(v, CONF_PCT)) if v.size > 5 else np.nan


# --------------------------------------------------------------------------- OBJECTIVE
def endocrine_objective(df) -> pd.DataFrame:
    """Per stratum (+ 'all'), over endocrine-typed cells: mean corrected keratin(3)+Vim, ratio to the donor
    Epithelial reference, and the impossible-co-expression rate (confident endocrine that are Pan_CK+)."""
    endo = _endocrine_typed(df)
    epi = (df["broad_lineage"] == "Epithelial").to_numpy()
    epi_ref = {m: (float(np.nanmean(df.loc[epi, m])) if epi.sum() else np.nan) for m in SPILL}
    hstr = _hormone_strength(df)
    thr = _conf_threshold(hstr, endo)
    conf = endo & (hstr >= thr)
    strat = df["stratum"].to_numpy()
    rows = []
    for s in STRATA + ["all"]:
        sel = np.ones(len(df), bool) if s == "all" else (strat == s)
        m, cm = endo & sel, conf & sel
        r = {"stratum": s, "n_endo": int(m.sum()), "n_conf": int(cm.sum())}
        for k in SPILL:
            r[f"endo_{k}"] = float(np.nanmean(df.loc[m, k])) if m.sum() else np.nan
            r[f"ratio_{k}"] = (r[f"endo_{k}"] / epi_ref[k]) if epi_ref[k] else np.nan
        r["conf_PanCK_pos_frac"] = float(df.loc[cm, "Pan_Cytokeratin_pos"].mean()) if cm.sum() else np.nan
        rows.append(r)
    out = pd.DataFrame(rows)
    for k in SPILL:
        out[f"epi_ref_{k}"] = epi_ref[k]
    return out


# --------------------------------------------------------------------------- GUARDRAILS
def guardrails(df) -> dict:
    """Donor-level quantities that must NOT regress under a lever."""
    epi = (df["broad_lineage"] == "Epithelial").to_numpy()
    di = df["dist_islet"].to_numpy(float)
    far = di >= np.nanpercentile(di, FAR_DIST_PCT)
    acinar = epi & (df["cell_region"].to_numpy() == "tissue") & far & df["Pan_Cytokeratin_pos"].to_numpy(bool)
    vessel = (df["CD31_pos"].to_numpy(bool) | df["SMA_pos"].to_numpy(bool))
    g = {
        "acinar_n": int(acinar.sum()),
        "acinar_PanCK": float(np.nanmean(df.loc[acinar, "Pan_Cytokeratin"])) if acinar.sum() else np.nan,
        "acinar_Ker8_18": float(np.nanmean(df.loc[acinar, "Ker8_18"])) if acinar.sum() else np.nan,
        "vessel_n": int(vessel.sum()),
        "vessel_Vimentin": float(np.nanmean(df.loc[vessel, "Vimentin"])) if vessel.sum() else np.nan,
        "vessel_Vim_pos_frac": float(df.loc[vessel, "Vimentin_pos"].mean()) if vessel.sum() else np.nan,
        "fibroblast_frac": float((df["broad_lineage"] == "Fibroblast").mean()),
        "INS_pos_frac": float(df["INS_pos"].mean()),      # β / insulin (expect ND>T1D)
        "CD3e_pos_frac": float(df["CD3e_pos"].mean()),    # T cells (expect T1D>ND)
        "n_cells": int(len(df)),
    }
    return g


# --------------------------------------------------------------------------- REAL vs FALSE endocrine
def realvsfalse(df) -> pd.DataFrame:
    """Per stratum: confident (hormone-high) vs marginal (near-threshold) endocrine — keratin load +
    morphology, so spillover false-positives / segmentation doublets are distinguishable from genuine cells."""
    endo = _endocrine_typed(df)
    hstr = _hormone_strength(df)
    thr = _conf_threshold(hstr, endo)
    strat = df["stratum"].to_numpy()
    def med(mask, col):
        v = pd.to_numeric(df.loc[mask, col], errors="coerce")
        return float(np.nanmedian(v)) if mask.sum() else np.nan
    rows = []
    for s in STRATA:
        base = endo & (strat == s)
        for grp, mask in [("confident", base & (hstr >= thr)), ("marginal", base & (hstr < thr))]:
            if mask.sum() == 0:
                continue
            rows.append({"stratum": s, "group": grp, "n": int(mask.sum()),
                         "hormone_norm_med": float(np.nanmedian(hstr[mask])),
                         "PanCK_corr_med": med(mask, "Pan_Cytokeratin"),
                         "PanCK_pos_frac": float(df.loc[mask, "Pan_Cytokeratin_pos"].mean()),
                         "solidity_med": med(mask, "cell_solidity"),
                         "maxdiam_med": med(mask, "cell_maxdiam"),
                         "nc_ratio_med": med(mask, "nc_area_ratio")})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- FALSE-ENDOCRINE burden
HORMONE_MARKERS = ["INS", "GCG", "SST"]


def hormone_false_burden(df) -> pd.DataFrame:
    """Per hormone, the false-positive burden + its discriminators — the yardstick a call-fix must move.
    `_norm` is threshold-relative (1.0 == exactly at the RESTORE threshold), so a call sitting at ~1 is a
    noise-floor call. Real endocrine sit well above (ND INS median ~8-25); false ones (acinar with a
    threshold-in-the-noise call) sit at ~1 AND carry acinar keratin AND scatter outside islets.
      strong_frac  = n_strong/n_pos  — signal-presence: high => a real separated population exists
      sep_p99_p90  = p99/p90 of corrected intensity — a separated bright population => >> 1 (ND ~80, β-loss ~2)
      pct_pos_in_islet — real endocrine cluster in islets; false scatter in tissue"""
    rows = []
    region = df["cell_region"].to_numpy()
    in_islet = (region == "core") | (region == "peri")
    for h in HORMONE_MARKERS:
        pos = df[f"{h}_pos"].to_numpy(bool)
        norm = df[f"{h}_norm"].to_numpy(float)
        val = df[h].to_numpy(float) if h in df.columns else np.full(len(df), np.nan)
        npos = int(pos.sum())
        weak = pos & (norm < WEAK_NORM)
        strong = pos & (norm >= STRONG_NORM)
        fin = np.isfinite(val)
        p90 = float(np.percentile(val[fin], 90)) if fin.any() else np.nan
        p99 = float(np.percentile(val[fin], 99)) if fin.any() else np.nan
        rows.append({
            "hormone": h, "n_pos": npos,
            "n_weak": int(weak.sum()), "n_strong": int(strong.sum()),
            "strong_frac": (strong.sum() / npos) if npos else np.nan,
            "ker_weak": float(np.nanmedian(df.loc[weak, "Pan_Cytokeratin"])) if weak.sum() else np.nan,
            "ker_strong": float(np.nanmedian(df.loc[strong, "Pan_Cytokeratin"])) if strong.sum() else np.nan,
            "sep_p99_p90": (p99 / p90) if (np.isfinite(p90) and p90 > 0) else np.nan,
            "pct_pos_in_islet": float(in_islet[pos].mean()) if npos else np.nan,
            "pct_strong_in_islet": float(in_islet[strong].mean()) if strong.sum() else np.nan,
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- driver
def compute_all(donors, cfg: PipelineConfig):
    smap = status_map(cfg)
    obj, gr, rvf, fbs = [], [], [], []
    for d in donors:
        df = load_donor(d, cfg)
        if df is None:
            print(f"[{d}] SKIP: missing inputs", flush=True)
            continue
        st = smap.get(d, "?")
        o = endocrine_objective(df); o.insert(0, "donor", d); o.insert(1, "status", st); obj.append(o)
        g = guardrails(df); g["donor"] = d; g["status"] = st; gr.append(g)
        r = realvsfalse(df); r.insert(0, "donor", d); r.insert(1, "status", st); rvf.append(r)
        fb = hormone_false_burden(df); fb.insert(0, "donor", d); fb.insert(1, "status", st); fbs.append(fb)
        ins = fb[fb.hormone == "INS"].iloc[0]
        print(f"[{d}] {st:3s} INS+ n={ins.n_pos:>6} | weak(false)={ins.n_weak:>6} strong(real)={ins.n_strong:>6} "
              f"strong_frac={ins.strong_frac:.2f} sep(p99/p90)={ins.sep_p99_p90:.1f} "
              f"in_islet={100*ins.pct_pos_in_islet:.0f}% | ker weak/strong={ins.ker_weak:.0f}/{ins.ker_strong:.0f}",
              flush=True)
    return (pd.concat(obj, ignore_index=True) if obj else pd.DataFrame(),
            pd.DataFrame(gr), pd.concat(rvf, ignore_index=True) if rvf else pd.DataFrame(),
            pd.concat(fbs, ignore_index=True) if fbs else pd.DataFrame())


def make_figure(obj, gr, path):
    fig, ax = plt.subplots(2, 2, figsize=(15, 11))
    sx = {s: i for i, s in enumerate(STRATA)}

    # A: endocrine corrected Pan_CK by stratum, one line per donor (the spatial gradient)
    for d, sub in obj[obj.stratum.isin(STRATA)].groupby("donor"):
        sub = sub.set_index("stratum").reindex(STRATA)
        ax[0, 0].plot([sx[s] for s in STRATA], sub["endo_Pan_Cytokeratin"], "o-", alpha=0.5, lw=1)
    ax[0, 0].set_xticks(range(len(STRATA))); ax[0, 0].set_xticklabels(STRATA, rotation=30, ha="right")
    ax[0, 0].set_ylabel("endocrine mean corrected Pan_Cytokeratin")
    ax[0, 0].set_title("Keratin on endocrine cells is SPATIAL\n(one line per donor)")

    # B: endocrine keratin/Vim ratio-to-Epithelial by stratum (pooled median across donors)
    piv = obj[obj.stratum.isin(STRATA)].groupby("stratum")[[f"ratio_{k}" for k in SPILL]].median().reindex(STRATA)
    w = 0.2
    for j, k in enumerate(SPILL):
        ax[0, 1].bar(np.arange(len(STRATA)) + (j - 1.5) * w, piv[f"ratio_{k}"], w, label=k)
    ax[0, 1].axhline(1.0, color="k", ls="--", lw=1)
    ax[0, 1].set_xticks(range(len(STRATA))); ax[0, 1].set_xticklabels(STRATA, rotation=30, ha="right")
    ax[0, 1].set_ylabel("endocrine / Epithelial-ref (median over donors)")
    ax[0, 1].set_title("Endocrine spillover relative to the epithelial reference"); ax[0, 1].legend(fontsize=8)

    # C: guardrail baselines to protect — acinar Pan_CK & vessel Vimentin per donor
    gs = gr.sort_values("status")
    x = np.arange(len(gs))
    ax[1, 0].bar(x - 0.2, gs["acinar_PanCK"], 0.4, label="acinar Pan_CK (protect)", color="#4477AA")
    ax[1, 0].bar(x + 0.2, gs["vessel_Vimentin"], 0.4, label="vessel Vimentin (protect)", color="#EE6677")
    ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels([f"{d}\n{s}" for d, s in zip(gs.donor, gs.status)], fontsize=7)
    ax[1, 0].set_ylabel("mean corrected MFI"); ax[1, 0].set_title("Guardrail baselines (must not drop)")
    ax[1, 0].legend(fontsize=8)

    # D: T1D hallmarks by status — MEDIAN (robust to the 6476 pancreatitis outlier, shown separately)
    for col, lab, mk in [("INS_pos_frac", "INS⁺ (β, expect ↓)", "o-"), ("CD3e_pos_frac", "CD3e⁺ (T, expect ↑)", "s-")]:
        y = [gr[gr.status == s][col].median() for s in STATUS_ORDER]
        ax[1, 1].plot(STATUS_ORDER, y, mk, lw=2, ms=9, label=lab)
    m6476 = gr[gr.donor.astype(str) == "6476"]
    if len(m6476):
        ax[1, 1].scatter([0], [m6476.CD3e_pos_frac.iloc[0]], marker="*", s=200, color="#CC3311", zorder=5)
        ax[1, 1].annotate("6476 pancreatitis\n(real CD3e, excluded\nfrom the trend)",
                          (0, m6476.CD3e_pos_frac.iloc[0]), fontsize=7, ha="left", va="top", color="#CC3311")
    ax[1, 1].set_ylabel("median fraction of cells")
    ax[1, 1].set_title("T1D hallmarks (median): β-loss + CD3e T-gain"); ax[1, 1].legend()

    fig.suptitle("REDSEA reassessment — baseline (α=1) endocrine keratin/Vimentin spillover + guardrails", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[saved] {path}")


def make_false_endocrine_figure(fb, path):
    """Per-donor INS false-positive burden + the three discriminators (the acceptance yardstick)."""
    from matplotlib.patches import Patch
    ins = fb[fb.hormone == "INS"].sort_values(["status", "strong_frac"])
    x = np.arange(len(ins))
    col = {"ND": "#4477AA", "AAB": "#CC6633", "T1D": "#228833"}
    colors = [col.get(s, "#888888") for s in ins.status]
    lbl = [f"{d}\n{s}" for d, s in zip(ins.donor, ins.status)]
    fig, ax = plt.subplots(1, 3, figsize=(21, 6))

    ax[0].bar(x, ins.n_weak, 0.8, color=colors, alpha=0.4)
    ax[0].bar(x, ins.n_strong, 0.8, bottom=ins.n_weak.to_numpy(), color=colors)
    ax[0].set_xticks(x); ax[0].set_xticklabels(lbl, fontsize=7); ax[0].set_ylabel("INS+ cell count")
    ax[0].set_title("INS+ calls: pale=weak (_norm<2, likely FALSE), solid=strong (_norm>=3, real β)")
    ax[0].legend(handles=[Patch(fc="grey", alpha=0.4, label="weak / noise-floor (FALSE)"),
                          Patch(fc="grey", label="strong / real β")], fontsize=8)

    ax[1].bar(x, 100 * ins.strong_frac, 0.8, color=colors)
    ax[1].plot(x, 100 * ins.pct_pos_in_islet, "k.-", lw=1, ms=8, label="% of INS+ inside an islet")
    ax[1].set_xticks(x); ax[1].set_xticklabels(lbl, fontsize=7); ax[1].set_ylabel("%")
    ax[1].set_title("Signal presence: strong_frac (bars) + % in islet (line)\nlow both ⇒ noise-floor false calls (β-loss)")
    ax[1].legend(fontsize=8)

    w = 0.4
    ax[2].bar(x - w / 2, ins.ker_weak, w, color="#CC3311", label="weak INS+ keratin (acinar)")
    ax[2].bar(x + w / 2, ins.ker_strong, w, color="#4477AA", label="strong INS+ keratin (real β)")
    ax[2].set_xticks(x); ax[2].set_xticklabels(lbl, fontsize=7); ax[2].set_ylabel("median corrected Pan_CK")
    ax[2].set_title("Weak INS+ carry acinar keratin; strong (real β) are keratin-clean")
    ax[2].legend(fontsize=8)

    fig.suptitle("False-endocrine diagnostic (INS): weak-_norm calls are acinar false-positives, prevalent in β-loss donors",
                 fontsize=13)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"[saved] {path}")


def run_reassess_diag(cfg: PipelineConfig, *, donors=None, figure=True, out_dir=None):
    """Read-only REDSEA reassessment for `donors` (default: all under ``cfg.restore_gated_dir``).

    Writes ``endocrine_by_stratum.csv``, ``guardrails.csv``, ``realvsfalse.csv``, ``false_endocrine.csv`` and
    (when `figure`) the two diagnostic PNGs under `out_dir` (default ``cfg.redsea_reassess_dir``). Modifies no
    pipeline data. Returns the four computed frames ``(obj, gr, rvf, fb)``.
    """
    donor_ids = donors or cfg.discover_donors(cfg.restore_gated_dir)
    out = Path(out_dir) if out_dir is not None else cfg.redsea_reassess_dir
    out.mkdir(parents=True, exist_ok=True)
    obj, gr, rvf, fb = compute_all(donor_ids, cfg)
    if obj.empty:
        raise SystemExit("no donors processed")
    obj.to_csv(out / "endocrine_by_stratum.csv", index=False)
    gr.to_csv(out / "guardrails.csv", index=False)
    rvf.to_csv(out / "realvsfalse.csv", index=False)
    fb.to_csv(out / "false_endocrine.csv", index=False)
    print(f"[saved] {out}/endocrine_by_stratum.csv, guardrails.csv, realvsfalse.csv, false_endocrine.csv")

    # false-endocrine burden: INS weak(false) vs strong(real) by status (median over donors)
    print("\n=== FALSE-ENDOCRINE burden — INS by status (median over donors) ===")
    ins = fb[fb.hormone == "INS"]
    print(ins.groupby("status")[["n_weak", "n_strong", "strong_frac", "sep_p99_p90",
                                 "pct_pos_in_islet", "ker_weak", "ker_strong"]].median().round(2).to_string())
    # hallmarks by status — MEDIAN (robust to the 6476 pancreatitis outlier that distorts the mean)
    print("\n=== T1D hallmarks by status (MEDIAN — robust to the 6476 pancreatitis outlier) ===")
    print(gr.groupby("status")[["INS_pos_frac", "CD3e_pos_frac"]].median().round(4).to_string())
    m6476 = gr[gr.donor.astype(str) == "6476"]
    if len(m6476):
        print(f"   (6476 ND pancreatitis CD3e+={m6476.CD3e_pos_frac.iloc[0]:.3f} — real, shown separately, not T1D immune gain)")
    if figure:
        make_figure(obj, gr, out / "redsea_reassess_baseline.png")
        make_false_endocrine_figure(fb, out / "false_endocrine_diag.png")
    return obj, gr, rvf, fb


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--donors", nargs="*", default=None, help="donor ids (default: all under restore_gated_redsea)")
    ap.add_argument("--out-dir", type=Path, default=None, help="output dir (default cfg.redsea_reassess_dir)")
    ap.add_argument("--no-figure", action="store_true")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    run_reassess_diag(cfg, donors=a.donors, figure=not a.no_figure, out_dir=a.out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
