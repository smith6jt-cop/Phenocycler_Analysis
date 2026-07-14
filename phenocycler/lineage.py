#!/usr/bin/env python3
"""
Broad phenotyping — 5-compartment ORDERED RESIDUAL partition (+ explicit "Other").

Every cell is typed by the first matching gate in ``COMPARTMENT_ORDER``, each gate running on the
RESIDUAL of the prior (so the priority order IS the multi-positive tie-break); cells failing every gate
land in ``Other`` (real panel gaps — mast/Schwann/adipocyte/quiescent-stellate — NOT force-assigned).
Priority high->low:

  1. Epithelial  <- E_cadherin_pos      (exocrine + endocrine; the only epithelial marker that stays +
                                         on endocrine). sub: Endocrine (hormone+: beta INS, alpha
                                         GCG, delta SST) vs Exocrine (hormone-; basal/ductal Keratin_5/TP63)
  2. Endothelial <- CD31_pos            (sub: lymphatic Podoplanin+, else vascular)
  3. Neural      <- B3TUBB_pos          (residual, after epithelial -> islet TUBB3 co-expression removed)
  4. Immune      <- ANY(CD3e/CD20/CD79a/CD68/CD163/CD206/Iba1/CD11b/CD11c/MPO)_pos
                     sub: T / B / Plasma / Neutrophil / Macrophage / DC / NK / Immune-other
  5. Mesenchymal <- SMA_pos (Muscle), then Vimentin_pos in the SMA-neg residual (Fibroblast)
  6. Other       <- failed every gate

`_pos` comes straight from RESTORE (no `_norm` floor). Reads ONE RESTORE-gated parquet per donor
(all markers, `{m}_pos`/`{m}_norm`).

Inputs : <restore_gated_dir>/donor_id=*/data_0.parquet   ({m}_pos, {m}_norm)
         <cells_dir>/donor_id=*/data_0.parquet           (cell_region, islet_num)
         <donor_metadata>                                 (disease.status; figure only)
Outputs: <phenotype_dir>/broad/donor_id=*/data_0.parquet  (object_id, donor_id, compartment, cell_type, ...)
         <phenotype_dir>/broad_lineage_composition.png

    python -m phenocycler.lineage --jobs 8
"""

from __future__ import annotations

import argparse
import functools
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
if not matplotlib.get_backend().startswith("module://"):  # keep a notebook's inline/widget backend; force Agg only headless
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import (PipelineConfig, load_config, COMPARTMENT_ORDER, COMPARTMENT_GATES,
                     ENDOCRINE_MARKERS, COMPARTMENT_COLORS, COMPARTMENT_ABBR, STATUS_ORDER, OTHER_LABEL)
from .parallel import map_donors

# 6 rows for composition (5 compartments + Other), in priority order.
CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL]


def status_map(cfg: PipelineConfig) -> dict:
    """donor_id -> disease.status from the donor metadata workbook (best-effort)."""
    try:
        m = pd.read_excel(cfg.donor_metadata)
        return dict(zip(m[cfg.metadata_donor_col].astype(str),
                        m[cfg.metadata_status_col].astype(str)))
    except Exception as e:  # metadata is optional for the assignment itself
        print(f"[warn] donor metadata unavailable ({e}); status will be '?'")
        return {}


def _pos_norm_getters(g: pd.DataFrame):
    """Return (pos, norm) closures giving `{m}_pos` (bool) / `{m}_norm` (float) for any marker, with
    all-False / 0.0 fallbacks for a marker absent from this donor's gated parquet (e.g. batch-specific)."""
    n = len(g)
    zeros_b = np.zeros(n, bool)
    zeros_f = np.zeros(n, np.float64)

    def pos(m: str) -> np.ndarray:
        c = f"{m}_pos"
        return g[c].fillna(False).to_numpy(bool) if c in g.columns else zeros_b

    def norm(m: str) -> np.ndarray:
        c = f"{m}_norm"
        return g[c].fillna(0.0).to_numpy(np.float64) if c in g.columns else zeros_f

    return pos, norm


def _cell_types(compartment: np.ndarray, pos, norm) -> np.ndarray:
    """Finer sub-split within each compartment (secondary annotation). Later writes win, so within a
    compartment the more-specific / higher-priority identity overrides a co-positive lower one."""
    n = len(compartment)
    ct = compartment.astype("<U24").copy()

    # --- Epithelial: endocrine (hormone+) vs exocrine (hormone-) ---
    epi = compartment == "Epithelial"
    horm = np.zeros(n, bool)
    for m in ENDOCRINE_MARKERS:
        horm |= pos(m)
    exo = epi & ~horm
    ct[exo] = "Exocrine"
    ct[exo & (pos("Keratin_5") | pos("TP63"))] = "Exocrine-ductal"
    endo = epi & horm
    # endocrine subtype by dominant hormone _norm: beta = INS, alpha = GCG, delta = SST
    # (IAPP removed 2026-07-10 — failed marker; beta was max(INS, IAPP))
    beta = norm("INS")
    sub = np.array(["Endocrine-beta", "Endocrine-alpha", "Endocrine-delta"])[
        np.argmax(np.vstack([beta, norm("GCG"), norm("SST")]).T, axis=1)]
    ct[endo] = sub[endo]

    # --- Endothelial: vascular vs lymphatic (Podoplanin+) ---
    eth = compartment == "Endothelial"
    ct[eth] = "Endothelial-vascular"
    ct[eth & pos("Podoplanin")] = "Endothelial-lymphatic"

    # --- Neural ---
    ct[compartment == "Neural"] = "Neural"

    # --- Immune (last-write-wins precedence: myeloid -> B/Plasma -> NK -> T, so lymphoid identity wins) ---
    imm = compartment == "Immune"
    ct[imm] = "Immune-other"
    ct[imm & (pos("CD11c") | pos("CD209"))] = "DC"
    ct[imm & (pos("CD68") | pos("CD163") | pos("CD206") | pos("Iba1"))] = "Macrophage"
    ct[imm & pos("MPO")] = "Neutrophil"
    b = imm & (pos("CD20") | pos("CD79a"))
    ct[b] = "B"
    ct[b & pos("CD38")] = "Plasma"                                   # CD38 is batch-1 only -> no Plasma in batch-2
    ct[imm & pos("CD56") & pos("CD57") & ~pos("CD3e")] = "NK"
    ct[imm & pos("CD3e")] = "T"

    # --- Mesenchymal: fibroblast (Vimentin) then muscle (SMA wins — claimed first in the ordered gate) ---
    mes = compartment == "Mesenchymal"
    ct[mes & pos("Vimentin")] = "Fibroblast"
    ct[mes & pos("SMA")] = "Muscle"

    return ct


def assign_donor(donor: str, gated_f: str, cells_dir):
    """Type every cell in one donor by the ordered residual gating tree (+ Other).

    Returns (out_df, mfrac) where mfrac is the {INS,GCG} positivity fraction (% of ALL cells) for the
    disease-trend panel (the Endocrine sub-branch lumps the hormones, masking beta-loss)."""
    g = pd.read_parquet(gated_f)
    g["object_id"] = g["object_id"].astype(str)
    pos, norm = _pos_norm_getters(g)
    n = len(g)

    # ordered residual gating: first matching compartment claims the cell; the rest -> Other.
    compartment = np.array([OTHER_LABEL] * n, dtype="<U16")   # "Endothelial"/"Mesenchymal" = 11 chars
    assigned = np.zeros(n, bool)
    for comp in COMPARTMENT_ORDER:
        hit = np.zeros(n, bool)
        for m in COMPARTMENT_GATES[comp]:
            hit |= pos(m)
        if comp == "Immune":
            # NK cells express NONE of the 10 canonical immune markers, so add the CD56+CD57+ gate to pull
            # them into Immune (they are NOT an "Other" gap per the user). Safe here: epithelial + neural
            # (the CD56 confounders) are already removed by this point in the residual. Sub-typed "NK".
            hit |= (pos("CD56") & pos("CD57"))
        take = (~assigned) & hit
        compartment[take] = comp
        assigned |= take

    cell_type = _cell_types(compartment, pos, norm)

    out = pd.DataFrame({"object_id": g["object_id"].astype(str), "donor_id": donor,
                        "compartment": compartment, "cell_type": cell_type})
    # attach cell_region / islet_num for downstream (spatial / drill-down)
    cf = sorted(glob.glob(str(Path(cells_dir) / f"donor_id={donor}" / "*.parquet")))
    if cf:
        ctx = pd.read_parquet(cf[0], columns=["object_id", "cell_region", "islet_num"])
        out = out.merge(ctx.astype({"object_id": str}), on="object_id", how="left")
    mfrac = {"INS": 100 * pos("INS").mean(), "GCG": 100 * pos("GCG").mean()}
    return out, mfrac


def _assign_and_write(donor: str, gated_dir: str, cells_dir: str, out_dir: str) -> dict:
    """Worker: assign one donor, write its partition, return composition + marker-fraction summary."""
    gf = sorted(glob.glob(str(Path(gated_dir) / f"donor_id={donor}" / "*.parquet")))
    if not gf:
        return {}
    out, mfrac = assign_donor(donor, gf[0], cells_dir)
    dst = Path(out_dir) / f"donor_id={donor}"
    dst.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dst / "data_0.parquet", index=False)
    vc = out["compartment"].value_counts()
    comp = (vc / len(out) * 100).reindex(CLASSES).fillna(0)
    n = len(out); other = int((out["compartment"] == OTHER_LABEL).sum())
    print(f"[{donor}] {n:,} cells | " + " ".join(f"{COMPARTMENT_ABBR[c]}={comp[c]:.0f}%" for c in CLASSES) +
          f" | Other={100*other/n:.1f}%", flush=True)
    return {"donor": donor, "comp": comp.to_dict(), "mfrac": mfrac, "n": n, "other": other}


def run_lineage(cfg: PipelineConfig, *, donors=None, n_jobs=None) -> pd.DataFrame:
    """Assign compartments for all donors (parallel), write partitions + composition figure."""
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    cfg.broad_dir.mkdir(parents=True, exist_ok=True)
    donor_ids = donors or cfg.discover_donors(cfg.restore_gated_dir)
    if not donor_ids:
        raise SystemExit(f"[err] no gated donors under {cfg.restore_gated_dir}")

    fn = functools.partial(_assign_and_write, gated_dir=str(cfg.restore_gated_dir),
                           cells_dir=str(cfg.cells_dir), out_dir=str(cfg.broad_dir))
    results = [r for r in map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True) if r]
    if not results:
        raise SystemExit("[err] no donors assigned")

    smap = status_map(cfg)
    comp = {r["donor"]: pd.Series(r["comp"]).reindex(CLASSES).fillna(0) for r in results}
    mfracs = {r["donor"]: r["mfrac"] for r in results}
    n_total = sum(r["n"] for r in results)
    other_total = sum(r["other"] for r in results)
    C = pd.DataFrame(comp).T
    C["status"] = [smap.get(d, "?") for d in C.index]
    M = pd.DataFrame(mfracs).T
    M["status"] = [smap.get(d, "?") for d in M.index]
    print(f"\n[total] {n_total:,} cells, {100*other_total/max(n_total,1):.1f}% Other "
          "(failed every gate — real panel gaps, NOT force-assigned)")

    _composition_figure(C, M, cfg.phenotype_dir / "broad_lineage_composition.png")
    print("\ncomposition by status (mean %):")
    print(C.groupby("status")[CLASSES].mean().reindex(STATUS_ORDER).round(1).to_string())
    print("\ndisease trend (mean % of cells positive, by status):")
    print(pd.concat([M.groupby("status")[["INS", "GCG"]].mean(),
                     C.groupby("status")["Immune"].mean()], axis=1)
            .reindex(STATUS_ORDER).round(2).to_string())
    return C


def _composition_figure(C: pd.DataFrame, M: pd.DataFrame, path: Path):
    """Per-donor stacked composition (broken y-axis) + per-compartment boxplots by disease status with
    non-parametric stats. Units of replication are donors, not cells."""
    plt.rcParams.update({"font.size": 15, "xtick.labelsize": 14, "ytick.labelsize": 14,
                         "axes.titlesize": 17, "axes.labelsize": 16, "legend.fontsize": 14})
    _rank = {s: i for i, s in enumerate(STATUS_ORDER)}
    Cs = C.sort_values("status", key=lambda col: col.map(lambda s: _rank.get(s, len(STATUS_ORDER))))
    SLABEL = {"ND": "ND", "AAB": "Aab+", "T1D": "T1D"}
    SCOL = {"ND": "#4477AA", "AAB": "#CC6633", "T1D": "#228833"}   # Paul-Tol status colours
    present_status = [s for s in STATUS_ORDER if (C["status"] == s).any()]
    gdata = {s: {L: C.loc[C["status"] == s, L].to_numpy(np.float64) for L in CLASSES} for s in present_status}
    PAIRS = [(present_status[i], present_status[j])
             for i in range(len(present_status)) for j in range(i + 1, len(present_status))]

    def _breakmarks(ax_up, ax_dn):   # diagonal cut marks at the shared broken-axis edge
        d = 0.008
        for ax_, ys in ((ax_up, (-d, +d)), (ax_dn, (1 - d, 1 + d))):
            kw = dict(transform=ax_.transAxes, color="k", clip_on=False, lw=1)
            ax_.plot((-d, +d), ys, **kw); ax_.plot((1 - d, 1 + d), ys, **kw)

    # ---- statistics: per-compartment Kruskal-Wallis omnibus (BH across compartments), then pairwise
    #      Mann-Whitney U (BH across all compartment x pair tests). Brackets only where omnibus survives.
    omni_p = {}
    for L in CLASSES:
        arrs = [gdata[s][L] for s in present_status if len(gdata[s][L]) > 0]
        try:
            omni_p[L] = float(kruskal(*arrs).pvalue) if len(arrs) >= 2 else np.nan
        except ValueError:      # all values identical across groups
            omni_p[L] = np.nan
    _Lok = [L for L in CLASSES if np.isfinite(omni_p[L])]
    omni_q = {L: np.nan for L in CLASSES}
    if _Lok:
        omni_q.update(dict(zip(_Lok, multipletests([omni_p[L] for L in _Lok], method="fdr_bh")[1])))

    pw = []   # (compartment, a, b, p)
    for L in CLASSES:
        for a, b in PAIRS:
            xa, xb = gdata[a][L], gdata[b][L]
            try:
                p = float(mannwhitneyu(xa, xb, alternative="two-sided").pvalue) \
                    if (len(xa) and len(xb) and (len(xa) + len(xb) >= 3)) else np.nan
            except ValueError:
                p = np.nan
            pw.append([L, a, b, p])
    pw = pd.DataFrame(pw, columns=["L", "a", "b", "p"])
    _fin = pw["p"].notna()
    pw["q"] = np.nan
    if _fin.any():
        pw.loc[_fin, "q"] = multipletests(pw.loc[_fin, "p"].to_numpy(), method="fdr_bh")[1]
    qlut = {(r.L, r.a, r.b): r.q for r in pw.itertuples()}

    def _stars(q):
        if not np.isfinite(q):
            return None
        return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else None

    fig = plt.figure(figsize=(23, 9))
    outer = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.95], wspace=0.13)
    gsA = outer[0].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.06)
    a_top, a_bot = fig.add_subplot(gsA[0]), fig.add_subplot(gsA[1])
    gsB = outer[1].subgridspec(2, 3, hspace=0.62, wspace=0.34)
    bax = {L: fig.add_subplot(gsB[i // 3, i % 3]) for i, L in enumerate(CLASSES)}

    # A: per-donor stacked composition — split y-axis (Epithelial dominates the low range; the split
    # expands the top so the minority compartments are legible). Ordered ND->Aab+->T1D.
    xA = np.arange(len(Cs)); bottom = np.zeros(len(Cs)); handles = []
    for L in CLASSES:
        h = a_top.bar(xA, Cs[L], bottom=bottom, color=COMPARTMENT_COLORS[L], width=0.85)
        a_bot.bar(xA, Cs[L], bottom=bottom, color=COMPARTMENT_COLORS[L], width=0.85)
        bottom += Cs[L].to_numpy(); handles.append(h)
    a_hi = max(15.0, float(np.floor(Cs["Epithelial"].min())) - 4)      # break below the lowest Epithelial
    a_top.set_ylim(a_hi, 100); a_bot.set_ylim(0, 8)
    a_top.spines["bottom"].set_visible(False); a_bot.spines["top"].set_visible(False)
    a_top.tick_params(bottom=False, labelbottom=False); _breakmarks(a_top, a_bot)
    a_bot.set_xticks(xA)   # rotate donor labels vertical so many large labels never overlap
    a_bot.set_xticklabels([f"{d}  {Cs.loc[d,'status']}" for d in Cs.index], rotation=90, fontsize=12)

    # B: per-compartment box plots by disease status — box = median/IQR, diamond = mean +/- SEM, dots =
    # individual donors (deterministic jitter). Significant pairwise differences bracketed above.
    for L in CLASSES:
        ax = bax[L]
        data = [gdata[s][L] for s in present_status]
        bp = ax.boxplot(data, positions=list(range(len(present_status))), widths=0.6,
                        patch_artist=True, showfliers=False, showmeans=False, zorder=1,
                        medianprops=dict(color="black", lw=1.3))
        for patch, s in zip(bp["boxes"], present_status):
            patch.set_facecolor(SCOL[s]); patch.set_alpha(0.30); patch.set_edgecolor(SCOL[s]); patch.set_linewidth(1.2)
        for grp in ("whiskers", "caps"):
            for ln in bp[grp]:
                ln.set_color("#555555")
        for i, s in enumerate(present_status):
            y = data[i]
            if len(y) == 0:
                continue
            jit = np.linspace(-0.17, 0.17, len(y)) if len(y) > 1 else np.array([0.0])
            ax.scatter(i + jit, y, s=24, color=SCOL[s], edgecolor="black", lw=0.4, alpha=0.9, zorder=3)
            m = float(np.mean(y)); sem = float(np.std(y, ddof=1) / np.sqrt(len(y))) if len(y) > 1 else 0.0
            ax.errorbar(i, m, yerr=sem, fmt="D", color="black", ms=6.5, mfc="white", mec="black",
                        capsize=4, elinewidth=1.4, zorder=4)
        ax.set_xticks(range(len(present_status)))
        ax.set_xticklabels([SLABEL[s] for s in present_status], fontsize=12)
        if CLASSES.index(L) % 3 == 0:
            ax.set_ylabel("% of cells", fontsize=13)
        ax.set_title(L, fontsize=14, color=COMPARTMENT_COLORS[L])

        allv = np.concatenate([d for d in data if len(d)]) if any(len(d) for d in data) else np.array([0.0, 1.0])
        lo, hi = float(allv.min()), float(allv.max())
        span = (hi - lo) if hi > lo else max(1.0, 0.1 * hi)
        pad = 0.15 * span
        step = 0.16 * span; base = hi + 0.10 * span; level = 0
        if np.isfinite(omni_q.get(L, np.nan)) and omni_q[L] < 0.05:
            for a, b in PAIRS:
                st = _stars(qlut.get((L, a, b), np.nan))
                if st is None:
                    continue
                x1, x2 = present_status.index(a), present_status.index(b)
                yb = base + level * step
                ax.plot([x1, x1, x2, x2], [yb, yb + 0.03 * span, yb + 0.03 * span, yb], color="black", lw=1.1)
                ax.text((x1 + x2) / 2, yb + 0.035 * span, st, ha="center", va="bottom", fontsize=13)
                level += 1
        ax.set_ylim(lo - pad, (base + (level + 0.4) * step) if level else hi + pad)

    ap_t, ap_b = a_top.get_position(), a_bot.get_position()
    fig.text(ap_t.x0 - 0.032, (ap_b.y0 + ap_t.y1) / 2, "% of cells",
             rotation="vertical", va="center", ha="center", fontsize=16)
    a_bot.legend(handles=handles, labels=CLASSES, loc="upper center", ncol=3,
                 fontsize=13, frameon=False, bbox_to_anchor=(0.5, -0.9))
    ap = a_top.get_position()
    bp0 = bax[CLASSES[0]].get_position(); bp2 = bax[CLASSES[2]].get_position()
    ty = ap.y1 + 0.05
    fig.text((ap.x0 + ap.x1) / 2, ty, "Compartment composition per donor",
             ha="center", va="bottom", fontsize=17)
    fig.text((bp0.x0 + bp2.x1) / 2, ty, "Composition by disease status (per-donor)",
             ha="center", va="bottom", fontsize=17)
    fig.text((bp0.x0 + bp2.x1) / 2, 0.045,
             "Box = median/IQR . whiskers = range . diamond = mean +/- SEM . dots = donors.",
             ha="center", va="center", fontsize=11)
    fig.suptitle("Cell Phenotyping: 5 compartments + Other (ordered residual gating)", fontsize=20, y=1.0)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {path}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=None, help="per-donor process pool size")
    ap.add_argument("--donors", nargs="*", default=None)
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    run_lineage(cfg, donors=a.donors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
