#!/usr/bin/env python3
"""
Broad phenotyping — 5-compartment ORDERED RESIDUAL partition over TRI-STATE RESTORE calls
(+ explicit ``Other`` and ``Unresolved``).

Every cell carries, per gate marker, a tri-state RESTORE call ``positive`` | ``negative`` |
``unavailable`` (emitted by ``restore_apply`` from the frozen ``ACCEPTED_PAIRS`` — one per gate marker).
Each cell is typed by the first POSITIVE gate in ``COMPARTMENT_ORDER`` (each gate on the RESIDUAL of the
prior, so the priority order IS the multi-positive tie-break). The first gate that is UNAVAILABLE
(indeterminate — a gate marker has no valid call for that cell, so the compartment can be neither
confirmed nor ruled out) BLOCKS assignment -> ``Unresolved``: no lower-priority gate may claim a cell we
could not rule out of a higher-priority compartment. A cell valid-and-NEGATIVE for every gate ->
``Other`` (a real panel gap — mast/Schwann/adipocyte/quiescent-stellate — NOT force-assigned).
Priority high->low:

  1. Immune      <- ANY(CD3e, CD68, CD11b) positive   (T / macrophage / myeloid — the 3 screened anchors)
  2. Epithelial  <- E_cadherin positive               (the only epithelial marker that stays + on endocrine)
  3. Endothelial <- CD31 positive
  4. Neural      <- B3TUBB positive                    (residual, after epithelial -> islet TUBB3 removed)
  5. Mesenchymal <- Vimentin positive
  6. Other       <- valid-and-negative for every gate
     Unresolved  <- blocked by an unavailable higher-priority gate (NOT a negative Other)

Gate marker calls come straight from RESTORE (``{m}_state``). Subtypes (endocrine/exocrine, immune/
mesenchymal fine types, NK) are DEFERRED to a validated level-2 pass (docs/level2_hierarchy.md); this
first deliverable is the 5 broad compartments only, so ``cell_type`` == ``compartment`` here.

Inputs : <restore_gated_dir>/donor_id=*/data_0.parquet   ({m}_pos, {m}_norm, {m}_state)
         <cells_dir>/donor_id=*/data_0.parquet           (cell_region, islet_num)
         <donor_metadata>                                 (disease.status; figure only)
Outputs: <phenotype_dir>/broad/donor_id=*/data_0.parquet  (object_id, donor_id, compartment, cell_type,
                                                            assignment_reason, blocked_by, {m}_state)
         <phenotype_dir>/broad_lineage_composition.png

    python -m phenocycler.lineage --jobs 8
"""

from __future__ import annotations

import argparse
import functools
import glob
import os
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
                     COMPARTMENT_COLORS, COMPARTMENT_ABBR, STATUS_ORDER, OTHER_LABEL, UNRESOLVED_LABEL)
from .parallel import map_donors

# 7 rows for composition (5 compartments + Other + Unresolved), in priority order.
CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL, UNRESOLVED_LABEL]

# tri-state {m}_state values (must match restore_apply.CELL_{POSITIVE,NEGATIVE,UNAVAILABLE})
CELL_POSITIVE = "positive"
CELL_UNAVAILABLE = "unavailable"

# assignment_reason vocabulary (kept beside the dtype widths so both stay in sync with the labels)
_REASONS = ("positive", "blocked_unavailable_gate", "all_negative")
# Fixed-width numpy string widths DERIVED from the actual labels, so a future longer compartment /
# reason / gate name can never silently truncate on assignment (numpy truncates without error).
_COMPARTMENT_W = max(len(c) for c in CLASSES)
_REASON_W = max(len(r) for r in _REASONS)
_BLOCKED_W = max(len(c) for c in COMPARTMENT_ORDER)


def status_map(cfg: PipelineConfig) -> dict:
    """donor_id -> disease.status from the donor metadata workbook (best-effort)."""
    try:
        m = pd.read_excel(cfg.donor_metadata)
        return dict(zip(m[cfg.metadata_donor_col].astype(str),
                        m[cfg.metadata_status_col].astype(str)))
    except Exception as e:  # metadata is optional for the assignment itself
        print(f"[warn] donor metadata unavailable ({e}); status will be '?'")
        return {}


def _required_gate_markers() -> list[str]:
    """The gate markers every donor's gated parquet must carry a ``{m}_state`` for (the accepted targets)."""
    return sorted({m for gate in COMPARTMENT_GATES.values() for m in gate})


def _state_getters(g: pd.DataFrame, required_markers):
    """Return (is_pos, is_unavail) closures over the tri-state ``{m}_state`` columns.

    Fail LOUD (Gate 4): a required gate marker with no ``{m}_state`` column is a real pipeline error
    (the gated parquet is incomplete), NOT an all-negative silent default. There is no ``.fillna``
    coercion here — every gate marker is one of the frozen ACCEPTED_PAIRS targets, so its column must
    exist for every eligible donor."""
    missing = [f"{m}_state" for m in required_markers if f"{m}_state" not in g.columns]
    if missing:
        raise KeyError(
            f"gated parquet missing required tri-state columns {missing}; "
            "restore_apply must emit a {marker}_state for every accepted target"
        )
    states = {m: g[f"{m}_state"].to_numpy() for m in required_markers}

    def is_pos(m: str) -> np.ndarray:
        return states[m] == CELL_POSITIVE

    def is_unavail(m: str) -> np.ndarray:
        return states[m] == CELL_UNAVAILABLE

    return is_pos, is_unavail


def _atomic_write_parquet(df: pd.DataFrame, dst: Path) -> None:
    """Write a parquet via a per-process temp file then ``os.replace`` (atomic on the same filesystem)."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f"{dst.name}.tmp.{os.getpid()}")
    try:
        df.to_parquet(tmp, index=False)
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()


def assign_donor(donor: str, gated_f: str, cells_dir) -> pd.DataFrame:
    """Type every cell in one donor by the tri-state ordered residual gating tree.

    Per gate the tri-state precedence is POSITIVE > UNAVAILABLE > NEGATIVE: a positive marker claims the
    compartment even if a sibling gate marker is unavailable; otherwise an unavailable gate marker makes
    the whole gate indeterminate. Walking ``COMPARTMENT_ORDER`` high->low, the first POSITIVE gate claims
    the cell; the first UNAVAILABLE gate (all higher gates having been valid-and-negative) BLOCKS it to
    ``Unresolved``; a cell that passes every gate valid-and-negative is ``Other``."""
    g = pd.read_parquet(gated_f)
    g["object_id"] = g["object_id"].astype(str)
    required = _required_gate_markers()
    is_pos, is_unavail = _state_getters(g, required)
    n = len(g)

    compartment = np.full(n, OTHER_LABEL, dtype=f"<U{_COMPARTMENT_W}")   # widths derived from the labels
    reason = np.full(n, "all_negative", dtype=f"<U{_REASON_W}")          # (never silently truncate)
    blocked_by = np.full(n, "", dtype=f"<U{_BLOCKED_W}")
    decided = np.zeros(n, bool)
    for comp in COMPARTMENT_ORDER:
        gate_pos = np.zeros(n, bool)
        gate_unavail = np.zeros(n, bool)
        for m in COMPARTMENT_GATES[comp]:
            gate_pos |= is_pos(m)
            gate_unavail |= is_unavail(m)
        gate_unavail &= ~gate_pos                                 # a positive marker dominates an unavailable sibling
        take_pos = (~decided) & gate_pos
        take_unavail = (~decided) & gate_unavail
        compartment[take_pos] = comp
        reason[take_pos] = "positive"
        compartment[take_unavail] = UNRESOLVED_LABEL
        reason[take_unavail] = "blocked_unavailable_gate"
        blocked_by[take_unavail] = comp
        decided |= take_pos | take_unavail
    # undecided cells are valid-and-negative for every gate -> Other/all_negative (as initialised).

    out = pd.DataFrame({
        "object_id": g["object_id"].astype(str),
        "donor_id": str(donor),
        "compartment": compartment,
        "cell_type": compartment,               # subtypes deferred to level 2 -> cell_type == compartment
        "assignment_reason": reason,
        "blocked_by": blocked_by,
    })
    for m in required:                           # per-marker tri-state passthrough (provenance)
        out[f"{m}_state"] = g[f"{m}_state"].to_numpy()
    # attach cell_region / islet_num for downstream (spatial / drill-down). Read EVERY partition file
    # (not just the first) and validate a one-to-one context key, so a duplicate object_id fails loud
    # instead of fanning the merge out into duplicate phenotype rows.
    cf = sorted(glob.glob(str(Path(cells_dir) / f"donor_id={donor}" / "*.parquet")))
    if cf:
        ctx = pd.concat(
            [pd.read_parquet(f, columns=["object_id", "cell_region", "islet_num"]) for f in cf],
            ignore_index=True,
        )
        out = out.merge(ctx.astype({"object_id": str}), on="object_id", how="left", validate="m:1")
    return out


def _assign_and_write(donor: str, gated_dir: str, cells_dir: str, out_dir: str) -> dict:
    """Worker: assign one donor, write its partition atomically, return the composition summary."""
    gf = sorted(glob.glob(str(Path(gated_dir) / f"donor_id={donor}" / "*.parquet")))
    if not gf:
        return {}
    out = assign_donor(donor, gf[0], cells_dir)
    dst = Path(out_dir) / f"donor_id={donor}"
    _atomic_write_parquet(out, dst / "data_0.parquet")
    vc = out["compartment"].value_counts()
    comp = (vc / len(out) * 100).reindex(CLASSES).fillna(0)
    n = len(out)
    other = int((out["compartment"] == OTHER_LABEL).sum())
    unresolved = int((out["compartment"] == UNRESOLVED_LABEL).sum())
    print(f"[{donor}] {n:,} cells | " + " ".join(f"{COMPARTMENT_ABBR[c]}={comp[c]:.0f}%" for c in CLASSES),
          flush=True)
    return {"donor": donor, "comp": comp.to_dict(), "n": n, "other": other, "unresolved": unresolved}


def run_lineage(cfg: PipelineConfig, *, donors=None, n_jobs=None) -> pd.DataFrame:
    """Assign compartments for all donors (parallel), write partitions + composition figure."""
    from .cohort import ensure_eligible_donors

    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    cfg.broad_dir.mkdir(parents=True, exist_ok=True)
    donor_ids = (
        list(ensure_eligible_donors(donors, context="lineage analysis"))
        if donors
        else cfg.discover_donors(cfg.restore_gated_dir)
    )
    if not donor_ids:
        raise SystemExit(f"[err] no gated donors under {cfg.restore_gated_dir}")

    fn = functools.partial(_assign_and_write, gated_dir=str(cfg.restore_gated_dir),
                           cells_dir=str(cfg.cells_dir), out_dir=str(cfg.broad_dir))
    results = [r for r in map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True) if r]
    if not results:
        raise SystemExit("[err] no donors assigned")

    smap = status_map(cfg)
    comp = {r["donor"]: pd.Series(r["comp"]).reindex(CLASSES).fillna(0) for r in results}
    n_total = sum(r["n"] for r in results)
    other_total = sum(r["other"] for r in results)
    unresolved_total = sum(r["unresolved"] for r in results)
    C = pd.DataFrame(comp).T
    C["status"] = [smap.get(d, "?") for d in C.index]
    print(f"\n[total] {n_total:,} cells | {100*other_total/max(n_total,1):.1f}% Other "
          "(valid-and-negative for every gate — real panel gaps, NOT force-assigned) | "
          f"{100*unresolved_total/max(n_total,1):.1f}% Unresolved (blocked by an unavailable gate)")

    _composition_figure(C, cfg.phenotype_dir / "broad_lineage_composition.png")
    print("\ncomposition by status (mean %):")
    print(C.groupby("status")[CLASSES].mean().reindex(STATUS_ORDER).round(1).to_string())
    return C


def _composition_figure(C: pd.DataFrame, path: Path):
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
    _ncol = 3                                       # rows grow with CLASSES (7: 5 compartments + Other + Unresolved)
    _nrow = int(np.ceil(len(CLASSES) / _ncol))
    gsB = outer[1].subgridspec(_nrow, _ncol, hspace=0.62, wspace=0.34)
    bax = {L: fig.add_subplot(gsB[i // _ncol, i % _ncol]) for i, L in enumerate(CLASSES)}

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
        if not present_status:          # no donor-status metadata -> the disease-comparison panel is moot
            ax.set_title(L, fontsize=14, color=COMPARTMENT_COLORS[L])
            ax.text(0.5, 0.5, "no disease-status\nmetadata", ha="center", va="center",
                    transform=ax.transAxes, fontsize=10, color="#888888")
            ax.set_xticks([]); ax.set_yticks([])
            continue
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
    fig.suptitle("Cell Phenotyping: 5 compartments + Other + Unresolved (tri-state ordered residual)",
                 fontsize=20, y=1.0)
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
