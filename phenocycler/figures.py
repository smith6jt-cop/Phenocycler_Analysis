#!/usr/bin/env python3
"""
Cell-type × marker QC figures for the broad lineages (dotplot + heatmap).

Faithful port of ``scripts/senior/plot_celltype_markers.py`` (Islet-Explorer-Senior). The
scientific logic is byte-for-byte the upstream "Step-F" identity QC generator; the only changes
are (a) hardcoded ``data/`` paths → a :class:`PipelineConfig`, and (b) an in-process
``run_figures(cfg)`` entry point + a ``main`` CLI so the pipeline/notebook can run it as a stage.

Reads the broad-lineage outputs and writes two PNGs under ``cfg.phenotype_dir`` (data/phenotype/):

  • celltype_marker_dotplot.png  — broad lineage × the RESTORE lineage markers.
        dot SIZE  = % of the lineage's cells `_pos` for the marker
        dot COLOR = mean RESTORE `_norm`, on a single shared ABSOLUTE scale (all markers directly
                    comparable), capped at a robust upper value (~97th pct) so a few very-bright
                    lineage markers saturate while background markers stay dark at their true level.
  • celltype_marker_heatmap.png  — broad lineage × full REDSEA-corrected TYPE-marker panel, z-scored
        per marker (|z| ≥ 1.2 annotated). Most markers do NOT separate the broad lineages (they are
        immune/endocrine SUB-type markers) so they sit near 0; the lineage-defining markers light up.

Inputs (per donor, joined by `object_id`):
  cfg.broad_dir              (data/phenotype/broad/donor_id=*)  — object_id, broad_lineage
  cfg.restore_gated_dir      (data/restore_gated_redsea/donor_id=*) — {m}_pos, {m}_norm (validated 10)
  cfg.restore_gated_extra_dir(data/restore_gated_redsea_extra/donor_id=*) — B3TUBB/CD99/MPO _pos/_norm
  cfg.cells_redsea_dir       (data/cells_redsea/donor_id=*)  — full REDSEA-corrected panel (heatmap)

NaN-safe aggregation: several `cells_redsea` donors carry a few NaN edge cells (`cell_area_px == 0`),
which default to Epithelial; a plain `np.sum` NaN-poisons the whole marker column (→ a blank
Epithelial row), so every accumulation uses `np.nansum` + a finite-value count. `restore_gated_redsea`
has no NaN, but its `_norm` is treated the same way for safety.

Run:  python -m phenocycler.figures [--config config.ini] [--donors 6539 6380 ...]
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.dataset as ds

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import (PipelineConfig, load_config, LINEAGES, LINEAGE_COLORS,  # noqa: F401
                     STATUS_ORDER, EXTRA_MARKERS)
from .marker_taxonomy import heatmap_markers, TYPE, PROCESS, EXCLUDED, CD99_BRIGHT  # noqa: F401
# TYPE/PROCESS/EXCLUDED are the taxonomy `heatmap_markers()` splits on (imported for API parity;
# the split itself is applied inside `heatmap_markers`). LINEAGE_COLORS/STATUS_ORDER are part of the
# figure API surface but this identity QC figure does not colour-by-lineage or split by disease status.

# 13 RESTORE lineage-gate markers, ordered by the lineage they define (readable dotplot). B3TUBB/CD99/MPO
# are gated separately (Pan_CK ref) into restore_gated_redsea_extra so the validated 10-marker gates stay
# untouched; read from there and merged by object_id below.
RMARK = ["Pan_Cytokeratin", "Vimentin", "SMA", "B3TUBB", "CD31",
         "INS", "GCG", "SST", "CD99", "CD3e", "CD20", "CD163", "MPO"]
# Heatmap row order (encodes the intended diagonal reading of the identity heatmap). These are the same
# 8 mutually-exclusive broad classes as config.LINEAGES; the literal order is preserved from upstream.
LIN_ORDER = ["Epithelial", "Fibroblast", "Muscle", "Neural", "Endothelial", "Endocrine", "Immune", "Neutrophil"]
assert set(LIN_ORDER) == set(LINEAGES), "figures LIN_ORDER drifted from config.LINEAGES (8 broad classes)"
# Non-marker columns to exclude from the full heatmap panel.
META = {"object_id", "donor_id", "image", "cell_area", "cell_area_px", "cell_region", "islet_num"}


def run_figures(cfg: PipelineConfig, *, donors=None) -> dict:
    """Build the cell-type × marker dotplot + heatmap and write them under ``cfg.phenotype_dir``.

    Streams one donor at a time (memory released each iteration), accumulating EXACT NaN-safe
    per-lineage sums. Returns the two output PNG paths. ``donors`` restricts the run to a subset
    (defaults to every donor discovered under ``cfg.broad_dir``).
    """
    donors = list(donors) if donors else cfg.discover_donors(cfg.broad_dir)
    if not donors:
        raise SystemExit(f"no Step-1 output found under {cfg.broad_dir} — run "
                         f"assign_broad_lineage (phenocycler.lineage) first")

    # Full marker panel = channels COMMON to every donor (for the heatmap). Batch 1 carries CD38 and
    # batch 2 carries b_Catenin1 — a physical panel change — so a single donor's schema is not uniform;
    # take the intersection (the common panel, matching build_singlecell_anndata.py). Order follows the
    # first donor's schema so the heatmap columns are stable.
    schemas = {d: set(ds.dataset(cfg.cells_redsea_dir / f"donor_id={d}", format="parquet").schema.names)
               for d in donors}
    common = set.intersection(*schemas.values())
    first = ds.dataset(cfg.cells_redsea_dir / f"donor_id={donors[0]}", format="parquet").schema.names
    avail = [c for c in first if c in common and c not in META and not c.startswith("__")]
    panel = heatmap_markers(avail)   # TYPE markers only — drops PROCESS/state markers + DAPI + IAPP
    dropped_batch = sorted((set.union(*schemas.values()) - common) - META)
    print(f"[panel] heatmap = {len(panel)} TYPE markers (from {len(avail)} common); "
          f"dropped {len(avail) - len(panel)} PROCESS/excluded + {len(dropped_batch)} batch-specific "
          f"{dropped_batch}", flush=True)

    # One streaming pass per donor: accumulate EXACT per-lineage sums (memory released each iteration)
    pos_sum = defaultdict(lambda: np.zeros(len(RMARK)))    # fraction positive  -> dot size
    norm_sum = defaultdict(lambda: np.zeros(len(RMARK)))   # mean RESTORE norm  -> dot color
    norm_cnt = defaultdict(lambda: np.zeros(len(RMARK)))   # finite _norm count (NaN-safe mean)
    expr_sum = defaultdict(lambda: np.zeros(len(panel)))   # mean REDSEA expr   -> heatmap
    expr_cnt = defaultdict(lambda: np.zeros(len(panel)))   # finite-value count per marker (NaN-safe mean)
    n = defaultdict(int)
    main_mark = [m for m in RMARK if m not in EXTRA_MARKERS]          # from restore_gated_redsea (validated 10)
    gcols_main = [f"{m}_pos" for m in main_mark] + [f"{m}_norm" for m in main_mark]
    gcols_extra = [f"{m}_pos" for m in EXTRA_MARKERS] + [f"{m}_norm" for m in EXTRA_MARKERS]  # from the extra dir
    for d in donors:
        bl = ds.dataset(cfg.broad_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id", "broad_lineage"]).to_pandas()
        g = ds.dataset(cfg.restore_gated_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id"] + gcols_main).to_pandas()
        gx = ds.dataset(cfg.restore_gated_extra_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id"] + gcols_extra).to_pandas()
        gx["CD99_pos"] = gx["CD99_norm"] >= CD99_BRIGHT   # bright-only CD99 (matches the Endocrine gate)
        e = ds.dataset(cfg.cells_redsea_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id"] + panel).to_pandas()
        df = bl.merge(g, on="object_id").merge(gx, on="object_id").merge(e, on="object_id")
        for lin, sub in df.groupby("broad_lineage"):
            pv = sub[[f"{m}_pos" for m in RMARK]].to_numpy(dtype=float)
            nv = sub[[f"{m}_norm" for m in RMARK]].to_numpy(dtype=float)
            ev = sub[panel].to_numpy(dtype=float)
            pos_sum[lin] += np.nansum(pv, axis=0)
            norm_sum[lin] += np.nansum(nv, axis=0); norm_cnt[lin] += np.isfinite(nv).sum(0)
            expr_sum[lin] += np.nansum(ev, axis=0); expr_cnt[lin] += np.isfinite(ev).sum(0)
            n[lin] += len(sub)
        print(f"[{d}] {len(df):,} cells aggregated", flush=True)

    lineages = [l for l in LIN_ORDER if l in n] + [l for l in n if l not in LIN_ORDER]
    frac = np.vstack([pos_sum[l] / n[l] for l in lineages])                          # L x 10 fraction positive
    cnorm = np.vstack([norm_sum[l] / np.maximum(norm_cnt[l], 1) for l in lineages])  # L x 10 NaN-safe mean norm
    expr = np.vstack([expr_sum[l] / np.maximum(expr_cnt[l], 1) for l in lineages])   # L x P  NaN-safe mean
    # Colour by ABSOLUTE mean RESTORE norm on a single shared scale (all markers directly comparable),
    # NOT the per-marker 0..1 rescale. A few lineage-defining markers are very bright (Muscle SMA ~20,
    # Endocrine INS ~6.6), so cap the colour scale at a robust upper value (~97th pct) — brighter cells
    # just saturate — while secondary/background markers render dark at their true low absolute level.
    cmax = float(np.percentile(cnorm, 97))
    print("cells per lineage:", {l: int(n[l]) for l in lineages})
    print(f"[dotplot] absolute colour scale: mean-norm range {cnorm.min():.2f}–{cnorm.max():.2f}, "
          f"saturating at {cmax:.1f}", flush=True)

    plt.rcParams.update({"font.size": 15, "xtick.labelsize": 14, "ytick.labelsize": 15,
                         "axes.titlesize": 17, "axes.labelsize": 16, "legend.fontsize": 14})

    # ---------- DOTPLOT: lineage x RESTORE marker ----------
    fig, ax = plt.subplots(figsize=(1.0 * len(RMARK) + 3, 0.62 * len(lineages) + 1.8))
    X, Y = np.meshgrid(np.arange(len(RMARK)), np.arange(len(lineages)))
    dots = ax.scatter(X.ravel(), Y.ravel(), s=(frac * 420 + 8).ravel(),
                      c=cnorm.ravel(), cmap="viridis", vmin=0, vmax=cmax, edgecolor="k", linewidth=0.3)
    ax.set_xticks(range(len(RMARK))); ax.set_xticklabels(RMARK, rotation=45, ha="right", fontsize=15)
    ax.set_yticks(range(len(lineages))); ax.set_yticklabels(lineages, fontsize=16)
    ax.set_ylim(len(lineages) - 0.5, -0.5); ax.set_xlim(-0.5, len(RMARK) - 0.5)
    ax.set_title("Cell type × marker — RESTORE positivity (size = % positive, color = absolute mean norm)", fontsize=17)
    cb = fig.colorbar(dots, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(f"mean RESTORE norm (absolute; ≥{cmax:.1f} sat.)", fontsize=15)
    cb.ax.tick_params(labelsize=13)
    for f in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=f * 420 + 8, c="lightgray", edgecolor="k", linewidth=0.3, label=f"{int(f*100)}%")
    ax.legend(title="% positive", bbox_to_anchor=(1.16, 1.0), loc="upper left", frameon=False,
              labelspacing=1.2, fontsize=14, title_fontsize=14)
    plt.tight_layout()
    fig.savefig(cfg.phenotype_dir / "celltype_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {cfg.phenotype_dir / 'celltype_marker_dotplot.png'}")

    # ---------- HEATMAP: lineage x full marker panel (z-scored per marker) ----------
    z = (expr - expr.mean(0)) / (expr.std(0) + 1e-9)
    fig, ax = plt.subplots(figsize=(0.46 * len(panel) + 4, 1.05 * len(lineages) + 3))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(panel))); ax.set_xticklabels(panel, rotation=90, fontsize=15)
    ax.set_yticks(range(len(lineages))); ax.set_yticklabels(lineages, fontsize=17)
    for i in range(len(lineages)):
        for j in range(len(panel)):
            if abs(z[i, j]) >= 1.2:
                ax.text(j, i, f"{z[i, j]:.1f}", ha="center", va="center", fontsize=11, color="black")
    ax.set_xticks(np.arange(-.5, len(panel), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(lineages), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.5); ax.tick_params(which="minor", length=0)
    ax.set_title("Cell type × marker — mean REDSEA-corrected expression (z-scored per marker)", fontsize=17)
    cb = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01); cb.set_label("z-score", fontsize=15)
    cb.ax.tick_params(labelsize=14)
    plt.tight_layout()
    fig.savefig(cfg.phenotype_dir / "celltype_marker_heatmap.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {cfg.phenotype_dir / 'celltype_marker_heatmap.png'}")

    return {"dotplot": cfg.phenotype_dir / "celltype_marker_dotplot.png",
            "heatmap": cfg.phenotype_dir / "celltype_marker_heatmap.png"}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--donors", nargs="*", default=None,
                    help="restrict to these donor ids (default: every donor under phenotype/broad)")
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    run_figures(cfg, donors=a.donors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
