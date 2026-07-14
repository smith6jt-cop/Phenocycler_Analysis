#!/usr/bin/env python3
"""
Cell-type × marker QC figures (dotplot + heatmap) for the 5-compartment ordered-residual phenotyping.

Reads the broad ``cell_type`` labels + the RESTORE `_pos`/`_norm` (ONE gated dir now — the extra pass is
folded in) + the full REDSEA-corrected panel, and writes two PNGs under ``cfg.phenotype_dir``:

  • celltype_marker_dotplot.png  — cell_type × the RESTORE gate markers.
        dot SIZE  = % of the cell_type's cells `_pos` for the marker
        dot COLOR = mean RESTORE `_norm`, on a single shared ABSOLUTE scale (capped ~97th pct so bright
                    markers saturate while background markers stay dark at their true level).
  • celltype_marker_heatmap.png  — cell_type × full REDSEA-corrected TYPE-marker panel, z-scored per
        marker (|z| ≥ 1.2 annotated). Each cell_type's defining markers light up.

Inputs (per donor, joined by `object_id`):
  cfg.broad_dir          (data/phenotype/broad/donor_id=*)      — object_id, cell_type
  cfg.restore_gated_dir  (data/restore_gated_redsea/donor_id=*) — {m}_pos, {m}_norm (all markers)
  cfg.cells_redsea_dir   (data/cells_redsea/donor_id=*)         — full REDSEA-corrected panel (heatmap)

NaN-safe aggregation: several `cells_redsea` donors carry a few NaN edge cells (`cell_area_px == 0`); a
plain `np.sum` NaN-poisons the whole marker column, so every accumulation uses `np.nansum` + a finite
count. `restore_gated_redsea` has no NaN but is treated the same for safety.

Run:  python -m phenocycler.figures [--config config.ini] [--donors 6539 6380 ...]
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.dataset as ds

import matplotlib
if not matplotlib.get_backend().startswith("module://"):  # keep a notebook's inline/widget backend; force Agg only headless
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import PipelineConfig, load_config
from .marker_taxonomy import heatmap_markers

# RESTORE gate markers for the dotplot, grouped by the compartment / sub-type they define (readable).
RMARK = ["E_cadherin", "Pan_Cytokeratin", "EpCAM", "Keratin_5", "TP63",   # epithelial (+ basal/ductal)
         "INS", "GCG", "SST",                                             # endocrine hormones (IAPP removed — failed)
         "CD31", "CD34", "Podoplanin",                                    # endothelial (+ lymphatic)
         "B3TUBB",                                                        # neural
         "CD3e", "CD20", "CD79a", "CD68", "CD163", "CD11c", "MPO", "CD56", # immune (lymphoid->myeloid->NK)
         "SMA", "Vimentin"]                                              # mesenchymal

# Row order (cell_type) — encodes the intended diagonal reading of the identity figures.
CELLTYPE_ORDER = ["Endocrine-beta", "Endocrine-alpha", "Endocrine-delta", "Exocrine", "Exocrine-ductal",
                  "Endothelial-vascular", "Endothelial-lymphatic", "Neural",
                  "T", "B", "Plasma", "NK", "Macrophage", "DC", "Neutrophil", "Immune-other",
                  "Muscle", "Fibroblast", "Other"]

# Non-marker columns to exclude from the full heatmap panel.
META = {"object_id", "donor_id", "image", "cell_area", "cell_area_px", "cell_region", "islet_num"}


def run_figures(cfg: PipelineConfig, *, donors=None) -> dict:
    """Build the cell_type × marker dotplot + heatmap and write them under ``cfg.phenotype_dir``.

    Streams one donor at a time (memory released each iteration), accumulating EXACT NaN-safe per-cell_type
    sums. Returns the two output PNG paths. ``donors`` restricts the run (defaults to every donor under
    ``cfg.broad_dir``)."""
    donors = list(donors) if donors else cfg.discover_donors(cfg.broad_dir)
    if not donors:
        raise SystemExit(f"no broad output under {cfg.broad_dir} — run phenocycler.lineage first")

    # Full marker panel = channels COMMON to every donor (for the heatmap). Batch 1 carries CD38 and
    # batch 2 carries b_Catenin1 (a physical panel change), so take the intersection; order follows the
    # first donor's schema so the heatmap columns are stable.
    schemas = {d: set(ds.dataset(cfg.cells_redsea_dir / f"donor_id={d}", format="parquet").schema.names)
               for d in donors}
    common = set.intersection(*schemas.values())
    first = ds.dataset(cfg.cells_redsea_dir / f"donor_id={donors[0]}", format="parquet").schema.names
    avail = [c for c in first if c in common and c not in META and not c.startswith("__")]
    panel = heatmap_markers(avail)   # TYPE markers only — drops PROCESS/state markers + DAPI
    dropped_batch = sorted((set.union(*schemas.values()) - common) - META)
    print(f"[panel] heatmap = {len(panel)} TYPE markers (from {len(avail)} common); "
          f"dropped {len(avail) - len(panel)} PROCESS/excluded + {len(dropped_batch)} batch-specific "
          f"{dropped_batch}", flush=True)

    # One streaming pass per donor: accumulate EXACT per-cell_type sums (memory released each iteration)
    pos_sum = defaultdict(lambda: np.zeros(len(RMARK)))    # fraction positive  -> dot size
    norm_sum = defaultdict(lambda: np.zeros(len(RMARK)))   # mean RESTORE norm  -> dot color
    norm_cnt = defaultdict(lambda: np.zeros(len(RMARK)))   # finite _norm count (NaN-safe mean)
    expr_sum = defaultdict(lambda: np.zeros(len(panel)))   # mean REDSEA expr   -> heatmap
    expr_cnt = defaultdict(lambda: np.zeros(len(panel)))   # finite-value count per marker (NaN-safe mean)
    n = defaultdict(int)
    gcols = [f"{m}_pos" for m in RMARK] + [f"{m}_norm" for m in RMARK]
    for d in donors:
        bl = ds.dataset(cfg.broad_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id", "cell_type"]).to_pandas()
        g = ds.dataset(cfg.restore_gated_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id"] + gcols).to_pandas()
        e = ds.dataset(cfg.cells_redsea_dir / f"donor_id={d}", format="parquet").to_table(
                columns=["object_id"] + panel).to_pandas()
        df = bl.merge(g, on="object_id").merge(e, on="object_id")
        for ct, sub in df.groupby("cell_type"):
            pv = sub[[f"{m}_pos" for m in RMARK]].to_numpy(dtype=float)
            nv = sub[[f"{m}_norm" for m in RMARK]].to_numpy(dtype=float)
            ev = sub[panel].to_numpy(dtype=float)
            pos_sum[ct] += np.nansum(pv, axis=0)
            norm_sum[ct] += np.nansum(nv, axis=0); norm_cnt[ct] += np.isfinite(nv).sum(0)
            expr_sum[ct] += np.nansum(ev, axis=0); expr_cnt[ct] += np.isfinite(ev).sum(0)
            n[ct] += len(sub)
        print(f"[{d}] {len(df):,} cells aggregated", flush=True)

    rows = [c for c in CELLTYPE_ORDER if c in n] + [c for c in n if c not in CELLTYPE_ORDER]
    frac = np.vstack([pos_sum[c] / n[c] for c in rows])                          # R x M fraction positive
    cnorm = np.vstack([norm_sum[c] / np.maximum(norm_cnt[c], 1) for c in rows])  # R x M NaN-safe mean norm
    expr = np.vstack([expr_sum[c] / np.maximum(expr_cnt[c], 1) for c in rows])   # R x P NaN-safe mean
    cmax = float(np.percentile(cnorm, 97))
    print("cells per cell_type:", {c: int(n[c]) for c in rows})
    print(f"[dotplot] absolute colour scale: mean-norm range {cnorm.min():.2f}–{cnorm.max():.2f}, "
          f"saturating at {cmax:.1f}", flush=True)

    plt.rcParams.update({"font.size": 15, "xtick.labelsize": 14, "ytick.labelsize": 14,
                         "axes.titlesize": 17, "axes.labelsize": 16, "legend.fontsize": 14})

    # ---------- DOTPLOT: cell_type × RESTORE marker ----------
    fig, ax = plt.subplots(figsize=(0.95 * len(RMARK) + 3, 0.55 * len(rows) + 1.8))
    X, Y = np.meshgrid(np.arange(len(RMARK)), np.arange(len(rows)))
    dots = ax.scatter(X.ravel(), Y.ravel(), s=(frac * 420 + 8).ravel(),
                      c=cnorm.ravel(), cmap="viridis", vmin=0, vmax=cmax, edgecolor="k", linewidth=0.3)
    ax.set_xticks(range(len(RMARK))); ax.set_xticklabels(RMARK, rotation=45, ha="right", fontsize=14)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=14)
    ax.set_ylim(len(rows) - 0.5, -0.5); ax.set_xlim(-0.5, len(RMARK) - 0.5)
    ax.set_title("Cell type × marker — RESTORE positivity (size = % positive, color = absolute mean norm)", fontsize=16)
    cb = fig.colorbar(dots, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label(f"mean RESTORE norm (absolute; ≥{cmax:.1f} sat.)", fontsize=14)
    cb.ax.tick_params(labelsize=12)
    for f in (0.25, 0.5, 1.0):
        ax.scatter([], [], s=f * 420 + 8, c="lightgray", edgecolor="k", linewidth=0.3, label=f"{int(f*100)}%")
    ax.legend(title="% positive", bbox_to_anchor=(1.16, 1.0), loc="upper left", frameon=False,
              labelspacing=1.2, fontsize=13, title_fontsize=13)
    plt.tight_layout()
    fig.savefig(cfg.phenotype_dir / "celltype_marker_dotplot.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {cfg.phenotype_dir / 'celltype_marker_dotplot.png'}")

    # ---------- HEATMAP: cell_type × full marker panel (z-scored per marker) ----------
    z = (expr - expr.mean(0)) / (expr.std(0) + 1e-9)
    fig, ax = plt.subplots(figsize=(0.46 * len(panel) + 4, 0.6 * len(rows) + 3))
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
    ax.set_xticks(range(len(panel))); ax.set_xticklabels(panel, rotation=90, fontsize=14)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=14)
    for i in range(len(rows)):
        for j in range(len(panel)):
            if abs(z[i, j]) >= 1.2:
                ax.text(j, i, f"{z[i, j]:.1f}", ha="center", va="center", fontsize=10, color="black")
    ax.set_xticks(np.arange(-.5, len(panel), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.5); ax.tick_params(which="minor", length=0)
    ax.set_title("Cell type × marker — mean REDSEA-corrected expression (z-scored per marker)", fontsize=16)
    cb = fig.colorbar(im, ax=ax, fraction=0.012, pad=0.01); cb.set_label("z-score", fontsize=14)
    cb.ax.tick_params(labelsize=13)
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
