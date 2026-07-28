"""
S10 (figures) — the plots you actually need to trust a registration.

Numbers in ``registration_qc.csv`` say whether the gates passed. An overlay says *how* it
failed when it did, which is the difference between "re-run with a different seed" and "these
sections are not from the same block". The registration overlay is the one figure worth
looking at for every donor.

Matplotlib with the Agg backend, matching the core pipeline's ``figures.py``, so everything
renders headless on a cluster node.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..config import LINEAGE_COLORS, STATUS_ORDER, PipelineConfig, load_config  # noqa: E402
from .contract import read_cell_table  # noqa: E402
from .manifest import load_manifest, paired_rows  # noqa: E402
from .rasterize import grid_for, tissue_mask  # noqa: E402
from .register import registration_dir  # noqa: E402
from .transform import Transform  # noqa: E402
from .vocab import COMMON_LINEAGES  # noqa: E402

#: Colours for the harmonised vocabulary. Reuses the core pipeline's lineage palette where
#: the names line up so a figure from either half of the project reads the same way.
COMMON_COLORS = {
    "Endocrine": LINEAGE_COLORS.get("Endocrine", "#CCBB44"),
    "Exocrine": LINEAGE_COLORS.get("Epithelial", "#4477AA"),
    "Stromal": LINEAGE_COLORS.get("Fibroblast", "#EE6677"),
    "Vascular": LINEAGE_COLORS.get("Endothelial", "#66CCEE"),
    "T_NK": "#228833",
    "B_Plasma": "#66AA55",
    "Myeloid": "#117733",
    "Granulocyte": LINEAGE_COLORS.get("Neutrophil", "#E69F00"),
    "Neural": LINEAGE_COLORS.get("Neural", "#B5838D"),
    "Unassigned": "#BBBBBB",
}


def registration_overlay(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    out_path: Optional[Path] = None,
) -> Optional[Path]:
    """Before/after overlay of the two tissue masks, plus the islet correspondence.

    Three panels: the two sections as they arrive, the same pair after registration, and the
    islet centroids with matched pairs joined by a line. A failed registration is obvious in
    panel 2 in a way it is not in a table of Dice values.
    """
    try:
        pheno = read_cell_table(cfg.cells_pheno_dir, donor, roi, modality="phenocycler")
        xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")
    except FileNotFoundError:
        return None
    reg = registration_dir(cfg, donor, roi)
    if not (reg / "transform.json").exists():
        return None
    tf = Transform.load(reg)

    grid = grid_for(pheno.xy(), 10.0, other=xen.xy(), margin_um=200)
    p_mask = tissue_mask(pheno.xy(), grid)
    x_mask = tissue_mask(xen.xy(), grid)
    x_warp = tf.warp_mask(x_mask, grid)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))

    def _rgb(a, b):
        img = np.zeros((*grid.shape, 3), dtype=float)
        img[..., 0] = a * 0.9          # PhenoCycler -> magenta
        img[..., 2] = a * 0.9
        img[..., 1] = b * 0.9          # Xenium -> green
        return img

    axes[0].imshow(_rgb(p_mask, x_mask), origin="lower")
    axes[0].set_title("before registration")
    axes[1].imshow(_rgb(p_mask, x_warp), origin="lower")
    from .register import dice

    axes[1].set_title(f"after — Dice {dice(x_warp, p_mask):.3f}")
    for ax in axes[:2]:
        ax.set_xticks([])
        ax.set_yticks([])

    ax = axes[2]
    from .structures import structures_path

    p_path = structures_path(cfg, "phenocycler", donor, roi, "islet")
    x_path = structures_path(cfg, "xenium", donor, roi, "islet")
    if p_path.exists() and x_path.exists():
        pi = pd.read_parquet(p_path)
        xi = pd.read_parquet(x_path)
        if len(xi):
            xy = tf.apply(xi.loc[:, ["x_um", "y_um"]].to_numpy(np.float64))
            ax.scatter(xy[:, 0], xy[:, 1], s=np.sqrt(xi["n_cells"]) * 4, c="tab:green",
                       alpha=0.6, label=f"Xenium islets (n={len(xi)})")
        if len(pi):
            ax.scatter(pi["x_um"], pi["y_um"], s=np.sqrt(pi["n_cells"]) * 4, c="magenta",
                       alpha=0.6, label=f"PhenoCycler islets (n={len(pi)})")

        pairs_path = cfg.paired_dir / "islets.parquet"
        if pairs_path.exists():
            pr = pd.read_parquet(pairs_path)
            pr = pr[pr["donor_id"].astype(str) == str(donor)]
            pmap = pi.set_index("structure_id") if len(pi) else None
            xmap = xi.set_index("structure_id") if len(xi) else None
            n_drawn = 0
            for _, r in pr.iterrows():
                if pmap is None or xmap is None:
                    break
                if r["fixed_id"] not in pmap.index or r["moving_id"] not in xmap.index:
                    continue
                f = pmap.loc[r["fixed_id"]]
                m = xmap.loc[r["moving_id"]]
                mx = tf.apply(np.array([[m["x_um"], m["y_um"]]]))[0]
                ax.plot([f["x_um"], mx[0]], [f["y_um"], mx[1]], "k-", lw=0.8, alpha=0.6)
                n_drawn += 1
            ax.set_title(f"islet correspondence ({n_drawn} matched)")
        else:
            ax.set_title("islets (registered frame)")
        ax.legend(fontsize=8, loc="best")
    ax.set_aspect("equal")
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")

    fig.suptitle(f"donor {donor} / {roi} — {cfg.integration_mode} mode "
                 f"(magenta = PhenoCycler, green = Xenium)")
    fig.tight_layout()

    out_path = out_path or (cfg.integration_figures_dir / f"registration_{donor}_{roi}.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def composition_figure(cfg: PipelineConfig, out_path: Optional[Path] = None) -> Optional[Path]:
    """Per-donor composition, both modalities side by side, evidence-positive fractions.

    Uses the evidence-positive fraction rather than the raw one, because PhenoCycler's
    Epithelial class is a default sink and Xenium's is not — plotting the raw fractions makes
    the two look more different than they are for a reason that is methodological, not
    biological.
    """
    path = cfg.integration_qc_dir / "composition.csv"
    if not path.exists():
        return None
    comp = pd.read_csv(path)
    if comp.empty:
        return None

    lineages = [c for c in COMMON_LINEAGES if c != "Unassigned"]
    donors = sorted(comp["donor_id"].astype(str).unique())
    fig, axes = plt.subplots(1, 2, figsize=(max(9, 1.1 * len(donors) + 4), 5), sharey=True)

    for ax, modality in zip(axes, ("phenocycler", "xenium")):
        sub = comp[comp["modality"] == modality]
        bottom = np.zeros(len(donors))
        for lineage in lineages:
            vals = np.array([
                sub.loc[(sub["donor_id"].astype(str) == d)
                        & (sub["lineage_common"] == lineage),
                        "frac_evidence_positive"].sum()
                for d in donors])
            ax.bar(donors, vals, bottom=bottom, label=lineage,
                   color=COMMON_COLORS.get(lineage, "#999999"), edgecolor="white", lw=0.3)
            bottom += vals
        ax.set_title(modality)
        ax.set_xlabel("donor")
        ax.tick_params(axis="x", rotation=90)
    axes[0].set_ylabel("fraction of evidence-positive cells")
    axes[-1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.suptitle("Harmonised lineage composition (evidence-positive cells only)")
    fig.tight_layout()

    out_path = out_path or (cfg.integration_figures_dir / "composition.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def islet_concordance_figure(cfg: PipelineConfig,
                             out_path: Optional[Path] = None) -> Optional[Path]:
    """Matched-islet agreement: the scientific acceptance test for the spatial pipeline.

    If registration and matching are real, the same islet measured two ways should agree on
    size, on beta fraction and on insulitis. A cloud with no structure here means the
    machinery ran but did not find the same islets.
    """
    path = cfg.paired_dir / "islets.parquet"
    if not path.exists():
        return None
    pairs = pd.read_parquet(path)
    if pairs.empty:
        return None

    mv = pairs["moving_modality"].iloc[0] if "moving_modality" in pairs.columns else "xenium"
    fx = pairs["fixed_modality"].iloc[0] if "fixed_modality" in pairs.columns else "phenocycler"

    specs = [("n_cells", "endocrine cells per islet"),
             ("frac_beta", "beta fraction"),
             ("immune_per_100_endo", "immune per 100 endocrine")]
    fig, axes = plt.subplots(1, len(specs), figsize=(5 * len(specs), 4.6))
    if len(specs) == 1:
        axes = [axes]

    from scipy import stats as sps

    for ax, (feat, label) in zip(axes, specs):
        a, b = f"{mv}_{feat}", f"{fx}_{feat}"
        if a not in pairs.columns or b not in pairs.columns:
            ax.set_visible(False)
            continue
        x = pd.to_numeric(pairs[a], errors="coerce")
        y = pd.to_numeric(pairs[b], errors="coerce")
        ok = x.notna() & y.notna()
        ax.scatter(x[ok], y[ok], s=22, alpha=0.7, c="tab:blue", edgecolor="none")
        if ok.sum() >= 3:
            rho, p = sps.spearmanr(x[ok], y[ok])
            lim = [min(x[ok].min(), y[ok].min()), max(x[ok].max(), y[ok].max())]
            ax.plot(lim, lim, "k--", lw=0.8, alpha=0.5)
            ax.set_title(f"{label}\nSpearman {rho:.3f} (p={p:.1e}, n={int(ok.sum())})")
        ax.set_xlabel(f"{mv}")
        ax.set_ylabel(f"{fx}")

    fig.suptitle("Matched islets — cross-modal concordance")
    fig.tight_layout()

    out_path = out_path or (cfg.integration_figures_dir / "islet_concordance.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def disease_trend_figure(cfg: PipelineConfig,
                         out_path: Optional[Path] = None) -> Optional[Path]:
    """Lineage fraction by ND/AAB/T1D, per modality.

    The cohort-level acceptance test. Two independent measurement systems reproducing the
    same ND -> AAB -> T1D gradient is much stronger evidence that the harmonisation is sound
    than agreement on any single donor.
    """
    path = cfg.integration_qc_dir / "disease_trend.csv"
    if not path.exists():
        return None
    trend = pd.read_csv(path)
    trend = trend[trend["status"].isin(STATUS_ORDER)]
    if trend.empty:
        return None

    lineages = [lin for lin in COMMON_LINEAGES
                if lin != "Unassigned" and lin in set(trend["lineage_common"])]
    ncol = min(5, max(1, len(lineages)))
    nrow = int(np.ceil(len(lineages) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.0 * nrow), squeeze=False)

    for i, lineage in enumerate(lineages):
        ax = axes[i // ncol][i % ncol]
        sub = trend[trend["lineage_common"] == lineage]
        for modality, marker in (("phenocycler", "o"), ("xenium", "s")):
            m = sub[sub["modality"] == modality].set_index("status").reindex(STATUS_ORDER)
            ax.errorbar(range(len(STATUS_ORDER)), m["mean_frac"], yerr=m["sd_frac"],
                        marker=marker, capsize=3, label=modality, lw=1.4)
        ax.set_xticks(range(len(STATUS_ORDER)))
        ax.set_xticklabels(STATUS_ORDER)
        ax.set_title(lineage, fontsize=10)
        ax.set_ylim(bottom=0)
    for j in range(len(lineages), nrow * ncol):
        axes[j // ncol][j % ncol].set_visible(False)
    axes[0][0].set_ylabel("fraction (evidence-positive)")
    axes[0][-1].legend(fontsize=8)
    fig.suptitle("Disease-status trend, measured independently by each modality")
    fig.tight_layout()

    out_path = out_path or (cfg.integration_figures_dir / "disease_trend.png")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def run_figures(cfg: PipelineConfig, donors: Optional[list[str]] = None,
                *, roi: Optional[str] = None, tissue: Optional[str] = None) -> list[Path]:
    """Stage entry point: render everything that has data behind it."""
    cfg.integration_figures_dir.mkdir(parents=True, exist_ok=True)
    made: list[Path] = []

    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    for _, r in man.iterrows():
        p = registration_overlay(cfg, r["donor_id"], r["roi"])
        if p:
            made.append(p)
            print(f"[figures] {p}", flush=True)

    for fn in (composition_figure, islet_concordance_figure, disease_trend_figure):
        p = fn(cfg)
        if p:
            made.append(p)
            print(f"[figures] {p}", flush=True)

    if not made:
        print("[figures] nothing to draw yet — run the earlier stages first", flush=True)
    return made


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Integration figures")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None,
                    help="one ROI (default: every ROI of --tissue, or all tissues)")
    ap.add_argument("--tissue", default=None, help="restrict to one tissue's ROIs")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    run_figures(cfg, args.donors, roi=args.roi, tissue=args.tissue)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
