#!/usr/bin/env python3
"""
Step 4 — broad lineage assignment (every cell typed, ZERO "Unassigned").

Faithful port of ``scripts/senior/assign_broad_lineage.py`` (Islet-Explorer-
Senior).  On the REDSEA-corrected, re-RESTORE'd data each cell is typed by a
strict hierarchy, guaranteeing no "Unassigned" bucket:

  1. Endocrine    if INS/GCG/SST `_pos`
  2. Immune       elif CD3e/CD20/CD163 `_pos`
  3. Endothelial  elif CD31 `_pos`
  4. else STRUCTURAL background → argmax of (Pan_Cytokeratin, Vimentin, SMA) `_norm`
        → Epithelial / Fibroblast / Muscle, with cells below ALL three structural
        thresholds defaulted to Epithelial (flagged `epi_default`).

Rationale (unchanged from upstream): RESTORE sets each threshold at mean+3σ of the
marker's NEGATIVE population — right for sparse specific markers, but it flags only
the brightest ~20% of any ABUNDANT marker, so a plain 6-way `_norm` argmax scatters
the dim acinar majority into low-threshold markers.  Hence sparse lineages use the
calibrated `_pos`; the structural background uses `_norm` argmax; and marker-negative
cells default to Epithelial (validated via Ker8_18/Pan_CK medians).

Inputs : <restore_gated_dir>/donor_id=*/data_0.parquet   ({m}_norm, {m}_pos)
         <cells_dir>/donor_id=*/data_0.parquet            (cell_region, islet_num)
         <donor_metadata>                                 (disease.status)
Outputs: <phenotype_dir>/broad/donor_id=*/data_0.parquet
         <phenotype_dir>/broad_lineage_composition.png

Additions: config-driven paths + optional per-donor parallel assignment.

    python -m phenocycler.lineage --jobs 4
"""

from __future__ import annotations

import argparse
import functools
import glob
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import (PipelineConfig, load_config, LINEAGES, STRUCT_LINEAGES,
                     LINEAGE_COLORS, STATUS_ORDER)
from .parallel import map_donors

LNAMES = list(LINEAGES)
STRUCT = list(STRUCT_LINEAGES)


def status_map(cfg: PipelineConfig) -> dict:
    """donor_id -> disease.status from the donor metadata workbook (best-effort)."""
    try:
        m = pd.read_excel(cfg.donor_metadata)
        return dict(zip(m[cfg.metadata_donor_col].astype(str),
                        m[cfg.metadata_status_col].astype(str)))
    except Exception as e:  # metadata is optional for the assignment itself
        print(f"[warn] donor metadata unavailable ({e}); status will be '?'")
        return {}


def assign_donor(donor: str, gated_f: str, cells_dir: Path) -> pd.DataFrame:
    """Type every cell in one donor by the deterministic hierarchy (zero Unassigned)."""
    g = pd.read_parquet(gated_f)
    norm = {m: g[f"{m}_norm"].to_numpy(np.float64) for ms in LINEAGES.values() for m in ms}
    pos = {m: g[f"{m}_pos"].to_numpy(bool) for ms in LINEAGES.values() for m in ms}
    scores = {L: np.max([norm[m] for m in LINEAGES[L]], axis=0) for L in LNAMES}

    # structural background: argmax over Pan_Cytokeratin / Vimentin / SMA _norm
    struct = np.vstack([scores[L] for L in STRUCT]).T            # N x 3
    sorder = np.argsort(-struct, axis=1)
    # dtype "<U16" so later longer labels ("Endothelial", 11 chars) are not silently
    # truncated by the narrow fixed width of np.array(STRUCT) (longest = "Epithelial", 10).
    struct_call = np.array(STRUCT)[sorder[:, 0]].astype("<U16")
    ss = np.take_along_axis(struct, sorder, axis=1)
    struct_margin = ss[:, 0] - ss[:, 1]
    # cells below every structural threshold → Epithelial (validated background acinar).
    struct_pos = pos["Pan_Cytokeratin"] | pos["Vimentin"] | pos["SMA"]
    struct_call[~struct_pos] = "Epithelial"

    # hierarchy: specific _pos lineages take precedence over the structural background
    endocrine = pos["INS"] | pos["GCG"] | pos["SST"]
    immune = pos["CD3e"] | pos["CD20"] | pos["CD163"]
    endoth = pos["CD31"]
    broad = struct_call.copy()
    margin = struct_margin.astype(np.float64)
    broad[endoth] = "Endothelial"; margin[endoth] = scores["Endothelial"][endoth] - 1.0
    broad[immune] = "Immune";      margin[immune] = scores["Immune"][immune] - 1.0
    broad[endocrine] = "Endocrine"; margin[endocrine] = scores["Endocrine"][endocrine] - 1.0
    specific = endocrine | immune | endoth
    epi_default = (~specific) & (~struct_pos)   # epithelial by background (validated)

    out = pd.DataFrame({"object_id": g["object_id"].astype(str), "donor_id": donor,
                        "broad_lineage": broad, "assign_margin": margin.astype(np.float32),
                        "epi_default": epi_default})
    for L in LNAMES:
        out[f"score_{L}"] = scores[L].astype(np.float32)
    # attach cell_region / islet_num for downstream
    cf = sorted(glob.glob(str(Path(cells_dir) / f"donor_id={donor}" / "*.parquet")))
    if cf:
        ctx = pd.read_parquet(cf[0], columns=["object_id", "cell_region", "islet_num"])
        out = out.merge(ctx.astype({"object_id": str}), on="object_id", how="left")
    return out


def _assign_and_write(donor: str, gated_dir: str, cells_dir: str, out_dir: str) -> dict:
    """Worker: assign one donor, write its partition, return composition summary."""
    gf = sorted(glob.glob(str(Path(gated_dir) / f"donor_id={donor}" / "*.parquet")))
    if not gf:
        return {}
    out = assign_donor(donor, gf[0], Path(cells_dir))
    dst = Path(out_dir) / f"donor_id={donor}"
    dst.mkdir(parents=True, exist_ok=True)
    out.to_parquet(dst / "data_0.parquet", index=False)
    vc = out["broad_lineage"].value_counts()
    comp = (vc / len(out) * 100).reindex(LNAMES).fillna(0)
    n = len(out); weak = int(out["epi_default"].sum())
    print(f"[{donor}] {n:,} cells | " + " ".join(f"{L[:3]}={comp[L]:.0f}%" for L in LNAMES) +
          f" | epi_default={100*weak/n:.1f}% | Unassigned=0", flush=True)
    return {"donor": donor, "comp": comp.to_dict(), "n": n, "weak": weak}


def run_lineage(cfg: PipelineConfig, *, donors=None, n_jobs=None) -> pd.DataFrame:
    """Assign broad lineage for all donors (parallel), write partitions + composition figure."""
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    cfg.broad_dir.mkdir(parents=True, exist_ok=True)
    donor_ids = donors or cfg.discover_donors(cfg.restore_gated_dir)
    if not donor_ids:
        raise SystemExit(f"[err] no gated donors under {cfg.restore_gated_dir}")

    fn = functools.partial(_assign_and_write, gated_dir=str(cfg.restore_gated_dir),
                           cells_dir=str(cfg.cells_dir), out_dir=str(cfg.broad_dir))
    results = [r for r in map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True) if r]

    smap = status_map(cfg)
    comp = {r["donor"]: pd.Series(r["comp"]).reindex(LNAMES).fillna(0) for r in results}
    n_total = sum(r["n"] for r in results)
    weak_total = sum(r["weak"] for r in results)
    C = pd.DataFrame(comp).T
    C["status"] = [smap.get(d, "?") for d in C.index]
    print(f"\n[total] {n_total:,} cells, 0 Unassigned, "
          f"{100*weak_total/max(n_total,1):.1f}% epithelial-by-background")

    _composition_figure(C, cfg.phenotype_dir / "broad_lineage_composition.png")
    print("\ncomposition by status (mean %):")
    print(C.groupby("status")[LNAMES].mean().round(1).to_string())
    return C


def _composition_figure(C: pd.DataFrame, path: Path):
    fig, ax = plt.subplots(1, 3, figsize=(22, 6))
    Cs = C.sort_values("status")
    bottom = np.zeros(len(Cs))
    for L in LNAMES:
        ax[0].bar(range(len(Cs)), Cs[L], bottom=bottom, color=LINEAGE_COLORS[L], label=L)
        bottom += Cs[L].to_numpy()
    ax[0].set_xticks(range(len(Cs)))
    ax[0].set_xticklabels([f"{d}\n{Cs.loc[d,'status']}" for d in Cs.index], fontsize=7)
    ax[0].set_ylabel("% of cells"); ax[0].set_title("Broad-lineage composition per donor (0 Unassigned)")
    ax[0].legend(fontsize=8, ncol=2, loc="upper right")
    xs = np.arange(len(LNAMES)); bw = 0.25
    for si, s in enumerate(STATUS_ORDER):
        m = C[C.status == s][LNAMES].mean()
        ax[1].bar(xs + (si - 1) * bw, m.values, bw, label=s)
    ax[1].set_xticks(xs); ax[1].set_xticklabels(LNAMES, rotation=45, ha="right", fontsize=9)
    ax[1].set_ylabel("mean % of cells"); ax[1].set_title("Composition by disease status"); ax[1].legend()
    endo = [C[C.status == s]["Endocrine"].mean() for s in STATUS_ORDER]
    imm = [C[C.status == s]["Immune"].mean() for s in STATUS_ORDER]
    ax[2].plot(STATUS_ORDER, endo, "o-", color=LINEAGE_COLORS["Endocrine"], lw=2, ms=9, label="Endocrine (expect ↓)")
    ax[2].plot(STATUS_ORDER, imm, "s-", color=LINEAGE_COLORS["Immune"], lw=2, ms=9, label="Immune (expect ↑)")
    ax[2].set_ylabel("mean % of cells")
    ax[2].set_title("Disease trend: endocrine loss + immune gain"); ax[2].legend()
    fig.suptitle("Step 1 — broad lineage (argmax RESTORE _norm on REDSEA-corrected data)", fontsize=14)
    fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight")
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
