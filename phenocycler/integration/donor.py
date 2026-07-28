"""
S8 — donor-level concordance. No registration, no matching, no assumptions.

This is the floor of the whole integration, and it is deliberately the least clever part.
Join on ``donor_id``, map both modalities into the common lineage vocabulary, and compare
composition. It works when registration fails, when the sections are too far apart to align,
when the morphology images are not mounted, and before any of the spatial machinery has run.
For a cohort where only some donors register cleanly, this is what keeps the other donors in
the analysis instead of silently dropping them.

The one thing it must get right is the **evidence-positive asymmetry**, and it is worth being
precise about why.

PhenoCycler's ``lineage.py`` assigns by a hierarchy of positive gates and then puts everything
that failed all of them into ``Epithelial``, flagged ``epi_default``. It is a *default sink*,
chosen because in pancreas the background parenchyma really is acinar — a defensible design,
and the flag is right there in the output. Xenium's broad tier has no such sink: cells with no
detected identity-core gene fall back to an absolute-score argmax and are flagged
``broad_lowsignal``, but they still land in whichever lineage scored highest rather than
piling into one class.

So a naive comparison of "fraction Epithelial vs fraction Exocrine" is measuring the
difference between a sink and a non-sink, not a difference between protein and RNA. Every
composition number here is therefore reported three ways — ``frac_all``,
``frac_evidence_positive`` and ``frac_fallback`` — and the headline concordance uses the
evidence-positive figure on both sides.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..config import STATUS_ORDER, PipelineConfig, load_config
from .contract import CellTable, discover_partitions, read_cell_table
from .manifest import load_manifest
from .vocab import COMMON_LINEAGES, RESOLVABLE

#: Lineages both modalities can resolve at the broad level without a caveat. Concordance is
#: reported for everything, but these are the ones a disagreement is actually informative
#: about — for the rest, a difference may just be one modality being unable to see the class.
COMPARABLE = tuple(
    lin for lin in COMMON_LINEAGES
    if lin != "Unassigned"
    and RESOLVABLE.get(lin, {}).get("phenocycler") in ("yes", "subclass")
    and RESOLVABLE.get(lin, {}).get("xenium") in ("yes", "surrogate")
)


def composition(table: CellTable,
                lineages: Sequence[str] = tuple(c for c in COMMON_LINEAGES
                                                if c != "Unassigned")) -> pd.DataFrame:
    """Per-lineage fractions for one section, reported three ways.

    ``frac_evidence_positive`` is the denominator-corrected figure: it counts only cells whose
    lineage came from positive evidence, on both the numerator and the denominator, so the two
    modalities' default/fallback behaviour cancels instead of biasing the comparison.
    """
    df = table.df
    total = len(df)
    ev = df["evidence_positive"].to_numpy(bool)
    total_ev = int(ev.sum())

    rows = []
    for lineage in lineages:
        sel = (df["lineage_common"] == lineage).to_numpy()
        n = int(sel.sum())
        n_ev = int((sel & ev).sum())
        rows.append({
            "modality": table.modality,
            "donor_id": table.donor_id,
            "roi": table.roi,
            "lineage_common": lineage,
            "n_cells": n,
            "n_evidence_positive": n_ev,
            "frac_all": n / total if total else np.nan,
            "frac_evidence_positive": n_ev / total_ev if total_ev else np.nan,
            "frac_fallback": (n - n_ev) / n if n else 0.0,
        })
    out = pd.DataFrame(rows)
    out.attrs["n_cells_total"] = total
    out.attrs["n_evidence_positive_total"] = total_ev
    return out


def marker_gene_concordance(
    pheno: CellTable,
    xen: CellTable,
    pairs: Sequence,
) -> pd.DataFrame:
    """Per-lineage mean protein vs mean gene-score, for the verified marker pairs.

    Cell-level correlation is impossible across serial sections — there are no shared cells.
    What *is* comparable is the per-lineage mean: if CD31 protein is high in the lineage
    PhenoCycler calls Vascular, PECAM1 expression should be high in the lineage Xenium calls
    Vascular. That is a real cross-modal check and needs no spatial correspondence at all.
    """
    rows = []
    p_feats = {c[len("feat__"):]: c for c in pheno.df.columns if c.startswith("feat__")}
    x_feats = {c[len("feat__"):]: c for c in xen.df.columns if c.startswith("feat__")}

    for lineage in COMPARABLE:
        p_sel = (pheno.df["lineage_common"] == lineage) & pheno.df["evidence_positive"]
        x_sel = (xen.df["lineage_common"] == lineage) & xen.df["evidence_positive"]
        if p_sel.sum() < 10 or x_sel.sum() < 10:
            continue

        score_col = x_feats.get(f"score_{_xenium_lineage_for(lineage)}")
        for mp in pairs:
            col = p_feats.get(mp.protein)
            if col is None:
                continue
            rows.append({
                "lineage_common": lineage,
                "protein": mp.protein,
                "genes": "|".join(mp.genes) if mp.genes else "",
                "on_panel": mp.on_panel,
                "pheno_mean_norm": float(pheno.df.loc[p_sel, col].mean()),
                "pheno_n": int(p_sel.sum()),
                "xen_lineage_score": (float(xen.df.loc[x_sel, score_col].mean())
                                      if score_col in xen.df.columns else np.nan),
                "xen_n": int(x_sel.sum()),
            })
    return pd.DataFrame(rows)


def _xenium_lineage_for(common: str) -> str:
    """Common lineage -> the Xenium broad name whose identity-core score represents it."""
    from .vocab import XENIUM_TO_COMMON

    for xen_name, com in XENIUM_TO_COMMON.items():
        if com == common:
            return xen_name
    return common


def compare_donor(pheno: CellTable, xen: CellTable) -> tuple[pd.DataFrame, dict]:
    """Composition comparison for one donor. Returns ``(long table, summary stats)``."""
    p = composition(pheno)
    x = composition(xen)

    merged = p.merge(x, on="lineage_common", suffixes=("_pheno", "_xen"))
    merged["delta_frac_all"] = merged["frac_all_pheno"] - merged["frac_all_xen"]
    merged["delta_frac_ev"] = (merged["frac_evidence_positive_pheno"]
                               - merged["frac_evidence_positive_xen"])
    merged["comparable"] = merged["lineage_common"].isin(COMPARABLE)

    from scipy import stats as sps

    sub = merged[merged["comparable"]]
    stats = {
        "donor_id": pheno.donor_id,
        "roi": pheno.roi,
        "n_cells_pheno": len(pheno.df),
        "n_cells_xen": len(xen.df),
        "frac_fallback_pheno": float(1 - pheno.df["evidence_positive"].mean()),
        "frac_fallback_xen": float(1 - xen.df["evidence_positive"].mean()),
        "n_comparable_lineages": int(len(sub)),
    }
    if len(sub) >= 3:
        for label, a, b in (
            ("all", "frac_all_pheno", "frac_all_xen"),
            ("ev", "frac_evidence_positive_pheno", "frac_evidence_positive_xen"),
        ):
            ok = sub[a].notna() & sub[b].notna()
            if ok.sum() >= 3:
                rho, pv = sps.spearmanr(sub.loc[ok, a], sub.loc[ok, b])
                stats[f"spearman_composition_{label}"] = float(rho)
                stats[f"p_composition_{label}"] = float(pv)
                stats[f"l1_composition_{label}"] = float(np.abs(sub.loc[ok, a]
                                                                - sub.loc[ok, b]).sum())
    return merged, stats


def disease_trend(compositions: pd.DataFrame, status_by_donor: dict) -> pd.DataFrame:
    """Per-lineage composition by disease status, per modality.

    The cross-modal acceptance test at cohort level: the endocrine fraction should fall
    ND -> AAB -> T1D and the immune fraction rise, *in both modalities independently*. Two
    measurement systems reproducing the same biological gradient is much stronger evidence
    that the harmonisation is sound than any single-donor agreement.
    """
    df = compositions.copy()
    df["status"] = df["donor_id"].map(status_by_donor).fillna("?")
    grp = (df.groupby(["modality", "lineage_common", "status"], observed=True)
             .agg(mean_frac=("frac_evidence_positive", "mean"),
                  sd_frac=("frac_evidence_positive", "std"),
                  n_donors=("donor_id", "nunique"))
             .reset_index())

    order = {s: i for i, s in enumerate(STATUS_ORDER)}
    grp["status_order"] = grp["status"].map(order).fillna(len(order))
    return grp.sort_values(["modality", "lineage_common", "status_order"]).reset_index(drop=True)


def run_donor(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
) -> pd.DataFrame:
    """Stage entry point. Runs on whatever is exported — pairing is a bonus, not a requirement."""
    man = load_manifest(cfg)
    status_by_donor = dict(zip(man["donor_id"], man["disease_status"]))

    pheno_parts = {d: r for d, r in discover_partitions(cfg.cells_pheno_dir)}
    xen_parts = {d: r for d, r in discover_partitions(cfg.cells_xen_dir)}
    want = set(donors) if donors else (set(pheno_parts) | set(xen_parts))

    all_comp, all_cmp, stats_rows, all_markers = [], [], [], []
    try:
        from .vocab import load_panel_taxonomy, protein_gene_pairs
        tax = load_panel_taxonomy(cfg)
        panel = sorted({g for gs in tax.core_genes.values() for g in gs})
        pairs = protein_gene_pairs(panel, tax=tax)
    except Exception:  # noqa: BLE001 - marker comparison is optional
        pairs = []

    for donor in sorted(want):
        p = x = None
        if donor in pheno_parts:
            try:
                p = read_cell_table(cfg.cells_pheno_dir, donor, pheno_parts[donor],
                                    modality="phenocycler")
                all_comp.append(composition(p))
            except FileNotFoundError:
                p = None
        if donor in xen_parts:
            try:
                x = read_cell_table(cfg.cells_xen_dir, donor, xen_parts[donor],
                                    modality="xenium")
                all_comp.append(composition(x))
            except FileNotFoundError:
                x = None

        if p is not None and x is not None:
            merged, stats = compare_donor(p, x)
            merged.insert(0, "donor_id", donor)
            merged["status"] = status_by_donor.get(donor, "")
            all_cmp.append(merged)
            stats["status"] = status_by_donor.get(donor, "")
            stats_rows.append(stats)
            print(f"[donor] {donor}: composition Spearman "
                  f"{stats.get('spearman_composition_ev', float('nan')):.3f} "
                  f"(evidence-positive), fallback {stats['frac_fallback_pheno']:.1%} pheno / "
                  f"{stats['frac_fallback_xen']:.1%} xenium", flush=True)
            if pairs:
                mk = marker_gene_concordance(p, x, pairs)
                if len(mk):
                    mk.insert(0, "donor_id", donor)
                    all_markers.append(mk)
        elif p is not None or x is not None:
            which = "PhenoCycler" if p is not None else "Xenium"
            print(f"[donor] {donor}: {which} only — composition recorded, no comparison",
                  flush=True)

    cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
    if all_comp:
        comp = pd.concat(all_comp, ignore_index=True)
        comp.to_csv(cfg.integration_qc_dir / "composition.csv", index=False)
        trend = disease_trend(comp, status_by_donor)
        trend.to_csv(cfg.integration_qc_dir / "disease_trend.csv", index=False)
        print(f"\n[donor] wrote composition.csv ({len(comp)} rows) and disease_trend.csv",
              flush=True)
    if all_cmp:
        pd.concat(all_cmp, ignore_index=True).to_csv(
            cfg.integration_qc_dir / "composition_comparison.csv", index=False)
    if all_markers:
        pd.concat(all_markers, ignore_index=True).to_csv(
            cfg.integration_qc_dir / "marker_gene_concordance.csv", index=False)

    stats_df = pd.DataFrame(stats_rows)
    if len(stats_df):
        stats_df.to_csv(cfg.integration_qc_dir / "donor_qc.csv", index=False)
    return stats_df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Donor-level composition concordance (no registration needed)")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_donor(cfg, args.donors, roi=args.roi)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
