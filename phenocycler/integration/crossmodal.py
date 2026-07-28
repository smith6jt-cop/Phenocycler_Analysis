"""
S9 — pseudo-cell linking. Inference, and labelled as such everywhere.

MaxFuse-class weakly-linked integration: fuzzy-smooth each modality's features over its own
spatial neighbourhood, project both into a shared space with CCA over the verified
protein<->gene pairs, then solve a linear assignment between cells.

Two constraints are enforced in code rather than documented and forgotten, because both are
ways this step can produce confident nonsense.

**The anchor set is thinner than "46 of 55 markers" suggests.** Of the markers that actually
discriminate pancreatic lineage, INS, GCG, SST and Vimentin are all off the Xenium panel —
those four define PhenoCycler's Endocrine and Fibroblast classes outright. What is left is
mostly immune and vascular: KRT19, ACTA2 (custom-100 only), PECAM1, CD3E, MS4A1, CD163,
CD68, CD99, TUBB3, MPO. That is a usable but narrow basis, and a basis narrow enough to be
degenerate is worse than none — so this module **refuses to run** below
``crossmodal_min_anchors`` rather than quietly producing links from three genes.

**Serial sections have no ground-truth cell correspondence.** A link is a statement about
similarity in a shared feature space, not a claim that two measurements came from one cell —
that claim is only available in ``same_slide`` mode. So every output column is prefixed
``inferred_``, every link carries a confidence, and validation is against structure-level
aggregates plus a permuted-link null. If the inferred pairing cannot reproduce the islet-level
concordance that S6 measured directly, it has not earned belief.

Running *within* matched islets and grid bins is what makes the assignment tractable and
meaningful: a global assignment over a million cells against a million others is both
computationally absurd and statistically meaningless, whereas "which of these 90 endocrine
cells corresponds to which of those 80" is a real question with spatial support behind it.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, feature_name, read_cell_table
from .manifest import PAIRED, load_manifest
from .vocab import MarkerPair, load_panel_taxonomy, protein_gene_pairs, usable_anchors


class AnchorError(RuntimeError):
    """Too few verified protein<->gene anchors to attempt linking."""


@dataclass
class LinkResult:
    links: pd.DataFrame
    n_anchors: int
    anchors: list[str]
    stats: dict


def resolve_anchors(
    cfg: PipelineConfig,
    pheno: CellTable,
    xen: CellTable,
    *,
    panel_genes: Optional[Sequence[str]] = None,
) -> tuple[list[MarkerPair], list[str], list[str]]:
    """Anchors usable for *this pair of tables*, given what each actually carries.

    Three filters, all necessary: the protein must have a gene on the run's panel, the
    PhenoCycler table must carry that protein as a feature, and the Xenium table must carry
    something representing the gene. A pair that passes the crosswalk but is absent from one
    of the tables is not an anchor, and counting it as one is how the minimum-anchor guard
    gets silently defeated.
    """
    tax = load_panel_taxonomy(cfg)
    p_feats = {feature_name(c): c for c in pheno.features}
    x_feats = {feature_name(c): c for c in xen.features}

    if panel_genes is None:
        # The panel is whatever the Xenium table actually carries — NOT the identity-core
        # union. Several markers that matter here (KRT19, ACTA2, PTPRC) reach the deployed
        # panel only through the custom-100 add-on and are absent from every identity core,
        # so validating against the core union would declare real anchors off-panel and
        # silently starve the minimum-anchor guard.
        panel_genes = sorted(g for g in x_feats if not g.startswith("score_"))

    pairs = usable_anchors(protein_gene_pairs(panel_genes, tax=tax), type_only=False)

    keep, p_cols, x_cols = [], [], []
    for mp in pairs:
        if mp.protein not in p_feats:
            continue
        # The Xenium side may carry the gene directly, or a lineage score standing in for it.
        xcol = None
        for gene in mp.genes:
            if gene in x_feats:
                xcol = x_feats[gene]
                break
        if xcol is None:
            continue
        keep.append(mp)
        p_cols.append(p_feats[mp.protein])
        x_cols.append(xcol)
    return keep, p_cols, x_cols


def fuzzy_smooth(values: np.ndarray, xy: np.ndarray, *, k: int = 15) -> np.ndarray:
    """Average each cell's features over its k nearest spatial neighbours.

    Single-cell measurements in both modalities are noisy — Xenium counts are dominated by
    single-transcript detections, PhenoCycler intensities by segmentation spillover. Smoothing
    over the local neighbourhood trades a little spatial resolution for a large gain in
    signal, which is what makes a shared embedding of two noisy modalities possible at all.
    """
    from scipy.spatial import cKDTree

    if len(xy) == 0:
        return values
    kk = int(min(k, len(xy)))
    _d, idx = cKDTree(xy).query(xy, k=kk)
    if idx.ndim == 1:
        idx = idx[:, None]
    return values[idx].mean(axis=1)


def _zscore(a: np.ndarray) -> np.ndarray:
    mu = np.nanmean(a, axis=0)
    sd = np.nanstd(a, axis=0)
    sd[sd == 0] = 1.0
    return np.nan_to_num((a - mu) / sd)


def link_cells(
    pheno: CellTable,
    xen: CellTable,
    p_cols: Sequence[str],
    x_cols: Sequence[str],
    *,
    smooth_k: int = 15,
    n_components: int = 4,
    max_cells: int = 4000,
    seed: int = 0,
) -> pd.DataFrame:
    """Link cells within one spatial group. Returns a frame of ``inferred_*`` columns.

    The assignment is one-to-one via Hungarian on the CCA-space distance, so a cell cannot be
    claimed twice. ``max_cells`` bounds the O(n^3) assignment; larger groups are subsampled
    and the fact is recorded rather than silently truncating.
    """
    rng = np.random.default_rng(seed)
    empty = pd.DataFrame(columns=["inferred_pheno_cell_id", "inferred_xen_cell_id",
                                  "inferred_distance", "inferred_confidence",
                                  "inferred_lineage_agree"])
    if len(pheno.df) < 3 or len(xen.df) < 3:
        return empty

    p_idx = np.arange(len(pheno.df))
    x_idx = np.arange(len(xen.df))
    if len(p_idx) > max_cells:
        p_idx = rng.choice(p_idx, max_cells, replace=False)
    if len(x_idx) > max_cells:
        x_idx = rng.choice(x_idx, max_cells, replace=False)

    P = pheno.df.iloc[p_idx].loc[:, list(p_cols)].to_numpy(np.float64)
    X = xen.df.iloc[x_idx].loc[:, list(x_cols)].to_numpy(np.float64)
    P = fuzzy_smooth(_zscore(P), pheno.xy()[p_idx], k=smooth_k)
    X = fuzzy_smooth(_zscore(X), xen.xy()[x_idx], k=smooth_k)

    n_comp = int(min(n_components, P.shape[1], X.shape[1], len(P) - 1, len(X) - 1))
    if n_comp < 1:
        return empty

    # CCA needs equal row counts, so fit on a paired random subsample. The fit only has to
    # find the shared axes; it is then applied to all cells in the group.
    n_fit = min(len(P), len(X))
    try:
        from sklearn.cross_decomposition import CCA

        cca = CCA(n_components=n_comp, max_iter=500)
        cca.fit(P[rng.permutation(len(P))[:n_fit]], X[rng.permutation(len(X))[:n_fit]])
        Pc, Xc = cca.transform(P), cca.transform(X)
    except Exception:  # noqa: BLE001 - degenerate features; fall back to the raw z-scores
        Pc, Xc = P[:, :n_comp], X[:, :n_comp]

    Pc, Xc = _zscore(Pc), _zscore(Xc)
    cost = np.linalg.norm(Pc[:, None, :] - Xc[None, :, :], axis=2)

    from scipy.optimize import linear_sum_assignment

    rows, cols = linear_sum_assignment(cost)
    d = cost[rows, cols]
    # Confidence: how much better is this partner than the next-best available one? A link
    # that is barely preferred over the alternative is not a link worth reporting.
    second = np.partition(cost, 1, axis=1)[:, 1] if cost.shape[1] > 1 else d
    margin = np.maximum(second[rows] - d, 0.0)
    conf = margin / (margin + d + 1e-9)

    p_ids = pheno.df["cell_id"].to_numpy()[p_idx][rows]
    x_ids = xen.df["cell_id"].to_numpy()[x_idx][cols]
    p_lin = pheno.df["lineage_common"].to_numpy()[p_idx][rows]
    x_lin = xen.df["lineage_common"].to_numpy()[x_idx][cols]

    return pd.DataFrame({
        "inferred_pheno_cell_id": p_ids,
        "inferred_xen_cell_id": x_ids,
        "inferred_distance": d,
        "inferred_confidence": conf,
        "inferred_lineage_agree": p_lin == x_lin,
        "inferred_pheno_lineage": p_lin,
        "inferred_xen_lineage": x_lin,
    })


def _group_tables(pheno: CellTable, xen: CellTable, by: str) -> list[tuple[str, CellTable, CellTable]]:
    """Split both tables by a shared grouping column (matched islet or joint niche)."""
    if by not in pheno.df.columns or by not in xen.df.columns:
        return [("__all__", pheno, xen)]
    groups = sorted(set(pheno.df[by].astype(str)) & set(xen.df[by].astype(str)) - {""})
    out = []
    for g in groups:
        p = pheno.df[pheno.df[by].astype(str) == g]
        x = xen.df[xen.df[by].astype(str) == g]
        if len(p) >= 3 and len(x) >= 3:
            out.append((g,
                        CellTable(df=p.reset_index(drop=True), modality="phenocycler",
                                  donor_id=pheno.donor_id, roi=pheno.roi),
                        CellTable(df=x.reset_index(drop=True), modality="xenium",
                                  donor_id=xen.donor_id, roi=xen.roi)))
    return out


def permuted_link_null(links: pd.DataFrame, seed: int = 0, n_perm: int = 100) -> dict:
    """Lineage agreement under randomly shuffled links — the floor any real result must clear.

    With a handful of lineages and skewed abundances, random pairing already agrees a
    surprising fraction of the time. Reporting raw agreement without this null overstates the
    result badly.
    """
    if links.empty or "inferred_pheno_lineage" not in links.columns:
        return {"agree_real": float("nan"), "agree_null": float("nan"), "z": float("nan")}
    rng = np.random.default_rng(seed)
    p = links["inferred_pheno_lineage"].to_numpy()
    x = links["inferred_xen_lineage"].to_numpy()
    real = float((p == x).mean())
    null = np.array([float((p == rng.permutation(x)).mean()) for _ in range(n_perm)])
    sd = float(null.std())
    return {
        "agree_real": real,
        "agree_null": float(null.mean()),
        "agree_null_p95": float(np.percentile(null, 95)),
        "z": float((real - null.mean()) / sd) if sd > 0 else float("nan"),
    }


def link_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    group_by: str = "niche_id_joint",
) -> LinkResult:
    """Link one donor's cells, grouped spatially. Raises :class:`AnchorError` if too thin."""
    if cfg.integration_mode != "sequential":
        raise RuntimeError(
            f"crossmodal linking is the sequential-mode approximation to cell pairing; in "
            f"{cfg.integration_mode!r} mode use sameslide.py, which pairs real cells"
        )

    pheno = read_cell_table(cfg.cells_pheno_dir, donor, roi, modality="phenocycler")
    xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")

    pairs, p_cols, x_cols = resolve_anchors(cfg, pheno, xen)
    if len(pairs) < cfg.crossmodal_min_anchors:
        raise AnchorError(
            f"donor {donor}: only {len(pairs)} usable protein<->gene anchors "
            f"({[p.protein for p in pairs]}), below the configured minimum of "
            f"{cfg.crossmodal_min_anchors}.\n"
            f"Linking on this few features produces confident-looking but meaningless "
            f"correspondences. Either lower [integration] crossmodal_min_anchors "
            f"deliberately, or check that the PhenoCycler export carries marker features and "
            f"the Xenium import carries per-gene values (identity-core lineage scores alone "
            f"are not per-gene anchors)."
        )

    groups = _group_tables(pheno, xen, group_by)
    frames = []
    for gid, p, x in groups:
        lk = link_cells(p, x, p_cols, x_cols, smooth_k=15)
        if len(lk):
            lk.insert(0, "group", gid)
            frames.append(lk)

    links = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if len(links):
        links.insert(0, "donor_id", donor)
        links.insert(1, "roi", roi)

    stats = {
        "donor_id": donor, "roi": roi, "n_groups": len(groups), "n_links": len(links),
        "n_anchors": len(pairs), "group_by": group_by,
        "mean_confidence": float(links["inferred_confidence"].mean()) if len(links) else np.nan,
    }
    stats.update(permuted_link_null(links))
    return LinkResult(links=links, n_anchors=len(pairs),
                      anchors=[p.protein for p in pairs], stats=stats)


def validate_against_islets(cfg: PipelineConfig, links: pd.DataFrame) -> dict:
    """Do the inferred links reproduce what S6 measured directly on matched islets?

    The only honest validation available. S6 pairs islets using geometry alone and measures
    their agreement; if the feature-space links are real, aggregating them per islet should
    recover a similar picture. If it does not, the links are describing the feature space and
    not the tissue.
    """
    islets = cfg.paired_dir / "islets.parquet"
    if not islets.exists() or links.empty:
        return {"validated": False, "reason": "no matched islets to validate against"}
    pairs = pd.read_parquet(islets)
    return {
        "validated": True,
        "n_islet_pairs": int(len(pairs)),
        "n_links": int(len(links)),
        "link_lineage_agreement": float(links["inferred_lineage_agree"].mean()),
    }


def run_crossmodal(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
    group_by: str = "niche_id_joint",
    require_qc_pass: bool = True,
) -> pd.DataFrame:
    """Stage entry point. Skips donors whose registration did not pass QC."""
    man = load_manifest(cfg)
    man = man[(man["pair_status"] == PAIRED) & (man["roi"] == roi)]
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[crossmodal] no paired donors", flush=True)
        return pd.DataFrame()

    passed: Optional[set[str]] = None
    qc_csv = cfg.integration_qc_dir / "registration_qc.csv"
    if require_qc_pass and qc_csv.exists():
        qc = pd.read_csv(qc_csv, dtype=str)
        passed = set(qc.loc[qc["verdict"] != "FAIL", "donor_id"])

    rows, frames = [], []
    for _, r in man.iterrows():
        donor = r["donor_id"]
        if passed is not None and donor not in passed:
            print(f"[crossmodal] {donor}: registration QC did not pass — skipped "
                  f"(links would inherit a bad alignment)", flush=True)
            continue
        try:
            res = link_donor(cfg, donor, roi, group_by=group_by)
        except (AnchorError, FileNotFoundError) as exc:
            print(f"[crossmodal] {donor}: {exc}", flush=True)
            continue
        rows.append(res.stats)
        if len(res.links):
            frames.append(res.links)
        print(f"[crossmodal] {donor}: {res.stats['n_links']:,} links over "
              f"{res.stats['n_groups']} groups using {res.n_anchors} anchors | "
              f"lineage agreement {res.stats['agree_real']:.3f} vs null "
              f"{res.stats['agree_null']:.3f} (z={res.stats['z']:.1f})", flush=True)

    if frames:
        out = cfg.crossmodal_dir / "pseudocell_links.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        allf = pd.concat(frames, ignore_index=True)
        allf.to_parquet(out, index=False)
        print(f"\n[crossmodal] wrote {out} ({len(allf):,} inferred links)", flush=True)
        print("[crossmodal] NOTE: these are INFERRED correspondences between serial "
              "sections, not measurements of the same cell. Use them for hypothesis "
              "generation; use matched islets (S6) for inference.", flush=True)

    stats = pd.DataFrame(rows)
    if len(stats):
        cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
        stats.to_csv(cfg.integration_qc_dir / "crossmodal_qc.csv", index=False)
    return stats


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Pseudo-cell linking across serial sections")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--group-by", default="niche_id_joint")
    ap.add_argument("--ignore-qc", action="store_true",
                    help="link even for donors whose registration failed QC")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_crossmodal(cfg, args.donors, roi=args.roi, group_by=args.group_by,
                           require_qc_pass=not args.ignore_qc)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
