"""
S7 — niches: two views of spatial organisation, one of which needs no registration at all.

**Niche view (primary).** Mirrors the Xenium pipeline's own niche definition exactly:
k-nearest-neighbour composition vector per cell (K=50), KMeans into 12 niches, then two
passes of majority-vote spatial smoothing (K=30). Running the identical algorithm on both
modalities, over the *same* harmonised lineage vocabulary, makes the two niche assignments
directly comparable.

The important property is that this comparison **requires no registration**. A niche is a
description of local tissue composition — "endocrine core", "immune-infiltrated islet
margin", "acinar parenchyma" — and those categories exist independently of where they sit on
the slide. So niche composition can be compared between two serial sections even when
registration fails, or when the sections are too far apart to align, or before any images are
mounted. For a sequential-slide design that is the most robust result available, and it is
the one this module leads with.

**Grid view (secondary).** A hex grid over the registered overlap, per-bin composition from
each modality, then per-bin similarity and a CCA for shared axes. This one *does* depend on
registration, which is precisely why it is useful: it is a spatially-resolved check on
whether the alignment is real.

Hex rather than square bins: hexagons have uniform neighbour distance, so a cell near a bin
edge is not arbitrarily closer to some neighbours than others, and the bin-to-bin comparisons
are not biased by the grid's own anisotropy.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, partition_path, read_cell_table, write_cell_table
from .manifest import PAIRED, load_manifest
from .register import registration_dir
from .transform import Transform
from .vocab import COMMON_LINEAGES

#: Lineages that form the composition vector. 'Unassigned' is excluded — it is an absence of
#: information, and including it would let two sections agree by both not knowing.
NICHE_LINEAGES = tuple(c for c in COMMON_LINEAGES if c != "Unassigned")


def composition_vectors(
    table: CellTable,
    *,
    k: int = 50,
    lineages: Sequence[str] = NICHE_LINEAGES,
    registered: bool = False,
) -> np.ndarray:
    """(n_cells, n_lineages) local composition from each cell's k nearest neighbours."""
    from scipy.spatial import cKDTree

    xy = table.xy(registered=registered)
    codes = pd.Categorical(table.df["lineage_common"], categories=list(lineages)).codes
    n_lin = len(lineages)
    if len(xy) == 0:
        return np.zeros((0, n_lin), dtype=np.float32)

    kk = int(min(k, len(xy)))
    tree = cKDTree(xy)
    _d, idx = tree.query(xy, k=kk)
    if idx.ndim == 1:
        idx = idx[:, None]

    neigh = codes[idx]                                   # (n, kk), -1 for unknown lineage
    out = np.zeros((len(xy), n_lin), dtype=np.float32)
    for j in range(n_lin):
        out[:, j] = (neigh == j).sum(axis=1)
    total = out.sum(axis=1, keepdims=True)
    np.divide(out, np.maximum(total, 1), out=out)
    return out


def majority_smooth(labels: np.ndarray, xy: np.ndarray, *, k: int = 30,
                    passes: int = 2) -> np.ndarray:
    """Replace each label with the majority among its k nearest neighbours, ``passes`` times.

    Cleans up the salt-and-pepper assignment KMeans produces on noisy composition vectors, and
    is what makes niches read as contiguous tissue regions rather than speckle.
    """
    from scipy.spatial import cKDTree

    if len(xy) == 0:
        return labels
    tree = cKDTree(xy)
    kk = int(min(k, len(xy)))
    _d, idx = tree.query(xy, k=kk)
    if idx.ndim == 1:
        idx = idx[:, None]

    out = labels.copy()
    n_lab = int(labels.max()) + 1 if len(labels) else 0
    for _ in range(passes):
        neigh = out[idx]
        counts = np.zeros((len(out), max(n_lab, 1)), dtype=np.int32)
        for j in range(max(n_lab, 1)):
            counts[:, j] = (neigh == j).sum(axis=1)
        out = counts.argmax(axis=1)
    return out


def call_niches(
    tables: Sequence[CellTable],
    *,
    k: int = 50,
    n_niches: int = 12,
    smooth_k: int = 30,
    seed: int = 0,
    registered: bool = False,
) -> dict[tuple[str, str, str], np.ndarray]:
    """Assign niches jointly across sections. Keyed by ``(modality, donor_id, roi)``.

    Fitting KMeans on the *pooled* composition vectors is what makes niche 7 mean the same
    thing in both modalities. Clustering each section separately would produce two unrelated
    labellings whose only correspondence would have to be guessed afterwards — which is
    exactly the ambiguity this step exists to remove.
    """
    from sklearn.cluster import KMeans

    comps, spans = [], []
    for t in tables:
        c = composition_vectors(t, k=k, registered=registered)
        spans.append((t.modality, t.donor_id, t.roi, len(c)))
        comps.append(c)
    if not comps:
        return {}

    pooled = np.vstack(comps)
    if len(pooled) < n_niches:
        n_niches = max(1, len(pooled))

    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=10)
    labels = km.fit_predict(pooled.astype(np.float32))

    out, cursor = {}, 0
    for (modality, donor, roi, n), table in zip(spans, tables):
        lab = labels[cursor:cursor + n]
        cursor += n
        lab = majority_smooth(lab, table.xy(registered=registered), k=smooth_k, passes=2)
        out[(modality, donor, roi)] = lab
    return out


def niche_profiles(table: CellTable, niches: np.ndarray,
                   lineages: Sequence[str] = NICHE_LINEAGES) -> pd.DataFrame:
    """Mean lineage composition of each niche — what makes a niche id interpretable."""
    df = pd.DataFrame({"niche": niches, "lineage": table.df["lineage_common"].to_numpy()})
    tab = pd.crosstab(df["niche"], df["lineage"])
    for lineage in lineages:
        if lineage not in tab.columns:
            tab[lineage] = 0
    tab = tab.loc[:, list(lineages)]
    prof = tab.div(tab.sum(axis=1).replace(0, 1), axis=0)
    prof.insert(0, "n_cells", tab.sum(axis=1))
    prof.insert(0, "modality", table.modality)
    prof.insert(0, "donor_id", table.donor_id)
    return prof.reset_index()


def niche_agreement(pheno_profiles: pd.DataFrame, xen_profiles: pd.DataFrame,
                    lineages: Sequence[str] = NICHE_LINEAGES) -> pd.DataFrame:
    """How similar is each niche's composition between the two modalities?

    Because niches were fitted jointly, niche ``j`` should describe the same kind of
    neighbourhood in both. A large divergence means either the lineage vocabularies are not
    as harmonised as intended, or one modality cannot see something the other can.
    """
    rows = []
    common = sorted(set(pheno_profiles["niche"]) & set(xen_profiles["niche"]))
    for niche in common:
        p = pheno_profiles.loc[pheno_profiles["niche"] == niche, list(lineages)].to_numpy().ravel()
        x = xen_profiles.loc[xen_profiles["niche"] == niche, list(lineages)].to_numpy().ravel()
        if p.size == 0 or x.size == 0:
            continue
        denom = np.linalg.norm(p) * np.linalg.norm(x)
        cos = float(p @ x / denom) if denom > 0 else float("nan")
        rows.append({
            "niche": niche,
            "cosine": cos,
            "l1_distance": float(np.abs(p - x).sum()),
            "pheno_n_cells": int(pheno_profiles.loc[pheno_profiles["niche"] == niche,
                                                    "n_cells"].iloc[0]),
            "xen_n_cells": int(xen_profiles.loc[xen_profiles["niche"] == niche,
                                                "n_cells"].iloc[0]),
            "dominant_lineage": list(lineages)[int(np.argmax(p + x))],
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# Hex grid (registration-dependent)
# --------------------------------------------------------------------------- #

def hex_bin(xy_um: np.ndarray, size_um: float) -> np.ndarray:
    """Axial hex-cell index (q, r) for each point. Flat-topped, ``size_um`` = circumradius."""
    xy = np.asarray(xy_um, dtype=np.float64)
    if len(xy) == 0:
        return np.zeros((0, 2), dtype=np.int64)
    q = (2.0 / 3.0 * xy[:, 0]) / size_um
    r = (-1.0 / 3.0 * xy[:, 0] + np.sqrt(3.0) / 3.0 * xy[:, 1]) / size_um

    # Cube rounding — plain rounding of axial coordinates lands points in the wrong hex near
    # the shared edges.
    x, z = q, r
    y = -x - z
    rx, ry, rz = np.round(x), np.round(y), np.round(z)
    dx, dy, dz = np.abs(rx - x), np.abs(ry - y), np.abs(rz - z)
    fix_x = (dx > dy) & (dx > dz)
    fix_z = (~fix_x) & (dz > dy)
    rx = np.where(fix_x, -ry - rz, rx)
    rz = np.where(fix_z, -rx - ry, rz)
    return np.column_stack([rx, rz]).astype(np.int64)


def grid_composition(
    table: CellTable,
    size_um: float,
    *,
    registered: bool = False,
    lineages: Sequence[str] = NICHE_LINEAGES,
) -> pd.DataFrame:
    """Per-hex-bin lineage composition. Index columns ``q``, ``r``; one column per lineage."""
    xy = table.xy(registered=registered)
    qr = hex_bin(xy, size_um)
    df = pd.DataFrame({"q": qr[:, 0], "r": qr[:, 1],
                       "lineage": table.df["lineage_common"].to_numpy()})
    tab = pd.crosstab([df["q"], df["r"]], df["lineage"])
    for lineage in lineages:
        if lineage not in tab.columns:
            tab[lineage] = 0
    tab = tab.loc[:, list(lineages)]
    n = tab.sum(axis=1)
    comp = tab.div(n.replace(0, 1), axis=0)
    comp["n_cells"] = n
    return comp.reset_index()


def compare_grids(
    pheno: pd.DataFrame,
    xen: pd.DataFrame,
    *,
    min_cells: int = 20,
    lineages: Sequence[str] = NICHE_LINEAGES,
) -> tuple[pd.DataFrame, dict]:
    """Join two grid compositions on (q, r) and score their agreement.

    ``min_cells`` drops sparse bins at the tissue edge, where a handful of cells produces a
    composition vector that is mostly sampling noise and would dominate the correlation.
    """
    lin = list(lineages)
    merged = pheno.merge(xen, on=["q", "r"], suffixes=("_pheno", "_xen"))
    merged = merged[(merged["n_cells_pheno"] >= min_cells)
                    & (merged["n_cells_xen"] >= min_cells)].reset_index(drop=True)
    if merged.empty:
        return merged, {"n_bins": 0}

    P = merged.loc[:, [f"{c}_pheno" for c in lin]].to_numpy(np.float64)
    X = merged.loc[:, [f"{c}_xen" for c in lin]].to_numpy(np.float64)

    num = (P * X).sum(axis=1)
    den = np.linalg.norm(P, axis=1) * np.linalg.norm(X, axis=1)
    merged["cosine"] = np.where(den > 0, num / np.maximum(den, 1e-12), np.nan)

    from scipy import stats as sps

    stats = {"n_bins": int(len(merged)),
             "median_cosine": float(np.nanmedian(merged["cosine"]))}
    for j, lineage in enumerate(lin):
        if P[:, j].std() == 0 or X[:, j].std() == 0:
            continue
        rho, p = sps.spearmanr(P[:, j], X[:, j])
        stats[f"spearman_{lineage}"] = float(rho)
        stats[f"p_{lineage}"] = float(p)

    # CCA: are there shared spatial axes even where individual lineages disagree? Guarded
    # because CCA needs more bins than components and fails noisily otherwise.
    n_comp = min(3, P.shape[1], X.shape[1], max(1, len(merged) - 1))
    if len(merged) > n_comp + 2:
        try:
            from sklearn.cross_decomposition import CCA

            cca = CCA(n_components=n_comp, max_iter=1000)
            U, V = cca.fit_transform(P, X)
            for c in range(n_comp):
                if U[:, c].std() > 0 and V[:, c].std() > 0:
                    stats[f"cca_{c + 1}"] = float(np.corrcoef(U[:, c], V[:, c])[0, 1])
        except Exception as exc:  # noqa: BLE001
            stats["cca_error"] = str(exc)
    return merged, stats


def run_grid(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
    use_registration: bool = True,
) -> pd.DataFrame:
    """Stage entry point: niches for every donor, plus the grid comparison where registered."""
    man = load_manifest(cfg)
    man = man[(man["pair_status"] == PAIRED) & (man["roi"] == roi)]
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[grid] no paired donors", flush=True)
        return pd.DataFrame()

    cfg.niches_dir.mkdir(parents=True, exist_ok=True)
    rows, all_profiles, all_bins = [], [], []

    for _, r in man.iterrows():
        donor = r["donor_id"]
        try:
            pheno = read_cell_table(cfg.cells_pheno_dir, donor, roi, modality="phenocycler")
            xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")
        except FileNotFoundError as exc:
            print(f"[grid] {donor}: {exc}", flush=True)
            continue

        # --- niche view: no registration needed ---
        niches = call_niches([pheno, xen], k=cfg.niche_k, n_niches=cfg.niche_n,
                             smooth_k=cfg.niche_smooth_k)
        p_lab = niches[("phenocycler", donor, roi)]
        x_lab = niches[("xenium", donor, roi)]

        pheno.df["niche_id_joint"] = p_lab.astype(str)
        xen.df["niche_id_joint"] = x_lab.astype(str)
        write_cell_table(pheno.df, partition_path(cfg.cells_pheno_dir, donor, roi))
        write_cell_table(xen.df, partition_path(cfg.cells_xen_dir, donor, roi))

        p_prof = niche_profiles(pheno, p_lab)
        x_prof = niche_profiles(xen, x_lab)
        all_profiles.extend([p_prof, x_prof])
        agree = niche_agreement(p_prof, x_prof)

        row = {"donor_id": donor, "roi": roi,
               "n_niches": int(max(p_lab.max(), x_lab.max()) + 1),
               "niche_median_cosine": float(agree["cosine"].median()) if len(agree) else np.nan}

        # --- grid view: needs registration ---
        reg = registration_dir(cfg, donor, roi)
        if use_registration and (reg / "transform.json").exists():
            tf = Transform.load(reg)
            xen_reg = xen.df.copy()
            tf.apply_to_frame(xen_reg)
            xen_r = CellTable(df=xen_reg, modality="xenium", donor_id=donor, roi=roi)
            write_cell_table(xen_reg, partition_path(cfg.cells_xen_dir, donor, roi))

            p_bins = grid_composition(pheno, cfg.grid_um)
            x_bins = grid_composition(xen_r, cfg.grid_um, registered=True)
            merged, stats = compare_grids(p_bins, x_bins)
            merged.insert(0, "donor_id", donor)
            all_bins.append(merged)
            row.update({k: v for k, v in stats.items() if not k.startswith("p_")})
        else:
            row["n_bins"] = 0
            row["grid_note"] = "no registration — niche view only"

        rows.append(row)
        print(f"[grid] {donor}: niche cosine median "
              f"{row.get('niche_median_cosine', float('nan')):.3f}, "
              f"{row.get('n_bins', 0)} comparable hex bins", flush=True)

    if all_profiles:
        pd.concat(all_profiles, ignore_index=True).to_csv(
            cfg.niches_dir / "niche_profiles.csv", index=False)
    if all_bins:
        allb = pd.concat(all_bins, ignore_index=True)
        allb.to_parquet(cfg.niches_dir / "grid_bins.parquet", index=False)
        print(f"[grid] wrote {cfg.niches_dir/'grid_bins.parquet'} ({len(allb)} bins)",
              flush=True)

    out = pd.DataFrame(rows)
    if len(out):
        out.to_csv(cfg.integration_qc_dir / "niche_qc.csv", index=False)
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Joint niches + registered-frame grid comparison")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--no-registration", action="store_true",
                    help="niche view only (works without any registration)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
    out = run_grid(cfg, args.donors, roi=args.roi,
                   use_registration=not args.no_registration)
    if len(out):
        print()
        print(out.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
