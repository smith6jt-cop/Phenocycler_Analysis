"""
S1b — Xenium outputs -> the common cell-table contract.

Three readers behind one interface, tried in order of how much they know:

``_from_h5ad``    ``{sample}_phenotyped.h5ad`` — the richest source. Carries the two-tier
                  identity-core labels (``celltype_broad`` / ``celltype`` /
                  ``celltype_fine``), ``broad_confidence``, ``broad_n_identity_genes``,
                  ``cycling``, ``spatial_niche`` and ``obsm['spatial']`` in microns.
``_from_zarr``    the SpatialData store — the table plus real cell polygons. Used when
                  phenotyping has not run, or when boundaries are needed.
``_from_bundle``  ``cells.parquet`` + ``cell_feature_matrix.h5`` straight from the
                  instrument — coordinates and morphology only, no labels.

Two things this stage fixes rather than inherits.

**The join key.** No Xenium artifact carries ``donor_id`` or ``roi``: everything is keyed by
the XETG slide serial, which encodes no donor and no region, and a single serial can carry
both a pancreas and a lymph-node section. The manifest resolves the mapping; this stage
*stamps* it into the data, which is where it needed to be all along.

**Evidence-positivity.** Xenium's broad tier is an expression-masked scANVI ``predict()``:
cells with no detected identity-core gene fall back to an absolute-score argmax and are
flagged ``broad_lowsignal``. Those calls are weaker than the classifier's, in the same way
PhenoCycler's ``epi_default`` calls are weaker than a positive gate. Both land in the
contract's ``evidence_positive`` so composition comparisons can be run on like against like.

Optional dependencies (``anndata``, ``spatialdata``, ``h5py``) are imported lazily and only
in the reader that needs them, so importing this module in the lean core environment works
and the error names the extra when a reader is actually used.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import as_feature, coerce_cell_table, partition_path, write_cell_table
from .manifest import PAIRED, load_manifest
from .vocab import (COMMON_LINEAGES, PanelTaxonomy, load_panel_taxonomy, map_xenium_series)

#: obs columns worth carrying when present. Everything else stays in the h5ad.
OBS_PASSTHROUGH = (
    "celltype_broad", "celltype", "celltype_fine", "celltype_lineage", "celltype_hier",
    "broad_confidence", "broad_n_identity_genes", "broad_lowsignal", "cycling",
    "spatial_niche", "donor_sex", "leiden_scvi",
)


class ImportError_(RuntimeError):
    """A Xenium sample cannot be imported. (Named to avoid shadowing the builtin.)"""


def _require(module: str, extra: str = "integration"):
    """Import an optional dependency with an actionable message."""
    try:
        return __import__(module)
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError_(
            f"reading this Xenium format needs `{module}`, which is not installed.\n"
            f"  pip install -e '.[{extra}]'\n"
            f"or use environment-integration.yml."
        ) from exc


def _obs_str(obs: pd.DataFrame, col: str, default: str = "") -> pd.Series:
    """One obs column as strings, with nulls replaced by ``default``.

    Order matters here. ``astype(str)`` converts NaN to the literal string ``"nan"``, after
    which ``fillna`` has nothing left to fill — so nulls must be replaced *first*. Getting
    this backwards has two consequences, and the second is worse than it looks:

    * ``"nan"`` leaks into ``lineage_native`` / ``celltype_fine`` and shows up as a phantom
      label in the "did not map to the common vocabulary" warning;
    * the ``celltype`` -> ``celltype_fine`` fallback in :func:`import_sample` tests
      ``fine.eq("").all()``, and an all-``"nan"`` series is not all-empty — so an h5ad with
      a null ``celltype`` but populated ``celltype_fine`` would silently lose its fine
      labels rather than falling back.

    ``where`` rather than ``fillna`` because obs columns are frequently Categorical, and
    ``fillna`` on a Categorical raises unless the fill value is already a category.
    """
    if col not in obs.columns:
        return pd.Series([default] * len(obs), index=obs.index, dtype="object")
    s = obs[col]
    if isinstance(s.dtype, pd.CategoricalDtype):
        s = s.astype(object)
    return s.where(s.notna(), default).astype(str)


def _obs_num(obs: pd.DataFrame, col: str) -> pd.Series:
    if col not in obs.columns:
        return pd.Series(np.nan, index=obs.index, dtype="float64")
    return pd.to_numeric(obs[col], errors="coerce")


def identity_core_scores(adata, tax: PanelTaxonomy, layer: str = "lognorm") -> pd.DataFrame:
    """Per-cell, per-lineage mean expression over each identity-core gene set.

    Deliberately the *absolute* mean over the curated core, matching how the Xenium
    phenotyper scores lineages. The upstream repo learned this the hard way: z-scoring per
    lineage forces a near-uniform argmax and lets rare classes vacuum up marker-negative
    cells, and scoring the full fat subpanel instead of the identity core means a
    proliferating acinar cell scores high on the B-cell panel because that panel is padded
    with cell-cycle genes. Absolute mean over the core avoids both.
    """
    var_names = list(adata.var_names)
    present = set(var_names)
    index = {g: i for i, g in enumerate(var_names)}

    X = adata.layers[layer] if (layer and layer in getattr(adata, "layers", {})) else adata.X

    out = {}
    lineages = sorted({v for v in tax.panel_to_lineage.values() if v})
    for lineage in lineages:
        genes = [g for g in tax.lineage_genes(lineage) if g in present]
        if not genes:
            continue
        cols = [index[g] for g in genes]
        sub = X[:, cols]
        mean = np.asarray(sub.mean(axis=1)).ravel() if hasattr(sub, "mean") else sub.mean(axis=1)
        out[lineage] = np.asarray(mean, dtype=np.float32)
    return pd.DataFrame(out, index=pd.RangeIndex(adata.n_obs))


def _from_h5ad(path: Path, tax: Optional[PanelTaxonomy]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read ``{sample}_phenotyped.h5ad`` -> (obs, features)."""
    ad = _require("anndata")
    adata = ad.read_h5ad(path)
    obs = adata.obs.copy().reset_index(drop=False)
    obs = obs.rename(columns={obs.columns[0]: "cell_id"})

    if "spatial" not in adata.obsm:
        raise ImportError_(
            f"{path} has no obsm['spatial'] — integration needs micron centroids "
            f"(the Xenium ingest writes them; a slimmed h5ad may have dropped them)"
        )
    xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    obs["x_um"] = xy[:, 0]
    obs["y_um"] = xy[:, 1]

    feats = pd.DataFrame(index=obs.index)
    if tax is not None:
        try:
            scores = identity_core_scores(adata, tax)
            for lineage in scores.columns:
                feats[as_feature(f"score_{lineage}")] = scores[lineage].to_numpy()
        except Exception as exc:  # noqa: BLE001 - scores are a bonus, not a requirement
            print(f"[import_xenium] identity-core scores unavailable for {path.name}: {exc}",
                  flush=True)
    return obs, feats


def _from_zarr(path: Path, tax: Optional[PanelTaxonomy]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a SpatialData zarr -> (obs, features)."""
    _require("spatialdata")
    import spatialdata as sd

    sdata = sd.read_zarr(path)
    tables = getattr(sdata, "tables", {})
    if not tables:
        raise ImportError_(f"{path} contains no tables")
    adata = tables.get("table") or next(iter(tables.values()))

    obs = adata.obs.copy().reset_index(drop=False)
    obs = obs.rename(columns={obs.columns[0]: "cell_id"})

    if "spatial" in adata.obsm:
        xy = np.asarray(adata.obsm["spatial"], dtype=np.float64)
    elif {"x_centroid", "y_centroid"} <= set(obs.columns):
        xy = obs.loc[:, ["x_centroid", "y_centroid"]].to_numpy(np.float64)
    else:
        raise ImportError_(f"{path}: no spatial coordinates in obsm or obs")
    obs["x_um"] = xy[:, 0]
    obs["y_um"] = xy[:, 1]

    feats = pd.DataFrame(index=obs.index)
    if tax is not None and adata.n_vars:
        try:
            scores = identity_core_scores(adata, tax, layer="")
            for lineage in scores.columns:
                feats[as_feature(f"score_{lineage}")] = scores[lineage].to_numpy()
        except Exception as exc:  # noqa: BLE001
            print(f"[import_xenium] identity-core scores unavailable for {path.name}: {exc}",
                  flush=True)
    return obs, feats


def _from_bundle(path: Path, tax: Optional[PanelTaxonomy]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a raw Xenium bundle directory -> (obs, features).

    Coordinates and morphology only — an un-phenotyped bundle has no labels, so the contract
    gets ``lineage_native = ''`` and everything downstream that needs labels will say so.
    """
    cells = Path(path) / "cells.parquet"
    if not cells.exists():
        alt = Path(path) / "cells.csv.gz"
        if not alt.exists():
            raise ImportError_(f"{path} has neither cells.parquet nor cells.csv.gz")
        obs = pd.read_csv(alt)
    else:
        obs = pd.read_parquet(cells)

    ren = {"x_centroid": "x_um", "y_centroid": "y_um"}
    obs = obs.rename(columns={k: v for k, v in ren.items() if k in obs.columns})
    if "cell_id" not in obs.columns:
        obs = obs.reset_index(drop=False).rename(columns={obs.columns[0]: "cell_id"})
    if not {"x_um", "y_um"} <= set(obs.columns):
        raise ImportError_(f"{path}/cells.parquet has no x_centroid/y_centroid columns")
    return obs, pd.DataFrame(index=obs.index)


def read_xenium(
    source: Path,
    tax: Optional[PanelTaxonomy] = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Dispatch on what ``source`` actually is."""
    source = Path(source)
    if not source.exists():
        raise ImportError_(f"Xenium source does not exist: {source}")
    if source.is_file() and source.suffix == ".h5ad":
        return _from_h5ad(source, tax)
    if source.suffix == ".zarr" or (source.is_dir() and (source / ".zgroup").exists()):
        return _from_zarr(source, tax)
    if source.is_dir():
        return _from_bundle(source, tax)
    raise ImportError_(f"unrecognised Xenium source: {source}")


def import_sample(
    source: Path,
    *,
    donor_id: str,
    roi: str,
    sample_key: str = "",
    tax: Optional[PanelTaxonomy] = None,
) -> pd.DataFrame:
    """Build the contract table for one Xenium section."""
    obs, feats = read_xenium(source, tax)

    broad = _obs_str(obs, "celltype_broad")
    fine = _obs_str(obs, "celltype")
    if fine.eq("").all():
        fine = _obs_str(obs, "celltype_fine")

    # Xenium's broad tier falls back to an absolute-score argmax when a cell has no detected
    # identity-core gene, flagging it `broad_lowsignal`. That is the same kind of weaker call
    # as PhenoCycler's `epi_default`, so it belongs in the same contract column.
    if "broad_lowsignal" in obs.columns:
        lowsignal = obs["broad_lowsignal"].fillna(False).astype(bool).to_numpy()
    else:
        lowsignal = np.zeros(len(obs), dtype=bool)

    area = pd.Series(np.nan, index=obs.index, dtype="float64")
    for col in ("cell_area", "cell_area_um2", "area"):
        if col in obs.columns:
            area = pd.to_numeric(obs[col], errors="coerce")
            break

    out = pd.DataFrame({
        "cell_id": obs["cell_id"].astype(str),
        "donor_id": str(donor_id),
        "roi": roi,
        "x_um": pd.to_numeric(obs["x_um"], errors="coerce"),
        "y_um": pd.to_numeric(obs["y_um"], errors="coerce"),
        "area_um2": area,
        "lineage_native": broad,
        "celltype_fine": fine,
        "confidence": _obs_num(obs, "broad_confidence"),
        "evidence_positive": ~lowsignal,
        "region": "",
        "structure_id": "",
    })
    out["lineage_common"] = map_xenium_series(broad, fine)

    # Provenance: the section identity that Xenium artifacts never carry.
    out["xenium_sample_key"] = sample_key
    for col in ("spatial_niche", "cycling", "broad_n_identity_genes", "donor_sex"):
        if col in obs.columns:
            out[col] = obs[col].astype(str).to_numpy()

    for c in feats.columns:
        out[c] = pd.to_numeric(feats[c], errors="coerce").to_numpy()

    if broad.eq("").all():
        print(f"[import_xenium] {source.name}: no cell-type labels found — importing "
              f"coordinates only. Run the Xenium phenotyping stage for a labelled table.",
              flush=True)
    else:
        unmapped = sorted(set(out.loc[out["lineage_common"] == "Unassigned", "lineage_native"])
                          - {""})
        if unmapped:
            print(f"[import_xenium] {source.name}: {unmapped} did not map to the common "
                  f"vocabulary (known: {[c for c in COMMON_LINEAGES if c != 'Unassigned']})",
                  flush=True)

    return coerce_cell_table(out, modality="xenium", donor_id=str(donor_id), roi=roi)


def _resolve_source(row: pd.Series) -> Optional[Path]:
    """Pick the richest available source for a manifest row."""
    for col in ("xenium_h5ad", "xenium_zarr", "xenium_bundle"):
        val = str(row.get(col, "") or "").strip()
        if val and Path(val).exists():
            return Path(val)
    return None


def run_import_xenium(
    cfg: PipelineConfig,
    *,
    donors: Optional[list[str]] = None,
    roi: str = "panc",
    only_paired: bool = True,
) -> pd.DataFrame:
    """Stage entry point: import every manifest row whose source is on disk."""
    man = load_manifest(cfg)
    if only_paired:
        man = man[man["pair_status"] == PAIRED]
    if roi:
        man = man[man["roi"] == roi]
    if donors:
        want = {str(d) for d in donors}
        man = man[man["donor_id"].isin(want)]

    if man.empty:
        print("[import_xenium] no manifest rows to import "
              "(check `python -m phenocycler.integration.manifest --report`)", flush=True)
        return pd.DataFrame(columns=["donor_id", "roi", "n_cells", "source"])

    try:
        tax = load_panel_taxonomy(cfg)
    except Exception as exc:  # noqa: BLE001 - scores are optional
        print(f"[import_xenium] panel taxonomy unavailable ({exc}); "
              f"importing without identity-core scores", flush=True)
        tax = None

    rows = []
    for _, r in man.iterrows():
        src = _resolve_source(r)
        if src is None:
            print(f"[import_xenium] donor {r['donor_id']} ({r['roi']}): no Xenium source on "
                  f"disk — skipped", flush=True)
            continue
        df = import_sample(src, donor_id=r["donor_id"], roi=r["roi"],
                           sample_key=r.get("xenium_sample_key", ""), tax=tax)
        path = partition_path(cfg.cells_xen_dir, r["donor_id"], r["roi"])
        write_cell_table(df, path)
        print(f"[import_xenium] {r['donor_id']}/{r['roi']}: {len(df):,} cells -> {path.parent}",
              flush=True)
        rows.append({"donor_id": r["donor_id"], "roi": r["roi"],
                     "n_cells": len(df), "source": str(src)})

    out = pd.DataFrame(rows)
    if len(out):
        print(f"\n[import_xenium] {len(out)} section(s), "
              f"{int(out['n_cells'].sum()):,} cells total", flush=True)
    return out


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Xenium -> integration cell-table contract")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--all", action="store_true",
                    help="import unpaired rows too (donor_unknown / xenium_only)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    summary = run_import_xenium(cfg, donors=args.donors, roi=args.roi,
                                only_paired=not args.all)
    if len(summary):
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
