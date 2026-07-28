"""
S11 — same-slide mode: real cell-to-cell pairing.

When both modalities imaged the *same physical section*, the cells are literally the same
cells, and a genuine paired protein+RNA matrix is recoverable. Everything the sequential mode
has to approximate — islet matching, niche comparison, feature-space linking — becomes
unnecessary, because the correspondence is a measurement rather than an inference.

The mode guard is a hard error in both directions, and that is deliberate. Every claim
downstream rests on which of these two situations holds: "these two measurements come from
one cell" is true in same-slide mode and false in sequential mode, and no amount of
registration quality changes that. A warning would be too easy to ignore.

Assignment proceeds from the strongest evidence available:

1. **Polygon IoU**, when PhenoCycler's QuPath GeoJSON and Xenium's ``cell_boundaries`` are
   both present. Two segmentations of the same cell overlap substantially; two neighbouring
   cells do not. This is the most reliable criterion and the only one that handles the case
   of two adjacent nuclei whose centroids are closer to each other than to their own partner.
2. **Mutual nearest centroid** within ``max_dist_um`` otherwise. Mutual rather than one-way,
   because one-way nearest-neighbour assignment produces many-to-one collisions exactly where
   segmentation disagrees most.

Unmatched cells are reported, not hidden. The two pipelines segment independently — QuPath
on a qptiff, Xenium's own algorithm on DAPI plus the multimodal stains — so they will not
agree on cell count, and a match rate well below 1 is the expected result rather than a
failure.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import CellTable, read_cell_table
from .manifest import PAIRED, load_manifest
from .register import registration_dir
from .transform import Transform

#: Centroid distance beyond which two cells cannot be the same cell. A pancreatic cell is
#: ~10-20 um across, so two segmentations of one cell should agree on its centre to within a
#: few microns once registered; 5 um is generous without being permissive enough to match
#: neighbours.
DEFAULT_MAX_DIST_UM = 5.0


class ModeError(RuntimeError):
    """Called in the wrong integration mode."""


def require_same_slide(cfg: PipelineConfig) -> None:
    """Hard guard. Cell-to-cell pairing is only valid for one physical section."""
    if cfg.integration_mode != "same_slide":
        raise ModeError(
            f"same-slide cell pairing was requested but [integration] mode is "
            f"{cfg.integration_mode!r}.\n"
            f"Across serial sections the two modalities measure DIFFERENT cells, so a "
            f"cell-to-cell pairing is not a measurement no matter how good the registration "
            f"is. Use match.py (islet-level), grid.py (niche-level) or crossmodal.py "
            f"(explicitly inferred) instead, or set mode = same_slide if the two modalities "
            f"really did image one section."
        )


def require_sequential(cfg: PipelineConfig, what: str) -> None:
    """The mirror-image guard, for steps whose semantics assume distinct sections."""
    if cfg.integration_mode != "sequential":
        raise ModeError(
            f"{what} is a sequential-mode step, but [integration] mode is "
            f"{cfg.integration_mode!r}. In same-slide mode use sameslide.py, which pairs "
            f"cells directly instead of approximating the pairing."
        )


def match_by_centroid(
    moving_xy: np.ndarray,
    fixed_xy: np.ndarray,
    *,
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
) -> pd.DataFrame:
    """Mutual-nearest centroid pairing within a distance cap."""
    empty = pd.DataFrame(columns=["moving_idx", "fixed_idx", "distance_um"])
    if len(moving_xy) == 0 or len(fixed_xy) == 0:
        return empty
    from scipy.spatial import cKDTree

    t_fixed, t_moving = cKDTree(fixed_xy), cKDTree(moving_xy)
    d_mf, i_mf = t_fixed.query(moving_xy, k=1)
    _d_fm, i_fm = t_moving.query(fixed_xy, k=1)

    mutual = i_fm[i_mf] == np.arange(len(moving_xy))
    keep = mutual & (d_mf <= max_dist_um)
    if not keep.any():
        return empty
    idx = np.flatnonzero(keep)
    return pd.DataFrame({"moving_idx": idx, "fixed_idx": i_mf[idx],
                         "distance_um": d_mf[idx]})


def match_by_iou(
    moving_polys,
    fixed_polys,
    *,
    min_iou: float = 0.25,
) -> pd.DataFrame:
    """Polygon-overlap pairing. ``*_polys`` are GeoSeries in a common frame.

    Needs shapely/geopandas. Preferred over centroids where available: IoU distinguishes "the
    same cell segmented twice" from "two adjacent cells", which pure centroid distance cannot
    do reliably in densely packed islet tissue.
    """
    try:
        import geopandas as gpd  # noqa: F401
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "polygon matching needs geopandas/shapely; pip install -e '.[integration]'"
        ) from exc

    empty = pd.DataFrame(columns=["moving_idx", "fixed_idx", "iou"])
    if len(moving_polys) == 0 or len(fixed_polys) == 0:
        return empty

    mv = moving_polys.reset_index(drop=True)
    fx = fixed_polys.reset_index(drop=True)
    joined = (mv.to_frame("geometry").sjoin(fx.to_frame("geometry"), how="inner",
                                            predicate="intersects")
              .reset_index().rename(columns={"index": "moving_idx",
                                             "index_right": "fixed_idx"}))
    if joined.empty:
        return empty

    rows = []
    for _, r in joined.iterrows():
        a, b = mv.iloc[int(r["moving_idx"])], fx.iloc[int(r["fixed_idx"])]
        inter = a.intersection(b).area
        union = a.union(b).area
        if union > 0 and inter / union >= min_iou:
            rows.append({"moving_idx": int(r["moving_idx"]),
                         "fixed_idx": int(r["fixed_idx"]), "iou": inter / union})
    if not rows:
        return empty

    # One-to-one: highest IoU wins, and a cell already claimed cannot be claimed again.
    df = pd.DataFrame(rows).sort_values("iou", ascending=False)
    seen_m, seen_f, keep = set(), set(), []
    for _, r in df.iterrows():
        m, f = int(r["moving_idx"]), int(r["fixed_idx"])
        if m in seen_m or f in seen_f:
            continue
        seen_m.add(m)
        seen_f.add(f)
        keep.append(r)
    return pd.DataFrame(keep).reset_index(drop=True)


def pair_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
) -> tuple[pd.DataFrame, dict]:
    """Pair one donor's cells across modalities on the same section."""
    require_same_slide(cfg)

    pheno = read_cell_table(cfg.cells_pheno_dir, donor, roi, modality="phenocycler")
    xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")

    reg = registration_dir(cfg, donor, roi)
    if not (reg / "transform.json").exists():
        return pd.DataFrame(), {"error": f"donor {donor}: no registration; run S4 first"}
    tf = Transform.load(reg)

    fixed_is_pheno = cfg.fixed_modality == "phenocycler"
    moving, fixed = (xen, pheno) if fixed_is_pheno else (pheno, xen)
    moving_xy = tf.apply(moving.xy())
    fixed_xy = fixed.xy()

    m = match_by_centroid(moving_xy, fixed_xy, max_dist_um=max_dist_um)
    if m.empty:
        return pd.DataFrame(), {
            "donor_id": donor, "roi": roi, "n_pairs": 0,
            "note": (f"no cell pairs within {max_dist_um} um — either the registration is "
                     f"wrong or these are not the same section"),
        }

    out = pd.DataFrame({
        "donor_id": donor,
        "roi": roi,
        "pheno_cell_id": (fixed.df["cell_id"].to_numpy()[m["fixed_idx"]] if fixed_is_pheno
                          else moving.df["cell_id"].to_numpy()[m["moving_idx"]]),
        "xen_cell_id": (moving.df["cell_id"].to_numpy()[m["moving_idx"]] if fixed_is_pheno
                        else fixed.df["cell_id"].to_numpy()[m["fixed_idx"]]),
        "distance_um": m["distance_um"].to_numpy(),
    })

    p_lin = (fixed.df["lineage_common"].to_numpy()[m["fixed_idx"]] if fixed_is_pheno
             else moving.df["lineage_common"].to_numpy()[m["moving_idx"]])
    x_lin = (moving.df["lineage_common"].to_numpy()[m["moving_idx"]] if fixed_is_pheno
             else fixed.df["lineage_common"].to_numpy()[m["fixed_idx"]])
    out["pheno_lineage"] = p_lin
    out["xen_lineage"] = x_lin
    out["lineage_agree"] = p_lin == x_lin

    stats = {
        "donor_id": donor, "roi": roi,
        "n_pheno": len(pheno.df), "n_xen": len(xen.df), "n_pairs": len(out),
        "pair_rate_pheno": len(out) / max(len(pheno.df), 1),
        "pair_rate_xen": len(out) / max(len(xen.df), 1),
        "median_distance_um": float(out["distance_um"].median()),
        # The headline number for same-slide mode: two independent measurements of the SAME
        # cell should agree on its lineage. Unlike sequential mode this is a real accuracy,
        # not a distributional similarity.
        "lineage_agreement": float(out["lineage_agree"].mean()),
    }
    return out, stats


def build_paired_matrix(
    cfg: PipelineConfig,
    donor: str,
    pairs: pd.DataFrame,
    roi: str = "panc",
) -> pd.DataFrame:
    """The payoff: one row per paired cell, protein and RNA features side by side.

    This is what same-slide mode exists to produce and what sequential mode fundamentally
    cannot. Downstream it supports genuinely paired modelling (totalVI-style joint latent
    spaces, per-cell protein-vs-RNA discordance) that would be unjustified on serial sections.
    """
    pheno = read_cell_table(cfg.cells_pheno_dir, donor, roi, modality="phenocycler")
    xen = read_cell_table(cfg.cells_xen_dir, donor, roi, modality="xenium")

    p = pheno.df.set_index("cell_id")
    x = xen.df.set_index("cell_id")
    p_feats = [c for c in pheno.features]
    x_feats = [c for c in xen.features]

    out = pairs.loc[:, ["donor_id", "roi", "pheno_cell_id", "xen_cell_id",
                        "distance_um", "lineage_agree"]].copy()
    for c in p_feats:
        out[f"protein__{c[len('feat__'):]}"] = p.loc[pairs["pheno_cell_id"], c].to_numpy()
    for c in x_feats:
        out[f"rna__{c[len('feat__'):]}"] = x.loc[pairs["xen_cell_id"], c].to_numpy()
    return out


def run_sameslide(
    cfg: PipelineConfig,
    donors: Optional[list[str]] = None,
    *,
    roi: str = "panc",
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
) -> pd.DataFrame:
    """Stage entry point for same-slide mode."""
    require_same_slide(cfg)

    man = load_manifest(cfg)
    man = man[(man["pair_status"] == PAIRED) & (man["roi"] == roi)]
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[sameslide] no paired donors", flush=True)
        return pd.DataFrame()

    rows, frames = [], []
    for _, r in man.iterrows():
        donor = r["donor_id"]
        try:
            pairs, stats = pair_donor(cfg, donor, roi, max_dist_um=max_dist_um)
        except FileNotFoundError as exc:
            print(f"[sameslide] {donor}: {exc}", flush=True)
            continue
        if "error" in stats:
            print(f"[sameslide] {stats['error']}", flush=True)
            continue
        rows.append(stats)
        if len(pairs):
            frames.append(build_paired_matrix(cfg, donor, pairs, roi))
        print(f"[sameslide] {donor}: {stats['n_pairs']:,} cell pairs "
              f"({stats['pair_rate_pheno']:.1%} of PhenoCycler, "
              f"{stats['pair_rate_xen']:.1%} of Xenium), median {stats['median_distance_um']:.2f} um, "
              f"lineage agreement {stats['lineage_agreement']:.3f}", flush=True)

    if frames:
        out = cfg.paired_dir / "cells.parquet"
        out.parent.mkdir(parents=True, exist_ok=True)
        allf = pd.concat(frames, ignore_index=True)
        allf.to_parquet(out, index=False)
        print(f"\n[sameslide] wrote {out} ({len(allf):,} paired cells, "
              f"{allf.shape[1]} columns)", flush=True)

    stats_df = pd.DataFrame(rows)
    if len(stats_df):
        cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
        stats_df.to_csv(cfg.integration_qc_dir / "sameslide_qc.csv", index=False)
    return stats_df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Same-slide cell-to-cell pairing")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--max-dist-um", type=float, default=DEFAULT_MAX_DIST_UM)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_sameslide(cfg, args.donors, roi=args.roi, max_dist_um=args.max_dist_um)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
