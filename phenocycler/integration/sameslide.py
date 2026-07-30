"""
S11 — same-slide mode: real cell-to-cell pairing.

When both modalities imaged the *same physical section*, the cells are literally the same
cells, and a genuine paired protein+RNA matrix is recoverable. Everything the sequential mode
has to approximate — islet matching, niche comparison, feature-space linking — becomes
unnecessary, because the correspondence is a measurement rather than an inference.

The mode guard is a hard error in both directions, and that is deliberate: same-slide mode
claims *whole-cell* correspondence for *every* cell, which serial sections cannot support.

That is narrower than what an earlier version of this docstring said. It claimed serial
sections never measure the same cell at all, which is simply false — pancreatic cells are
10-15 um and sections are 5 um, so most cells straddle the cut. What serial sections cannot
give is correspondence for every cell: only ~28% are detectable in both, the recoverable
population is heavily biased toward large cells (0% of 8 um lymphocytes), and each side sees
a different *part* of the cell. ``cellmatch.py`` (S6b) recovers the fraction that is real —
endocrine cells inside matched islets, under a verified local registration residual — and
refuses the rest. This module remains for the case where one section was imaged twice.

**How ``pair_donor`` matches.** Two criteria, selected by ``method``:

``iou`` (preferred)
    Polygon overlap. Two segmentations of one cell overlap substantially; two neighbouring
    cells do not. That distinction is what centroid distance cannot make reliably in packed
    tissue, where an adjacent nucleus is often nearer than the partner outline's centroid.
    Polygons come from :mod:`phenocycler.integration.boundaries`, which reads QuPath GeoJSON
    and Xenium ``cell_boundaries.parquet``, reconciles their units, and warps the moving side
    through the transform — non-rigid field included, since applying only the affine would
    leave exactly the local error the B-spline stage exists to remove.

``centroid``
    Mutual-nearest within ``max_dist_um``. Mutual rather than one-way, because one-way
    nearest-neighbour produces many-to-one collisions exactly where segmentation disagrees
    most.

``auto`` (default) tries IoU and falls back to centroids when polygons are unavailable,
**recording which ran** in ``stats['method']``. A silent fallback would make a polygon-loading
failure look like a successful centroid run, and the two have different error characteristics.

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
from .contract import read_cell_table
from .manifest import load_manifest, paired_rows, rows_for_scope
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


def _match(cfg, donor, roi, moving, fixed, tf, fixed_is_pheno, moving_xy, fixed_xy, *,
           method, max_dist_um, min_iou, pheno_base=None):
    """Run the requested matcher. Returns ``(pairs, method_used)``.

    ``auto`` prefers IoU and falls back to centroids, recording which ran — a silent fallback
    would make a polygon-loading failure look like a successful centroid run, and the two have
    different error characteristics. ``iou`` re-raises instead, so an explicit request cannot
    be quietly downgraded.
    """
    if method not in ("auto", "iou", "centroid"):
        raise ValueError(f"method must be auto|iou|centroid, got {method!r}")

    if method in ("auto", "iou"):
        try:
            mv_poly, fx_poly = _load_polygons(cfg, donor, roi, moving, fixed, tf,
                                              fixed_is_pheno, pheno_base=pheno_base)
            pairs = match_tables_by_iou(mv_poly, fx_poly, moving, fixed, min_iou=min_iou)
            return pairs, "iou"
        except Exception as exc:  # noqa: BLE001 - any loader failure falls back, loudly
            if method == "iou":
                raise
            print(f"[sameslide] {donor}/{roi}: polygon matching unavailable ({exc}); "
                  f"falling back to centroids", flush=True)

    return match_by_centroid(moving_xy, fixed_xy, max_dist_um=max_dist_um), "centroid"


def match_tables_by_iou(mv_poly, fx_poly, moving, fixed, *, min_iou):
    """IoU-match polygons, then re-index the result onto the cell tables.

    :func:`match_by_iou` returns positions into the GeoSeries it was handed, which are not the
    positions in the cell tables the caller indexes with. Both loaders restrict to ids the
    tables carry, so the lookup is total; the explicit filter is there because a KeyError here
    would be swallowed by ``auto``'s fallback and read as "no polygons available".
    """
    pairs = match_by_iou(mv_poly, fx_poly, min_iou=min_iou)
    if pairs.empty:
        return pairs
    mv_pos = {c: i for i, c in enumerate(moving.df["cell_id"])}
    fx_pos = {c: i for i, c in enumerate(fixed.df["cell_id"])}
    mv = [mv_pos.get(mv_poly.index[int(i)]) for i in pairs["moving_idx"]]
    fx = [fx_pos.get(fx_poly.index[int(i)]) for i in pairs["fixed_idx"]]
    keep = [j for j, (a, b) in enumerate(zip(mv, fx)) if a is not None and b is not None]
    pairs = pairs.iloc[keep].reset_index(drop=True)
    return pairs.assign(moving_idx=[mv[j] for j in keep], fixed_idx=[fx[j] for j in keep])


def _load_polygons(cfg, donor, roi, moving, fixed, tf, fixed_is_pheno, *, pheno_base=None):
    """Both sides' polygons in the FIXED frame, restricted to cells in the tables."""
    from .boundaries import load_pheno_polygons, load_xenium_polygons, warp_polygons

    xen = moving if fixed_is_pheno else fixed
    pheno_poly = load_pheno_polygons(cfg, donor, roi, base=pheno_base)

    man = rows_for_scope(load_manifest(cfg), roi=roi)
    row = man[man["donor_id"].astype(str) == str(donor)]
    bundle = str(row["xenium_bundle"].iloc[0]) if len(row) else ""
    if not bundle:
        raise RuntimeError(f"no xenium_bundle recorded for {donor}/{roi} in the manifest")
    xen_poly = load_xenium_polygons(Path(bundle), cell_ids=set(xen.df["cell_id"]))

    # Only the moving side is warped; the fixed frame is by definition already correct.
    if fixed_is_pheno:
        return warp_polygons(xen_poly, tf), pheno_poly
    return warp_polygons(pheno_poly, tf), xen_poly


def pair_donor(
    cfg: PipelineConfig,
    donor: str,
    roi: str = "panc",
    *,
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
    method: str = "auto",
    min_iou: float = 0.25,
    pheno_base=None,
) -> tuple[pd.DataFrame, dict]:
    """Pair one donor's cells across modalities on the same section.

    ``method`` is ``'iou'`` (polygon overlap), ``'centroid'``, or ``'auto'`` — try IoU and
    fall back to centroids when polygons are unavailable, recording which was used. IoU is
    preferred where it can run: two segmentations of one cell overlap substantially and two
    neighbours do not, which is the distinction centroid distance cannot make reliably in
    packed tissue, where an adjacent nucleus is often nearer than the partner outline.
    """
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

    m, used = _match(cfg, donor, roi, moving, fixed, tf, fixed_is_pheno,
                     moving_xy, fixed_xy, method=method, max_dist_um=max_dist_um,
                     min_iou=min_iou, pheno_base=pheno_base)
    if m.empty:
        crit = (f"overlapping by IoU >= {min_iou}" if used == "iou"
                else f"within {max_dist_um} um")
        return pd.DataFrame(), {
            "donor_id": donor, "roi": roi, "n_pairs": 0, "method": used,
            "note": (f"no cell pairs {crit} — either the registration is wrong or these "
                     f"are not the same section"),
        }

    out = pd.DataFrame({
        "donor_id": donor,
        "roi": roi,
        "pheno_cell_id": (fixed.df["cell_id"].to_numpy()[m["fixed_idx"]] if fixed_is_pheno
                          else moving.df["cell_id"].to_numpy()[m["moving_idx"]]),
        "xen_cell_id": (moving.df["cell_id"].to_numpy()[m["moving_idx"]] if fixed_is_pheno
                        else fixed.df["cell_id"].to_numpy()[m["fixed_idx"]]),
    })
    for col in ("distance_um", "iou"):
        if col in m.columns:
            out[col] = m[col].to_numpy()

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
        "method": used,
        "median_distance_um": (float(out["distance_um"].median())
                               if "distance_um" in out.columns else float("nan")),
        "median_iou": float(out["iou"].median()) if "iou" in out.columns else float("nan"),
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
    roi: Optional[str] = None,
    tissue: Optional[str] = None,
    max_dist_um: float = DEFAULT_MAX_DIST_UM,
) -> pd.DataFrame:
    """Stage entry point for same-slide mode. One call covers the whole selection."""
    require_same_slide(cfg)

    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)
    if donors:
        man = man[man["donor_id"].isin({str(d) for d in donors})]
    if man.empty:
        print("[sameslide] no paired donors", flush=True)
        return pd.DataFrame()

    rows, frames = [], []
    for _, r in man.iterrows():
        donor, roi_name = r["donor_id"], r["roi"]
        try:
            pairs, stats = pair_donor(cfg, donor, roi_name, max_dist_um=max_dist_um)
        except FileNotFoundError as exc:
            print(f"[sameslide] {donor}/{roi_name}: {exc}", flush=True)
            continue
        if "error" in stats:
            print(f"[sameslide] {stats['error']}", flush=True)
            continue
        rows.append(stats)
        if len(pairs):
            frames.append(build_paired_matrix(cfg, donor, pairs, roi_name))
        print(f"[sameslide] {donor}/{roi_name}: {stats['n_pairs']:,} cell pairs "
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
    ap.add_argument("--roi", default=None,
                    help="one ROI (default: every ROI of --tissue, or all tissues)")
    ap.add_argument("--tissue", default=None, help="restrict to one tissue's ROIs")
    ap.add_argument("--max-dist-um", type=float, default=DEFAULT_MAX_DIST_UM)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    stats = run_sameslide(cfg, args.donors, roi=args.roi, tissue=args.tissue,
                          max_dist_um=args.max_dist_um)
    if len(stats):
        print()
        print(stats.to_string(index=False))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
