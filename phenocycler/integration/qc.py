"""
S10 — QC gates and the integration report.

The job of this module is to decide, per donor, which levels of integration are believable,
and to make the decision visible rather than implicit.

Registration either worked or it did not, and the difference is not always obvious by eye.
A donor whose sections came from different blocks, or are too far apart, or whose tissue tore
during sectioning, will still produce *a* transform — every optimiser returns something. The
gates here are what stop that transform from silently propagating into matched islets,
grid comparisons and inferred cell links.

Three verdicts, and the graded response matters: a FAIL donor is not dropped from the study,
it is dropped from the *spatial* analyses. It still contributes donor-level composition and
niche comparison, both of which are registration-free. Losing a donor entirely because its
sections would not align would be throwing away good data for a bad reason.

    PASS   tissue Dice >= threshold AND islet RMSE <= threshold
           -> everything, including matched islets and inferred links
    WARN   one gate met, or metrics missing
           -> structure and grid results are produced but flagged
    FAIL   neither gate met
           -> registration-dependent steps skipped; donor-level and niche results stand

The islet-RMSE gate applies only to tissues that *have* a structural unit. A pancreatic
lymph node has no islets, so there is no landmark residual to compute and the metric is
undefined rather than missing. Grading on "one of two gates met" there would cap every
lymph-node donor at WARN forever and permanently withhold ``crossmodal`` from half the
cohort — for a reason that is an artefact of the gate design, not of the registration. So
for those tissues the verdict rests on tissue Dice alone, which is the metric that actually
exists. See ``phenocycler.integration.tissues``.
"""

from __future__ import annotations

import argparse
import json
from typing import Optional

import numpy as np
import pandas as pd

from ..config import PipelineConfig, load_config
from .manifest import load_manifest, paired_rows
from .register import registration_dir
from .tissues import DISPLAY_NAME, has_structures

PASS, WARN, FAIL = "PASS", "WARN", "FAIL"

#: What each verdict permits. Consumed by the orchestrator and by crossmodal.py.
ALLOWED_BY_VERDICT = {
    PASS: ("donor", "niche", "structure", "grid", "crossmodal"),
    WARN: ("donor", "niche", "structure", "grid"),
    FAIL: ("donor", "niche"),
}


def verdict_for(tissue_dice: float, islet_rmse_um: float, *,
                dice_min: float, rmse_max: float,
                require_structures: bool = True) -> tuple[str, list[str]]:
    """Grade one registration. Returns ``(verdict, reasons)``.

    ``require_structures=False`` for a tissue with no structural unit (lymph node): the
    landmark residual is undefined there, so it is dropped from the decision entirely rather
    than counted as a failed gate.
    """
    reasons = []
    dice_ok = np.isfinite(tissue_dice) and tissue_dice >= dice_min
    rmse_ok = np.isfinite(islet_rmse_um) and islet_rmse_um <= rmse_max

    if not np.isfinite(tissue_dice):
        reasons.append("tissue Dice not computed")
    elif not dice_ok:
        reasons.append(f"tissue Dice {tissue_dice:.3f} < {dice_min:.2f}")

    if not require_structures:
        # One gate, honestly labelled. Do not report the missing RMSE as a shortcoming of
        # this registration — no landmark exists in this tissue to measure it against.
        reasons.append("islet RMSE n/a — tissue has no structural unit; graded on Dice alone")
        return (PASS if dice_ok else FAIL), reasons

    if not np.isfinite(islet_rmse_um):
        reasons.append("islet RMSE not computed (too few structures to check)")
    elif not rmse_ok:
        reasons.append(f"islet RMSE {islet_rmse_um:.0f} um > {rmse_max:.0f} um")

    if dice_ok and rmse_ok:
        return PASS, reasons
    if dice_ok or rmse_ok:
        return WARN, reasons
    return FAIL, reasons


def collect_registration_qc(cfg: PipelineConfig, roi: Optional[str] = None,
                            tissue: Optional[str] = None) -> pd.DataFrame:
    """Read every stored ``qc.json`` and grade it.

    Scope follows ``paired_rows``: one ROI, one tissue's ROIs, or — with ``roi=None`` and no
    tissue — every paired section in the cohort, which is what a combined run grades.
    """
    man = paired_rows(load_manifest(cfg), roi=roi, tissue=tissue)

    rows = []
    for _, r in man.iterrows():
        donor, roi_name = r["donor_id"], r["roi"]
        roi_tissue = r.get("tissue", "") or cfg.tissue_for_roi(roi_name)
        graded_on_dice = not has_structures(roi_tissue) if roi_tissue else False
        path = registration_dir(cfg, donor, roi_name) / "qc.json"
        if not path.exists():
            rows.append({"donor_id": donor, "roi": roi_name, "tissue": roi_tissue,
                         "verdict": FAIL, "reasons": "no registration produced", "track": "",
                         "tissue_dice": np.nan, "islet_rmse_um": np.nan})
            continue
        d = json.loads(path.read_text())
        v, reasons = verdict_for(
            float(d.get("tissue_dice", np.nan)), float(d.get("rmse_um", np.nan)),
            dice_min=cfg.qc_tissue_dice_min, rmse_max=cfg.qc_islet_rmse_max_um,
            require_structures=not graded_on_dice)
        tf = d.get("transform", {})
        rows.append({
            "donor_id": donor, "roi": roi_name, "tissue": roi_tissue, "verdict": v,
            "reasons": "; ".join(reasons),
            "track": d.get("track", ""),
            "tissue_dice": d.get("tissue_dice", np.nan),
            "islet_dice": d.get("islet_dice", np.nan),
            "islet_rmse_um": d.get("rmse_um", np.nan),
            "n_inliers": d.get("n_inliers", 0),
            "mirrored": tf.get("mirrored", False),
            "rotation_deg": tf.get("rotation_deg", np.nan),
            "scale": tf.get("scale", np.nan),
            "max_disp_um": tf.get("max_displacement_um", 0.0),
            "stages": "|".join(d.get("stages", [])),
            "allows": ",".join(ALLOWED_BY_VERDICT[v]),
        })
    return pd.DataFrame(rows)


def _scope_label(cfg: PipelineConfig, roi: Optional[str], tissue: Optional[str]) -> str:
    """What this report covers, so a per-tissue run is never mistaken for the whole cohort."""
    if tissue:
        return f"{DISPLAY_NAME.get(tissue, tissue)} ({'/'.join(cfg.rois_for_tissue(tissue))})"
    if roi:
        return f"roi {roi} ({DISPLAY_NAME.get(cfg.tissue_for_roi(roi), '?')})"
    return "all tissues — " + ", ".join(DISPLAY_NAME.get(t, t) for t in cfg.tissue_list)


def sanity_flags(qc: pd.DataFrame, cfg: PipelineConfig) -> pd.DataFrame:
    """Flag transforms that are numerically fine but physically implausible.

    A registration can clear the Dice and RMSE gates while still being wrong in a way those
    metrics do not see. Both instruments report calibrated microns, so a scale far from 1 is a
    unit error rather than a biological observation; and a displacement field pinned at the
    cap means the non-rigid stage wanted to deform further than allowed, which usually means
    it was compensating for a bad affine rather than for real tissue distortion.
    """
    if qc.empty:
        return qc
    out = qc.copy()
    flags = []
    for _, r in out.iterrows():
        f = []
        scale = float(r.get("scale", np.nan))
        if np.isfinite(scale) and not (0.9 <= scale <= 1.11):
            f.append(f"scale {scale:.3f} far from 1 — check units, not biology")
        disp = float(r.get("max_disp_um", 0.0) or 0.0)
        if disp >= cfg.reg_max_disp_um * 0.99:
            f.append(f"displacement pinned at the {cfg.reg_max_disp_um:.0f} um cap")
        if bool(r.get("mirrored", False)):
            f.append("section mounted mirrored (expected sometimes; noted for the record)")
        flags.append("; ".join(f))
    out["sanity_flags"] = flags
    return out


def summarize(cfg: PipelineConfig, roi: Optional[str] = None,
              tissue: Optional[str] = None) -> tuple[pd.DataFrame, str]:
    """Grade every donor and render the report."""
    qc = sanity_flags(collect_registration_qc(cfg, roi, tissue), cfg)

    lines = [
        "PhenoCycler <-> Xenium integration QC",
        "=" * 68,
        f"mode           : {cfg.integration_mode}",
        f"fixed frame    : {cfg.fixed_modality}",
        f"scope          : {_scope_label(cfg, roi, tissue)}",
        f"gates          : tissue Dice >= {cfg.qc_tissue_dice_min:.2f}, "
        f"islet RMSE <= {cfg.qc_islet_rmse_max_um:.0f} um "
        f"(RMSE gate applies only where a structural unit exists)",
        "",
    ]
    if qc.empty:
        lines.append("no registered donors yet.")
        return qc, "\n".join(lines)

    counts = qc["verdict"].value_counts().to_dict()
    lines.append("verdicts: " + ", ".join(f"{k}={counts.get(k, 0)}"
                                          for k in (PASS, WARN, FAIL)))
    if qc["tissue"].nunique() > 1:
        by_tissue = qc.groupby("tissue")["verdict"].value_counts().unstack(fill_value=0)
        lines.append("")
        for t, row in by_tissue.iterrows():
            lines.append(f"  {DISPLAY_NAME.get(t, t):<22} " + ", ".join(
                f"{k}={int(row.get(k, 0))}" for k in (PASS, WARN, FAIL)))
    lines.append("")
    lines.append(f"  {'donor':<8} {'roi':<7} {'verdict':<7} {'track':<10} "
                 f"{'dice':>6} {'rmse_um':>8}  notes")
    for _, r in qc.sort_values(["verdict", "donor_id", "roi"]).iterrows():
        note = "; ".join(x for x in (r.get("reasons", ""), r.get("sanity_flags", "")) if x)
        lines.append(f"  {r['donor_id']:<8} {str(r['roi']):<7} {r['verdict']:<7} "
                     f"{str(r['track']):<10} {float(r['tissue_dice']):>6.3f} "
                     f"{float(r['islet_rmse_um']):>8.1f}  {note}")

    n_fail = counts.get(FAIL, 0)
    if n_fail:
        lines.append("")
        lines.append(f"{n_fail} donor(s) failed registration. They are excluded from islet "
                     f"matching, grid comparison and pseudo-cell linking, but STILL "
                     f"contribute donor-level composition and niche comparison — both of "
                     f"which need no registration.")
    return qc, "\n".join(lines)


def gather_results(cfg: PipelineConfig) -> dict[str, pd.DataFrame]:
    """Load whichever result tables exist, for the summary section of the report."""
    out = {}
    candidates = {
        "registration": cfg.integration_qc_dir / "registration_qc.csv",
        "match": cfg.integration_qc_dir / "match_qc.csv",
        "niche": cfg.integration_qc_dir / "niche_qc.csv",
        "donor": cfg.integration_qc_dir / "donor_qc.csv",
        "crossmodal": cfg.integration_qc_dir / "crossmodal_qc.csv",
        "sameslide": cfg.integration_qc_dir / "sameslide_qc.csv",
        "composition": cfg.integration_qc_dir / "composition.csv",
        "disease_trend": cfg.integration_qc_dir / "disease_trend.csv",
    }
    for name, path in candidates.items():
        if path.exists():
            try:
                out[name] = pd.read_csv(path)
            except Exception:  # noqa: BLE001
                pass
    return out


def concordance_summary(results: dict[str, pd.DataFrame]) -> str:
    """The scientific bottom line: did the two modalities agree, and above what null?"""
    lines = ["", "Cross-modal concordance", "-" * 68]

    d = results.get("donor")
    if d is not None and len(d) and "spearman_composition_ev" in d.columns:
        vals = pd.to_numeric(d["spearman_composition_ev"], errors="coerce").dropna()
        if len(vals):
            lines.append(f"donor-level composition (evidence-positive): "
                         f"median Spearman {vals.median():.3f} over {len(vals)} donor(s)")

    m = results.get("match")
    if m is not None and len(m):
        if "n_matched" in m.columns:
            lines.append(f"matched islets: {int(pd.to_numeric(m['n_matched'], errors='coerce').sum())} "
                         f"pairs over {len(m)} donor(s)")
        if "z" in m.columns:
            z = pd.to_numeric(m["z"], errors="coerce").dropna()
            if len(z):
                lines.append(f"  vs permuted null: median z = {z.median():.1f} "
                             f"({'above chance' if z.median() > 3 else 'NOT clearly above chance'})")
        for col, label in (("spearman_frac_beta", "per-islet beta fraction"),
                           ("insulitis_agreement", "insulitis grade agreement")):
            if col in m.columns:
                v = pd.to_numeric(m[col], errors="coerce").dropna()
                if len(v):
                    lines.append(f"  {label}: median {v.median():.3f}")

    n = results.get("niche")
    if n is not None and len(n) and "niche_median_cosine" in n.columns:
        v = pd.to_numeric(n["niche_median_cosine"], errors="coerce").dropna()
        if len(v):
            lines.append(f"niche composition cosine (registration-free): median {v.median():.3f}")

    c = results.get("crossmodal")
    if c is not None and len(c) and "agree_real" in c.columns:
        real = pd.to_numeric(c["agree_real"], errors="coerce").dropna()
        null = pd.to_numeric(c.get("agree_null", pd.Series(dtype=float)),
                             errors="coerce").dropna()
        if len(real):
            extra = f" vs null {null.median():.3f}" if len(null) else ""
            lines.append(f"inferred cell links: lineage agreement {real.median():.3f}{extra} "
                         f"[INFERENCE, not measurement]")

    s = results.get("sameslide")
    if s is not None and len(s) and "lineage_agreement" in s.columns:
        v = pd.to_numeric(s["lineage_agreement"], errors="coerce").dropna()
        if len(v):
            lines.append(f"same-slide paired cells: lineage agreement {v.median():.3f} "
                         f"[MEASUREMENT — same physical cells]")

    if len(lines) == 3:
        lines.append("(nothing computed yet)")
    return "\n".join(lines)


def run_qc(cfg: PipelineConfig, *, roi: Optional[str] = None,
           tissue: Optional[str] = None, write: bool = True) -> pd.DataFrame:
    """Stage entry point: grade registrations, write the gate table and the report."""
    qc, report = summarize(cfg, roi, tissue)
    results = gather_results(cfg)
    full = report + "\n" + concordance_summary(results)
    print(full, flush=True)

    if write:
        cfg.integration_qc_dir.mkdir(parents=True, exist_ok=True)
        if len(qc):
            qc.to_csv(cfg.integration_qc_dir / "registration_qc.csv", index=False)
        (cfg.integration_qc_dir / "integration_report.txt").write_text(full)
        print(f"\n[qc] wrote {cfg.integration_qc_dir/'integration_report.txt'}", flush=True)
    return qc


def passing_donors(cfg: PipelineConfig, level: str, roi: Optional[str] = None,
                   tissue: Optional[str] = None) -> list[str]:
    """Donors whose verdict permits ``level`` ('structure', 'grid', 'crossmodal', ...)."""
    qc = collect_registration_qc(cfg, roi, tissue)
    if qc.empty:
        return []
    ok = qc["verdict"].map(lambda v: level in ALLOWED_BY_VERDICT.get(v, ()))
    return sorted(set(qc.loc[ok, "donor_id"].tolist()))


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Integration QC gates and report")
    ap.add_argument("--config", default=None)
    ap.add_argument("--roi", default=None,
                    help="grade one ROI (default: every ROI of --tissue, or all tissues)")
    ap.add_argument("--tissue", default=None,
                    help="grade one tissue's ROIs, e.g. pancreatic_lymph_node")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    run_qc(cfg, roi=args.roi, tissue=args.tissue, write=not args.no_write)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
