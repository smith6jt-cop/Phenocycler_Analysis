"""Gate-2 reference-selection dossier for the RESTORE broad-lineage rebuild.

The locked screen (``restore_validation.evaluate_locked_cohort``) evaluates every ``(target,
reference)`` candidate pair **independently** and emits a per-donor candidate divisor -- it never
combines references.  Choosing *one* production reference per target is the pending Gate-2 human
decision.  This module turns the screen into the decision-support the maintainer needs, under the
policy fixed 2026-07-24:

  1. one predeclared reference per target (manuscript-faithful); the other references stay comparators;
  2. selection condition = the reference-positive negative control must be at least as large as the
     target-positive population -- a **pure proportion** ``reference_n / target_n >= 1`` with no
     absolute cell-count floor (images vary widely in total and per-marker counts).

It emits, from an existing screen (no NNMF rerun): a ``reference_matrix.csv`` (one row per candidate
``(target, reference)`` with the frequency check, divisor consistency, SBR, and a transparent
recommendation), one diagnostic figure per target, and a manifest.  It writes **no** divisor, no
normalized value, and no phenotype call -- it is a review artifact, not a production step.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .restore_validation import (
    BASELINE_PAIRS,
    CANDIDATE_PAIRS,
    METHOD_VERSION,
    PairValidationConfig,
)

# Targets whose single candidate reference is structural and predeclared; the immune targets carry
# several candidate references and are where the frequency rule does its work.
#
# Pan_Cytokeratin / Ker8_18 / EpCAM were added 2026-07-27. They have NO predeclared baseline -- they
# are candidates for promotion out of the epithelial comparator, not existing production targets whose
# baseline is being re-examined (see restore_validation.EPITHELIAL_TARGET_COMPARATOR_PAIRS and
# docs/ANATOMY_CHECK_2026-07-27.md). `recommend_references` handles the missing baseline explicitly.
_TARGET_ORDER = (
    "E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM",
    "CD31", "Vimentin", "B3TUBB", "CD3e", "CD68", "CD11b",
)
# Targets with no entry in BASELINE_PAIRS: nothing was predeclared for them, so "the baseline fails"
# is not a meaningful statement about them.
_NO_BASELINE = ("Pan_Cytokeratin", "Ker8_18", "EpCAM")

# Frozen Gate-2 maintainer decisions (2026-07-24) that OVERRIDE the deterministic "keep the passing
# baseline" recommendation. Any target not listed here uses its predeclared baseline reference where it
# passes the frequency rule.
#   CD11b -> EpCAM: the baseline E_cadherin passes only marginally (undersized on 30% of donors); EpCAM
#     is the stronger passing reference (median ratio 3.66, undersized on 15%) with identical divisor
#     consistency, so it is the production reference.
#   E_cadherin -> CD31: accepted DESPITE failing the pure-proportion rule (0.19x). Endothelium is the
#     only biologically exclusive negative control for epithelium in this panel; it is sparse but its
#     absolute control is large and stable. For a structural target the frequency rule is advisory, not
#     a gate -- the maintainer accepts CD31 with that caveat recorded.
ACCEPTED_REFERENCE_OVERRIDES: dict[str, str] = {
    "CD11b": "EpCAM",
    "E_cadherin": "CD31",
}


@dataclass(frozen=True)
class ReferenceDecision:
    """One target's recommended reference plus the disposition of every candidate reference."""

    target: str
    recommended_reference: str | None
    status: str  # "accepted" | "open"
    rationale: str


def _candidate_references() -> dict[str, list[str]]:
    refs: dict[str, list[str]] = {}
    for target, reference in CANDIDATE_PAIRS:
        refs.setdefault(target, [])
        if reference not in refs[target]:
            refs[target].append(reference)
    return refs


def _baseline_reference() -> dict[str, str]:
    return {target: reference for target, reference in BASELINE_PAIRS}


def build_reference_matrix(
    pair_evaluations: pd.DataFrame,
    *,
    config: PairValidationConfig | None = None,
    overrides: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Aggregate the per-donor screen into a per-``(target, reference)`` decision matrix.

    The frequency rule is applied ``in general`` -- a reference passes when its **median** ratio across
    donors is >= ``min_reference_to_target_ratio`` (a pure proportion). ``undersized_frac`` reports the
    per-donor fraction that individually fail, so a "passes in general but marginal" case (e.g. CD11b
    <- E_cadherin) is visible rather than hidden.  Divisor consistency is computed against the
    per-donor consensus of the rule-passing references, so a sparse outlier reference (CD3e for CD11b)
    shows up as both undersized AND a low-divisor outlier -- an independent cross-check.
    """
    config = config or PairValidationConfig()
    required = {
        "donor",
        "target",
        "reference",
        "reference_n",
        "target_n",
        "reference_to_target_ratio",
        "reference_undersized",
        "candidate_full_maximum",
        "reference_fold",
        "threshold_positive_n",
        "n_assigned",
        "double_high_fraction",
    }
    missing = required.difference(pair_evaluations.columns)
    if missing:
        raise ValueError(
            f"screen is missing columns {sorted(missing)} -- regenerate with METHOD_VERSION "
            f">= v8 (needs reference_to_target_ratio / reference_undersized)"
        )
    pe = pair_evaluations.copy()
    pe["called_frac"] = pe["threshold_positive_n"] / pe["n_assigned"].where(
        pe["n_assigned"] > 0
    )

    # Per-donor consensus divisor over the rule-passing references of each target (used only to flag
    # divisor outliers, never to combine references into a threshold).
    #
    # Fallback (2026-07-27): for a MAJORITY-POSITIVE target no reference can pass the pure-proportion
    # rule -- the target-positive population is simply larger than any non-target population in the
    # tissue. Measured on the 20-donor QC-clean screen, every epithelial target (E_cadherin,
    # Pan_Cytokeratin, Ker8_18, EpCAM) is undersized against every reference on every donor, ratios
    # 0.03-0.18. Restricting the consensus to passing references then left it empty, so
    # divisor_vs_consensus was NaN for exactly the targets whose reference choice most needs an
    # independent cross-check. Fall back to all of that target's references and say so in
    # `consensus_over_all_references`, rather than emitting nothing.
    undersized = pe["reference_undersized"].astype(bool)
    passing = pe[~undersized]
    consensus = (
        passing.groupby(["donor", "target"])["candidate_full_maximum"]
        .median()
        .rename("consensus_divisor")
        .reset_index()
    )
    fallback = (
        pe.groupby(["donor", "target"])["candidate_full_maximum"]
        .median()
        .rename("consensus_divisor_all")
        .reset_index()
    )
    pe = pe.merge(consensus, on=["donor", "target"], how="left")
    pe = pe.merge(fallback, on=["donor", "target"], how="left")
    pe["consensus_over_all_references"] = pe["consensus_divisor"].isna()
    pe["consensus_divisor"] = pe["consensus_divisor"].fillna(pe["consensus_divisor_all"])
    pe["divisor_vs_consensus"] = pe["candidate_full_maximum"] / pe["consensus_divisor"]

    rows = []
    for (target, reference), grp in pe.groupby(["target", "reference"]):
        ratio = grp["reference_to_target_ratio"].astype(float)
        double = grp["double_high_fraction"].astype(float)
        rows.append(
            {
                "target": target,
                "reference": reference,
                "n_donors": int(grp["donor"].nunique()),
                "ratio_median": float(ratio.median()),
                "ratio_min": float(ratio.min()),
                "ratio_max": float(ratio.max()),
                # Mutual exclusivity: the manuscript's premise, measured. The per-donor MAX is reported
                # beside the median because three pairs are clean at the median and catastrophic on a
                # single donor (INS <- Pan_CK max 56.6%, SST <- Pan_CK 46.1%, B3TUBB <- EpCAM 26.5%),
                # which a cohort median hides entirely.
                "double_positive_median": float(double.median()),
                "double_positive_max": float(double.max()),
                "double_positive_donors_over_limit": float(
                    (double > config.max_double_positive_fraction).mean()
                ),
                "undersized_frac": float(grp["reference_undersized"].astype(bool).mean()),
                "divisor_median": float(grp["candidate_full_maximum"].median()),
                "divisor_vs_consensus_median": (
                    float(grp["divisor_vs_consensus"].dropna().median())
                    if grp["divisor_vs_consensus"].notna().any()
                    else np.nan
                ),
                "reference_fold_median": float(grp["reference_fold"].median()),
                "called_frac_median": float(grp["called_frac"].median()),
                "consensus_over_all_references": bool(
                    grp["consensus_over_all_references"].all()
                ),
                "passes_frequency_rule": bool(
                    ratio.median() >= config.min_reference_to_target_ratio
                ),
                # HARD criterion, unlike the frequency rule. A reference that is not mutually exclusive
                # with its target cannot serve as a negative control at all -- the premise fails, so no
                # divisor derived from it means what the method says it means. The frequency rule is
                # advisory (a sparse-but-valid control is still a control); this one is not.
                "passes_exclusivity_rule": bool(
                    double.median() <= config.max_double_positive_fraction
                ),
                "is_baseline_reference": (target, reference) in set(BASELINE_PAIRS),
            }
        )
    matrix = pd.DataFrame(rows)

    decisions = recommend_references(matrix, overrides=overrides)
    # annotate each row with the per-target recommendation
    rec_ref = {d.target: d.recommended_reference for d in decisions}
    matrix["recommendation"] = matrix.apply(
        lambda r: _row_recommendation(r, rec_ref.get(r["target"])), axis=1
    )
    order = {t: i for i, t in enumerate(_TARGET_ORDER)}
    matrix = matrix.sort_values(
        by=["target", "ratio_median"],
        key=lambda s: s.map(order) if s.name == "target" else -s,
    ).reset_index(drop=True)
    return matrix


def _row_recommendation(row: pd.Series, recommended: str | None) -> str:
    # Mutual exclusivity is checked FIRST and outranks the recommendation: a non-exclusive reference is
    # inadmissible however it was chosen, so it must never render as "accepted".
    if not row["passes_exclusivity_rule"]:
        return "rejected (not exclusive)"
    # The maintainer-accepted reference then wins even if it fails the frequency rule
    # (E_cadherin <- CD31), so check the recommendation before that advisory rule.
    if row["reference"] == recommended:
        return "accepted"
    if not row["passes_frequency_rule"]:
        return "rejected"
    return "comparator"


def recommend_references(
    matrix: pd.DataFrame, *, overrides: dict[str, str] | None = None
) -> list[ReferenceDecision]:
    """Pick one reference per target.

    A frozen maintainer ``overrides`` entry wins outright (and may accept a reference that fails the
    frequency rule, with the caveat recorded). Otherwise: keep the predeclared baseline where it
    passes, else the highest-ratio passing reference; ``open`` when no candidate passes and there is no
    override.
    """
    overrides = overrides or {}
    baseline = _baseline_reference()
    decisions: list[ReferenceDecision] = []
    for target, grp in matrix.groupby("target"):
        # MUTUAL EXCLUSIVITY IS A HARD VETO and is applied BEFORE the override branch, because a
        # maintainer override cannot make a non-exclusive pair into a valid negative control -- the
        # manuscript's premise either holds or it does not. `CD11b <- EpCAM` was frozen by exactly such
        # an override on 2026-07-24 while carrying a 17.3% double-positive rate that nothing tested.
        exclusive = grp[grp["passes_exclusivity_rule"]]
        if exclusive.empty:
            worst = float(grp["double_positive_median"].min())
            decisions.append(
                ReferenceDecision(
                    target=target,
                    recommended_reference=None,
                    status="open",
                    rationale=(
                        f"NO admissible reference: every candidate violates mutual exclusivity "
                        f"(best median double-positive {worst:.1%}). The manuscript partitions a pair "
                        "into reference-positive / target-positive / negative-for-both, so a shared "
                        "population means no reference-positive cell can be assumed target-negative "
                        "and no divisor from it is interpretable. Defer this target or find an "
                        "exclusive reference; do NOT patch it in the lineage gate."
                    ),
                )
            )
            continue
        if target in overrides and overrides[target] not in set(exclusive["reference"]):
            row = grp[grp["reference"] == overrides[target]]
            dp = float(row.iloc[0]["double_positive_median"]) if not row.empty else float("nan")
            best = exclusive.sort_values("ratio_median", ascending=False).iloc[0]
            decisions.append(
                ReferenceDecision(
                    target=target,
                    recommended_reference=str(best["reference"]),
                    status="open",
                    rationale=(
                        f"OVERRIDE VOIDED: the frozen choice {overrides[target]} violates mutual "
                        f"exclusivity ({dp:.1%} double-positive). {best['reference']} is exclusive "
                        f"({float(best['double_positive_median']):.2%}) but fails the frequency rule "
                        f"(median ratio {float(best['ratio_median']):.2f}) -- the two criteria "
                        "disagree, so this is a maintainer decision, not a computation."
                        if not bool(best["passes_frequency_rule"])
                        else
                        f"OVERRIDE VOIDED: the frozen choice {overrides[target]} violates mutual "
                        f"exclusivity ({dp:.1%} double-positive); {best['reference']} is exclusive "
                        f"({float(best['double_positive_median']):.2%}) and passes the frequency rule."
                    ),
                )
            )
            continue
        grp = exclusive
        if target in overrides:
            chosen = overrides[target]
            chosen_row = grp[grp["reference"] == chosen]
            passes = (
                bool(chosen_row.iloc[0]["passes_frequency_rule"])
                if not chosen_row.empty
                else False
            )
            note = f"maintainer decision 2026-07-24: {chosen}"
            if not passes:
                ratio = (
                    float(chosen_row.iloc[0]["ratio_median"])
                    if not chosen_row.empty
                    else float("nan")
                )
                note += (
                    f"; accepted DESPITE reference<target proportion (median ratio {ratio:.2f}) -- "
                    "the frequency rule is advisory for this structural target, whose only "
                    "biologically exclusive reference (endothelium) is sparse but has a large "
                    "absolute control"
                )
            decisions.append(
                ReferenceDecision(
                    target=target,
                    recommended_reference=chosen,
                    status="accepted",
                    rationale=note,
                )
            )
            continue
        passing = grp[grp["passes_frequency_rule"]].sort_values(
            "ratio_median", ascending=False
        )
        base = baseline.get(target)
        base_row = grp[grp["reference"] == base]
        if passing.empty:
            best = float(grp["ratio_median"].max())
            note = (
                f"no candidate reference satisfies reference>=target frequency "
                f"(best median ratio {best:.2f}); reference is a Gate-2 design decision"
            )
            # EVERY candidate undersized on EVERY donor is a statement about the target, not about its
            # references: a marker positive on a large share of cells cannot have a larger negative
            # control unless some other marker is positive on more, and in pancreas none is. Saying
            # only "no reference is good enough" would send the reviewer looking for a better
            # reference that cannot exist. Keyed on the measured pattern rather than a cutoff on the
            # ratio, so it fires exactly when the rule is inadmissible rather than merely unmet.
            if float(grp["undersized_frac"].min()) >= 1.0:
                note += (
                    f". Every candidate is undersized on every donor, which for a majority-positive "
                    f"structural target is inevitable rather than a property of the references -- the "
                    f"pure-proportion rule has no admissible answer here. Compare "
                    f"reference_fold_median (arm separation) and the absolute control size instead, "
                    f"and note the same exception is already recorded for E_cadherin <- CD31"
                )
            decisions.append(
                ReferenceDecision(
                    target=target,
                    recommended_reference=None,
                    status="open",
                    rationale=note,
                )
            )
            continue
        if not base_row.empty and bool(base_row.iloc[0]["passes_frequency_rule"]):
            chosen = base
            marginal = base_row.iloc[0]["undersized_frac"] > 0.1
            stronger = passing.iloc[0]["reference"]
            note = (
                f"baseline reference passes (median ratio "
                f"{base_row.iloc[0]['ratio_median']:.2f})"
            )
            if marginal and stronger != base:
                note += (
                    f"; MARGINAL -- undersized on "
                    f"{base_row.iloc[0]['undersized_frac']*100:.0f}% of donors, "
                    f"stronger alternative {stronger} "
                    f"({passing.iloc[0]['ratio_median']:.2f})"
                )
        else:
            chosen = passing.iloc[0]["reference"]
            if base is None:
                # A promotion candidate, not a production target with a failing baseline. Say so, or
                # the dossier reads as though something was tried and rejected.
                note = (
                    f"no predeclared baseline (promotion candidate); strongest passing reference is "
                    f"{chosen} (median ratio {passing.iloc[0]['ratio_median']:.2f}). Promotion to a "
                    "production target still requires its own Gate-1->3 validation."
                )
            else:
                note = (
                    f"baseline {base} fails the frequency rule; "
                    f"strongest passing reference is {chosen} "
                    f"(median ratio {passing.iloc[0]['ratio_median']:.2f})"
                )
        decisions.append(
            ReferenceDecision(
                target=target,
                recommended_reference=chosen,
                status="accepted",
                rationale=note,
            )
        )
    return decisions


def plot_target_diagnostic(
    pair_evaluations: pd.DataFrame,
    target: str,
    matrix: pd.DataFrame,
    out_png: Path,
    *,
    config: PairValidationConfig | None = None,
) -> None:
    """Three aligned panels for one target: reference-frequency ratio, divisor consistency across
    references, and the fraction of cells each reference calls positive."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config = config or PairValidationConfig()
    pe = pair_evaluations[pair_evaluations["target"] == target].copy()
    pe["called_frac"] = pe["threshold_positive_n"] / pe["n_assigned"].where(
        pe["n_assigned"] > 0
    )
    refs = [r for r in _candidate_references().get(target, []) if r in set(pe["reference"])]
    tmatrix = matrix[matrix["target"] == target].set_index("reference")
    rec = tmatrix.index[tmatrix["recommendation"] == "accepted"].tolist()
    recommended = rec[0] if rec else None

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
    fig.suptitle(
        f"Gate-2 reference selection — target {target}   "
        f"(recommended reference: {recommended or 'OPEN'})",
        fontsize=16,
        fontweight="bold",
    )
    xs = np.arange(len(refs))

    def _color(r: str) -> str:
        rr = tmatrix.loc[r, "recommendation"] if r in tmatrix.index else "comparator"
        return {"accepted": "#228833", "rejected": "#CC3311", "rejected (not exclusive)": "#882255",
         "comparator": "#4477AA"}.get(
            rr, "#4477AA"
        )

    # Panel 1: reference/target frequency ratio (the rule) -- per-donor points + median.
    ax = axes[0]
    for i, r in enumerate(refs):
        vals = pe[pe["reference"] == r]["reference_to_target_ratio"].astype(float)
        ax.scatter(
            np.full(len(vals), i) + np.linspace(-0.15, 0.15, len(vals)),
            vals,
            s=18,
            color=_color(r),
            alpha=0.6,
        )
        ax.hlines(vals.median(), i - 0.3, i + 0.3, color=_color(r), lw=3)
    ax.axhline(
        config.min_reference_to_target_ratio,
        color="black",
        ls="--",
        lw=1.5,
        label=f"rule: ratio ≥ {config.min_reference_to_target_ratio:g}",
    )
    ax.set_yscale("log")
    ax.set_xticks(xs)
    ax.set_xticklabels(refs, rotation=40, ha="right", fontsize=11)
    ax.set_ylabel("reference_n / target_n  (per donor)", fontsize=12)
    ax.set_title("Reference ≥ target frequency (pure proportion)", fontsize=13)
    ax.legend(fontsize=10)

    # Panel 2: divisor by reference relative to per-donor consensus (outlier detection).
    ax = axes[1]
    passing_refs = [
        r for r in refs if bool(tmatrix.loc[r, "passes_frequency_rule"])
    ] if len(tmatrix) else refs
    if not passing_refs:  # no reference passes (E_cadherin): consensus falls back to all references
        passing_refs = refs
    for i, r in enumerate(refs):
        sub = pe[pe["reference"] == r]
        cons = (
            pe[pe["reference"].isin(passing_refs)]
            .groupby("donor")["candidate_full_maximum"]
            .median()
        )
        rel = sub.set_index("donor")["candidate_full_maximum"] / cons
        rel = rel.dropna()
        ax.scatter(
            np.full(len(rel), i) + np.linspace(-0.15, 0.15, len(rel)),
            rel,
            s=18,
            color=_color(r),
            alpha=0.6,
        )
        ax.hlines(rel.median(), i - 0.3, i + 0.3, color=_color(r), lw=3)
    ax.axhline(1.0, color="black", ls="--", lw=1.5, label="passing-reference consensus")
    ax.set_xticks(xs)
    ax.set_xticklabels(refs, rotation=40, ha="right", fontsize=11)
    ax.set_ylabel("divisor ÷ per-donor consensus", fontsize=12)
    ax.set_title("Divisor consistency across references", fontsize=13)
    ax.legend(fontsize=10)

    # Panel 3: fraction of cells called positive per reference (the overcall, made visible).
    ax = axes[2]
    for i, r in enumerate(refs):
        vals = pe[pe["reference"] == r]["called_frac"].astype(float) * 100
        ax.scatter(
            np.full(len(vals), i) + np.linspace(-0.15, 0.15, len(vals)),
            vals,
            s=18,
            color=_color(r),
            alpha=0.6,
        )
        ax.hlines(vals.median(), i - 0.3, i + 0.3, color=_color(r), lw=3)
    ax.set_xticks(xs)
    ax.set_xticklabels(refs, rotation=40, ha="right", fontsize=11)
    ax.set_ylabel("% of assigned cells called positive", fontsize=12)
    ax.set_title("Called fraction per reference", fontsize=13)

    handles = [
        plt.Line2D([0], [0], marker="o", ls="", color="#228833", label="accepted"),
        plt.Line2D([0], [0], marker="o", ls="", color="#4477AA", label="comparator (passes rule)"),
        plt.Line2D([0], [0], marker="o", ls="", color="#CC3311", label="rejected (fails rule)"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=11, frameon=False)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def build_bundle(
    results_dir: Path,
    output_dir: Path,
    *,
    config: PairValidationConfig | None = None,
    overrides: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Write the reference-selection dossier (matrix + per-target figures + manifest) atomically."""
    config = config or PairValidationConfig()
    if overrides is None:
        overrides = dict(ACCEPTED_REFERENCE_OVERRIDES)
    results_dir = Path(results_dir)
    pe = pd.read_csv(results_dir / "pair_evaluations.csv")
    matrix = build_reference_matrix(pe, config=config, overrides=overrides)
    decisions = recommend_references(matrix, overrides=overrides)

    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    (temp_dir / "figures").mkdir(parents=True)
    try:
        matrix.to_csv(temp_dir / "reference_matrix.csv", index=False)
        # the frozen Gate-2 exit artifact: one row per target with the accepted reference
        pd.DataFrame(
            [
                {
                    "target": d.target,
                    "accepted_reference": d.recommended_reference,
                    "status": d.status,
                    "rationale": d.rationale,
                }
                for d in decisions
            ]
        ).to_csv(temp_dir / "frozen_reference_matrix.csv", index=False)
        targets = [t for t in _TARGET_ORDER if t in set(matrix["target"])]
        for target in targets:
            plot_target_diagnostic(
                pe, target, matrix, temp_dir / "figures" / f"{target}.png", config=config
            )
        manifest = {
            "method_version": METHOD_VERSION,
            "screen": str(results_dir),
            "policy": (
                "one predeclared reference per target; selection condition = pure proportion "
                "reference_n/target_n >= min_reference_to_target_ratio (no absolute floor)"
            ),
            "min_reference_to_target_ratio": config.min_reference_to_target_ratio,
            "decisions": [
                {
                    "target": d.target,
                    "recommended_reference": d.recommended_reference,
                    "status": d.status,
                    "rationale": d.rationale,
                }
                for d in decisions
            ],
        }
        (temp_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
        _write_summary(temp_dir / "SUMMARY.md", matrix, decisions, config)
        temp_dir.rename(output_dir)
    finally:
        if temp_dir.exists():
            import shutil

            shutil.rmtree(temp_dir, ignore_errors=True)
    return matrix


def _write_summary(
    path: Path,
    matrix: pd.DataFrame,
    decisions: list[ReferenceDecision],
    config: PairValidationConfig,
) -> None:
    lines = [
        "# Gate-2 reference-selection dossier",
        "",
        f"Method: `{METHOD_VERSION}`. Policy: one reference per target; selection condition = "
        f"pure proportion `reference_n / target_n >= {config.min_reference_to_target_ratio:g}` "
        "(no absolute cell-count floor).",
        "",
        "## Recommended reference per target",
        "",
        "| target | recommended | status | rationale |",
        "|---|---|---|---|",
    ]
    for d in decisions:
        lines.append(
            f"| {d.target} | {d.recommended_reference or '—'} | {d.status} | {d.rationale} |"
        )
    lines += [
        "",
        "## Rejected candidate references (fail the frequency rule)",
        "",
        "| target | reference | median ratio | undersized on % donors |",
        "|---|---|---|---|",
    ]
    for _, r in matrix[matrix["recommendation"] == "rejected"].iterrows():
        lines.append(
            f"| {r['target']} | {r['reference']} | {r['ratio_median']:.2f} | "
            f"{r['undersized_frac']*100:.0f}% |"
        )
    path.write_text("\n".join(lines) + "\n")


def _main(argv: list[str] | None = None) -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Build the Gate-2 reference-selection dossier from a locked screen. Emits no divisor, "
            "no normalized value, and no phenotype call."
        )
    )
    parser.add_argument("--results", required=True, help="a locked screen output directory")
    parser.add_argument("--output", required=True, help="dossier output directory")
    parser.add_argument(
        "--min-reference-to-target-ratio",
        type=float,
        default=PairValidationConfig().min_reference_to_target_ratio,
    )
    args = parser.parse_args(argv)
    config = PairValidationConfig(
        min_reference_to_target_ratio=args.min_reference_to_target_ratio
    )
    matrix = build_bundle(Path(args.results), Path(args.output), config=config)
    print(matrix.to_string(index=False))


if __name__ == "__main__":
    _main()
