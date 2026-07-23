#!/usr/bin/env python3
"""Create manuscript-style figures for the RESTORE pair-validation pilot."""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path
from typing import Sequence

import matplotlib

if not matplotlib.get_backend().startswith("module://"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.ndimage import gaussian_filter1d  # noqa: E402

from .restore_validation import (
    BASELINE_PAIRS,
    METHOD_VERSION,
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    evaluate_locked_pair,
    load_expression_reviews,
    load_pair_reviews,
    load_validation_sample,
    locked_pair_metrics,
)


DONOR_ORDER = ("6539", "6476")
FEATURE_METHOD = "nnmf_membrane_cell"
TARGET_COLOR = "#D55E00"
REFERENCE_COLOR = "#0072B2"
INVALID_COLOR = "#CC79A7"
OTHER_COLOR = "#777777"
MAX_COLOR = "#009E73"
FULL_MAX_COLOR = "#005A42"
SENSITIVITY_COLOR = "#7A9A01"
STATUS_COLORS = {
    PairState.VALID.value: "#228833",
    PairState.NO_TARGET_EXPRESSION_CONFIRMED.value: "#6A3D9A",
    PairState.LOW_SIGNAL_UNRESOLVED.value: "#9467BD",
    PairState.NO_REFERENCE_POPULATION.value: "#CC8800",
    PairState.INVALID_PAIR.value: "#CC3311",
    PairState.MODEL_UNSTABLE.value: "#E69F00",
    PairState.PAIR_REVIEW_PENDING.value: "#56B4E9",
    PairState.TECHNICAL_FAILURE.value: "#333333",
}
STATUS_LABELS = {
    PairState.VALID.value: "VALID",
    PairState.NO_TARGET_EXPRESSION_CONFIRMED.value: "NO TARGET",
    PairState.LOW_SIGNAL_UNRESOLVED.value: "LOW / UNRESOLVED",
    PairState.NO_REFERENCE_POPULATION.value: "NO REFERENCE",
    PairState.INVALID_PAIR.value: "INVALID PAIR",
    PairState.MODEL_UNSTABLE.value: "UNSTABLE",
    PairState.PAIR_REVIEW_PENDING.value: "PAIR REVIEW PENDING",
    PairState.TECHNICAL_FAILURE.value: "TECHNICAL FAILURE",
}
DEFAULT_CONFIG = PairValidationConfig()
MIN_TARGET_ARM_FOLD = DEFAULT_CONFIG.min_target_arm_fold


def marker_label(marker: str) -> str:
    return "E-cadherin" if marker == "E_cadherin" else marker


def pair_label(target: str, reference: str) -> str:
    return f"{marker_label(target)} \u2190 {marker_label(reference)}"


def _display_indices(
    labels: np.ndarray, cap_per_group: int, seed: int
) -> np.ndarray:
    """Return an equal-N component sample while preserving actual counts elsewhere."""
    labels = np.asarray(labels)
    groups = np.unique(labels)
    group_indices = [np.flatnonzero(labels == group) for group in groups]
    sample_n = min(cap_per_group, *(len(indices) for indices in group_indices))
    if sample_n < 1:
        raise ValueError("every displayed component must contain at least one cell")
    rng = np.random.default_rng(seed)
    selected = [
        (
            indices
            if len(indices) == sample_n
            else rng.choice(indices, size=sample_n, replace=False)
        )
        for indices in group_indices
    ]
    return np.sort(np.concatenate(selected))


def _density(values: np.ndarray, limits: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    histogram, edges = np.histogram(values, bins=100, range=limits, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, gaussian_filter1d(histogram, sigma=1.25)


def _manuscript_coordinates(
    target_raw: np.ndarray,
    reference_raw: np.ndarray,
    *,
    target_bounds: tuple[float, float],
    reference_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Return robust min-max-scaled linear coordinates in Figure 7 orientation."""
    target_lower, target_upper = target_bounds
    reference_lower, reference_upper = reference_bounds
    target = np.clip(
        (np.asarray(target_raw, dtype=float) - target_lower)
        / (target_upper - target_lower),
        0.0,
        1.0,
    )
    reference = np.clip(
        (np.asarray(reference_raw, dtype=float) - reference_lower)
        / (reference_upper - reference_lower),
        0.0,
        1.0,
    )
    return reference, target


def _prepare_panel(
    donor_df: pd.DataFrame,
    donor: str,
    target: str,
    reference: str,
    *,
    config: PairValidationConfig,
    expression_review: ExpressionReview = ExpressionReview.UNREVIEWED,
    pair_review: PairReview = PairReview.UNREVIEWED,
) -> dict:
    evaluation = evaluate_locked_pair(
        donor_df,
        donor,
        target,
        reference,
        config=config,
        expression_review=expression_review,
        pair_review=pair_review,
    )
    plot_mask = evaluation["qc_retained"]
    target_plot_bounds = evaluation["display_target_bounds"]
    reference_plot_bounds = evaluation["display_reference_bounds"]
    x, y = _manuscript_coordinates(
        evaluation["target_redsea"][plot_mask],
        evaluation["reference_redsea"][plot_mask],
        target_bounds=target_plot_bounds,
        reference_bounds=reference_plot_bounds,
    )
    return {
        **evaluation,
        "status": evaluation["state"],
        "x": x,
        "y": y,
        "labels": evaluation["labels"][plot_mask],
        "target_plot_bounds": target_plot_bounds,
        "reference_plot_bounds": reference_plot_bounds,
    }


def _group_style(panel: dict, group: int) -> tuple[str, str]:
    if panel["distinct"]:
        if group == panel["reference_group"]:
            return REFERENCE_COLOR, "Reference-positive / target-negative"
        if group == panel["target_group"]:
            return TARGET_COLOR, "Target-positive arm"
    if group == panel["reference_group"] == panel["target_group"]:
        return INVALID_COLOR, "Same component dominates both markers"
    return OTHER_COLOR, f"Other NNMF component {group}"


def _scale_plot_value(value: float, bounds: tuple[float, float]) -> float:
    lower, upper = bounds
    return float(np.clip((value - lower) / (upper - lower), 0.0, 1.0))


def _panel_limits(panels: Sequence[dict]) -> tuple[tuple[float, float], tuple[float, float]]:
    if not panels:
        raise ValueError("at least one panel is required")
    return (0.0, 1.02), (0.0, 1.02)


def _visible_scaled_interval(
    lower: float, upper: float, bounds: tuple[float, float]
) -> tuple[float, float] | None:
    """Return a visible scaled interval, or None when it collapses after clipping."""
    scaled_lower = _scale_plot_value(lower, bounds)
    scaled_upper = _scale_plot_value(upper, bounds)
    if scaled_upper - scaled_lower <= np.finfo(float).eps:
        return None
    return scaled_lower, scaled_upper


def _draw_thresholds(ax: plt.Axes, panel: dict) -> None:
    if not panel["distinct"]:
        return
    target_stats = panel["target_background"]
    reference_stats = panel["reference_background"]
    target_full = panel["target_full_stats"]
    reference_full = panel["reference_full_stats"]
    if any(stats is None for stats in (target_stats, reference_stats, target_full, reference_full)):
        return
    reference_interval = _visible_scaled_interval(
        reference_stats["maximum_q05"],
        reference_stats["maximum_q95"],
        panel["reference_plot_bounds"],
    )
    target_interval = _visible_scaled_interval(
        target_stats["maximum_q05"],
        target_stats["maximum_q95"],
        panel["target_plot_bounds"],
    )
    if reference_interval is not None:
        ax.axvspan(
            *reference_interval,
            color=MAX_COLOR,
            alpha=0.12,
            zorder=4,
        )
    if target_interval is not None:
        ax.axhspan(
            *target_interval,
            color=MAX_COLOR,
            alpha=0.12,
            zorder=4,
        )
    ax.axvline(
        _scale_plot_value(
            reference_stats["maximum"], panel["reference_plot_bounds"]
        ),
        color=MAX_COLOR,
        lw=1.6,
        ls="--",
        zorder=5,
    )
    ax.axhline(
        _scale_plot_value(
            target_stats["maximum"], panel["target_plot_bounds"]
        ),
        color=MAX_COLOR,
        lw=1.6,
        ls="--",
        zorder=5,
    )
    ax.axvline(
        _scale_plot_value(
            reference_full["maximum"], panel["reference_plot_bounds"]
        ),
        color=FULL_MAX_COLOR,
        lw=1.8,
        ls="--",
        zorder=6,
    )
    ax.axhline(
        _scale_plot_value(target_full["maximum"], panel["target_plot_bounds"]),
        color=FULL_MAX_COLOR,
        lw=1.8,
        ls="-",
        zorder=6,
    )
    ax.axvline(
        _scale_plot_value(
            reference_full["mean_plus_3sd"], panel["reference_plot_bounds"]
        ),
        color=SENSITIVITY_COLOR,
        lw=1.25,
        ls=":",
        zorder=5,
    )
    ax.axhline(
        _scale_plot_value(
            target_full["mean_plus_3sd"], panel["target_plot_bounds"]
        ),
        color=SENSITIVITY_COLOR,
        lw=1.25,
        ls=":",
        zorder=5,
    )


def _draw_scatter(
    ax: plt.Axes,
    panel: dict,
    limits: tuple[tuple[float, float], tuple[float, float]],
    *,
    cap_per_group: int,
    seed: int,
) -> None:
    selected = _display_indices(panel["labels"], cap_per_group, seed)
    for group in np.unique(panel["labels"]):
        group_idx = selected[panel["labels"][selected] == group]
        color, _ = _group_style(panel, int(group))
        ax.scatter(
            panel["x"][group_idx],
            panel["y"][group_idx],
            s=3.0,
            color=color,
            alpha=0.28,
            edgecolors="none",
            rasterized=True,
        )
    _draw_thresholds(ax, panel)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.grid(color="#DDDDDD", lw=0.45, alpha=0.65)
    ax.tick_params(labelsize=8)
    status = panel["status"]
    ax.text(
        0.98,
        0.98,
        STATUS_LABELS[status],
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=7,
        fontweight="bold",
        color="white",
        bbox=dict(boxstyle="round,pad=0.22", fc=STATUS_COLORS[status], ec="none"),
    )


def _draw_manuscript_panel(
    fig: plt.Figure,
    subplot_spec,
    panel: dict,
    limits: tuple[tuple[float, float], tuple[float, float]],
    *,
    seed: int,
) -> plt.Axes:
    grid = subplot_spec.subgridspec(
        2,
        2,
        width_ratios=(1.25, 5.0),
        height_ratios=(5.0, 1.25),
        wspace=0.04,
        hspace=0.04,
    )
    main = fig.add_subplot(grid[0, 1])
    left = fig.add_subplot(grid[0, 0], sharey=main)
    bottom = fig.add_subplot(grid[1, 1], sharex=main)
    fig.add_subplot(grid[1, 0]).axis("off")

    _draw_scatter(main, panel, limits, cap_per_group=4_000, seed=seed)
    for group in np.unique(panel["labels"]):
        color, _ = _group_style(panel, int(group))
        group_mask = panel["labels"] == group
        x_centers, x_density = _density(panel["x"][group_mask], limits[0])
        y_centers, y_density = _density(panel["y"][group_mask], limits[1])
        bottom.plot(x_centers, x_density, color=color, lw=1.35)
        bottom.fill_between(x_centers, 0, x_density, color=color, alpha=0.12)
        left.plot(y_density, y_centers, color=color, lw=1.35)
        left.fill_betweenx(y_centers, 0, y_density, color=color, alpha=0.12)

    left.invert_xaxis()
    left.axis("off")
    bottom.axis("off")
    main.set_title(
        f"Donor {panel['donor']}",
        fontsize=11,
        fontweight="bold",
        color=STATUS_COLORS[panel["status"]],
    )
    main.text(
        0.02,
        0.98,
        (
            f"groups R/T={panel['reference_n']}/{panel['target_n']}\n"
            f"fit/arm={panel['fit_per_arm']}, threshold n={panel['threshold_sample_n']}\n"
            f"target arm fold={panel['target_fold']:.2f}\n"
            f"seed Jaccard min={panel['min_control_jaccard']:.2f}\n"
            f"retained={100 * panel['qc_retained_fraction']:.1f}%\n"
            f"double-high={100 * panel['double_high_fraction']:.1f}%"
        ),
        transform=main.transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#AAAAAA", alpha=0.9),
    )
    main.set_xlabel(
        f"{marker_label(panel['reference'])} robust min-max REDSEA MFI", fontsize=9
    )
    main.set_ylabel(
        f"{marker_label(panel['target'])} robust min-max REDSEA MFI", fontsize=9
    )
    return main


def _legend_handles() -> list:
    return [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=REFERENCE_COLOR,
               markersize=6, label="Reference-positive / target-negative component"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=TARGET_COLOR,
               markersize=6, label="Target-positive component"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=INVALID_COLOR,
               markersize=6, label="No distinct directional arms"),
        Line2D([0], [0], color=FULL_MAX_COLOR, lw=1.8, ls="--",
               label="Step 1 vertical: reference maximum in target anchor"),
        Line2D([0], [0], color=FULL_MAX_COLOR, lw=1.8, ls="-",
               label="Step 2 horizontal: target divisor, not retained-target exclusion"),
        Line2D([0], [0], color=MAX_COLOR, lw=1.8, ls="--",
               label="Median equal-N sampled maximum (axis-specific; QC only)"),
        Patch(facecolor=MAX_COLOR, edgecolor="none", alpha=0.12,
              label="5-95% interval; absent when zero-width after clipping"),
        Line2D([0], [0], color=SENSITIVITY_COLOR, lw=1.5, ls=":",
               label="mean + 3SD (released-code sensitivity only)"),
    ]


def _add_bounded_footer(
    fig: plt.Figure,
    axes: Sequence[plt.Axes],
    text: str,
    *,
    bottom: float,
    height: float,
    fontsize: float,
) -> None:
    """Add wrapped explanatory text constrained to the subplot lateral span."""
    if not axes:
        raise ValueError("at least one subplot is required for a bounded footer")
    fig.canvas.draw()
    left = min(axis.get_position().x0 for axis in axes)
    right = max(axis.get_position().x1 for axis in axes)
    width_inches = fig.get_size_inches()[0] * (right - left)
    wrapped = textwrap.fill(text, width=max(50, int(width_inches * 13)))
    footer = fig.add_axes([left, bottom, right - left, height])
    footer.axis("off")
    footer.text(
        0.5,
        0.5,
        wrapped,
        ha="center",
        va="center",
        fontsize=fontsize,
        transform=footer.transAxes,
        clip_on=True,
    )


def create_pair_cloud_figures(
    donor_data: dict[str, pd.DataFrame],
    output_dir: Path,
    *,
    config: PairValidationConfig,
    expression_reviews: dict[tuple[str, str], ExpressionReview],
    pair_reviews: dict[tuple[str, str], PairReview],
    display_seed: int,
) -> tuple[dict[tuple[str, str], dict[str, dict]], pd.DataFrame]:
    """Create one Figure-7-style page per pair and return the recomputed panels."""
    pair_dir = output_dir / "pair_clouds"
    pair_dir.mkdir(parents=True, exist_ok=True)
    all_pdf = output_dir / "figureS1_all_pair_clouds.pdf"
    panels_by_pair: dict[tuple[str, str], dict[str, dict]] = {}
    metrics = []
    with PdfPages(all_pdf) as pdf:
        for target, reference in BASELINE_PAIRS:
            panels = {}
            for donor, donor_df in donor_data.items():
                panel = _prepare_panel(
                    donor_df,
                    donor,
                    target,
                    reference,
                    config=config,
                    expression_review=expression_reviews.get(
                        (donor, target), ExpressionReview.UNREVIEWED
                    ),
                    pair_review=pair_reviews.get(
                        (target, reference), PairReview.UNREVIEWED
                    ),
                )
                panels[donor] = panel
                metrics.append(locked_pair_metrics(panel))
            panels_by_pair[(target, reference)] = panels
            limits = _panel_limits(list(panels.values()))
            fig = plt.figure(figsize=(15.5, 6.8))
            outer = fig.add_gridspec(1, len(DONOR_ORDER), wspace=0.25)
            main_axes = []
            for column, donor in enumerate(DONOR_ORDER):
                main_axes.append(
                    _draw_manuscript_panel(
                        fig,
                        outer[0, column],
                        panels[donor],
                        limits,
                        seed=display_seed + column,
                    )
                )
            fig.suptitle(
                f"{pair_label(target, reference)}: locked RESTORE pair evaluation",
                fontsize=16,
                fontweight="bold",
                y=0.995,
            )
            footer_text = (
                "NNMF uses equal-arm training after robust per-feature min-max scaling "
                "(0.5-99.5%); axes show scaled REDSEA whole-cell MFI without a log transform. "
                "Points are equal-N component samples, while annotations report actual "
                "group counts across all assigned cells. Solid lines are the manuscript "
                "full-control maximum candidates. Green intervals are fixed-N QC and are "
                "calculated independently for x and y; "
                "a band is absent when its interval is zero-width because the arm is "
                "sampled in full on every repeat, or when it collapses at a clipped 0/1 "
                "display boundary. A candidate is accepted only when the panel state is VALID."
            )
            fig.subplots_adjust(
                top=0.88, bottom=0.26, left=0.04, right=0.99
            )
            fig.legend(
                handles=_legend_handles(),
                loc="lower center",
                bbox_to_anchor=(0.5, 0.105),
                ncol=3,
                frameon=False,
                fontsize=7.5,
                title="NNMF arms and background diagnostics",
                title_fontsize=8,
            )
            _add_bounded_footer(
                fig,
                main_axes,
                footer_text,
                bottom=0.015,
                height=0.07,
                fontsize=8,
            )
            stem = f"{target}_from_{reference}".lower()
            fig.savefig(pair_dir / f"{stem}.png", dpi=220, bbox_inches="tight")
            pdf.savefig(fig, dpi=220, bbox_inches="tight")
            plt.close(fig)
    return panels_by_pair, pd.DataFrame(metrics)


def create_example_figure(
    panels_by_pair: dict[tuple[str, str], dict[str, dict]], output_dir: Path, *, seed: int
) -> None:
    # CD20 <- E_cadherin was the immune exemplar here until CD20 was excluded from RESTORE in both
    # roles (2026-07-23); CD68 <- EpCAM is the shortlisted myeloid pair that replaces it.
    selected_pairs = (
        ("E_cadherin", "CD31"),
        ("B3TUBB", "EpCAM"),
        ("CD3e", "E_cadherin"),
        ("CD68", "EpCAM"),
    )
    fig, axes = plt.subplots(
        len(selected_pairs), len(DONOR_ORDER), figsize=(14.5, 16.5), squeeze=False
    )
    for row, pair in enumerate(selected_pairs):
        panels = panels_by_pair[pair]
        limits = _panel_limits(list(panels.values()))
        for column, donor in enumerate(DONOR_ORDER):
            ax = axes[row, column]
            panel = panels[donor]
            _draw_scatter(
                ax,
                panel,
                limits,
                cap_per_group=3_000,
                seed=seed + row * 10 + column,
            )
            if row == 0:
                ax.set_title(f"Donor {donor}", fontsize=12, fontweight="bold")
            if column == 0:
                ax.set_ylabel(
                    f"{pair_label(*pair)}\n{marker_label(pair[0])} robust min-max REDSEA MFI",
                    fontsize=10,
                )
            else:
                ax.set_ylabel(
                    f"{marker_label(pair[0])} robust min-max REDSEA MFI", fontsize=9
                )
            ax.set_xlabel(
                f"{marker_label(pair[1])} robust min-max REDSEA MFI", fontsize=9
            )
            ax.text(
                0.02,
                0.98,
                (
                    f"groups R/T={panel['reference_n']}/{panel['target_n']}\n"
                    f"target fold={panel['target_fold']:.2f}\n"
                    f"retained={100 * panel['qc_retained_fraction']:.1f}%\n"
                    f"double-high={100 * panel['double_high_fraction']:.1f}%"
                ),
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=7.3,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#AAAAAA", alpha=0.88),
            )
    fig.suptitle(
        "RESTORE pair-validation pilot: locked REDSEA divisor evaluation",
        fontsize=17,
        fontweight="bold",
        y=0.995,
    )
    fig.legend(
        handles=_legend_handles(),
        loc="lower center",
        bbox_to_anchor=(0.5, 0.07),
        ncol=3,
        frameon=False,
        fontsize=8,
        title="NNMF arms and background diagnostics",
        title_fontsize=8.5,
    )
    footer_text = (
        "Axes show robust min-max REDSEA MFI without log transformation (reference on x, "
        "target on y). Points are equal-N component samples; annotations report actual "
        "assigned group counts. Solid lines are full-control maximum divisor candidates. "
        "Green fixed-N 5-95% intervals are calculated independently for x and y, "
        "so a band is absent when it has zero visible width after full-arm sampling or "
        "0/1 clipping. LOW / UNRESOLVED marks target-arm separation below "
        f"{MIN_TARGET_ARM_FOLD:.1f}-fold without confirming biological absence. "
        "Only VALID is eligible for a RESTORE call."
    )
    fig.tight_layout(rect=(0.02, 0.16, 1, 0.97))
    _add_bounded_footer(
        fig,
        list(axes.ravel()),
        footer_text,
        bottom=0.012,
        height=0.045,
        fontsize=8.5,
    )
    for suffix in ("png", "pdf"):
        fig.savefig(
            output_dir / f"figure1_pair_cloud_examples.{suffix}",
            dpi=240,
            bbox_inches="tight",
        )
    plt.close(fig)


def _annotated_heatmap(
    ax: plt.Axes,
    values: np.ndarray,
    row_labels: Sequence[str],
    column_labels: Sequence[str],
    *,
    cmap,
    norm=None,
    value_format: str,
    title: str,
    text_values: np.ndarray | None = None,
) -> None:
    image = ax.imshow(values, aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(np.arange(len(column_labels)))
    ax.set_xticklabels(column_labels, rotation=35, ha="right", fontsize=9)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=9)
    annotations = values if text_values is None else text_values
    for row in range(values.shape[0]):
        for column in range(values.shape[1]):
            ax.text(
                column,
                row,
                value_format.format(annotations[row, column]),
                ha="center",
                va="center",
                fontsize=8,
                color="black",
            )
    ax.set_xticks(np.arange(-0.5, len(column_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.2)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title(title, fontsize=12, fontweight="bold")
    return image


def create_summary_figure(panel_metrics: pd.DataFrame, output_dir: Path) -> None:
    """Summarize locked pair states, fitting QC, and all-cell assignments."""
    pair_labels = [pair_label(*pair) for pair in BASELINE_PAIRS]
    state_order = [
        PairState.TECHNICAL_FAILURE.value,
        PairState.NO_REFERENCE_POPULATION.value,
        PairState.INVALID_PAIR.value,
        PairState.MODEL_UNSTABLE.value,
        PairState.LOW_SIGNAL_UNRESOLVED.value,
        PairState.NO_TARGET_EXPRESSION_CONFIRMED.value,
        PairState.PAIR_REVIEW_PENDING.value,
        PairState.VALID.value,
    ]
    status_code = {state: index for index, state in enumerate(state_order)}
    status_matrix = np.zeros((len(BASELINE_PAIRS), len(DONOR_ORDER)), dtype=int)
    status_text = np.empty_like(status_matrix, dtype=object)
    retained = np.zeros_like(status_matrix, dtype=float)
    minority_arm = np.zeros_like(status_matrix, dtype=float)
    double_high = np.zeros_like(status_matrix, dtype=float)
    for row, pair in enumerate(BASELINE_PAIRS):
        for column, donor in enumerate(DONOR_ORDER):
            selected = panel_metrics[
                (panel_metrics["target"] == pair[0])
                & (panel_metrics["reference"] == pair[1])
                & (panel_metrics["donor"].astype(str) == donor)
            ]
            if len(selected) != 1:
                raise ValueError(
                    f"expected one panel metric for {donor} {pair_label(*pair)}"
                )
            status = selected.iloc[0]["state"]
            status_matrix[row, column] = status_code[status]
            status_text[row, column] = STATUS_LABELS[status]
            retained[row, column] = float(
                selected.iloc[0]["qc_retained_fraction"]
            )
            reference_n = float(selected.iloc[0]["reference_n"])
            target_n = float(selected.iloc[0]["target_n"])
            minority_arm[row, column] = min(reference_n, target_n) / (
                reference_n + target_n
            )
            double_high[row, column] = float(
                selected.iloc[0]["double_high_fraction"]
            )

    fig = plt.figure(figsize=(15.5, 17.5))
    grid = fig.add_gridspec(4, 1, hspace=0.55)
    status_cmap = ListedColormap(
        [STATUS_COLORS[state] for state in state_order]
    )
    status_norm = BoundaryNorm(
        np.arange(-0.5, len(state_order) + 0.5), status_cmap.N
    )
    ax1 = fig.add_subplot(grid[0])
    _annotated_heatmap(
        ax1,
        status_matrix,
        pair_labels,
        DONOR_ORDER,
        cmap=status_cmap,
        norm=status_norm,
        value_format="{}",
        title="A. Auditable donor/pair state (only VALID permits calls)",
        text_values=status_text,
    )

    ax2 = fig.add_subplot(grid[1])
    _annotated_heatmap(
        ax2,
        retained,
        pair_labels,
        DONOR_ORDER,
        cmap="Blues",
        norm=plt.Normalize(0, 1),
        value_format="{:.0%}",
        title="B. Cells eligible for fitting after raw-linear double-low QC",
    )

    ax3 = fig.add_subplot(grid[2])
    _annotated_heatmap(
        ax3,
        minority_arm,
        pair_labels,
        DONOR_ORDER,
        cmap="Purples",
        norm=plt.Normalize(0, 0.5),
        value_format="{:.0%}",
        title="C. Minority NNMF arm fraction (all assigned cells; fitting was equal-N)",
    )

    ax4 = fig.add_subplot(grid[3])
    _annotated_heatmap(
        ax4,
        double_high,
        pair_labels,
        DONOR_ORDER,
        cmap="Oranges",
        norm=plt.Normalize(0, 1),
        value_format="{:.0%}",
        title="D. Retained cells above both marker floors",
    )
    fig.suptitle(
        "Locked RESTORE input QC, NNMF stability, and review-state summary",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )
    fig.text(
        0.5,
        0.015,
        (
            "Donor-local raw-MFI arm floors are explicitly marker-specific: CD3e "
            "uses max(Otsu, Triangle), while the other screened markers use Otsu. Cells "
            "are retained when either marker clears its floor for fitting only; every "
            "complete canonical cell is then assigned. NNMF uses equal-arm fitting, robust "
            "per-feature scaling, and three locked seeds. Full REDSEA control maxima are "
            "divisor candidates; sampled maxima are QC only. Weak signal remains unresolved "
            "unless independent image review confirms absence."
        ),
        ha="center",
        fontsize=9,
    )
    fig.subplots_adjust(top=0.95, bottom=0.06, left=0.2, right=0.98)
    for suffix in ("png", "pdf"):
        fig.savefig(
            output_dir / f"figure2_input_qc_summary.{suffix}",
            dpi=240,
            bbox_inches="tight",
        )
    plt.close(fig)


def _config_from_specification(specification: dict) -> PairValidationConfig:
    parameters = dict(specification["parameters"])
    for key in ("feature_compartments", "seeds", "threshold_sample_sizes"):
        parameters[key] = tuple(parameters[key])
    return PairValidationConfig(**parameters)


def _assert_metrics_match(expected: pd.DataFrame, observed: pd.DataFrame) -> None:
    """Prove figure recomputation used exactly the locked table implementation."""
    keys = ["donor", "target", "reference"]
    expected = expected.copy()
    observed = observed.copy()
    expected["donor"] = expected["donor"].astype(str)
    observed["donor"] = observed["donor"].astype(str)
    expected = expected.sort_values(keys).reset_index(drop=True)
    observed = observed.sort_values(keys).reset_index(drop=True)
    if expected[keys].to_dict("records") != observed[keys].to_dict("records"):
        raise ValueError("figure and locked-result donor/pair keys differ")
    missing = set(observed.columns).difference(expected.columns)
    if missing:
        raise ValueError(f"locked result is missing figure metric columns: {sorted(missing)}")
    for column in observed.columns:
        first = expected[column]
        second = observed[column]
        if pd.api.types.is_numeric_dtype(first) and pd.api.types.is_numeric_dtype(second):
            if not np.allclose(
                first.to_numpy(float),
                second.to_numpy(float),
                rtol=1e-10,
                atol=1e-12,
                equal_nan=True,
            ):
                raise ValueError(f"figure recomputation differs in metric {column}")
        else:
            if not first.fillna("").astype(str).equals(
                second.fillna("").astype(str)
            ):
                raise ValueError(f"figure recomputation differs in metric {column}")


def run_figures(
    sample_root: Path,
    raw_cells: Path,
    redsea_cells: Path,
    results_dir: Path,
    output_dir: Path,
    *,
    donors: Sequence[str] = DONOR_ORDER,
) -> dict[str, Path]:
    """Recompute locked evaluations, verify their tables, and render figures."""
    donors = tuple(str(donor) for donor in donors)
    if donors != DONOR_ORDER:
        raise ValueError(f"pilot figures require donor order {DONOR_ORDER}")
    results_dir = Path(results_dir).resolve()
    method_spec_path = results_dir / "method_spec.json"
    pair_results_path = results_dir / "pair_evaluations.csv"
    review_path = results_dir / "expression_reviews.csv"
    pair_review_path = results_dir / "pair_reviews.csv"
    for required in (
        method_spec_path,
        pair_results_path,
        review_path,
        pair_review_path,
    ):
        if not required.is_file():
            raise FileNotFoundError(f"locked result artifact not found: {required}")
    method_spec = json.loads(method_spec_path.read_text())
    config = _config_from_specification(method_spec)
    if method_spec != config.specification():
        raise ValueError("locked method specification is not canonical")
    review_map, _ = load_expression_reviews(review_path)
    pair_review_map, _ = load_pair_reviews(pair_review_path)

    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    temp_dir.mkdir(parents=True)

    try:
        donor_data = {}
        for donor in donors:
            donor_data[donor], _ = load_validation_sample(
                Path(sample_root), Path(raw_cells), Path(redsea_cells), donor
            )
        panels, panel_metrics = create_pair_cloud_figures(
            donor_data,
            temp_dir,
            config=config,
            expression_reviews=review_map,
            pair_reviews=pair_review_map,
            display_seed=config.seeds[0],
        )
        expected_metrics = pd.read_csv(pair_results_path, dtype={"donor": str})
        _assert_metrics_match(expected_metrics, panel_metrics)
        panel_metrics.to_csv(temp_dir / "figure_panel_metrics.csv", index=False)
        create_example_figure(panels, temp_dir, seed=config.seeds[0])
        create_summary_figure(panel_metrics, temp_dir)
        manifest = {
            "method_version": METHOD_VERSION,
            "sample_root": str(Path(sample_root).resolve()),
            "raw_cells": str(Path(raw_cells).resolve()),
            "redsea_cells": str(Path(redsea_cells).resolve()),
            "results_dir": str(results_dir),
            "donors": list(donors),
            "method": FEATURE_METHOD,
            "method_specification": method_spec,
            "table_verification": (
                "Every figure metric was recomputed through evaluate_locked_pair and "
                "matched pair_evaluations.csv."
            ),
            "axis_convention": (
                "RESTORE Figure 7 orientation: robust min-max-scaled linear REDSEA "
                "reference whole-cell MFI on x and target whole-cell MFI on y; no log."
            ),
            "threshold_display": (
                "Solid lines are full inferred-control maxima (manuscript divisor "
                "candidates). Dashed lines and bands are fixed-N sampled-maximum QC. "
                "Dotted lines are full-control mean+3SD sensitivity values."
            ),
            "output_policy": (
                "These are validation artifacts. No RESTORE-normalized values or "
                "phenotype calls were written."
            ),
        }
        (temp_dir / "figure_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        temp_dir.rename(output_dir)
    except Exception:
        import shutil

        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise
    return {
        "examples": output_dir / "figure1_pair_cloud_examples.png",
        "summary": output_dir / "figure2_input_qc_summary.png",
        "all_pairs": output_dir / "figureS1_all_pair_clouds.pdf",
    }


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--raw-cells", type=Path, required=True)
    parser.add_argument("--redsea-cells", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    outputs = run_figures(
        args.sample,
        args.raw_cells,
        args.redsea_cells,
        args.results,
        args.output,
    )
    for name, path in outputs.items():
        print(f"[figure] {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
