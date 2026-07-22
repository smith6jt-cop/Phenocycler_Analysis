#!/usr/bin/env python3
"""Build deterministic spatial-review artifacts for shortlisted RESTORE pairs."""

from __future__ import annotations

import argparse
import json
import re
import shutil
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib

if not matplotlib.get_backend().startswith("module://"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import tifffile  # noqa: E402
from matplotlib.colors import to_rgba  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from .restore_validation import (  # noqa: E402
    ExpressionReview,
    MARKER_INPUT_FLOOR_METHODS,
    METHOD_VERSION,
    PairReview,
    PairValidationConfig,
    evaluate_locked_pair,
    load_expression_reviews,
    load_pair_reviews,
    load_validation_sample,
)
from .cohort import DONOR_EXCLUSIONS, ensure_eligible_donors
from .spatial import read_ome_region


SHORTLIST_PAIRS: tuple[tuple[str, str], ...] = (
    ("E_cadherin", "CD31"),
    ("CD31", "E_cadherin"),
    ("B3TUBB", "EpCAM"),
    ("Vimentin", "E_cadherin"),
    ("CD3e", "E_cadherin"),
    ("CD3e", "CD68"),
    ("CD3e", "CD163"),
    ("CD20", "CD68"),
    ("CD20", "CD163"),
    ("CD68", "EpCAM"),
    ("CD68", "CD20"),
    ("CD11b", "Pan_Cytokeratin"),
    ("CD11b", "CD20"),
)

PAIR_RATIONALE: dict[tuple[str, str], str] = {
    ("E_cadherin", "CD31"): "reciprocal epithelial/endothelial structural anchor",
    ("CD31", "E_cadherin"): "reciprocal endothelial/epithelial structural anchor",
    ("B3TUBB", "EpCAM"): "neural candidate; must exclude endocrine TUBB3",
    ("Vimentin", "E_cadherin"): "mesenchymal candidate; must exclude genuine epithelial Vimentin",
    ("CD3e", "E_cadherin"): "best fixed-reference availability; epithelial adjacency risk",
    ("CD3e", "CD68"): "literature-supported T/myeloid exclusion; weak target separation in many donors",
    ("CD3e", "CD163"): "literature-supported T/macrophage exclusion; reference absent in some donors",
    ("CD20", "CD68"): "strongest fixed B-cell/myeloid candidate in the expanded screen",
    ("CD20", "CD163"): "B-cell/macrophage alternative; reference absent in some donors",
    ("CD68", "EpCAM"): "strongest computational CD68 candidate; epithelial adjacency risk",
    ("CD68", "CD20"): "biologically exclusive comparator; CD20 reference absent in some donors",
    ("CD11b", "Pan_Cytokeratin"): "best fixed-reference availability; epithelial adjacency risk",
    ("CD11b", "CD20"): "least-confounded lymphoid reference; often lacks a supported reference arm",
}

CATEGORY_CODE = {
    "Reference control arm": 1,
    "High-confidence target anchor (sets Step 1)": 2,
    "Lower-confidence target component (reference-negative; retained)": 5,
    "Double-high": 3,
    "Other retained (not a divisor control)": 4,
    "Double-low (fit excluded)": 0,
    "Unavailable": -1,
}
CATEGORY_COLOR = {
    "Reference control arm": "#0072B2",
    "High-confidence target anchor (sets Step 1)": "#D55E00",
    "Lower-confidence target component (reference-negative; retained)": "#6A3D9A",
    "Double-high": "#CC79A7",
    "Other retained (not a divisor control)": "#7F7F7F",
    "Double-low (fit excluded)": "#D9D9D9",
    "Unavailable": "#222222",
}
CATEGORY_ORDER = tuple(CATEGORY_CODE)
CATEGORY_LEGEND = {
    "Reference control arm": "Reference controls: set Step 2 target divisor",
    "High-confidence target anchor (sets Step 1)": (
        "High-confidence target anchor: sets Step 1 vertical separator"
    ),
    "Lower-confidence target component (reference-negative; retained)": (
        "Input-QC-retained lower-confidence target component left of Step 1: retained even below Step 2"
    ),
    "Double-high": "Double-high: excluded from both maxima",
    "Other retained (not a divisor control)": "Other retained: excluded from both maxima",
    "Double-low (fit excluded)": "Double-low: excluded from fitting and maxima",
    "Unavailable": "Unavailable required measurements",
}
THRESHOLD_COLOR = "#005A42"
PDF_SUPTITLE_FONT_SIZE = 20
PDF_ROW_TITLE_FONT_SIZE = 12
PDF_AXIS_LABEL_FONT_SIZE = 14
PDF_TICK_FONT_SIZE = 11
PDF_IMAGE_TITLE_FONT_SIZE = 10
PDF_IMAGE_LABEL_FONT_SIZE = 9
PDF_IMAGE_TICK_FONT_SIZE = 8
PDF_LEGEND_FONT_SIZE = 11
PDF_LEGEND_TITLE_FONT_SIZE = 14
PDF_ROW_HEIGHT_IN = 5.6
PDF_TITLE_HEIGHT_IN = 2.25
PDF_LEGEND_HEIGHT_IN = 3.2
MARKER_IMAGE_MAX_WIDTH_PX = 1_000
MARKER_IMAGE_MAX_HEIGHT_PX = 600
REPRESENTATIVE_ZOOM_SPAN_FRACTION = 0.10
MARKER_DISPLAY_LOWER_QUANTILE = 0.01
MARKER_DISPLAY_UPPER_QUANTILE = 0.998
MARKER_DISPLAY_GAMMA = 0.65
TARGET_MARKER_RGB = (1.0, 0.12, 0.52)
REFERENCE_MARKER_RGB = (0.0, 0.78, 1.0)
QUPATH_IMPORTER_NAME = "import_restore_pair_review.groovy"
QUPATH_REVIEW_DIR_TOKEN = "__RESTORE_REVIEW_CSV_DIR__"
REVIEW_WORKFLOW_NAME = "START_HERE_RESTORE_REVIEW.txt"
REVIEW_CHECKLIST_NAME = "review_decisions.csv"


@dataclass(frozen=True)
class MarkerZoom:
    """Raw QPTIFF target/reference pixels aligned to a physical crop."""

    target_pixels: np.ndarray
    reference_pixels: np.ndarray
    x_bounds_um: tuple[float, float]
    y_bounds_um: tuple[float, float]
    target_channel: str
    reference_channel: str
    source_image: str
    source_path: Path
    pyramid_level: int
    downsample_x: float
    downsample_y: float
    pixel_size_x_um: float
    pixel_size_y_um: float

    @property
    def extent(self) -> tuple[float, float, float, float]:
        return (
            self.x_bounds_um[0],
            self.x_bounds_um[1],
            self.y_bounds_um[1],
            self.y_bounds_um[0],
        )


def _config_from_specification(specification: dict) -> PairValidationConfig:
    parameters = dict(specification["parameters"])
    for key in ("feature_compartments", "seeds", "threshold_sample_sizes"):
        parameters[key] = tuple(parameters[key])
    return PairValidationConfig(**parameters)


def summarize_pairs(evaluations: pd.DataFrame) -> pd.DataFrame:
    """Summarize donor-local states without turning the summary into acceptance."""
    evaluations = evaluations.copy()
    evaluations["sample_to_full"] = (
        evaluations["target_sample_maximum"]
        / evaluations["candidate_full_maximum"]
    )
    states = (
        evaluations.groupby(["target", "reference", "state"])
        .size()
        .unstack(fill_value=0)
    )
    metrics = (
        evaluations.groupby(["target", "reference"])
        .agg(
            donor_n=("donor", "nunique"),
            target_fold_min=("target_fold", "min"),
            target_fold_median=("target_fold", "median"),
            reference_fold_min=("reference_fold", "min"),
            reference_fold_median=("reference_fold", "median"),
            seed_jaccard_min=("min_control_jaccard", "min"),
            seed_jaccard_median=("min_control_jaccard", "median"),
            double_high_median=("double_high_fraction", "median"),
            double_high_max=("double_high_fraction", "max"),
            sample_to_full_min=("sample_to_full", "min"),
            sample_to_full_median=("sample_to_full", "median"),
            unavailable_cell_n=("n_unavailable", "sum"),
            target_anchor_n=("target_anchor_n", "sum"),
            target_supported_n=("target_supported_n", "sum"),
            target_supported_additional_n=(
                "target_supported_additional_n",
                "sum",
            ),
            target_supported_below_divisor_n=(
                "target_supported_below_divisor_n",
                "sum",
            ),
            target_supported_additional_below_divisor_n=(
                "target_supported_additional_below_divisor_n",
                "sum",
            ),
            target_supported_rescued_fraction_median=(
                "target_supported_rescued_fraction",
                "median",
            ),
        )
        .join(states)
        .reset_index()
    )
    metrics["shortlisted"] = [
        (target, reference) in SHORTLIST_PAIRS
        for target, reference in metrics[["target", "reference"]].itertuples(
            index=False, name=None
        )
    ]
    metrics["review_rationale"] = [
        PAIR_RATIONALE.get((target, reference), "expanded-screen comparator")
        for target, reference in metrics[["target", "reference"]].itertuples(
            index=False, name=None
        )
    ]
    return metrics


def _extreme_donor(
    pair_rows: pd.DataFrame, column: str, *, largest: bool
) -> str:
    ordered = pair_rows.sort_values(
        [column, "donor"],
        ascending=[not largest, True],
        na_position="last",
    )
    return str(ordered.iloc[0]["donor"])


def select_review_queue(
    evaluations: pd.DataFrame,
    *,
    pairs: Sequence[tuple[str, str]] = SHORTLIST_PAIRS,
) -> pd.DataFrame:
    """Select technical extremes and controls deterministically for image review."""
    evaluations = evaluations.copy()
    evaluations["donor"] = evaluations["donor"].astype(str)
    evaluations["sample_to_full"] = (
        evaluations["target_sample_maximum"]
        / evaluations["candidate_full_maximum"]
    )
    rows: list[dict] = []
    for target, reference in pairs:
        pair_rows = evaluations[
            (evaluations["target"] == target)
            & (evaluations["reference"] == reference)
        ].copy()
        if pair_rows.empty:
            raise ValueError(f"missing expanded-screen pair: {target} <- {reference}")
        reasons: dict[str, set[str]] = {}

        def add(donor: str, reason: str) -> None:
            if donor in set(pair_rows["donor"]):
                reasons.setdefault(donor, set()).add(reason)

        add(
            _extreme_donor(pair_rows, "target_fold", largest=False),
            "lowest target-arm fold",
        )
        add(
            _extreme_donor(pair_rows, "reference_fold", largest=False),
            "lowest reference-arm fold",
        )
        add(
            _extreme_donor(pair_rows, "min_control_jaccard", largest=False),
            "lowest seed stability",
        )
        add(
            _extreme_donor(pair_rows, "double_high_fraction", largest=True),
            "highest double-high fraction",
        )
        add(
            _extreme_donor(pair_rows, "sample_to_full", largest=False),
            "largest sampled/full-maximum divergence",
        )
        median_fold = float(pair_rows["target_fold"].median())
        median_donor = (
            pair_rows.assign(
                _median_distance=(pair_rows["target_fold"] - median_fold).abs()
            )
            .sort_values(
                ["_median_distance", "donor"],
                ascending=[True, True],
                na_position="last",
            )
            .iloc[0]["donor"]
        )
        add(str(median_donor), "median target-arm fold")
        if target == "CD20":
            add("6539", "maintainer-confirmed target absence")
        if target == "B3TUBB":
            arm_fraction = pair_rows["target_n"] / pair_rows["n_assigned"]
            add(
                str(pair_rows.loc[arm_fraction.idxmin(), "donor"]),
                "sparsest inferred target arm",
            )

        for donor in sorted(reasons):
            selected = pair_rows[pair_rows["donor"] == donor]
            if len(selected) != 1:
                raise ValueError(
                    f"expected one screen row for {donor} {target} <- {reference}"
                )
            row = selected.iloc[0].to_dict()
            row["selection_reasons"] = "; ".join(sorted(reasons[donor]))
            row["review_rationale"] = PAIR_RATIONALE[(target, reference)]
            rows.append(row)
    return pd.DataFrame(rows).sort_values(
        ["target", "reference", "donor"]
    ).reset_index(drop=True)


def review_categories(evaluation: dict, n_rows: int) -> np.ndarray:
    """Map every evaluated cell to an explicit review category."""
    categories = np.full(n_rows, "Unavailable", dtype=object)
    valid_idx = np.asarray(evaluation["valid_idx"], dtype=int)
    valid_categories = np.full(
        len(valid_idx), "Double-low (fit excluded)", dtype=object
    )
    retained = np.asarray(evaluation["qc_retained"], dtype=bool)
    valid_categories[retained] = "Other retained (not a divisor control)"
    target_high = evaluation["target_raw"] > evaluation["target_input_floor"]
    reference_high = (
        evaluation["reference_raw"] > evaluation["reference_input_floor"]
    )
    valid_categories[retained & target_high & reference_high] = "Double-high"
    target_anchor = np.asarray(evaluation["target_population"], dtype=bool)
    target_supported = np.asarray(evaluation["target_supported"], dtype=bool)
    valid_categories[
        retained & target_supported & ~target_anchor
    ] = "Lower-confidence target component (reference-negative; retained)"
    valid_categories[target_anchor] = "High-confidence target anchor (sets Step 1)"
    valid_categories[
        np.asarray(evaluation["reference_control"], dtype=bool)
    ] = "Reference control arm"
    categories[valid_idx] = valid_categories
    return categories


def project_review_crop(
    full_donor_df: pd.DataFrame,
    evaluation: dict,
    x_bounds_um: tuple[float, float],
    y_bounds_um: tuple[float, float],
) -> pd.DataFrame:
    """Project the locked sample-fitted model onto every canonical crop cell."""
    x = full_donor_df["x_canonical"].to_numpy(float)
    y = full_donor_df["y_canonical"].to_numpy(float)
    inside = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= x_bounds_um[0])
        & (x <= x_bounds_um[1])
        & (y >= y_bounds_um[0])
        & (y <= y_bounds_um[1])
    )
    crop = full_donor_df.loc[inside].copy()
    if crop.empty:
        raise ValueError("representative marker crop contains no canonical cells")
    if not crop["object_id"].is_unique:
        raise ValueError("representative marker crop contains duplicate object IDs")

    feature_columns = list(evaluation["feature_columns"])
    target = evaluation["target"]
    reference = evaluation["reference"]
    required = {
        "object_id",
        "image",
        "x_canonical",
        "y_canonical",
        f"raw__{target}",
        f"raw__{reference}",
        f"redsea__{target}",
        f"redsea__{reference}",
        *feature_columns,
    }
    missing = required.difference(crop.columns)
    if missing:
        raise ValueError(
            f"full-cell review projection is missing columns: {sorted(missing)}"
        )

    features_all = crop[feature_columns].to_numpy(float)
    target_raw_all = crop[f"raw__{target}"].to_numpy(float)
    reference_raw_all = crop[f"raw__{reference}"].to_numpy(float)
    target_redsea_all = crop[f"redsea__{target}"].to_numpy(float)
    reference_redsea_all = crop[f"redsea__{reference}"].to_numpy(float)
    valid = (
        np.isfinite(features_all).all(axis=1)
        & (features_all >= 0).all(axis=1)
        & np.isfinite(target_raw_all)
        & (target_raw_all >= 0)
        & np.isfinite(reference_raw_all)
        & (reference_raw_all >= 0)
        & np.isfinite(target_redsea_all)
        & (target_redsea_all >= 0)
        & np.isfinite(reference_redsea_all)
        & (reference_redsea_all >= 0)
    )
    valid_idx = np.flatnonzero(valid)
    if not len(valid_idx):
        raise ValueError("representative marker crop has no projectable cells")

    scale_lower = np.asarray(evaluation["scale_lower"], dtype=float)
    scale_upper = np.asarray(evaluation["scale_upper"], dtype=float)
    components = np.asarray(evaluation["components"], dtype=float)
    if (
        scale_lower.shape != (len(feature_columns),)
        or scale_upper.shape != (len(feature_columns),)
        or components.shape[1] != len(feature_columns)
    ):
        raise ValueError("locked NNMF projection parameters do not match features")
    scaled = np.clip(
        (features_all[valid_idx] - scale_lower) / (scale_upper - scale_lower),
        0.0,
        1.0,
    )
    from sklearn.decomposition import non_negative_factorization

    weights, _, _ = non_negative_factorization(
        scaled,
        H=components.copy(),
        n_components=len(components),
        init="custom",
        update_H=False,
        solver="cd",
        max_iter=2_000,
        tol=1e-5,
    )
    labels = np.argmax(weights, axis=1)
    target_raw = target_raw_all[valid_idx]
    reference_raw = reference_raw_all[valid_idx]
    target_redsea = target_redsea_all[valid_idx]
    reference_redsea = reference_redsea_all[valid_idx]
    target_floor = float(evaluation["target_input_floor"])
    reference_floor = float(evaluation["reference_input_floor"])
    qc_retained = (target_raw > target_floor) | (
        reference_raw > reference_floor
    )
    target_population = (
        qc_retained
        & (labels == evaluation["target_group"])
        & (target_raw > target_floor)
        & (reference_raw <= reference_floor)
    )
    reference_candidates = (
        qc_retained
        & (labels == evaluation["reference_group"])
        & (target_raw <= target_floor)
        & (reference_raw > reference_floor)
    )
    separator = evaluation["reference_separator"]
    if separator is None:
        target_supported = np.zeros(len(valid_idx), dtype=bool)
        reference_control = np.zeros(len(valid_idx), dtype=bool)
    else:
        target_supported = (
            qc_retained
            & (labels == evaluation["target_group"])
            & (reference_redsea <= float(separator))
        )
        reference_control = reference_candidates & (
            reference_redsea > float(separator)
        )
    candidate_divisor = evaluation["candidate_full_maximum"]
    target_supported_below_divisor = np.zeros(len(valid_idx), dtype=bool)
    if candidate_divisor is not None:
        target_supported_below_divisor = target_supported & (
            target_redsea <= float(candidate_divisor)
        )
    projected_evaluation = {
        "valid_idx": valid_idx,
        "qc_retained": qc_retained,
        "reference_control": reference_control,
        "target_population": target_population,
        "target_supported": target_supported,
        "target_raw": target_raw,
        "reference_raw": reference_raw,
        "target_input_floor": target_floor,
        "reference_input_floor": reference_floor,
    }
    categories = review_categories(projected_evaluation, len(crop))
    full_target_population = np.zeros(len(crop), dtype=bool)
    full_target_population[valid_idx] = target_population
    full_qc_retained = np.zeros(len(crop), dtype=bool)
    full_qc_retained[valid_idx] = qc_retained
    full_target_supported = np.zeros(len(crop), dtype=bool)
    full_target_supported[valid_idx] = target_supported
    full_below_divisor = np.zeros(len(crop), dtype=bool)
    full_below_divisor[valid_idx] = target_supported_below_divisor

    projected = crop[
        ["object_id", "image", "x_canonical", "y_canonical"]
    ].copy()
    projected["review_category"] = categories
    projected["input_qc_retained"] = full_qc_retained
    projected["target_anchor"] = full_target_population
    projected["target_component_supported"] = full_target_supported
    projected["retained_below_divisor"] = full_below_divisor
    return projected.reset_index(drop=True)


def _validate_projection_overlap(
    sample_df: pd.DataFrame,
    sample_categories: np.ndarray,
    projected_crop: pd.DataFrame,
    x_bounds_um: tuple[float, float],
    y_bounds_um: tuple[float, float],
) -> int:
    x = sample_df["x_canonical"].to_numpy(float)
    y = sample_df["y_canonical"].to_numpy(float)
    inside = (
        np.isfinite(x)
        & np.isfinite(y)
        & (x >= x_bounds_um[0])
        & (x <= x_bounds_um[1])
        & (y >= y_bounds_um[0])
        & (y <= y_bounds_um[1])
    )
    expected = sample_df.loc[inside, ["object_id"]].copy()
    expected["expected_category"] = sample_categories[inside]
    observed = projected_crop[["object_id", "review_category"]]
    comparison = expected.merge(
        observed,
        on="object_id",
        how="left",
        validate="one_to_one",
    )
    if comparison["review_category"].isna().any():
        raise ValueError("full-cell crop projection omitted modeling-sample cells")
    mismatched = (
        comparison["expected_category"] != comparison["review_category"]
    )
    if mismatched.any():
        raise ValueError(
            "full-cell crop projection does not reproduce sample classifications"
        )
    return int(len(comparison))


def _scaled(values: np.ndarray, bounds: tuple[float, float]) -> np.ndarray:
    lower, upper = bounds
    return np.clip((np.asarray(values) - lower) / (upper - lower), 0.0, 1.0)


def _point_colors(
    categories: np.ndarray,
    *,
    emphasis_alpha: float,
    other_alpha: float,
    double_low_alpha: float,
) -> np.ndarray:
    emphasized = {
        "Reference control arm",
        "High-confidence target anchor (sets Step 1)",
        "Lower-confidence target component (reference-negative; retained)",
        "Double-high",
    }
    categories = np.asarray(categories, dtype=object)
    colors = np.empty((len(categories), 4), dtype=float)
    for category in CATEGORY_ORDER:
        if category == "Double-low (fit excluded)":
            alpha = double_low_alpha
        elif category in emphasized:
            alpha = emphasis_alpha
        else:
            alpha = other_alpha
        colors[categories == category] = to_rgba(
            CATEGORY_COLOR[category], alpha
        )
    return colors


def representative_zoom_bounds(
    donor_df: pd.DataFrame,
    categories: np.ndarray,
    *,
    span_fraction: float = REPRESENTATIVE_ZOOM_SPAN_FRACTION,
    display_aspect: float = 2.2,
    bins: int = 24,
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Choose a deterministic target-focused crop with undistorted display geometry."""
    if not 0 < span_fraction <= 1:
        raise ValueError("span_fraction must be in (0, 1]")
    if not np.isfinite(display_aspect) or display_aspect <= 0:
        raise ValueError("display_aspect must be finite and positive")
    if bins < 2:
        raise ValueError("bins must be >= 2")
    x = donor_df["x_canonical"].to_numpy(float)
    y = donor_df["y_canonical"].to_numpy(float)
    finite = np.isfinite(x) & np.isfinite(y)
    if not finite.any():
        raise ValueError("spatial review requires finite coordinates")
    full_x = (float(x[finite].min()), float(x[finite].max()))
    full_y = (float(y[finite].min()), float(y[finite].max()))
    if full_x[0] == full_x[1] or full_y[0] == full_y[1]:
        raise ValueError("spatial review coordinates must span both axes")

    focus = finite & np.isin(
        categories,
        (
            "High-confidence target anchor (sets Step 1)",
            "Lower-confidence target component (reference-negative; retained)",
        ),
    )
    if not focus.any():
        for category in (
            "Double-high",
            "Reference control arm",
            "Other retained (not a divisor control)",
        ):
            focus = finite & (categories == category)
            if focus.any():
                break
    if not focus.any():
        focus = finite
    hist, x_edges, y_edges = np.histogram2d(
        x[focus], y[focus], bins=bins, range=(full_x, full_y)
    )
    x_bin, y_bin = np.unravel_index(np.argmax(hist), hist.shape)
    in_peak = (
        focus
        & (x >= x_edges[x_bin])
        & (x <= x_edges[x_bin + 1])
        & (y >= y_edges[y_bin])
        & (y <= y_edges[y_bin + 1])
    )
    center_x = float(np.median(x[in_peak]))
    center_y = float(np.median(y[in_peak]))
    full_width = full_x[1] - full_x[0]
    full_height = full_y[1] - full_y[0]
    crop_area = full_width * full_height * span_fraction**2
    width = min(float(np.sqrt(crop_area * display_aspect)), full_width)
    height = width / display_aspect
    if height > full_height:
        height = full_height
        width = height * display_aspect

    def clipped(center: float, size: float, limits: tuple[float, float]):
        lower = max(limits[0], center - size / 2)
        upper = min(limits[1], lower + size)
        lower = max(limits[0], upper - size)
        return float(lower), float(upper)

    return clipped(center_x, width, full_x), clipped(center_y, height, full_y)


def _normalized_marker_name(marker: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(marker).lower())


def _qptiff_channel_names(series) -> list[str]:
    names = []
    for index, page in enumerate(series.levels[0].pages):
        description = page.description or ""
        match = re.search(r"<Biomarker>([^<]+)</Biomarker>", description)
        if match is None:
            match = re.search(r"<Name>([^<]+)</Name>", description)
        names.append(match.group(1).strip() if match else f"ch{index}")
    return names


def _marker_channel_index(
    channel_names: Sequence[str], marker: str, source: Path
) -> int:
    normalized = _normalized_marker_name(marker)
    matches = [
        index
        for index, channel_name in enumerate(channel_names)
        if _normalized_marker_name(channel_name) == normalized
    ]
    if len(matches) != 1:
        raise ValueError(
            f"{source.name}: expected exactly one QPTIFF channel for {marker}; "
            f"found {[channel_names[index] for index in matches]}"
        )
    return matches[0]


def _microns_per_pixel(page, resolution_tag: str) -> float:
    if resolution_tag not in page.tags or "ResolutionUnit" not in page.tags:
        raise ValueError(f"QPTIFF is missing {resolution_tag}/ResolutionUnit")
    numerator, denominator = page.tags[resolution_tag].value
    pixels_per_unit = float(numerator) / float(denominator)
    unit = int(page.tags["ResolutionUnit"].value)
    if unit == 2:
        microns_per_unit = 25_400.0
    elif unit == 3:
        microns_per_unit = 10_000.0
    else:
        raise ValueError(f"unsupported QPTIFF ResolutionUnit: {unit}")
    pixel_size = microns_per_unit / pixels_per_unit
    if not np.isfinite(pixel_size) or pixel_size <= 0:
        raise ValueError(f"invalid QPTIFF {resolution_tag}: {pixel_size}")
    return pixel_size


def _marker_pyramid_level(
    series,
    zoom_x: tuple[float, float],
    zoom_y: tuple[float, float],
    pixel_size_x_um: float,
    pixel_size_y_um: float,
) -> int:
    full_height, full_width = series.levels[0].shape[-2:]
    tiled_levels = []
    for index, level in enumerate(series.levels):
        if not level.pages[0].is_tiled:
            continue
        tiled_levels.append(index)
        level_height, level_width = level.shape[-2:]
        downsample_x = full_width / level_width
        downsample_y = full_height / level_height
        crop_width = (
            (zoom_x[1] - zoom_x[0]) / pixel_size_x_um / downsample_x
        )
        crop_height = (
            (zoom_y[1] - zoom_y[0]) / pixel_size_y_um / downsample_y
        )
        if (
            crop_width <= MARKER_IMAGE_MAX_WIDTH_PX
            and crop_height <= MARKER_IMAGE_MAX_HEIGHT_PX
        ):
            return index
    if not tiled_levels:
        raise ValueError("QPTIFF has no tiled pyramid level for regional reads")
    return tiled_levels[-1]


def load_marker_zoom(
    images_dir: Path,
    donor_df: pd.DataFrame,
    target: str,
    reference: str,
    zoom_x: tuple[float, float],
    zoom_y: tuple[float, float],
) -> MarkerZoom:
    """Read only the target/reference QPTIFF tiles needed for one review crop."""
    images = donor_df["image"].dropna().astype(str).unique()
    if len(images) != 1:
        raise ValueError(
            f"marker zoom requires exactly one donor image; found {images}"
        )
    source_image = images[0]
    image_name = Path(source_image.split(" - ", 1)[0].strip()).name
    source_path = Path(images_dir).resolve() / image_name
    if not source_path.is_file():
        raise FileNotFoundError(f"QPTIFF image not found: {source_path}")

    with tifffile.TiffFile(source_path) as image:
        series = image.series[0]
        channel_names = _qptiff_channel_names(series)
        target_index = _marker_channel_index(channel_names, target, source_path)
        reference_index = _marker_channel_index(
            channel_names, reference, source_path
        )
        full_page = series.levels[0].pages[0]
        pixel_size_x_um = _microns_per_pixel(full_page, "XResolution")
        pixel_size_y_um = _microns_per_pixel(full_page, "YResolution")
        level_index = _marker_pyramid_level(
            series,
            zoom_x,
            zoom_y,
            pixel_size_x_um,
            pixel_size_y_um,
        )
        full_height, full_width = series.levels[0].shape[-2:]
        level_height, level_width = series.levels[level_index].shape[-2:]
        downsample_x = full_width / level_width
        downsample_y = full_height / level_height
        x0 = max(
            0,
            int(np.floor(zoom_x[0] / pixel_size_x_um / downsample_x)),
        )
        x1 = min(
            level_width,
            int(np.ceil(zoom_x[1] / pixel_size_x_um / downsample_x)),
        )
        y0 = max(
            0,
            int(np.floor(zoom_y[0] / pixel_size_y_um / downsample_y)),
        )
        y1 = min(
            level_height,
            int(np.ceil(zoom_y[1] / pixel_size_y_um / downsample_y)),
        )
        if x1 <= x0 or y1 <= y0:
            raise ValueError(
                f"{source_path.name}: marker crop is outside the QPTIFF bounds"
            )
        target_pixels = read_ome_region(
            image, level_index, target_index, y0, y1, x0, x1
        )
        reference_pixels = read_ome_region(
            image, level_index, reference_index, y0, y1, x0, x1
        )

    if target_pixels.shape != reference_pixels.shape or target_pixels.ndim != 2:
        raise ValueError(
            f"{source_path.name}: target/reference marker crops are not aligned"
        )
    return MarkerZoom(
        target_pixels=target_pixels,
        reference_pixels=reference_pixels,
        x_bounds_um=(
            x0 * downsample_x * pixel_size_x_um,
            x1 * downsample_x * pixel_size_x_um,
        ),
        y_bounds_um=(
            y0 * downsample_y * pixel_size_y_um,
            y1 * downsample_y * pixel_size_y_um,
        ),
        target_channel=channel_names[target_index],
        reference_channel=channel_names[reference_index],
        source_image=source_image,
        source_path=source_path,
        pyramid_level=level_index,
        downsample_x=float(downsample_x),
        downsample_y=float(downsample_y),
        pixel_size_x_um=float(pixel_size_x_um),
        pixel_size_y_um=float(pixel_size_y_um),
    )


def _display_marker_channel(
    pixels: np.ndarray,
) -> tuple[np.ndarray, tuple[float, float]]:
    values = np.asarray(pixels, dtype=np.float32)
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        return np.zeros(values.shape, dtype=np.float32), (0.0, 0.0)
    lower, upper = np.quantile(
        positive,
        [MARKER_DISPLAY_LOWER_QUANTILE, MARKER_DISPLAY_UPPER_QUANTILE],
    )
    if upper <= lower:
        lower = 0.0
        upper = float(positive.max())
    if upper <= lower:
        return np.zeros(values.shape, dtype=np.float32), (
            float(lower),
            float(upper),
        )
    scaled = np.clip((values - lower) / (upper - lower), 0.0, 1.0)
    scaled[~np.isfinite(scaled)] = 0.0
    return np.power(scaled, MARKER_DISPLAY_GAMMA), (
        float(lower),
        float(upper),
    )


def _tinted_marker(channel: np.ndarray, color: Sequence[float]) -> np.ndarray:
    return np.clip(channel[..., np.newaxis] * np.asarray(color), 0.0, 1.0)


def _marker_panel_axes(
    spatial: plt.Axes,
) -> tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]:
    spatial.set_axis_off()
    positions = (
        [0.01, 0.53, 0.48, 0.40],
        [0.51, 0.53, 0.48, 0.40],
        [0.01, 0.04, 0.48, 0.40],
        [0.51, 0.04, 0.48, 0.40],
    )
    return tuple(
        spatial.inset_axes(position, facecolor="#080808")
        for position in positions
    )


def _marker_zoom_metadata(marker_zoom: MarkerZoom) -> dict:
    _, target_window = _display_marker_channel(marker_zoom.target_pixels)
    _, reference_window = _display_marker_channel(marker_zoom.reference_pixels)
    return {
        "source_qptiff": str(marker_zoom.source_path),
        "pyramid_level": marker_zoom.pyramid_level,
        "downsample_x": marker_zoom.downsample_x,
        "downsample_y": marker_zoom.downsample_y,
        "crop_shape_yx": list(marker_zoom.target_pixels.shape),
        "crop_x_um": list(marker_zoom.x_bounds_um),
        "crop_y_um": list(marker_zoom.y_bounds_um),
        "target_channel": marker_zoom.target_channel,
        "reference_channel": marker_zoom.reference_channel,
        "target_display_window": list(target_window),
        "reference_display_window": list(reference_window),
    }


def _draw_review_row(
    axes: Sequence[plt.Axes],
    donor_df: pd.DataFrame,
    evaluation: dict,
    categories: np.ndarray,
    reasons: str,
    marker_zoom: MarkerZoom,
    projected_crop: pd.DataFrame,
) -> None:
    if len(axes) != 2:
        raise ValueError("review rows require expression and composite spatial axes")
    if (
        marker_zoom.target_pixels.shape != marker_zoom.reference_pixels.shape
        or marker_zoom.target_pixels.ndim != 2
        or marker_zoom.target_pixels.size == 0
    ):
        raise ValueError("marker zoom requires aligned nonempty 2D channels")
    projected_required = {
        "x_canonical",
        "y_canonical",
        "review_category",
    }
    projected_missing = projected_required.difference(projected_crop.columns)
    if projected_missing or projected_crop.empty:
        raise ValueError(
            "full-cell projected crop is empty or missing columns: "
            f"{sorted(projected_missing)}"
        )
    cloud, spatial = axes
    target_axis, reference_axis, overlay_axis, merged_axis = _marker_panel_axes(
        spatial
    )
    whole_tissue = overlay_axis.inset_axes(
        [0.72, 0.59, 0.26, 0.37],
        facecolor="#080808",
        zorder=10,
    )
    selected = np.arange(len(donor_df))
    valid = evaluation["full_labels"] >= 0
    target_scaled = np.full(len(donor_df), np.nan)
    reference_scaled = np.full(len(donor_df), np.nan)
    target_scaled[valid] = _scaled(
        evaluation["target_redsea"], evaluation["display_target_bounds"]
    )
    reference_scaled[valid] = _scaled(
        evaluation["reference_redsea"], evaluation["display_reference_bounds"]
    )
    zoom_x = marker_zoom.x_bounds_um
    zoom_y = marker_zoom.y_bounds_um
    x = donor_df["x_canonical"].to_numpy(float)
    y = donor_df["y_canonical"].to_numpy(float)
    crop_x = projected_crop["x_canonical"].to_numpy(float)
    crop_y = projected_crop["y_canonical"].to_numpy(float)
    crop_categories = projected_crop["review_category"].to_numpy(object)
    cloud_visible = selected[valid[selected]]
    if len(cloud_visible):
        cloud.scatter(
            reference_scaled[cloud_visible],
            target_scaled[cloud_visible],
            s=1.5,
            color=_point_colors(
                categories[cloud_visible],
                emphasis_alpha=0.65,
                other_alpha=0.15,
                double_low_alpha=0.06,
            ),
            edgecolors="none",
            rasterized=True,
        )
    whole_tissue.scatter(
        donor_df.iloc[selected]["x_canonical"],
        donor_df.iloc[selected]["y_canonical"],
        s=0.3,
        color=_point_colors(
            categories[selected],
            emphasis_alpha=0.75,
            other_alpha=0.2,
            double_low_alpha=0.07,
        ),
        edgecolors="none",
        rasterized=True,
    )
    target_display, _ = _display_marker_channel(marker_zoom.target_pixels)
    reference_display, _ = _display_marker_channel(marker_zoom.reference_pixels)
    target_rgb = _tinted_marker(target_display, TARGET_MARKER_RGB)
    reference_rgb = _tinted_marker(reference_display, REFERENCE_MARKER_RGB)
    merged_rgb = np.clip(target_rgb + reference_rgb, 0.0, 1.0)
    target_axis.imshow(
        target_rgb,
        origin="upper",
        extent=marker_zoom.extent,
        interpolation="nearest",
    )
    reference_axis.imshow(
        reference_rgb,
        origin="upper",
        extent=marker_zoom.extent,
        interpolation="nearest",
    )
    merged_axis.imshow(
        merged_rgb,
        origin="upper",
        extent=marker_zoom.extent,
        interpolation="nearest",
    )
    overlay_colors = _point_colors(
        crop_categories,
        emphasis_alpha=0.9,
        other_alpha=0.45,
        double_low_alpha=0.22,
    )
    overlay_axis.scatter(
        crop_x,
        crop_y,
        s=2.0,
        color=overlay_colors,
        edgecolors="none",
        rasterized=True,
    )
    merged_axis.scatter(
        crop_x,
        crop_y,
        s=5.0,
        facecolors="none",
        edgecolors=overlay_colors,
        linewidths=0.35,
        rasterized=True,
    )
    if evaluation["target_full_stats"]:
        cloud.axhline(
            _scaled(
                [evaluation["target_full_stats"]["maximum"]],
                evaluation["display_target_bounds"],
            )[0],
            color=THRESHOLD_COLOR,
            lw=2.2,
            ls="-",
        )
    if evaluation["reference_full_stats"]:
        cloud.axvline(
            _scaled(
                [evaluation["reference_full_stats"]["maximum"]],
                evaluation["display_reference_bounds"],
            )[0],
            color=THRESHOLD_COLOR,
            lw=2.2,
            ls="--",
        )
    cloud.set_xlim(0, 1.02)
    cloud.set_ylim(0, 1.02)
    cloud.grid(color="#DDDDDD", lw=0.4)
    cloud.set_ylabel(
        f"Donor {evaluation['donor']}\n{evaluation['target']} scaled REDSEA MFI",
        fontsize=PDF_AXIS_LABEL_FONT_SIZE,
        labelpad=10,
    )
    cloud.set_title(
        (
            f"{evaluation['state']}\nTarget/reference fold "
            f"{evaluation['target_fold']:.2f}/{evaluation['reference_fold']:.2f} | "
            f"control Jaccard {evaluation['min_control_jaccard']:.2f}\n"
            f"High-confidence anchor n={evaluation['target_anchor_n']:,} | "
            "input-QC-retained lower-confidence target "
            f"n={evaluation['target_supported_additional_n']:,}\n"
            "Lower-confidence at/below Step 2 "
            f"n={evaluation['target_supported_additional_below_divisor_n']:,}"
        ),
        fontsize=PDF_ROW_TITLE_FONT_SIZE,
        pad=10,
        linespacing=1.3,
    )
    finite = np.isfinite(x) & np.isfinite(y)
    x_padding = max((x[finite].max() - x[finite].min()) * 0.01, 1.0)
    y_padding = max((y[finite].max() - y[finite].min()) * 0.01, 1.0)
    whole_tissue.set_xlim(x[finite].min() - x_padding, x[finite].max() + x_padding)
    whole_tissue.set_ylim(y[finite].max() + y_padding, y[finite].min() - y_padding)
    whole_tissue.add_patch(
        Rectangle(
            (zoom_x[0], zoom_y[0]),
            zoom_x[1] - zoom_x[0],
            zoom_y[1] - zoom_y[0],
            fill=False,
            edgecolor="white",
            linewidth=1.4,
        )
    )
    whole_tissue.set_aspect("equal", adjustable="box")
    whole_tissue.set_xticks([])
    whole_tissue.set_yticks([])
    for spine in whole_tissue.spines.values():
        spine.set_color("white")
        spine.set_linewidth(0.7)
    for axis in (target_axis, reference_axis, overlay_axis, merged_axis):
        axis.set_xlim(*zoom_x)
        axis.set_ylim(zoom_y[1], zoom_y[0])
        axis.set_aspect("equal", adjustable="box")
        axis.tick_params(
            axis="both",
            which="major",
            labelsize=PDF_IMAGE_TICK_FONT_SIZE,
        )
    target_axis.set_title(
        f"Raw target marker pixels — {evaluation['target']}",
        fontsize=PDF_IMAGE_TITLE_FONT_SIZE,
        pad=3,
    )
    reference_axis.set_title(
        f"Raw reference marker pixels — {evaluation['reference']}",
        fontsize=PDF_IMAGE_TITLE_FONT_SIZE,
        pad=3,
    )
    overlay_axis.set_title(
        "RESTORE review overlay",
        fontsize=PDF_IMAGE_TITLE_FONT_SIZE,
        pad=3,
    )
    merged_axis.set_title(
        "Merged markers + classification rings",
        fontsize=PDF_IMAGE_TITLE_FONT_SIZE,
        pad=3,
    )
    for axis in (target_axis, reference_axis):
        axis.set_xticks([])
        axis.set_yticks([])
    overlay_axis.set_xlabel(
        "X position (µm)", fontsize=PDF_IMAGE_LABEL_FONT_SIZE, labelpad=2
    )
    overlay_axis.set_ylabel(
        "Y position (µm)", fontsize=PDF_IMAGE_LABEL_FONT_SIZE, labelpad=2
    )
    merged_axis.set_xlabel(
        "X position (µm)", fontsize=PDF_IMAGE_LABEL_FONT_SIZE, labelpad=2
    )
    merged_axis.tick_params(labelleft=False)
    spatial.set_title(
        "Representative QPTIFF close-up (overview inset shows whole tissue)\n"
        f"All {len(projected_crop):,} canonical crop cells labeled; "
        "no display subsampling\n"
        + textwrap.fill(reasons, width=88),
        fontsize=PDF_ROW_TITLE_FONT_SIZE,
        pad=10,
    )
    cloud.tick_params(
        axis="both",
        which="major",
        labelsize=PDF_TICK_FONT_SIZE,
    )


def _create_review_figure(
    n_rows: int,
) -> tuple[plt.Figure, np.ndarray, plt.Axes]:
    """Create data rows plus a physically separate legend band."""
    if n_rows < 1:
        raise ValueError("review figure requires at least one row")
    height = (
        PDF_TITLE_HEIGHT_IN
        + n_rows * PDF_ROW_HEIGHT_IN
        + PDF_LEGEND_HEIGHT_IN
        + 0.35
    )
    fig = plt.figure(figsize=(20.0, height))
    grid = fig.add_gridspec(
        n_rows + 1,
        2,
        height_ratios=[1.0] * n_rows
        + [PDF_LEGEND_HEIGHT_IN / PDF_ROW_HEIGHT_IN],
        width_ratios=(1.0, 1.8),
        left=0.07,
        right=0.985,
        bottom=0.25 / height,
        top=1.0 - PDF_TITLE_HEIGHT_IN / height,
        wspace=0.12,
        hspace=0.44,
    )
    axes = np.empty((n_rows, 2), dtype=object)
    for row in range(n_rows):
        axes[row, 0] = fig.add_subplot(grid[row, 0])
        axes[row, 1] = fig.add_subplot(grid[row, 1])
    legend_axis = fig.add_subplot(grid[-1, :])
    legend_axis.set_axis_off()
    return fig, axes, legend_axis


def _review_legend_handles() -> list[Line2D]:
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=CATEGORY_COLOR[category],
            label=CATEGORY_LEGEND[category],
            markersize=9,
        )
        for category in CATEGORY_ORDER
    ]
    handles.extend(
        [
            Line2D(
                [0],
                [0],
                marker="s",
                color="none",
                markerfacecolor=TARGET_MARKER_RGB,
                label="Raw target-marker pixels (magenta)",
                markersize=9,
            ),
            Line2D(
                [0],
                [0],
                marker="s",
                color="none",
                markerfacecolor=REFERENCE_MARKER_RGB,
                label="Raw reference-marker pixels (cyan)",
                markersize=9,
            ),
            Line2D(
                [0],
                [0],
                color=THRESHOLD_COLOR,
                lw=2.2,
                ls="--",
                label=(
                    "Step 1 vertical: max reference MFI in high-confidence "
                    "target anchor"
                ),
            ),
            Line2D(
                [0],
                [0],
                color=THRESHOLD_COLOR,
                lw=2.2,
                ls="-",
                label=(
                    "Step 2 horizontal: RESTORE target divisor; not a "
                    "lower-confidence-target exclusion"
                ),
            ),
        ]
    )
    return handles


def _add_review_legend(legend_axis: plt.Axes):
    """Place the complete legend inside its dedicated non-data axis."""
    legend = legend_axis.legend(
        handles=_review_legend_handles(),
        loc="center",
        ncol=3,
        fontsize=PDF_LEGEND_FONT_SIZE,
        title="Figure legend",
        title_fontsize=PDF_LEGEND_TITLE_FONT_SIZE,
        frameon=True,
        facecolor="white",
        edgecolor="#999999",
        columnspacing=1.8,
        handletextpad=0.8,
        labelspacing=0.9,
        borderpad=0.9,
    )
    return legend


def _add_review_suptitle(
    fig: plt.Figure,
    target: str,
    reference: str,
):
    """Place the pair heading inside the reserved page-title band."""
    target_method = MARKER_INPUT_FLOOR_METHODS[target].replace("_", " ")
    reference_method = MARKER_INPUT_FLOOR_METHODS[reference].replace("_", " ")
    return fig.suptitle(
        f"RESTORE review: {target} <- {reference}\n"
        f"{PAIR_RATIONALE[(target, reference)]}\n"
        f"Raw arm floors: {target} = {target_method}; "
        f"{reference} = {reference_method}",
        fontsize=PDF_SUPTITLE_FONT_SIZE,
        fontweight="bold",
        y=1.0 - 0.12 / fig.get_figheight(),
        va="top",
        linespacing=1.25,
    )


def _write_qupath_review(
    projected_crop: pd.DataFrame,
    evaluation: dict,
    output: Path,
) -> dict:
    categories = projected_crop["review_category"].to_numpy(object)
    exported = projected_crop[
        [
            "object_id",
            "image",
            "x_canonical",
            "y_canonical",
            "input_qc_retained",
            "target_anchor",
            "target_component_supported",
            "retained_below_divisor",
        ]
    ].copy()
    exported["review_code"] = [
        CATEGORY_CODE[category] for category in categories
    ]
    exported["review_label"] = categories
    exported["target"] = evaluation["target"]
    exported["reference"] = evaluation["reference"]
    exported["pair_state"] = evaluation["state"]
    exported["target_input_floor_method"] = evaluation["target_input_floor_method"]
    exported["target_input_floor"] = evaluation["target_input_floor"]
    exported["target_input_floor_linear_otsu"] = evaluation[
        "target_input_floor_linear_otsu"
    ]
    exported["target_input_floor_linear_triangle"] = evaluation[
        "target_input_floor_linear_triangle"
    ]
    exported["reference_input_floor_method"] = evaluation[
        "reference_input_floor_method"
    ]
    exported["reference_input_floor"] = evaluation["reference_input_floor"]
    exported["reference_input_floor_linear_otsu"] = evaluation[
        "reference_input_floor_linear_otsu"
    ]
    exported["reference_input_floor_linear_triangle"] = evaluation[
        "reference_input_floor_linear_triangle"
    ]
    exported.to_csv(output, index=False)
    images = exported["image"].dropna().astype(str).unique()
    if len(images) != 1:
        raise ValueError(
            f"QuPath review export must contain exactly one image; found {images}"
        )
    return {
        "file": output.name,
        "donor": evaluation["donor"],
        "image": images[0],
        "target": evaluation["target"],
        "reference": evaluation["reference"],
        "target_input_floor_method": evaluation["target_input_floor_method"],
        "target_input_floor": evaluation["target_input_floor"],
        "target_input_floor_linear_otsu": evaluation[
            "target_input_floor_linear_otsu"
        ],
        "target_input_floor_linear_triangle": evaluation[
            "target_input_floor_linear_triangle"
        ],
        "reference_input_floor_method": evaluation["reference_input_floor_method"],
        "reference_input_floor": evaluation["reference_input_floor"],
        "reference_input_floor_linear_otsu": evaluation[
            "reference_input_floor_linear_otsu"
        ],
        "reference_input_floor_linear_triangle": evaluation[
            "reference_input_floor_linear_triangle"
        ],
        "model_sample_cell_n": int(evaluation["n_input"]),
        "canonical_crop_cell_n": int(len(projected_crop)),
        "export_cell_n": int(len(exported)),
        "all_canonical_crop_cells_exported": True,
        "category_counts": {
            category: int((categories == category).sum())
            for category in CATEGORY_ORDER
        },
    }


def _write_qupath_importer(output: Path, review_csv_dir: Path) -> None:
    """Render a ready-to-run importer into the review bundle."""
    source = (
        Path(__file__).resolve().parents[1]
        / "scripts"
        / "groovy"
        / QUPATH_IMPORTER_NAME
    )
    if not source.is_file():
        raise FileNotFoundError(f"QuPath importer template not found: {source}")
    template = source.read_text()
    if template.count(QUPATH_REVIEW_DIR_TOKEN) != 1:
        raise ValueError(
            "QuPath importer template must contain the review-directory token once"
        )
    escaped_dir = str(review_csv_dir.resolve()).replace("\\", "\\\\").replace('"', '\\"')
    output.write_text(template.replace(QUPATH_REVIEW_DIR_TOKEN, escaped_dir))


def _write_review_workflow(output: Path) -> None:
    """Copy the authoritative operating workflow into a review bundle."""
    source = (
        Path(__file__).resolve().parents[1]
        / "docs"
        / "RESTORE_PAIR_REVIEW_WORKFLOW.txt"
    )
    if not source.is_file():
        raise FileNotFoundError(f"RESTORE review workflow not found: {source}")
    output.write_text(source.read_text())


def build_review_checklist(
    queue: pd.DataFrame, export_manifest: Sequence[dict]
) -> pd.DataFrame:
    """Create the minimal human decision sheet from the technical review queue."""
    required = {
        "donor",
        "target",
        "reference",
        "selection_reasons",
        "state",
        "state_reason",
    }
    missing = required.difference(queue.columns)
    if missing:
        raise ValueError(f"review queue is missing columns: {sorted(missing)}")
    export_index = pd.DataFrame(export_manifest)
    export_required = {"donor", "target", "reference", "image", "file"}
    export_missing = export_required.difference(export_index.columns)
    if export_missing:
        raise ValueError(
            f"QuPath export manifest is missing columns: {sorted(export_missing)}"
        )

    checklist = queue[list(required)].copy()
    checklist["donor"] = checklist["donor"].astype(str)
    export_index = export_index[list(export_required)].copy()
    export_index["donor"] = export_index["donor"].astype(str)
    checklist = checklist.merge(
        export_index,
        on=["donor", "target", "reference"],
        how="left",
        validate="one_to_one",
    )
    if checklist[["image", "file"]].isna().any().any():
        raise ValueError("review checklist is missing a matching QuPath export")

    pair_order = {pair: index + 1 for index, pair in enumerate(SHORTLIST_PAIRS)}
    checklist["pair_order"] = [
        pair_order[(target, reference)]
        for target, reference in checklist[
            ["target", "reference"]
        ].itertuples(index=False, name=None)
    ]
    checklist = checklist.sort_values(
        ["pair_order", "donor"], kind="stable"
    ).reset_index(drop=True)
    checklist["pdf_row"] = (
        checklist.groupby("pair_order", sort=False).cumcount() + 1
    )
    checklist["review_order"] = np.arange(1, len(checklist) + 1)
    checklist["pair"] = (
        checklist["target"].astype(str)
        + " <- "
        + checklist["reference"].astype(str)
    )
    checklist["pdf_file"] = [
        f"spatial_review/{target}_from_{reference}.pdf"
        for target, reference in checklist[
            ["target", "reference"]
        ].itertuples(index=False, name=None)
    ]
    checklist = checklist.rename(
        columns={
            "file": "qupath_overlay",
            "state": "computational_state",
            "state_reason": "computational_reason",
            "selection_reasons": "why_selected",
        }
    )
    for column in (
        "step1_target_separation",
        "step2_reference_controls",
        "spatial_coherence",
        "donor_disposition",
        "notes",
    ):
        checklist[column] = ""
    return checklist[
        [
            "review_order",
            "pair_order",
            "pdf_row",
            "pair",
            "donor",
            "image",
            "why_selected",
            "computational_state",
            "computational_reason",
            "pdf_file",
            "qupath_overlay",
            "step1_target_separation",
            "step2_reference_controls",
            "spatial_coherence",
            "donor_disposition",
            "notes",
        ]
    ]


def build_review_bundle(
    sample_root: Path,
    full_compartment_cells: Path,
    raw_cells: Path,
    redsea_cells: Path,
    images_dir: Path,
    results_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    """Build pair summaries, review queue, spatial PDFs, and QuPath CSVs atomically."""
    results_dir = Path(results_dir).resolve()
    full_compartment_cells = Path(full_compartment_cells).resolve()
    images_dir = Path(images_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if not full_compartment_cells.is_dir():
        raise FileNotFoundError(
            "full-cell compartment directory not found: "
            f"{full_compartment_cells}"
        )
    if not images_dir.is_dir():
        raise FileNotFoundError(f"QPTIFF image directory not found: {images_dir}")
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    required = [
        results_dir / "pair_evaluations.csv",
        results_dir / "method_spec.json",
        results_dir / "expression_reviews.csv",
        results_dir / "pair_reviews.csv",
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(f"expanded-screen artifact not found: {path}")

    evaluations = pd.read_csv(
        results_dir / "pair_evaluations.csv", dtype={"donor": str}
    )
    ensure_eligible_donors(
        evaluations["donor"].unique(),
        context="RESTORE review bundle",
    )
    specification = json.loads((results_dir / "method_spec.json").read_text())
    if specification.get("method_version") != METHOD_VERSION:
        raise ValueError(
            "expanded-screen method version does not match the current review method"
        )
    evaluation_versions = set(evaluations["method_version"].dropna())
    if evaluation_versions != {METHOD_VERSION}:
        raise ValueError(
            f"pair-evaluation versions do not match {METHOD_VERSION}: "
            f"{sorted(evaluation_versions)}"
        )
    config = _config_from_specification(specification)
    expression_reviews, _ = load_expression_reviews(
        results_dir / "expression_reviews.csv"
    )
    pair_reviews, _ = load_pair_reviews(results_dir / "pair_reviews.csv")
    pair_summary = summarize_pairs(evaluations)
    queue = select_review_queue(evaluations)

    temp_dir.mkdir(parents=True)
    figure_dir = temp_dir / "spatial_review"
    qupath_dir = temp_dir / "qupath_review"
    figure_dir.mkdir()
    qupath_dir.mkdir()
    _write_qupath_importer(
        temp_dir / QUPATH_IMPORTER_NAME,
        output_dir / "qupath_review",
    )
    _write_review_workflow(temp_dir / REVIEW_WORKFLOW_NAME)
    donor_cache: dict[str, pd.DataFrame] = {}
    full_donor_cache: dict[str, pd.DataFrame] = {}
    export_manifest = []
    try:
        pair_summary.to_csv(temp_dir / "pair_summary.csv", index=False)
        queue.to_csv(temp_dir / "review_queue.csv", index=False)
        for target, reference in SHORTLIST_PAIRS:
            pair_queue = queue[
                (queue["target"] == target)
                & (queue["reference"] == reference)
            ]
            fig, axes, legend_axis = _create_review_figure(len(pair_queue))
            for row_index, row in enumerate(pair_queue.itertuples(index=False)):
                donor = str(row.donor)
                if donor not in donor_cache:
                    donor_cache[donor], _ = load_validation_sample(
                        Path(sample_root),
                        Path(raw_cells),
                        Path(redsea_cells),
                        donor,
                    )
                donor_df = donor_cache[donor]
                evaluation = evaluate_locked_pair(
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
                categories = review_categories(evaluation, len(donor_df))
                zoom_x, zoom_y = representative_zoom_bounds(
                    donor_df, categories
                )
                marker_zoom = load_marker_zoom(
                    images_dir,
                    donor_df,
                    target,
                    reference,
                    zoom_x,
                    zoom_y,
                )
                if donor not in full_donor_cache:
                    full_donor_df, full_audit = load_validation_sample(
                        full_compartment_cells,
                        Path(raw_cells),
                        Path(redsea_cells),
                        donor,
                    )
                    if full_audit["n"] != full_audit["canonical_n"]:
                        raise ValueError(
                            f"donor {donor}: full-cell review source contains "
                            f"{full_audit['n']:,} of "
                            f"{full_audit['canonical_n']:,} canonical cells"
                        )
                    full_donor_cache[donor] = full_donor_df
                projected_crop = project_review_crop(
                    full_donor_cache[donor],
                    evaluation,
                    marker_zoom.x_bounds_um,
                    marker_zoom.y_bounds_um,
                )
                sample_overlap_n = _validate_projection_overlap(
                    donor_df,
                    categories,
                    projected_crop,
                    marker_zoom.x_bounds_um,
                    marker_zoom.y_bounds_um,
                )
                _draw_review_row(
                    axes[row_index],
                    donor_df,
                    evaluation,
                    categories,
                    row.selection_reasons,
                    marker_zoom,
                    projected_crop,
                )
                stem = f"{target}_from_{reference}__{donor}".lower()
                export_record = _write_qupath_review(
                    projected_crop,
                    evaluation,
                    qupath_dir / f"{stem}.csv",
                )
                export_record["marker_zoom"] = _marker_zoom_metadata(marker_zoom)
                export_record["model_sample_overlap_n"] = sample_overlap_n
                export_record["model_sample_category_match"] = True
                export_manifest.append(export_record)
            axes[-1, 0].set_xlabel(
                f"{reference} scaled REDSEA MFI",
                fontsize=PDF_AXIS_LABEL_FONT_SIZE,
                labelpad=9,
            )
            _add_review_legend(legend_axis)
            _add_review_suptitle(fig, target, reference)
            fig.savefig(
                figure_dir / f"{target}_from_{reference}.pdf",
                dpi=160,
            )
            plt.close(fig)
        review_checklist = build_review_checklist(queue, export_manifest)
        review_checklist.to_csv(temp_dir / REVIEW_CHECKLIST_NAME, index=False)
        manifest = {
            "method_version": METHOD_VERSION,
            "cohort_exclusions": DONOR_EXCLUSIONS,
            "source_results": str(results_dir),
            "sample_root": str(Path(sample_root).resolve()),
            "full_compartment_cells": str(full_compartment_cells),
            "images_dir": str(images_dir),
            "donor_status_policy": (
                "Donor 6579 is the pancreatitis donor and remains excluded by "
                "maintainer decision. Donor 6476 is not a pancreatitis donor and "
                "receives no fixed stress-case selection."
            ),
            "marker_input_floor_methods": MARKER_INPUT_FLOOR_METHODS,
            "input_floor_review_policy": (
                "CD3e and CD20 use the donor-local maximum of raw-linear Otsu and "
                "raw-linear Triangle. This requires both criteria to classify an arm "
                "cell as high. Other screened markers use their explicitly configured "
                "raw-linear Otsu policy. Component statistics and final floors are "
                "recorded in the review queue and QuPath export manifest."
            ),
            "selection_policy": (
                "For each fixed shortlisted pair: lowest target/reference fold, lowest "
                "seed stability, highest double-high fraction, largest sampled/full-"
                "maximum divergence, median target fold, plus confirmed-absence and "
                "sparse-arm cases where applicable. Selection is for stress review only "
                "and is not a cell-frequency acceptance rule."
            ),
            "threshold_display_policy": (
                "The dashed vertical green line is Step 1: the full maximum REDSEA "
                "reference intensity in the high-confidence target anchor. Every input-"
                "QC-retained NNMF target-component cell at or left of that separator is "
                "retained as target lineage evidence, even if it lies below Step 2. "
                "Cells below both raw floors remain double-low. Only reference-"
                "component cells that are raw "
                "reference-high/target-low and strictly right of Step 1 may set Step 2, "
                "the solid horizontal green RESTORE target divisor. Step 2 is not an "
                "exclusion boundary for reference-negative target-component cells."
            ),
            "spatial_zoom_policy": (
                "Each row uses one deterministic representative QPTIFF crop in four "
                "aligned panels: raw target-marker pixels, raw reference-marker pixels, "
                "the RESTORE classification overlay, and both marker channels merged "
                "with classification rings. The overlay panel includes the whole tissue "
                "and crop rectangle as an untitled compact inset. The crop preserves the "
                "area of a 10%-per-axis tissue window but uses a 2.2:1 landscape aspect so equal "
                "X/Y spatial scaling fills each panel without distortion. The zoom is "
                "centered on the densest combined target-"
                "anchor/retained-target-component bin, falling back in order to double-"
                "high, reference-control, other-retained, then any finite-coordinate "
                "cells. Every canonical cell whose centroid lies in the resulting crop "
                "is labeled in the overlay and merged panels; there is no display "
                "subsampling. Crop selection itself remains a display rule only."
            ),
            "category_calculation_order": (
                "Required finite measurements are checked first. Valid cells below both "
                "raw marker floors are assigned Double-low before component support is "
                "considered. Only input-QC-retained cells may then become Double-high, "
                "lower-confidence target, high-confidence target anchor, or reference "
                "control; later categories cannot overwrite Double-low."
            ),
            "full_cell_projection_policy": (
                "RESTORE fitting, arm floors, separators, and displayed expression-cloud "
                "metrics remain locked to the deterministic modeling sample. For visual "
                "review only, its fixed robust feature bounds and NNMF components are "
                "projected onto every canonical cell in the representative crop using "
                "full QuPath membrane/cell measurements joined by object_id to canonical "
                "raw and REDSEA data. Every modeling-sample cell in each crop must retain "
                "exactly the same category after projection or bundle generation fails."
            ),
            "marker_image_policy": (
                "Only QPTIFF tiles overlapping the physical crop are decoded. The first "
                "tiled pyramid level producing at most 1000 x 600 crop pixels is used. "
                "If no tiled level reaches that target, the deepest tiled level is used. "
                "Each channel is independently windowed from the 1st to 99.8th percentile "
                "of positive crop pixels and displayed with gamma 0.65. Target pixels are "
                "magenta and reference pixels cyan; this contrast transform affects only "
                "display, never RESTORE fitting or calls. Exact source, level, crop, and "
                "display windows are recorded per row in qupath_exports.marker_zoom."
            ),
            "figure_layout": (
                "expression_cloud_plus_all_cell_target_reference_overlay_and_merged_closeup"
            ),
            "pdf_typography": {
                "suptitle_pt": PDF_SUPTITLE_FONT_SIZE,
                "row_title_pt": PDF_ROW_TITLE_FONT_SIZE,
                "axis_label_pt": PDF_AXIS_LABEL_FONT_SIZE,
                "tick_pt": PDF_TICK_FONT_SIZE,
                "image_title_pt": PDF_IMAGE_TITLE_FONT_SIZE,
                "image_label_pt": PDF_IMAGE_LABEL_FONT_SIZE,
                "image_tick_pt": PDF_IMAGE_TICK_FONT_SIZE,
                "row_summary_pt": PDF_ROW_TITLE_FONT_SIZE,
                "row_summary_placement": (
                    "scatter title area above the data axis; no in-data text box"
                ),
                "title_band_in": PDF_TITLE_HEIGHT_IN,
                "legend_pt": PDF_LEGEND_FONT_SIZE,
                "legend_placement": "dedicated non-data axis below all plot rows",
            },
            "target_review_colors": {
                "high_confidence_anchor": CATEGORY_COLOR[
                    "High-confidence target anchor (sets Step 1)"
                ],
                "lower_confidence_component": CATEGORY_COLOR[
                    "Lower-confidence target component (reference-negative; retained)"
                ],
            },
            "marker_pixel_colors": {
                "target_rgb": TARGET_MARKER_RGB,
                "reference_rgb": REFERENCE_MARKER_RGB,
            },
            "workflow": REVIEW_WORKFLOW_NAME,
            "decision_checklist": REVIEW_CHECKLIST_NAME,
            "operating_workflow": (
                "Review one pair PDF at a time and record only in "
                f"{REVIEW_CHECKLIST_NAME}. Use QuPath only for rows that need closer "
                "inspection; the importer presents pairs valid for the open image."
            ),
            "shortlist": [
                {
                    "target": target,
                    "reference": reference,
                    "rationale": PAIR_RATIONALE[(target, reference)],
                }
                for target, reference in SHORTLIST_PAIRS
            ],
            "qupath_exports": export_manifest,
            "qupath_importer": QUPATH_IMPORTER_NAME,
            "qupath_import_policy": (
                "The importer detects the open image and presents valid target <- "
                "reference overlays in a dialog; no filename editing is required. "
                "Before applying temporary classes, it stores every original PathClass "
                "in per-object metadata. Every canonical cell inside the representative "
                "crop is included with no category cap; detections outside the crop keep "
                "their existing classes. Rerunning the importer may replace prior "
                "RESTORE review classes on the imported crop. The Clear dialog action "
                "restores every saved original PathClass exactly, including Unclassified."
            ),
            "acceptance_note": (
                "These are review labels, not RESTORE calls or phenotypes. No pair was "
                "accepted by generating this bundle."
            ),
        }
        (temp_dir / "review_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n"
        )
        temp_dir.rename(output_dir)
    except Exception:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise
    return {
        "summary": output_dir / "pair_summary.csv",
        "queue": output_dir / "review_queue.csv",
        "workflow": output_dir / REVIEW_WORKFLOW_NAME,
        "checklist": output_dir / REVIEW_CHECKLIST_NAME,
        "figures": output_dir / "spatial_review",
        "qupath": output_dir / "qupath_review",
    }


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--full-compartment-cells", type=Path, required=True)
    parser.add_argument("--raw-cells", type=Path, required=True)
    parser.add_argument("--redsea-cells", type=Path, required=True)
    parser.add_argument("--images-dir", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    outputs = build_review_bundle(
        args.sample,
        args.full_compartment_cells,
        args.raw_cells,
        args.redsea_cells,
        args.images_dir,
        args.results,
        args.output,
    )
    for name, path in outputs.items():
        print(f"[review] {name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
