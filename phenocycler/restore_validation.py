#!/usr/bin/env python3
"""Manuscript-grounded validation helpers for RESTORE marker pairs.

This module does not write normalized values or phenotype calls. It extracts a
deterministic sample of raw QuPath compartment means and evaluates candidate
cell-group separators before a RESTORE implementation is accepted.
"""

from __future__ import annotations

import argparse
import json
import shutil
import warnings
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .cohort import DONOR_EXCLUSIONS, ensure_eligible_donors


# v7: the positive CALL is divisor-gated only (target > Step 2 divisor). The former target-component
# rescue -- which also called every NNMF-target-component cell left of Step 1 even below the divisor --
# is removed from the call: negligible for clean markers (<=0.3 pp) but it doubled a high-background
# marker like CD11b by retaining background cells the divisor correctly excludes (6521 CD11b<-PanCK
# 14% -> 28% supported). Below-divisor component cells stay tracked (target_supported_below_divisor)
# and are shown in the bundle as a distinct "retained, not called" category. Adds threshold_positive_n.
# v6: robust min-max tails adapt when the balanced fit sample cannot span the configured pair, so a
# single bright cell in a sparse arm no longer sets the feature scale (22 of 600 screened rows, all with
# immune targets). Adds anchor_recovery, arm-saturation/stability_informative reporting, and a
# saturated-arm probe Jaccard. v5 artifacts were produced under the degenerate bound rule.
METHOD_VERSION = "restore-pair-validation-v8"
# NB: ``BroadLineageRevised.md`` lists 8 first-pass pairs.  ``CD20 <- E_cadherin`` was the eighth and
# is deliberately absent — CD20 is deferred out of RESTORE entirely (see RESTORE_EXCLUDED_MARKERS).
BASELINE_PAIRS: tuple[tuple[str, str], ...] = (
    ("E_cadherin", "CD31"),
    ("CD31", "E_cadherin"),
    ("CD3e", "E_cadherin"),
    ("CD68", "E_cadherin"),
    ("CD11b", "E_cadherin"),
    ("B3TUBB", "EpCAM"),
    ("Vimentin", "E_cadherin"),
)
# Pooled-immune comparator marker set: the minimum union that nets the immune populations present with
# no CD45 in the panel. Used ONLY by evaluate_immune_joint, never by the locked pairwise path.
# CD20 is retained HERE ONLY, as a diagnostic co-factor of the joint fit — dropping it would change the
# fit for the other three and make the comparator non-comparable with the run already on record. The
# comparator emits no divisor and no call, and its finding that pooling lifts CD20's worst-case anchor
# recovery (0.217 -> 0.968) does NOT license reintroducing CD20: it is excluded from RESTORE by
# maintainer decision on separate grounds (RESTORE_EXCLUDED_MARKERS), which the guard above enforces.
IMMUNE_JOINT_MARKERS: tuple[str, ...] = ("CD3e", "CD20", "CD68", "CD11b")
# Markers held out of the RESTORE pair web in BOTH roles — neither target nor reference.  This is
# stronger than "cannot be a reference": such a marker gets no threshold, no divisor and no positivity
# call from this method, and is deferred to a validated second pass once the broad cell types are
# settled.  Its compartment means are still extracted (DEFERRED_EXTRACTION_MARKERS below) so that pass
# needs no re-extraction, and its input-floor policy entry is retained for the same reason.
RESTORE_EXCLUDED_MARKERS: dict[str, str] = {
    "CD20": (
        "B cells are sparse in pancreas even with inflammation or disease, and CD20 is the weakest "
        "channel in the panel after CD163 (display ceiling median ~503 counts vs ~6,800 for "
        "E-cadherin). As a REFERENCE its positive cells would BE the negative control, and its "
        "exclusive arm ran 41-2,026 cells per donor against 2,206-50,000 for the structural markers, "
        "producing inverted fits (6591 CD11b <- CD20: 50 reference cells thresholding 30,507 target "
        "cells). As a TARGET its Step 1 anchor bottomed out at 14 cells with anchor recovery down to "
        "0.217, and reviewer inspection of the CD20 review PDFs found overcalling throughout with "
        "both qptiff channels stretched from near-background into block-like tile artifacts. "
        "Maintainer decision 2026-07-23: CD20 takes no part in RESTORE; add it after the broad cell "
        "types are worked out."
    ),
}
IMMUNE_REFERENCE_SCREEN_PAIRS: tuple[tuple[str, str], ...] = (
    ("CD3e", "CD68"),
    ("CD3e", "CD163"),
    ("CD68", "CD3e"),
    ("CD11b", "CD3e"),
)
EPITHELIAL_REFERENCE_COMPARATOR_PAIRS: tuple[tuple[str, str], ...] = tuple(
    (target, reference)
    for target in ("CD3e", "CD68", "CD11b")
    for reference in ("Pan_Cytokeratin", "Ker8_18", "EpCAM")
)
CANDIDATE_PAIRS: tuple[tuple[str, str], ...] = tuple(
    dict.fromkeys(
        BASELINE_PAIRS
        + IMMUNE_REFERENCE_SCREEN_PAIRS
        + EPITHELIAL_REFERENCE_COMPARATOR_PAIRS
    )
)
for _target, _reference in CANDIDATE_PAIRS:  # fail at import, not mid-screen
    for _marker, _role in ((_target, "target"), (_reference, "reference")):
        if _marker in RESTORE_EXCLUDED_MARKERS:
            raise ValueError(
                f"{_target} <- {_reference}: {_marker} is excluded from RESTORE in both roles "
                f"(here as the {_role}). {RESTORE_EXCLUDED_MARKERS[_marker]}"
            )
    if _target == _reference:
        raise ValueError(f"{_target} cannot be its own reference")
del _target, _reference, _marker, _role
# The Gate-2-FROZEN production pairs (maintainer sign-off 2026-07-24): one predeclared reference per
# target under the reference>=target frequency rule. Identical to BASELINE_PAIRS except CD11b takes
# EpCAM (a maintainer override of the marginal E_cadherin baseline; see
# reference_selection.ACCEPTED_REFERENCE_OVERRIDES). E_cadherin<-CD31 is retained despite failing the
# proportion because endothelium is the only biologically exclusive epithelial reference (advisory for
# a structural target). Used by the Gate-3 full-cell pilot and Gate-5 all-donor validation. It is NOT
# yet wired into lineage.py/config.py -- that is Gate 4.
ACCEPTED_PAIRS: tuple[tuple[str, str], ...] = (
    ("E_cadherin", "CD31"),
    ("CD31", "E_cadherin"),
    ("CD3e", "E_cadherin"),
    ("CD68", "E_cadherin"),
    ("CD11b", "EpCAM"),
    ("B3TUBB", "EpCAM"),
    ("Vimentin", "E_cadherin"),
)
for _t, _r in ACCEPTED_PAIRS:  # fail at import if the frozen set drifts from the candidate web
    if (_t, _r) not in set(CANDIDATE_PAIRS):
        raise ValueError(f"accepted pair {_t} <- {_r} is not in CANDIDATE_PAIRS")
if len({_t for _t, _ in ACCEPTED_PAIRS}) != len(ACCEPTED_PAIRS):
    raise ValueError("ACCEPTED_PAIRS must carry exactly one reference per target")
del _t, _r
PAIR_SETS: dict[str, tuple[tuple[str, str], ...]] = {
    "baseline": BASELINE_PAIRS,
    "immune-screen": IMMUNE_REFERENCE_SCREEN_PAIRS,
    "accepted": ACCEPTED_PAIRS,
    "all": CANDIDATE_PAIRS,
}
COMPARTMENTS: tuple[str, ...] = ("Nucleus", "Cytoplasm", "Membrane", "Cell")
LOCKED_FEATURE_COMPARTMENTS: tuple[str, ...] = ("Membrane", "Cell")
QUPATH_MARKER_ALIASES = {
    "E_cadherin": "E-cadherin",
    "Pan_Cytokeratin": "Pan-Cytokeratin",
    "Ker8_18": "Ker8-18",
}
# Markers extracted from QuPath but NOT screened: they take no part in the pair web, yet their
# compartment means ride along so the deferred second pass (and any re-extraction) needs no rerun of
# the expensive QuPath export.  Keeping CD20 here also keeps the extracted sample schema-identical to
# the files already under qupath_compartment_*_20donor/.
DEFERRED_EXTRACTION_MARKERS: tuple[str, ...] = tuple(RESTORE_EXCLUDED_MARKERS)
QUPATH_MARKERS: dict[str, str] = {
    marker: QUPATH_MARKER_ALIASES.get(marker, marker)
    for marker in [m for pair in CANDIDATE_PAIRS for m in pair] + list(DEFERRED_EXTRACTION_MARKERS)
}
LINEAR_OTSU_FLOOR = "linear_otsu"
MAX_LINEAR_OTSU_TRIANGLE_FLOOR = "max_linear_otsu_triangle"
MARKER_INPUT_FLOOR_METHODS: dict[str, str] = {
    "E_cadherin": LINEAR_OTSU_FLOOR,
    "CD31": LINEAR_OTSU_FLOOR,
    "CD3e": MAX_LINEAR_OTSU_TRIANGLE_FLOOR,
    "CD20": MAX_LINEAR_OTSU_TRIANGLE_FLOOR,
    "CD68": LINEAR_OTSU_FLOOR,
    "CD11b": LINEAR_OTSU_FLOOR,
    "B3TUBB": LINEAR_OTSU_FLOOR,
    "EpCAM": LINEAR_OTSU_FLOOR,
    "Vimentin": LINEAR_OTSU_FLOOR,
    "CD163": LINEAR_OTSU_FLOOR,
    "Pan_Cytokeratin": LINEAR_OTSU_FLOOR,
    "Ker8_18": LINEAR_OTSU_FLOOR,
}


class PairState(str, Enum):
    """Auditable donor/pair disposition; only VALID permits RESTORE calls."""

    VALID = "VALID"
    NO_TARGET_EXPRESSION_CONFIRMED = "NO_TARGET_EXPRESSION_CONFIRMED"
    LOW_SIGNAL_UNRESOLVED = "LOW_SIGNAL_UNRESOLVED"
    NO_REFERENCE_POPULATION = "NO_REFERENCE_POPULATION"
    INVALID_PAIR = "INVALID_PAIR"
    MODEL_UNSTABLE = "MODEL_UNSTABLE"
    PAIR_REVIEW_PENDING = "PAIR_REVIEW_PENDING"
    TECHNICAL_FAILURE = "TECHNICAL_FAILURE"


class ExpressionReview(str, Enum):
    """Independent image-review evidence used to interpret weak target signal."""

    UNREVIEWED = "UNREVIEWED"
    CONFIRMED_PRESENT = "CONFIRMED_PRESENT"
    CONFIRMED_ABSENT = "CONFIRMED_ABSENT"


class PairReview(str, Enum):
    """Cohort-level biological review of a fixed target/reference rule."""

    UNREVIEWED = "UNREVIEWED"
    ACCEPTED = "ACCEPTED"
    REJECTED = "REJECTED"


@dataclass(frozen=True)
class PairValidationConfig:
    """Frozen parameters for the locked manuscript-grounded evaluation path."""

    feature_compartments: tuple[str, ...] = LOCKED_FEATURE_COMPARTMENTS
    n_components: int = 2
    seeds: tuple[int, ...] = (0, 1, 2)
    fit_cap: int = 50_000
    input_qc_bins: int = 512
    input_qc_upper_quantile: float = 0.995
    input_qc_triangle_bins: int = 256
    scale_lower_quantile: float = 0.005
    scale_upper_quantile: float = 0.995
    threshold_sample_sizes: tuple[int, ...] = (50, 100, 250, 500, 1_000, 2_000)
    threshold_resamples: int = 200
    min_target_arm_fold: float = 2.0
    min_reference_arm_fold: float = 2.0
    min_arm_n: int = 2
    min_control_jaccard: float = 0.90
    # Reference-frequency selection rule (2026-07-24): a target's negative control (the reference-
    # positive, target-negative arm) must be at least as large as the target-positive population, or
    # the background maximum is estimated from too few cells to threshold many (the inverted-fit failure
    # mode, e.g. CD11b <- CD3e). Reported per pair as reference_to_target_ratio = reference_n / target_n
    # and the boolean reference_undersized = ratio < this value. PURE PROPORTION -- no absolute cell-
    # count floor, because images vary widely in total and per-marker counts so an absolute threshold is
    # not comparable across donors. DIAGNOSTIC ONLY: it never gates PairState or alters a divisor; the
    # reviewer uses it and records any reject as a human decision (so E_cadherin <- CD31 stays
    # reviewable rather than auto-failed).
    min_reference_to_target_ratio: float = 1.0
    # Minimum cells each robust min-max tail must span for the balanced fit sample to estimate the
    # scaling bounds. Below it the quantile degenerates into an order statistic of a handful of cells
    # (a sparse-marker arm), so bounds fall back to the input-QC-retained population instead.
    min_quantile_support_cells: int = 1
    # When the configured tails are degenerate, widen them on the same balanced sample until each spans
    # this many cells (capped at the quartile). Chosen from the 22 affected pairs: 10 gave median anchor
    # recovery 0.976 (min 0.787) against 0.938/0.730 at 5, while preserving the inverted pairs at 1.00.
    degenerate_tail_support_cells: int = 10
    # Arm shrinkage used ONLY by the saturated-arm stability probe (reported, never gated).
    stability_probe_fraction: float = 0.8

    def __post_init__(self) -> None:
        if not self.feature_compartments:
            raise ValueError("at least one feature compartment is required")
        unknown = set(self.feature_compartments).difference(COMPARTMENTS)
        if unknown:
            raise ValueError(f"unknown feature compartments: {sorted(unknown)}")
        if self.n_components != 2:
            raise ValueError("the locked pair evaluation requires two NNMF components")
        if not self.seeds or len(set(self.seeds)) != len(self.seeds):
            raise ValueError("seeds must be a nonempty unique sequence")
        if self.fit_cap < 2 * self.min_arm_n:
            raise ValueError("fit_cap is too small for balanced arm fitting")
        if self.input_qc_bins < 2:
            raise ValueError("input_qc_bins must be >= 2")
        if self.input_qc_triangle_bins < 2:
            raise ValueError("input_qc_triangle_bins must be >= 2")
        if not 0 < self.input_qc_upper_quantile <= 1:
            raise ValueError("input_qc_upper_quantile must be in (0, 1]")
        if not (
            0
            <= self.scale_lower_quantile
            < self.scale_upper_quantile
            <= 1
        ):
            raise ValueError("invalid robust min-max quantiles")
        if (
            not self.threshold_sample_sizes
            or any(size < 1 for size in self.threshold_sample_sizes)
            or tuple(sorted(set(self.threshold_sample_sizes)))
            != self.threshold_sample_sizes
        ):
            raise ValueError(
                "threshold_sample_sizes must be unique, positive, and increasing"
            )
        if self.threshold_resamples < 1:
            raise ValueError("threshold_resamples must be positive")
        if self.min_target_arm_fold <= 1:
            raise ValueError("min_target_arm_fold must be > 1")
        if self.min_reference_arm_fold <= 1:
            raise ValueError("min_reference_arm_fold must be > 1")
        if self.min_arm_n < 2:
            raise ValueError("min_arm_n must be >= 2")
        if not 0 <= self.min_control_jaccard <= 1:
            raise ValueError("min_control_jaccard must be in [0, 1]")
        if self.min_reference_to_target_ratio <= 0:
            raise ValueError("min_reference_to_target_ratio must be > 0")
        if self.min_quantile_support_cells < 1:
            raise ValueError("min_quantile_support_cells must be >= 1")
        if self.degenerate_tail_support_cells < self.min_quantile_support_cells:
            raise ValueError(
                "degenerate_tail_support_cells must be >= min_quantile_support_cells"
            )
        if not 0 < self.stability_probe_fraction <= 1:
            raise ValueError("stability_probe_fraction must be in (0, 1]")

    def specification(self) -> dict:
        """Return the complete machine-readable method and data-role contract."""
        parameters = asdict(self)
        for key in ("feature_compartments", "seeds", "threshold_sample_sizes"):
            parameters[key] = list(parameters[key])
        return {
            "method_version": METHOD_VERSION,
            "parameters": parameters,
            "data_roles": {
                "group_inference": (
                    "raw linear QuPath compartment means from canonical object IDs"
                ),
                "divisor_and_calls": (
                    "matched nonnegative REDSEA whole-cell target intensities"
                ),
                "raw_whole_cell": (
                    "input-QC arm discovery, matching, and raw-to-REDSEA diagnostics"
                ),
            },
            "marker_input_floor_methods": MARKER_INPUT_FLOOR_METHODS,
            "input_floor_policy": (
                "Raw arm floors are donor-local and explicitly marker-specific. CD3e "
                "and CD20 use max(raw-linear Otsu, raw-linear Triangle), requiring both "
                "statistics to classify a cell as high. All other currently screened "
                "markers use raw-linear Otsu pending their own biological review. These "
                "floors define provisional fitting arms, the Step 1 anchor, and the "
                "input-QC-retained set eligible for target-component support. Cells "
                "below both floors remain double-low. These floors are not RESTORE "
                "divisors or production positivity thresholds."
            ),
            "fitting_policy": (
                "Explicit marker-specific raw-linear floors define provisional mutually "
                "exclusive arms and remove only double-low cells from fitting. NNMF is fit "
                "to equal-size provisional target-only and reference-only arms after "
                "robust per-feature min-max scaling, then every valid canonical cell is "
                "assigned."
            ),
            "ordered_control_policy": (
                "Control inference is directional and sequential. First, intersect the "
                "NNMF target component with the raw target-high/reference-low arm; this "
                "high-confidence anchor alone sets the vertical reference separator from "
                "its full maximum REDSEA reference intensity. Every input-QC-retained "
                "NNMF target-component cell at or left of that separator is then retained "
                "as target lineage evidence without requiring anchor membership. Cells "
                "below both raw floors remain double-low and cannot be rescued by NNMF "
                "assignment. Second, select NNMF reference-component "
                "cells that are raw reference-high/target-low and strictly right of the "
                "vertical separator. Only this refined reference-positive/target-negative "
                "population sets the horizontal target divisor. Double-high, double-low, "
                "and other assigned cells never set either maximum."
            ),
            "divisor_policy": (
                "The canonical RESTORE divisor is the full maximum REDSEA target "
                "intensity in the ordered reference-positive/target-negative control "
                "population selected after the vertical reference baseline. Fixed-N "
                "sampled maxima and mean+3SD are diagnostics only."
            ),
            "calling_policy": (
                "normalized = REDSEA value / divisor. Lineage evidence is positive when "
                "normalized > 1 OR the cell is input-QC-retained, is in the NNMF target "
                "component, and its matched REDSEA reference intensity is at or left of "
                "the Step 1 separator. The latter override retains reference-negative "
                "target-component cells even below the horizontal divisor; double-low "
                "cells are never rescued. The horizontal line remains the RESTORE "
                "normalization divisor, not an exclusion boundary."
            ),
            "cell_availability_policy": (
                "A canonical cell missing any required finite nonnegative measurement is "
                "recorded as unavailable and receives no marker call. Cell-level "
                "unavailability does not invalidate an otherwise estimable image-level pair."
            ),
            "expression_support_policy": (
                "Configured target- and reference-arm folds are provisional donor-local "
                "SBR diagnostics. Weak target support cannot confirm biological absence "
                "without independent image review; weak reference support makes the pair "
                "unusable for that image."
            ),
            "pair_acceptance_policy": (
                "Computational separation cannot make a pair VALID. The fixed pair must "
                "also have an explicit ACCEPTED biological review after cohort validation."
            ),
            "fallback_policy": (
                "No cohort threshold, adaptive counterpart, or success-shaped fallback."
            ),
        }


@dataclass(frozen=True)
class PairAssessment:
    state: PairState
    reason: str


def assess_pair_state(
    *,
    reference_n: int,
    target_n: int,
    distinct: bool,
    directional_mean: bool,
    directional_median: bool,
    converged: bool,
    stable: bool,
    target_arm_fold: float,
    reference_arm_fold: float = np.inf,
    expression_review: ExpressionReview = ExpressionReview.UNREVIEWED,
    pair_review: PairReview = PairReview.UNREVIEWED,
    config: PairValidationConfig | None = None,
    technical_error: str | None = None,
) -> PairAssessment:
    """Classify one donor/pair without turning uncertainty into negativity."""
    config = config or PairValidationConfig()
    expression_review = ExpressionReview(expression_review)
    pair_review = PairReview(pair_review)
    if technical_error:
        return PairAssessment(PairState.TECHNICAL_FAILURE, technical_error)
    if target_n < config.min_arm_n:
        if expression_review is ExpressionReview.CONFIRMED_ABSENT:
            return PairAssessment(
                PairState.NO_TARGET_EXPRESSION_CONFIRMED,
                "target arm is absent and independent image review confirms absence",
            )
        if expression_review is ExpressionReview.CONFIRMED_PRESENT:
            return PairAssessment(
                PairState.INVALID_PAIR,
                "the pair failed to recover independently confirmed target expression",
            )
        return PairAssessment(
            PairState.LOW_SIGNAL_UNRESOLVED,
            "target arm is too small to establish donor-local expression",
        )
    if not distinct or not directional_mean or not directional_median:
        return PairAssessment(
            PairState.INVALID_PAIR,
            "NNMF groups do not form distinct reference and target arms by mean and median",
        )
    if pair_review is PairReview.REJECTED:
        return PairAssessment(
            PairState.INVALID_PAIR,
            "the fixed target/reference rule was rejected by biological review",
        )
    if reference_n < config.min_arm_n:
        return PairAssessment(
            PairState.NO_REFERENCE_POPULATION,
            (
                f"ordered reference control has {reference_n} cells after the vertical "
                f"separator; requires {config.min_arm_n}"
            ),
        )
    if (
        not np.isfinite(reference_arm_fold)
        or reference_arm_fold < config.min_reference_arm_fold
    ):
        return PairAssessment(
            PairState.NO_REFERENCE_POPULATION,
            (
                f"reference-control fold {reference_arm_fold:.3g} is below the provisional "
                f"{config.min_reference_arm_fold:g} donor-local SBR criterion"
            ),
        )
    if not np.isfinite(target_arm_fold):
        return PairAssessment(
            PairState.TECHNICAL_FAILURE,
            "target-arm fold is nonfinite",
        )
    if target_arm_fold < config.min_target_arm_fold:
        if expression_review is ExpressionReview.CONFIRMED_ABSENT:
            return PairAssessment(
                PairState.NO_TARGET_EXPRESSION_CONFIRMED,
                (
                    f"target-arm fold {target_arm_fold:.3g} is below "
                    f"{config.min_target_arm_fold:g} and image review confirms absence"
                ),
            )
        if expression_review is ExpressionReview.CONFIRMED_PRESENT:
            return PairAssessment(
                PairState.INVALID_PAIR,
                (
                    f"target-arm fold {target_arm_fold:.3g} is below "
                    f"{config.min_target_arm_fold:g} despite confirmed target expression"
                ),
            )
        return PairAssessment(
            PairState.LOW_SIGNAL_UNRESOLVED,
            (
                f"target-arm fold {target_arm_fold:.3g} is below the provisional "
                f"{config.min_target_arm_fold:g} donor-local SBR criterion"
            ),
        )
    if expression_review is ExpressionReview.CONFIRMED_ABSENT:
        return PairAssessment(
            PairState.INVALID_PAIR,
            "NNMF target signal contradicts independently confirmed target absence",
        )
    if not converged or not stable:
        return PairAssessment(
            PairState.MODEL_UNSTABLE,
            "NNMF did not converge reproducibly across the locked seeds",
        )
    if pair_review is PairReview.UNREVIEWED:
        return PairAssessment(
            PairState.PAIR_REVIEW_PENDING,
            "computational checks passed but biological pair review is pending",
        )
    return PairAssessment(
        PairState.VALID,
        "directional arms, target support, convergence, and seed stability passed",
    )


def compartment_column(compartment: str, marker: str) -> str:
    """Return the compact parquet column for one QuPath compartment mean."""
    if compartment not in COMPARTMENTS:
        raise ValueError(f"unknown compartment: {compartment}")
    if marker not in QUPATH_MARKERS:
        raise ValueError(f"unknown candidate marker: {marker}")
    return f"{compartment}__{marker}"


def qupath_measurement(compartment: str, marker: str) -> str:
    """Return the raw QuPath measurement name for one compartment mean."""
    return f"{compartment}: {QUPATH_MARKERS[marker]}: Mean"


def _sql_string(value: str | Path) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def _sql_identifier(value: str) -> str:
    return '"' + value.replace('"', '""') + '"'


def _source_select(csv_path: Path, sample_modulus: int, donors: Sequence[str]) -> str:
    donor_sql = ", ".join(_sql_string(str(d)) for d in donors)
    measurements = []
    for compartment in COMPARTMENTS:
        for marker in QUPATH_MARKERS:
            source = _sql_identifier(qupath_measurement(compartment, marker))
            alias = _sql_identifier(compartment_column(compartment, marker))
            measurements.append(f"TRY_CAST({source} AS FLOAT) AS {alias}")
    measurement_sql = ",\n        ".join(measurements)
    source = (
        f"read_csv({_sql_string(csv_path)}, header=true, all_varchar=true, "
        "parallel=true)"
    )
    return f"""
    SELECT
        regexp_extract("Image", '[0-9]+') AS donor_id,
        CAST("Image" AS VARCHAR) AS image,
        CAST("Object ID" AS VARCHAR) AS object_id,
        TRY_CAST("Cell: Area µm^2" AS FLOAT) AS cell_area,
        TRY_CAST("Centroid X µm" AS FLOAT) AS x_um,
        TRY_CAST("Centroid Y µm" AS FLOAT) AS y_um,
        {_sql_string(csv_path.name)} AS source_csv,
        {measurement_sql}
    FROM {source}
    WHERE regexp_extract("Image", '[0-9]+') IN ({donor_sql})
      AND hash(CAST("Object ID" AS VARCHAR)) % {sample_modulus} = 0
    """


def build_extraction_sql(
    csv_paths: Sequence[Path],
    canonical_cells: Path,
    output_dir: Path,
    donors: Sequence[str],
    sample_modulus: int,
) -> str:
    """Build the one-pass DuckDB query for a canonical raw-compartment sample."""
    if not csv_paths:
        raise ValueError("at least one QuPath CSV is required")
    if not donors:
        raise ValueError("at least one donor is required")
    if sample_modulus < 1:
        raise ValueError("sample_modulus must be >= 1")

    source_sql = "\nUNION ALL\n".join(
        _source_select(Path(path), sample_modulus, donors) for path in csv_paths
    )
    canonical_glob = canonical_cells / "donor_id=*" / "*.parquet"
    canonical = f"""
    SELECT
        CAST(donor_id AS VARCHAR) AS donor_id,
        CAST(object_id AS VARCHAR) AS object_id
    FROM read_parquet({_sql_string(canonical_glob)}, hive_partitioning=true)
    WHERE hash(CAST(object_id AS VARCHAR)) % {sample_modulus} = 0
    """
    return f"""
    COPY (
        WITH source AS (
            {source_sql}
        ),
        canonical AS (
            {canonical}
        )
        SELECT source.*
        FROM source
        INNER JOIN canonical USING (donor_id, object_id)
    ) TO {_sql_string(output_dir)}
      (FORMAT PARQUET, PARTITION_BY (donor_id), COMPRESSION ZSTD)
    """


def extract_compartment_sample(
    csv_paths: Sequence[Path],
    canonical_cells: Path,
    output_dir: Path,
    *,
    donors: Sequence[str],
    sample_modulus: int = 20,
    threads: int = 8,
) -> pd.DataFrame:
    """Extract and verify a deterministic sample against canonical raw object IDs."""
    import duckdb

    donors = ensure_eligible_donors(
        donors,
        context="RESTORE validation sample extraction",
    )
    csv_paths = [Path(path).resolve() for path in csv_paths]
    canonical_cells = Path(canonical_cells).resolve()
    output_dir = Path(output_dir).resolve()
    missing = [str(path) for path in csv_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing QuPath CSVs: {missing}")
    if not canonical_cells.is_dir():
        raise FileNotFoundError(f"canonical cells directory not found: {canonical_cells}")
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")

    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    temp_dir.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={int(threads)}")
    con.execute("SET enable_progress_bar=false")
    try:
        sql = build_extraction_sql(
            csv_paths, canonical_cells, temp_dir, donors, sample_modulus
        )
        con.execute(sql)
        sample_glob = temp_dir / "donor_id=*" / "*.parquet"
        canonical_glob = canonical_cells / "donor_id=*" / "*.parquet"
        audit = con.execute(
            f"""
            WITH expected AS (
                SELECT CAST(donor_id AS VARCHAR) AS donor_id,
                       count(*) AS expected_n,
                       count(DISTINCT CAST(object_id AS VARCHAR)) AS expected_unique
                FROM read_parquet({_sql_string(canonical_glob)}, hive_partitioning=true)
                WHERE hash(CAST(object_id AS VARCHAR)) % {sample_modulus} = 0
                  AND CAST(donor_id AS VARCHAR) IN (
                      {", ".join(_sql_string(str(d)) for d in donors)}
                  )
                GROUP BY donor_id
            ),
            observed AS (
                SELECT CAST(donor_id AS VARCHAR) AS donor_id,
                       count(*) AS observed_n,
                       count(DISTINCT object_id) AS observed_unique,
                       count(DISTINCT image) AS image_n
                FROM read_parquet({_sql_string(sample_glob)}, hive_partitioning=true)
                GROUP BY donor_id
            )
            SELECT expected.*, observed.observed_n, observed.observed_unique,
                   observed.image_n
            FROM expected
            LEFT JOIN observed USING (donor_id)
            ORDER BY donor_id
            """
        ).fetchdf()
        if set(audit["donor_id"]) != {str(d) for d in donors}:
            raise ValueError("canonical donor set does not match requested donors")
        mismatch = audit[
            audit["observed_n"].isna()
            | (audit["expected_n"] != audit["observed_n"])
            | (audit["observed_n"] != audit["observed_unique"])
            | (audit["expected_n"] != audit["expected_unique"])
            | (audit["image_n"] != 1)
        ]
        if not mismatch.empty:
            raise ValueError(
                "QuPath sample does not exactly cover canonical sampled IDs:\n"
                + mismatch.to_string(index=False)
            )
        temp_dir.rename(output_dir)
        return audit
    except Exception:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise
    finally:
        con.close()


@dataclass(frozen=True)
class DirectionalGroups:
    reference_group: int
    target_group: int
    directional_mean: bool
    directional_median: bool
    groups: pd.DataFrame


def directional_groups(
    target: np.ndarray, reference: np.ndarray, labels: np.ndarray
) -> DirectionalGroups:
    """Identify the reference-positive/target-negative arm directionally.

    No cluster number or variance heuristic has biological meaning. The reference
    arm must have the highest reference center, while a distinct target arm must
    have the highest target center.
    """
    target = np.asarray(target, dtype=float)
    reference = np.asarray(reference, dtype=float)
    labels = np.asarray(labels)
    if not (target.shape == reference.shape == labels.shape):
        raise ValueError("target, reference, and labels must have matching shapes")
    finite = np.isfinite(target) & np.isfinite(reference)
    rows = []
    for group in np.unique(labels[finite]):
        mask = finite & (labels == group)
        rows.append(
            {
                "group": int(group),
                "n": int(mask.sum()),
                "target_mean": float(target[mask].mean()),
                "target_median": float(np.median(target[mask])),
                "reference_mean": float(reference[mask].mean()),
                "reference_median": float(np.median(reference[mask])),
            }
        )
    groups = pd.DataFrame(rows).sort_values("group").reset_index(drop=True)
    if len(groups) < 2:
        raise ValueError("at least two nonempty groups are required")

    reference_group = int(groups.loc[groups["reference_mean"].idxmax(), "group"])
    target_group = int(groups.loc[groups["target_mean"].idxmax(), "group"])
    by_group = groups.set_index("group")
    distinct = reference_group != target_group
    directional_mean = bool(
        distinct
        and by_group.loc[reference_group, "target_mean"]
        < by_group.loc[target_group, "target_mean"]
        and by_group.loc[reference_group, "reference_mean"]
        > by_group.loc[target_group, "reference_mean"]
    )

    reference_group_median = int(
        groups.loc[groups["reference_median"].idxmax(), "group"]
    )
    target_group_median = int(groups.loc[groups["target_median"].idxmax(), "group"])
    directional_median = bool(
        reference_group_median != target_group_median
        and by_group.loc[reference_group_median, "target_median"]
        < by_group.loc[target_group_median, "target_median"]
        and by_group.loc[reference_group_median, "reference_median"]
        > by_group.loc[target_group_median, "reference_median"]
    )
    return DirectionalGroups(
        reference_group=reference_group,
        target_group=target_group,
        directional_mean=directional_mean,
        directional_median=directional_median,
        groups=groups,
    )


def negative_control_statistics(values: np.ndarray) -> dict[str, float | int]:
    """Compute the manuscript divisor and released-code sensitivity statistic."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        raise ValueError("negative-control population is empty")
    if np.any(values < 0):
        raise ValueError("RESTORE requires nonnegative linear intensities")
    divisor = float(values.max())
    if divisor <= 0:
        raise ValueError("negative-control maximum must be positive")
    normalized = values / divisor
    if np.any(normalized > 1 + 1e-12):
        raise AssertionError("negative controls must normalize to <= 1")
    return {
        "n": int(values.size),
        "maximum": divisor,
        "mean_plus_3sd": float(values.mean() + 3.0 * values.std(ddof=0)),
        "median": float(np.median(values)),
        "p95": float(np.quantile(values, 0.95)),
    }


def normalize_and_call(
    values: np.ndarray,
    divisor: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize by the manuscript divisor; positivity is strictly ``value > divisor``.

    The call is divisor-gated only (2026-07-23). The former ``| target_component_supported`` rescue
    -- which additionally called every NNMF-target-component cell left of the Step 1 separator, even
    below the Step 2 divisor -- is removed: it is negligible for clean markers (<=0.3 pp) but doubles
    a high-background marker like CD11b by retaining background cells the divisor correctly excludes.
    Those below-divisor component cells remain tracked (``target_supported`` /
    ``target_supported_below_divisor``) and are shown in the review bundle as a distinct
    "retained, not called" category, so any marker they would help is visible for a per-marker
    reviewer decision -- they are simply not counted as positive here.
    """
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("RESTORE values must be finite nonnegative linear intensities")
    if not np.isfinite(divisor) or divisor <= 0:
        raise ValueError("divisor must be finite and positive")
    normalized = values / divisor
    return normalized, normalized > 1.0


def normalize_for_assessment(
    values: np.ndarray,
    assessment: PairAssessment,
    divisor: float | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply calls only for a valid pair or an independently confirmed absence."""
    values = np.asarray(values, dtype=float)
    if assessment.state is PairState.VALID:
        if divisor is None:
            raise ValueError("a VALID pair requires a canonical divisor")
        return normalize_and_call(values, divisor)
    if assessment.state is PairState.NO_TARGET_EXPRESSION_CONFIRMED:
        if divisor is not None:
            raise ValueError("confirmed target absence must not fabricate a divisor")
        return np.full(values.shape, np.nan), np.zeros(values.shape, dtype=bool)
    raise ValueError(
        f"cannot produce biological calls for pair state {assessment.state.value}"
    )


def raw_otsu_floor(
    values: np.ndarray,
    *,
    bins: int = 512,
    upper_quantile: float = 0.995,
) -> float:
    """Estimate a raw-linear low/high fitting floor without a frequency target."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values >= 0)]
    if values.size < 2:
        raise ValueError("at least two finite nonnegative intensities are required")
    if bins < 2:
        raise ValueError("input-QC histogram must have at least two bins")
    if not 0 < upper_quantile <= 1:
        raise ValueError("input-QC upper quantile must be in (0, 1]")

    upper = float(np.quantile(values, upper_quantile))
    if upper <= 0:
        raise ValueError("input-QC intensities must contain a positive value")
    clipped = np.clip(values, 0, upper)
    histogram, edges = np.histogram(clipped, bins=bins, range=(0, upper))
    centers = 0.5 * (edges[:-1] + edges[1:])
    probability = histogram / histogram.sum()
    cumulative_probability = np.cumsum(probability)
    cumulative_mean = np.cumsum(probability * centers)
    total_mean = cumulative_mean[-1]
    denominator = cumulative_probability * (1.0 - cumulative_probability)
    between_class_variance = np.divide(
        (total_mean * cumulative_probability - cumulative_mean) ** 2,
        denominator,
        out=np.full_like(denominator, -np.inf),
        where=denominator > 0,
    )
    return float(centers[int(np.argmax(between_class_variance))])


def raw_triangle_floor(
    values: np.ndarray,
    *,
    bins: int = 256,
) -> float:
    """Estimate the bright-tail onset for a skewed raw-linear marker histogram."""
    from skimage.filters import threshold_triangle

    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values >= 0)]
    if values.size < 2:
        raise ValueError("at least two finite nonnegative intensities are required")
    if bins < 2:
        raise ValueError("Triangle histogram must have at least two bins")
    if values.min() == values.max():
        raise ValueError("Triangle intensities must not be constant")
    return float(threshold_triangle(values, nbins=bins))


@dataclass(frozen=True)
class MarkerInputFloor:
    """One donor-local marker floor and the statistics that produced it."""

    marker: str
    method: str
    value: float
    linear_otsu: float
    linear_triangle: float | None


def estimate_marker_input_floor(
    values: np.ndarray,
    marker: str,
    *,
    otsu_bins: int = 512,
    otsu_upper_quantile: float = 0.995,
    triangle_bins: int = 256,
) -> MarkerInputFloor:
    """Apply the explicit arm-floor statistic assigned to one marker."""
    if marker not in MARKER_INPUT_FLOOR_METHODS:
        raise ValueError(f"no input-floor method is configured for marker {marker}")
    method = MARKER_INPUT_FLOOR_METHODS[marker]
    otsu = raw_otsu_floor(
        values,
        bins=otsu_bins,
        upper_quantile=otsu_upper_quantile,
    )
    triangle = None
    if method == LINEAR_OTSU_FLOOR:
        value = otsu
    elif method == MAX_LINEAR_OTSU_TRIANGLE_FLOOR:
        triangle = raw_triangle_floor(values, bins=triangle_bins)
        value = max(otsu, triangle)
    else:
        raise ValueError(f"unknown input-floor method for {marker}: {method}")
    return MarkerInputFloor(
        marker=marker,
        method=method,
        value=float(value),
        linear_otsu=float(otsu),
        linear_triangle=None if triangle is None else float(triangle),
    )


def pair_input_qc(
    target_raw: np.ndarray,
    reference_raw: np.ndarray,
    *,
    target_marker: str,
    reference_marker: str,
    bins: int = 512,
    upper_quantile: float = 0.995,
    triangle_bins: int = 256,
) -> tuple[np.ndarray, MarkerInputFloor, MarkerInputFloor]:
    """Identify fitting-eligible cells while preserving both exclusive arms."""
    target_raw = np.asarray(target_raw, dtype=float)
    reference_raw = np.asarray(reference_raw, dtype=float)
    if target_raw.shape != reference_raw.shape:
        raise ValueError("target and reference intensities must have matching shapes")
    finite = (
        np.isfinite(target_raw)
        & np.isfinite(reference_raw)
        & (target_raw >= 0)
        & (reference_raw >= 0)
    )
    target_floor = estimate_marker_input_floor(
        target_raw[finite],
        target_marker,
        otsu_bins=bins,
        otsu_upper_quantile=upper_quantile,
        triangle_bins=triangle_bins,
    )
    reference_floor = estimate_marker_input_floor(
        reference_raw[finite],
        reference_marker,
        otsu_bins=bins,
        otsu_upper_quantile=upper_quantile,
        triangle_bins=triangle_bins,
    )
    retained = finite & (
        (target_raw > target_floor.value) | (reference_raw > reference_floor.value)
    )
    if retained.sum() < 2:
        raise ValueError("input QC retained fewer than two cells")
    return retained, target_floor, reference_floor


@dataclass(frozen=True)
class OrderedControlGroups:
    """High-confidence controls inferred in the Figure 7 axis order."""

    target_population: np.ndarray
    target_supported: np.ndarray
    reference_candidates: np.ndarray
    reference_control: np.ndarray
    reference_separator: float | None


def ordered_control_groups(
    target_raw: np.ndarray,
    reference_raw: np.ndarray,
    target_redsea: np.ndarray,
    reference_redsea: np.ndarray,
    labels: np.ndarray,
    qc_retained: np.ndarray,
    *,
    target_group: int,
    reference_group: int,
    target_floor: float,
    reference_floor: float,
) -> OrderedControlGroups:
    """Infer the vertical separator before selecting target-negative controls.

    A high-confidence target anchor sets the vertical reference separator. Every
    input-QC-retained target-component cell at or left of that separator remains
    target evidence; double-low cells remain excluded, and reference controls must
    lie strictly to the separator's right to set the target divisor.
    """
    arrays = [
        np.asarray(values)
        for values in (
            target_raw,
            reference_raw,
            target_redsea,
            reference_redsea,
            labels,
            qc_retained,
        )
    ]
    if len({values.shape for values in arrays}) != 1:
        raise ValueError("ordered-control inputs must have matching shapes")
    target_raw, reference_raw, target_redsea, reference_redsea = [
        np.asarray(values, dtype=float) for values in arrays[:4]
    ]
    labels = arrays[4]
    qc_retained = np.asarray(arrays[5], dtype=bool)
    finite = (
        np.isfinite(target_raw)
        & np.isfinite(reference_raw)
        & np.isfinite(target_redsea)
        & np.isfinite(reference_redsea)
        & (target_raw >= 0)
        & (reference_raw >= 0)
        & (target_redsea >= 0)
        & (reference_redsea >= 0)
    )
    target_population = (
        finite
        & qc_retained
        & (labels == target_group)
        & (target_raw > target_floor)
        & (reference_raw <= reference_floor)
    )
    reference_candidates = (
        finite
        & qc_retained
        & (labels == reference_group)
        & (target_raw <= target_floor)
        & (reference_raw > reference_floor)
    )
    if not target_population.any():
        return OrderedControlGroups(
            target_population=target_population,
            target_supported=np.zeros_like(target_population),
            reference_candidates=reference_candidates,
            reference_control=np.zeros_like(target_population),
            reference_separator=None,
        )

    reference_separator = float(reference_redsea[target_population].max())
    target_supported = (
        finite
        & qc_retained
        & (labels == target_group)
        & (reference_redsea <= reference_separator)
    )
    reference_control = reference_candidates & (
        reference_redsea > reference_separator
    )
    if np.any(reference_redsea[target_population] > reference_separator):
        raise AssertionError("vertical baseline does not bound the target population")
    if np.any(reference_redsea[reference_control] <= reference_separator):
        raise AssertionError("reference controls must lie right of the vertical baseline")
    if np.any(target_supported & ~qc_retained):
        raise AssertionError("double-low cells cannot be target-component supported")
    if np.any(target_population & ~target_supported):
        raise AssertionError("the target anchor must be included in target support")
    return OrderedControlGroups(
        target_population=target_population,
        target_supported=target_supported,
        reference_candidates=reference_candidates,
        reference_control=reference_control,
        reference_separator=reference_separator,
    )


def balanced_arm_fit_indices(
    target_raw: np.ndarray,
    reference_raw: np.ndarray,
    target_floor: float,
    reference_floor: float,
    *,
    fit_cap: int,
    seed: int,
    min_arm_n: int = 2,
    probe_fraction: float | None = None,
) -> tuple[np.ndarray, dict[str, int]]:
    """Sample provisional exclusive arms equally for separator fitting.

    ``probe_fraction`` (stability probing only, never the reported fit) shrinks each arm to that
    fraction of the balanced size so a SATURATED arm -- one where the balanced draw would otherwise take
    every cell it has, making re-seeding a no-op -- is actually perturbed between seeds.
    """
    target_raw = np.asarray(target_raw, dtype=float)
    reference_raw = np.asarray(reference_raw, dtype=float)
    if target_raw.shape != reference_raw.shape:
        raise ValueError("target and reference intensities must have matching shapes")
    if fit_cap < 2 * min_arm_n:
        raise ValueError("fit_cap is too small for the required cells per arm")
    target_only = np.flatnonzero(
        (target_raw > target_floor) & (reference_raw <= reference_floor)
    )
    reference_only = np.flatnonzero(
        (target_raw <= target_floor) & (reference_raw > reference_floor)
    )
    double_high = np.flatnonzero(
        (target_raw > target_floor) & (reference_raw > reference_floor)
    )
    per_arm = min(len(target_only), len(reference_only), fit_cap // 2)
    if per_arm < min_arm_n:
        raise ValueError(
            f"input QC retained fewer than {min_arm_n} cells in an exclusive arm"
        )
    saturated_target = per_arm == len(target_only)
    saturated_reference = per_arm == len(reference_only)
    if probe_fraction is not None:
        if not 0 < probe_fraction <= 1:
            raise ValueError("probe_fraction must be in (0, 1]")
        per_arm = max(min_arm_n, int(round(probe_fraction * per_arm)))
    rng = np.random.default_rng(seed)
    target_fit = rng.choice(target_only, size=per_arm, replace=False)
    reference_fit = rng.choice(reference_only, size=per_arm, replace=False)
    return np.sort(np.concatenate([target_fit, reference_fit])), {
        "candidate_target_n": int(len(target_only)),
        "candidate_reference_n": int(len(reference_only)),
        "candidate_double_high_n": int(len(double_high)),
        "fit_per_arm": int(per_arm),
        # An arm is SATURATED when the balanced sample takes every cell it has: the without-replacement
        # draw is then a permutation, so re-seeding cannot perturb it. Seed agreement across a saturated
        # arm is arithmetically guaranteed and is NOT evidence of robustness (see stability_informative).
        # Reported for the UNPROBED balanced size, so a probe run still describes the real fit.
        "target_arm_saturated": bool(saturated_target),
        "reference_arm_saturated": bool(saturated_reference),
    }


def robust_minmax_scale(
    features: np.ndarray,
    fit_indices: np.ndarray,
    *,
    lower_quantile: float = 0.005,
    upper_quantile: float = 0.995,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fit robust bounds on balanced training rows and scale every supplied row."""
    features = np.asarray(features, dtype=float)
    fit_indices = np.asarray(fit_indices)
    if features.ndim != 2 or not np.isfinite(features).all():
        raise ValueError("features must be a finite two-dimensional matrix")
    if np.any(features < 0):
        raise ValueError("features must be nonnegative")
    if (
        fit_indices.ndim != 1
        or not np.issubdtype(fit_indices.dtype, np.integer)
        or len(fit_indices) == 0
        or np.any(fit_indices < 0)
        or np.any(fit_indices >= len(features))
    ):
        raise ValueError("fit_indices must contain in-range integer rows")
    if not 0 <= lower_quantile < upper_quantile <= 1:
        raise ValueError("invalid robust min-max quantiles")
    fit_features = features[fit_indices.astype(int, copy=False)]
    lower = np.quantile(fit_features, lower_quantile, axis=0)
    upper = np.quantile(fit_features, upper_quantile, axis=0)
    if np.any(upper <= lower):
        raise ValueError("robust min-max bounds must have positive ranges")
    scaled = np.clip((features - lower) / (upper - lower), 0.0, 1.0)
    return scaled, lower, upper


def quantile_support_cells(
    n: int, *, lower_quantile: float, upper_quantile: float
) -> float:
    """Cells spanned by the narrower robust min-max tail for a sample of ``n`` rows."""
    if n < 0:
        raise ValueError("n must be nonnegative")
    if not 0 <= lower_quantile < upper_quantile <= 1:
        raise ValueError("invalid robust min-max quantiles")
    return float(min(lower_quantile * n, (1.0 - upper_quantile) * n))


def adaptive_bound_quantiles(
    n: int,
    *,
    lower_quantile: float,
    upper_quantile: float,
    min_support: int = 1,
    repair_support: int = 10,
    max_tail: float = 0.25,
) -> tuple[float, float, bool, float]:
    """Widen the robust min-max tails when ``n`` rows cannot estimate the configured pair.

    A percentile is only robust if the sample spans it. A sparse marker's exclusive arm can hold a few
    dozen cells, and ``per_arm`` caps the balanced sample at the SMALLER arm, so the configured 0.5%
    tail can work out to well under one cell -- ``np.quantile`` then interpolates between the two most
    extreme values and a single bright cell sets the scale for every feature.

    The repair widens the tails on the SAME balanced sample rather than estimating them from a larger
    but unbalanced population. Balance is what keeps either marker from dominating the Frobenius loss,
    and dropping it simply moves the domination to the other marker: measured on the 22 affected donor
    pairs, re-estimating bounds from the input-QC-retained population raised anchor recovery to 1.00 on
    the 18 sparse-target rows but collapsed it (1.00 -> 0.10) on the 4 inverted rows, where the target
    arm is abundant and the REFERENCE arm is the small one.

    ``max_tail`` stops the widening from becoming meaningless on a very small sample: a tail can never
    exceed the quartile. Returns ``(lower, upper, adapted, support_cells)``.
    """
    if n < 0:
        raise ValueError("n must be nonnegative")
    if not 0 <= lower_quantile < upper_quantile <= 1:
        raise ValueError("invalid robust min-max quantiles")
    if min_support < 1 or repair_support < min_support:
        raise ValueError("require 1 <= min_support <= repair_support")
    if not 0 < max_tail < 0.5:
        raise ValueError("max_tail must be in (0, 0.5)")

    support = quantile_support_cells(
        n, lower_quantile=lower_quantile, upper_quantile=upper_quantile
    )
    if support >= min_support or n == 0:
        return lower_quantile, upper_quantile, False, support

    tail = min(repair_support / n, max_tail)
    lower = max(lower_quantile, tail)
    upper = min(upper_quantile, 1.0 - tail)
    if lower >= upper:  # degenerate sample -- keep the configured pair rather than invert it
        return lower_quantile, upper_quantile, False, support
    return (
        lower,
        upper,
        True,
        quantile_support_cells(n, lower_quantile=lower, upper_quantile=upper),
    )


def resampled_negative_control_statistics(
    values: np.ndarray,
    *,
    sample_n: int,
    repeats: int,
    seed: int,
) -> dict[str, float | int]:
    """Summarize fixed-N sampled maxima without replacing the full maximum."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if np.any(values < 0):
        raise ValueError("negative-control intensities must be nonnegative")
    if sample_n < 1 or sample_n > len(values):
        raise ValueError("sample_n must be within the negative-control population")
    if repeats < 1:
        raise ValueError("at least one threshold resample is required")

    rng = np.random.default_rng(seed)
    maxima = np.empty(repeats)
    mean_plus_3sd = np.empty(repeats)
    all_indices = np.arange(len(values))
    for repeat in range(repeats):
        selected = (
            all_indices
            if sample_n == len(values)
            else rng.choice(all_indices, size=sample_n, replace=False)
        )
        stats = negative_control_statistics(values[selected])
        maxima[repeat] = stats["maximum"]
        mean_plus_3sd[repeat] = stats["mean_plus_3sd"]
    full = negative_control_statistics(values)
    return {
        "n_available": int(len(values)),
        "sample_n": int(sample_n),
        "repeats": int(repeats),
        "full_maximum": float(full["maximum"]),
        "maximum": float(np.median(maxima)),
        "maximum_q05": float(np.quantile(maxima, 0.05)),
        "maximum_q95": float(np.quantile(maxima, 0.95)),
        "mean_plus_3sd": float(np.median(mean_plus_3sd)),
        "full_mean_plus_3sd": float(full["mean_plus_3sd"]),
    }


def sampled_maximum_grid(
    values: np.ndarray,
    *,
    sample_sizes: Sequence[int],
    repeats: int,
    seed: int,
    max_sample_n: int | None = None,
) -> pd.DataFrame:
    """Evaluate order-statistic sensitivity over a predeclared fixed-N grid."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        raise ValueError("negative-control population is empty")
    limit = len(values) if max_sample_n is None else min(len(values), max_sample_n)
    if limit < 1:
        raise ValueError("max_sample_n leaves no negative-control cells")
    sizes = sorted({int(size) for size in sample_sizes if 1 <= int(size) <= limit})
    if not sizes:
        sizes = [limit]
    rows = []
    for offset, sample_n in enumerate(sizes):
        stats = resampled_negative_control_statistics(
            values,
            sample_n=sample_n,
            repeats=repeats,
            seed=seed + offset,
        )
        rows.append(
            {
                **stats,
                "maximum_to_full": float(
                    stats["maximum"] / stats["full_maximum"]
                ),
            }
        )
    return pd.DataFrame(rows)


def fit_nnmf(
    features: np.ndarray,
    *,
    n_components: int,
    seed: int,
    max_iter: int = 2000,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | bool]]:
    """Fit manuscript NNMF and return scale-invariant dominant-component labels."""
    from sklearn.decomposition import NMF
    from sklearn.exceptions import ConvergenceWarning

    features = np.asarray(features, dtype=float)
    if features.ndim != 2:
        raise ValueError("features must be a two-dimensional matrix")
    if not np.isfinite(features).all():
        raise ValueError("NNMF features must be finite")
    if np.any(features < 0):
        raise ValueError("NNMF features must be nonnegative")
    if n_components < 2 or n_components > min(features.shape):
        raise ValueError("invalid NNMF component count")

    model = NMF(
        n_components=n_components,
        init="nndsvdar",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=max_iter,
        tol=1e-5,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        weights = model.fit_transform(features)
    component_norm = np.linalg.norm(model.components_, axis=1)
    if np.any(component_norm <= 0):
        raise ValueError("NNMF produced an empty component")
    labels = np.argmax(weights * component_norm, axis=1)
    unit_components = model.components_ / component_norm[:, None]
    converged = not any(issubclass(w.category, ConvergenceWarning) for w in caught)
    metadata: dict[str, float | int | bool] = {
        "converged": converged,
        "n_iter": int(model.n_iter_),
        "reconstruction_error": float(model.reconstruction_err_),
    }
    return labels, unit_components, metadata


def _fit_indices(n: int, cap: int | None, seed: int) -> np.ndarray:
    if not cap or n <= cap:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, size=cap, replace=False))


def fit_nnmf_predict(
    features: np.ndarray,
    *,
    n_components: int,
    seed: int,
    fit_cap: int | None,
    fit_indices: np.ndarray | None = None,
    max_iter: int = 2000,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | bool]]:
    """Fit NNMF on a deterministic subset and assign every supplied cell."""
    from sklearn.decomposition import NMF
    from sklearn.exceptions import ConvergenceWarning

    features = np.asarray(features, dtype=float)
    if features.ndim != 2 or not np.isfinite(features).all():
        raise ValueError("NNMF features must be a finite two-dimensional matrix")
    if np.any(features < 0):
        raise ValueError("NNMF features must be nonnegative")
    if fit_indices is None:
        fit_idx = _fit_indices(len(features), fit_cap, seed)
    else:
        if fit_cap is not None:
            raise ValueError("fit_cap and fit_indices cannot both be specified")
        fit_idx = np.asarray(fit_indices)
        if fit_idx.ndim != 1 or not np.issubdtype(fit_idx.dtype, np.integer):
            raise ValueError("fit_indices must be a one-dimensional integer array")
        if (
            len(fit_idx) == 0
            or np.any(fit_idx < 0)
            or np.any(fit_idx >= len(features))
            or len(np.unique(fit_idx)) != len(fit_idx)
        ):
            raise ValueError("fit_indices must be unique in-range row indices")
        fit_idx = np.sort(fit_idx.astype(int, copy=False))
    if n_components < 2 or n_components > min(len(fit_idx), features.shape[1]):
        raise ValueError("invalid NNMF component count")

    model = NMF(
        n_components=n_components,
        init="nndsvdar",
        solver="cd",
        beta_loss="frobenius",
        random_state=seed,
        max_iter=max_iter,
        tol=1e-5,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        model.fit(features[fit_idx])
        weights = model.transform(features)
    component_norm = np.linalg.norm(model.components_, axis=1)
    if np.any(component_norm <= 0):
        raise ValueError("NNMF produced an empty component")
    labels = np.argmax(weights * component_norm, axis=1)
    unit_components = model.components_ / component_norm[:, None]
    converged = not any(issubclass(w.category, ConvergenceWarning) for w in caught)
    metadata: dict[str, float | int | bool] = {
        "converged": converged,
        "n_iter": int(model.n_iter_),
        "fit_n": int(len(fit_idx)),
        "reconstruction_error": float(model.reconstruction_err_),
    }
    return labels, unit_components, metadata


def fit_gmm_predict(
    features: np.ndarray,
    *,
    n_components: int,
    seed: int,
    fit_cap: int | None,
) -> tuple[np.ndarray, np.ndarray, dict[str, float | int | bool]]:
    """Fit a manuscript-permitted GMM comparator and assign every cell."""
    from sklearn.mixture import GaussianMixture

    features = np.asarray(features, dtype=float)
    if features.ndim != 2 or not np.isfinite(features).all():
        raise ValueError("GMM features must be a finite two-dimensional matrix")
    fit_idx = _fit_indices(len(features), fit_cap, seed)
    model = GaussianMixture(
        n_components=n_components,
        covariance_type="full",
        reg_covar=1e-6,
        n_init=1,
        random_state=seed,
        max_iter=300,
    ).fit(features[fit_idx])
    labels = model.predict(features)
    metadata: dict[str, float | int | bool] = {
        "converged": bool(model.converged_),
        "n_iter": int(model.n_iter_),
        "fit_n": int(len(fit_idx)),
        "lower_bound": float(model.lower_bound_),
    }
    return labels, model.means_, metadata


def _sample_file(root: Path, donor: str) -> Path:
    files = sorted((root / f"donor_id={donor}").glob("*.parquet"))
    if len(files) != 1:
        raise ValueError(
            f"expected one compartment sample for donor {donor}, found {len(files)}"
        )
    return files[0]


def _canonical_file(root: Path, donor: str) -> Path:
    files = sorted((root / f"donor_id={donor}").glob("*.parquet"))
    if len(files) != 1:
        raise ValueError(
            f"expected one canonical parquet for donor {donor}, found {len(files)}"
        )
    return files[0]


def load_validation_sample(
    sample_root: Path,
    raw_cells: Path,
    redsea_cells: Path,
    donor: str,
) -> tuple[pd.DataFrame, dict[str, float | int | str]]:
    """Load one sampled donor and prove it matches raw and REDSEA key lineages."""
    donor = ensure_eligible_donors(
        [donor],
        context="RESTORE validation sample loading",
    )[0]
    markers = list(QUPATH_MARKERS)
    sample = pd.read_parquet(_sample_file(sample_root, donor))
    raw = pd.read_parquet(
        _canonical_file(raw_cells, donor),
        columns=["object_id", "X_centroid", "Y_centroid"] + markers,
    )
    redsea = pd.read_parquet(
        _canonical_file(redsea_cells, donor),
        columns=["object_id"] + markers,
    )
    for name, frame in (("sample", sample), ("raw", raw), ("REDSEA", redsea)):
        if "object_id" not in frame or frame["object_id"].isna().any():
            raise ValueError(f"donor {donor}: {name} object_id is missing or null")
    sample["object_id"] = sample["object_id"].astype(str)
    raw["object_id"] = raw["object_id"].astype(str)
    redsea["object_id"] = redsea["object_id"].astype(str)
    if not sample["object_id"].is_unique:
        raise ValueError(f"donor {donor}: duplicate sampled object_id")
    if not raw["object_id"].is_unique or not redsea["object_id"].is_unique:
        raise ValueError(f"donor {donor}: canonical object_id is not unique")
    raw_keys = raw["object_id"].sort_values(ignore_index=True)
    redsea_keys = redsea["object_id"].sort_values(ignore_index=True)
    if not raw_keys.equals(redsea_keys):
        raw_only = raw_keys[~raw_keys.isin(redsea_keys)].head(5).tolist()
        redsea_only = redsea_keys[~redsea_keys.isin(raw_keys)].head(5).tolist()
        raise ValueError(
            f"donor {donor}: raw/REDSEA canonical key mismatch "
            f"(raw={len(raw_keys):,}, REDSEA={len(redsea_keys):,}, "
            f"raw-only examples={raw_only}, REDSEA-only examples={redsea_only})"
        )
    if "donor_id" in sample and set(sample["donor_id"].astype(str)) != {str(donor)}:
        raise ValueError(f"donor {donor}: sampled partition contains another donor")
    if sample["image"].nunique(dropna=False) != 1:
        raise ValueError(f"donor {donor}: expected exactly one image")

    raw = raw.rename(
        columns={
            "X_centroid": "x_canonical",
            "Y_centroid": "y_canonical",
            **{marker: f"raw__{marker}" for marker in markers},
        }
    )
    raw["_raw_present"] = True
    redsea = redsea.rename(
        columns={marker: f"redsea__{marker}" for marker in markers}
    )
    redsea["_redsea_present"] = True
    merged = sample.merge(raw, on="object_id", how="left", validate="one_to_one")
    merged = merged.merge(redsea, on="object_id", how="left", validate="one_to_one")
    if merged["_raw_present"].isna().any():
        raise ValueError(f"donor {donor}: sampled ID absent from canonical raw cells")
    if merged["_redsea_present"].isna().any():
        raise ValueError(f"donor {donor}: sampled ID absent from REDSEA cells")
    merged = merged.drop(columns=["_raw_present", "_redsea_present"])

    max_delta = 0.0
    for marker in markers:
        qupath = merged[compartment_column("Cell", marker)].to_numpy(float)
        canonical = merged[f"raw__{marker}"].to_numpy(float)
        finite = np.isfinite(qupath) & np.isfinite(canonical)
        if not finite.all():
            raise ValueError(f"donor {donor}: nonfinite whole-cell raw value for {marker}")
        delta = float(np.max(np.abs(qupath - canonical), initial=0.0))
        max_delta = max(max_delta, delta)
        if not np.allclose(qupath, canonical, rtol=1e-5, atol=1e-3):
            raise ValueError(
                f"donor {donor}: QuPath/canonical raw mismatch for {marker}; "
                f"max delta={delta}"
            )

    if not np.allclose(
        merged["x_um"], merged["x_canonical"], rtol=1e-5, atol=1e-3
    ) or not np.allclose(
        merged["y_um"], merged["y_canonical"], rtol=1e-5, atol=1e-3
    ):
        raise ValueError(f"donor {donor}: QuPath/canonical coordinate mismatch")

    audit: dict[str, float | int | str] = {
        "donor": donor,
        "n": int(len(merged)),
        "canonical_n": int(len(raw)),
        "raw_redsea_key_match": True,
        "qupath_raw_max_abs_delta": max_delta,
        "image_n": int(merged["image"].nunique()),
        "image": str(merged["image"].iloc[0]),
    }
    for name, compartments in {
        "membrane_cell": ("Membrane", "Cell"),
        "nonuclear": ("Cytoplasm", "Membrane", "Cell"),
        "all_compartments": COMPARTMENTS,
    }.items():
        cols = [
            compartment_column(compartment, marker)
            for compartment in compartments
            for marker in markers
        ]
        audit[f"{name}_complete_fraction"] = float(
            merged[cols].notna().all(axis=1).mean()
        )
    return merged, audit


def _pair_metrics(values: pd.DataFrame, target: str, reference: str) -> dict[str, float]:
    from scipy.stats import spearmanr

    result: dict[str, float] = {}
    for source in ("raw", "redsea"):
        target_values = values[f"{source}__{target}"].to_numpy(float)
        reference_values = values[f"{source}__{reference}"].to_numpy(float)
        finite = np.isfinite(target_values) & np.isfinite(reference_values)
        matrix = np.column_stack(
            [target_values[finite], reference_values[finite]]
        )
        singular = np.linalg.svd(matrix, compute_uv=False)
        prefix = f"{source}_"
        result[prefix + "svd_ratio"] = (
            float(singular[1] / singular[0]) if singular[0] > 0 else np.nan
        )
        result[prefix + "pearson"] = float(
            np.corrcoef(matrix[:, 0], matrix[:, 1])[0, 1]
        )
        result[prefix + "spearman"] = float(
            spearmanr(matrix[:, 0], matrix[:, 1]).statistic
        )
        result[prefix + "double_zero_fraction"] = float(
            ((matrix[:, 0] == 0) & (matrix[:, 1] == 0)).mean()
        )
    return result


def _spatial_agreement(
    coordinates: np.ndarray,
    control: np.ndarray,
    *,
    cap: int = 20_000,
    seed: int = 0,
) -> tuple[float, float]:
    from sklearn.neighbors import NearestNeighbors

    finite = np.isfinite(coordinates).all(axis=1)
    coordinates = coordinates[finite]
    control = np.asarray(control, dtype=bool)[finite]
    if len(coordinates) < 3 or control.all() or (~control).all():
        return np.nan, np.nan
    idx = _fit_indices(len(coordinates), cap, seed)
    coordinates = coordinates[idx]
    control = control[idx]
    neighbors = NearestNeighbors(n_neighbors=2).fit(coordinates)
    nearest = neighbors.kneighbors(
        coordinates, n_neighbors=2, return_distance=False
    )[:, 1]
    agreement = float((control == control[nearest]).mean())
    prevalence = float(control.mean())
    chance = prevalence**2 + (1.0 - prevalence) ** 2
    excess = (
        float((agreement - chance) / (1.0 - chance))
        if chance < 1.0
        else np.nan
    )
    return agreement, excess


def _safe_control_stats(values: np.ndarray) -> dict[str, float | int]:
    try:
        return negative_control_statistics(values)
    except ValueError:
        return {
            "n": int(np.isfinite(values).sum()),
            "maximum": np.nan,
            "mean_plus_3sd": np.nan,
            "median": np.nan,
            "p95": np.nan,
        }


def _summarize_assignment(
    donor_df: pd.DataFrame,
    valid_idx: np.ndarray,
    labels: np.ndarray,
    *,
    donor: str,
    target: str,
    reference: str,
    method: str,
    source: str,
    n_components: int,
    seed: int,
    metadata: dict[str, float | int | bool],
    components: np.ndarray,
) -> tuple[dict, list[dict], np.ndarray]:
    target_values = donor_df[f"{source}__{target}"].to_numpy(float)[valid_idx]
    reference_values = donor_df[f"{source}__{reference}"].to_numpy(float)[valid_idx]
    direction = directional_groups(target_values, reference_values, labels)
    control = labels == direction.reference_group
    raw_target = donor_df[f"raw__{target}"].to_numpy(float)[valid_idx]
    redsea_target = donor_df[f"redsea__{target}"].to_numpy(float)[valid_idx]
    raw_stats = _safe_control_stats(raw_target[control])
    redsea_stats = _safe_control_stats(redsea_target[control])
    coords = donor_df[["x_canonical", "y_canonical"]].to_numpy(float)[valid_idx]
    spatial, spatial_excess = _spatial_agreement(coords, control)

    by_group = direction.groups.set_index("group")
    ref_group = direction.reference_group
    target_group = direction.target_group
    run = {
        "donor": donor,
        "target": target,
        "reference": reference,
        "method": method,
        "source": source,
        "n_components": n_components,
        "seed": seed,
        "n_input": int(len(donor_df)),
        "n_evaluated": int(len(valid_idx)),
        "coverage_fraction": float(len(valid_idx) / len(donor_df)),
        "reference_group": ref_group,
        "target_group": target_group,
        "directional_mean": direction.directional_mean,
        "directional_median": direction.directional_median,
        "control_n": int(control.sum()),
        "control_fraction": float(control.mean()),
        "target_mean_fold": float(
            by_group.loc[target_group, "target_mean"]
            / max(by_group.loc[ref_group, "target_mean"], np.finfo(float).tiny)
        ),
        "reference_mean_fold": float(
            by_group.loc[ref_group, "reference_mean"]
            / max(by_group.loc[target_group, "reference_mean"], np.finfo(float).tiny)
        ),
        "sample_raw_control_max": raw_stats["maximum"],
        "sample_raw_control_mean_plus_3sd": raw_stats["mean_plus_3sd"],
        "sample_redsea_control_max": redsea_stats["maximum"],
        "sample_redsea_control_mean_plus_3sd": redsea_stats["mean_plus_3sd"],
        "redsea_control_zero_fraction": float((redsea_target[control] == 0).mean()),
        "spatial_neighbor_agreement": spatial,
        "spatial_neighbor_excess": spatial_excess,
        "components": json.dumps(np.asarray(components).tolist(), separators=(",", ":")),
        **metadata,
    }
    group_rows = []
    for row in direction.groups.to_dict("records"):
        group_rows.append(
            {
                "donor": donor,
                "target": target,
                "reference": reference,
                "method": method,
                "source": source,
                "n_components": n_components,
                "seed": seed,
                **row,
            }
        )
    full_control = np.full(len(donor_df), -1, dtype=np.int8)
    full_control[valid_idx] = control.astype(np.int8)
    return run, group_rows, full_control


def _binary_agreement(first: np.ndarray, second: np.ndarray) -> dict[str, float | int]:
    from sklearn.metrics import adjusted_rand_score

    common = (first >= 0) & (second >= 0)
    first = first[common].astype(bool)
    second = second[common].astype(bool)
    union = (first | second).sum()
    return {
        "n_common": int(common.sum()),
        "control_jaccard": float((first & second).sum() / union) if union else 1.0,
        "binary_ari": float(adjusted_rand_score(first, second)),
    }


def load_expression_reviews(
    path: Path | None,
) -> tuple[dict[tuple[str, str], ExpressionReview], pd.DataFrame]:
    """Load explicit donor/target image reviews; absence is never inferred here."""
    columns = ["donor", "target", "review", "evidence"]
    if path is None:
        return {}, pd.DataFrame(columns=columns)
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"expression-review file not found: {path}")
    reviews = pd.read_csv(path, dtype={"donor": str, "target": str})
    missing = set(columns).difference(reviews.columns)
    if missing:
        raise ValueError(f"expression-review file is missing columns: {sorted(missing)}")
    reviews = reviews[columns].copy()
    if reviews[["donor", "target"]].duplicated().any():
        raise ValueError("expression-review donor/target keys must be unique")
    reviews["review"] = reviews["review"].map(lambda value: ExpressionReview(value).value)
    reviewed = reviews["review"] != ExpressionReview.UNREVIEWED.value
    if reviews.loc[reviewed, "evidence"].fillna("").str.strip().eq("").any():
        raise ValueError("reviewed expression states require nonempty evidence")
    mapping = {
        (str(row.donor), str(row.target)): ExpressionReview(row.review)
        for row in reviews.itertuples(index=False)
    }
    return mapping, reviews


def load_pair_reviews(
    path: Path | None,
) -> tuple[dict[tuple[str, str], PairReview], pd.DataFrame]:
    """Load explicit cohort-level acceptance of fixed target/reference rules."""
    columns = ["target", "reference", "review", "evidence"]
    if path is None:
        return {}, pd.DataFrame(columns=columns)
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"pair-review file not found: {path}")
    reviews = pd.read_csv(path, dtype={"target": str, "reference": str})
    missing = set(columns).difference(reviews.columns)
    if missing:
        raise ValueError(f"pair-review file is missing columns: {sorted(missing)}")
    reviews = reviews[columns].copy()
    if reviews[["target", "reference"]].duplicated().any():
        raise ValueError("pair-review target/reference keys must be unique")
    reviews["review"] = reviews["review"].map(lambda value: PairReview(value).value)
    decided = reviews["review"] != PairReview.UNREVIEWED.value
    if reviews.loc[decided, "evidence"].fillna("").str.strip().eq("").any():
        raise ValueError("accepted or rejected pair reviews require nonempty evidence")
    mapping = {
        (str(row.target), str(row.reference)): PairReview(row.review)
        for row in reviews.itertuples(index=False)
    }
    return mapping, reviews


def _locked_fit_once(
    features: np.ndarray,
    target_raw: np.ndarray,
    reference_raw: np.ndarray,
    qc_retained: np.ndarray,
    target_floor: float,
    reference_floor: float,
    *,
    config: PairValidationConfig,
    seed: int,
    probe_fraction: float | None = None,
) -> dict:
    fit_indices, arm_counts = balanced_arm_fit_indices(
        target_raw,
        reference_raw,
        target_floor,
        reference_floor,
        fit_cap=config.fit_cap,
        seed=seed,
        min_arm_n=config.min_arm_n,
        probe_fraction=probe_fraction,
    )
    bound_lower_q, bound_upper_q, bound_adapted, bound_support = adaptive_bound_quantiles(
        len(fit_indices),
        lower_quantile=config.scale_lower_quantile,
        upper_quantile=config.scale_upper_quantile,
        min_support=config.min_quantile_support_cells,
        repair_support=config.degenerate_tail_support_cells,
    )
    scaled, scale_lower, scale_upper = robust_minmax_scale(
        features,
        fit_indices,
        lower_quantile=bound_lower_q,
        upper_quantile=bound_upper_q,
    )
    labels, components, metadata = fit_nnmf_predict(
        scaled,
        n_components=config.n_components,
        seed=seed,
        fit_cap=None,
        fit_indices=fit_indices,
    )
    direction = directional_groups(
        target_raw[qc_retained],
        reference_raw[qc_retained],
        labels[qc_retained],
    )
    return {
        "seed": seed,
        "fit_indices": fit_indices,
        "arm_counts": arm_counts,
        "labels": labels,
        "components": components,
        "metadata": metadata,
        "direction": direction,
        "scale_lower": scale_lower,
        "scale_upper": scale_upper,
        "scale_bound_adapted": bound_adapted,
        "scale_bound_lower_quantile": bound_lower_q,
        "scale_bound_upper_quantile": bound_upper_q,
        "scale_bound_n": int(len(fit_indices)),
        "scale_bound_support_cells": bound_support,
    }


def evaluate_locked_pair(
    donor_df: pd.DataFrame,
    donor: str,
    target: str,
    reference: str,
    *,
    config: PairValidationConfig | None = None,
    expression_review: ExpressionReview = ExpressionReview.UNREVIEWED,
    pair_review: PairReview = PairReview.UNREVIEWED,
) -> dict:
    """Run the single locked pair workflow used by tables, figures, and pilots."""
    donor = ensure_eligible_donors(
        [donor],
        context="RESTORE pair evaluation",
    )[0]
    config = config or PairValidationConfig()
    if (target, reference) not in CANDIDATE_PAIRS:
        raise ValueError(f"pair is not in the candidate specification: {target} <- {reference}")
    feature_columns = [
        compartment_column(compartment, marker)
        for compartment in config.feature_compartments
        for marker in (target, reference)
    ]
    required = [
        f"raw__{target}",
        f"raw__{reference}",
        f"redsea__{target}",
        f"redsea__{reference}",
        *feature_columns,
    ]
    missing = set(required).difference(donor_df.columns)
    if missing:
        raise ValueError(
            f"donor {donor} {target} <- {reference}: missing columns {sorted(missing)}"
        )

    features_all = donor_df[feature_columns].to_numpy(float)
    target_raw_all = donor_df[f"raw__{target}"].to_numpy(float)
    reference_raw_all = donor_df[f"raw__{reference}"].to_numpy(float)
    target_redsea_all = donor_df[f"redsea__{target}"].to_numpy(float)
    reference_redsea_all = donor_df[f"redsea__{reference}"].to_numpy(float)
    canonical_valid = (
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
    valid_idx = np.flatnonzero(canonical_valid)
    if len(valid_idx) < 2 * config.min_arm_n:
        raise ValueError(
            f"donor {donor} {target} <- {reference}: insufficient complete canonical cells"
        )
    unavailable_n = int(len(donor_df) - len(valid_idx))
    unavailable_reason = (
        f"{unavailable_n:,} canonical cells lack finite nonnegative required measurements"
        if unavailable_n
        else ""
    )

    features = features_all[valid_idx]
    target_raw = target_raw_all[valid_idx]
    reference_raw = reference_raw_all[valid_idx]
    target_redsea = target_redsea_all[valid_idx]
    reference_redsea = reference_redsea_all[valid_idx]
    qc_retained, target_floor_result, reference_floor_result = pair_input_qc(
        target_raw,
        reference_raw,
        target_marker=target,
        reference_marker=reference,
        bins=config.input_qc_bins,
        upper_quantile=config.input_qc_upper_quantile,
        triangle_bins=config.input_qc_triangle_bins,
    )
    target_floor = target_floor_result.value
    reference_floor = reference_floor_result.value

    fits = [
        _locked_fit_once(
            features,
            target_raw,
            reference_raw,
            qc_retained,
            target_floor,
            reference_floor,
            config=config,
            seed=seed,
        )
        for seed in config.seeds
    ]
    primary = fits[0]
    direction: DirectionalGroups = primary["direction"]
    reference_group = direction.reference_group
    target_group = direction.target_group
    for fit in fits:
        fit_direction: DirectionalGroups = fit["direction"]
        fit["ordered_controls"] = ordered_control_groups(
            target_raw,
            reference_raw,
            target_redsea,
            reference_redsea,
            fit["labels"],
            qc_retained,
            target_group=fit_direction.target_group,
            reference_group=fit_direction.reference_group,
            target_floor=target_floor,
            reference_floor=reference_floor,
        )
    controls: OrderedControlGroups = primary["ordered_controls"]
    target_population = controls.target_population
    target_supported = controls.target_supported
    reference_control = controls.reference_control
    reference_component = primary["labels"] == reference_group
    target_component = primary["labels"] == target_group

    stability_rows = []
    primary_component = reference_component.astype(np.int8)
    primary_control = reference_control.astype(np.int8)
    for fit in fits:
        fit_direction: DirectionalGroups = fit["direction"]
        fit_component = (
            fit["labels"] == fit_direction.reference_group
        ).astype(np.int8)
        fit_controls: OrderedControlGroups = fit["ordered_controls"]
        component_agreement = _binary_agreement(
            primary_component, fit_component
        )
        control_agreement = _binary_agreement(
            primary_control, fit_controls.reference_control.astype(np.int8)
        )
        stability_rows.append(
            {
                "seed": int(fit["seed"]),
                "directional_mean": fit_direction.directional_mean,
                "directional_median": fit_direction.directional_median,
                "converged": bool(fit["metadata"]["converged"]),
                "component_jaccard": component_agreement["control_jaccard"],
                "component_binary_ari": component_agreement["binary_ari"],
                "scale_bound_adapted": fit["scale_bound_adapted"],
                **control_agreement,
            }
        )
    # Re-seeding only perturbs an arm that was SUBSAMPLED. Because per_arm = min(|target|, |reference|,
    # cap), the SMALLER arm is taken whole unless the cap binds -- so on almost every pair the seed sweep
    # leaves one arm identical across seeds and min_control_jaccard structurally under-tests. It is
    # "informative" only when neither arm is saturated; otherwise read probe_control_jaccard, which
    # shrinks both arms and therefore actually perturbs the one the seed sweep cannot.
    target_arm_saturated = bool(primary["arm_counts"]["target_arm_saturated"])
    reference_arm_saturated = bool(primary["arm_counts"]["reference_arm_saturated"])
    stability_informative = not (target_arm_saturated or reference_arm_saturated)

    # When either arm is saturated the seed sweep cannot perturb it, so run an extra probe that shrinks
    # both arms and re-measures agreement against the primary control set. Reported, never gated: acting
    # on it would move rows into MODEL_UNSTABLE, which is an acceptance decision for the pair review.
    # (Sub-sampling, not resampling with replacement -- fit_nnmf_predict requires unique fit rows.)
    probe_control_jaccard = None
    if target_arm_saturated or reference_arm_saturated:
        probe_values = []
        for seed in config.seeds[1:]:
            try:
                probe = _locked_fit_once(
                    features,
                    target_raw,
                    reference_raw,
                    qc_retained,
                    target_floor,
                    reference_floor,
                    config=config,
                    seed=seed,
                    probe_fraction=config.stability_probe_fraction,
                )
                probe_direction: DirectionalGroups = probe["direction"]
                probe_controls = ordered_control_groups(
                    target_raw,
                    reference_raw,
                    target_redsea,
                    reference_redsea,
                    probe["labels"],
                    qc_retained,
                    target_group=probe_direction.target_group,
                    reference_group=probe_direction.reference_group,
                    target_floor=target_floor,
                    reference_floor=reference_floor,
                )
            except (ValueError, AssertionError):
                probe_values.append(0.0)
                continue
            probe_values.append(
                _binary_agreement(
                    primary_control,
                    probe_controls.reference_control.astype(np.int8),
                )["control_jaccard"]
            )
        if probe_values:
            probe_control_jaccard = float(min(probe_values))
    min_jaccard = float(
        min(row["control_jaccard"] for row in stability_rows)
    )
    converged = all(row["converged"] for row in stability_rows)
    stable = bool(
        all(
            row["directional_mean"] and row["directional_median"]
            for row in stability_rows
        )
        and min_jaccard >= config.min_control_jaccard
    )

    by_seed_target_folds = []
    for fit in fits:
        fit_controls: OrderedControlGroups = fit["ordered_controls"]
        if (
            fit_controls.target_population.any()
            and fit_controls.reference_control.any()
        ):
            by_seed_target_folds.append(
                float(
                    target_redsea[fit_controls.target_population].mean()
                    / max(
                        target_redsea[fit_controls.reference_control].mean(),
                        np.finfo(float).tiny,
                    )
                )
            )
        else:
            by_seed_target_folds.append(np.nan)
    finite_target_folds = np.asarray(by_seed_target_folds, dtype=float)
    finite_target_folds = finite_target_folds[np.isfinite(finite_target_folds)]
    target_fold = (
        float(np.median(finite_target_folds))
        if len(finite_target_folds)
        else np.nan
    )
    reference_fold = (
        float(
            reference_redsea[reference_control].mean()
            / max(
                reference_redsea[target_population].mean(),
                np.finfo(float).tiny,
            )
        )
        if reference_control.any() and target_population.any()
        else np.nan
    )

    target_full_stats = None
    reference_full_stats = None
    target_background = None
    reference_background = None
    target_grid = pd.DataFrame()
    reference_grid = pd.DataFrame()
    background_error = None
    threshold_sample_n = 0
    if reference_group != target_group and target_population.any():
        try:
            reference_full_stats = negative_control_statistics(
                reference_redsea[target_population]
            )
            if reference_control.any():
                target_full_stats = negative_control_statistics(
                    target_redsea[reference_control]
                )
                threshold_sample_n = min(
                    int(reference_control.sum()),
                    int(target_population.sum()),
                    max(config.threshold_sample_sizes),
                )
                target_background = resampled_negative_control_statistics(
                    target_redsea[reference_control],
                    sample_n=threshold_sample_n,
                    repeats=config.threshold_resamples,
                    seed=config.seeds[0] + 101,
                )
                reference_background = resampled_negative_control_statistics(
                    reference_redsea[target_population],
                    sample_n=threshold_sample_n,
                    repeats=config.threshold_resamples,
                    seed=config.seeds[0] + 202,
                )
                common_arm_n = min(
                    int(reference_control.sum()),
                    int(target_population.sum()),
                )
                target_grid = sampled_maximum_grid(
                    target_redsea[reference_control],
                    sample_sizes=config.threshold_sample_sizes,
                    repeats=config.threshold_resamples,
                    seed=config.seeds[0] + 301,
                    max_sample_n=common_arm_n,
                )
                reference_grid = sampled_maximum_grid(
                    reference_redsea[target_population],
                    sample_sizes=config.threshold_sample_sizes,
                    repeats=config.threshold_resamples,
                    seed=config.seeds[0] + 401,
                    max_sample_n=common_arm_n,
                )
        except ValueError as exc:
            threshold_sample_n = 0
            background_error = f"background estimation failed: {exc}"

    assessment = assess_pair_state(
        reference_n=int(reference_control.sum()),
        target_n=int(target_population.sum()),
        distinct=reference_group != target_group,
        directional_mean=direction.directional_mean,
        directional_median=direction.directional_median,
        converged=converged,
        stable=stable,
        target_arm_fold=target_fold,
        reference_arm_fold=reference_fold,
        expression_review=expression_review,
        pair_review=pair_review,
        config=config,
        technical_error=background_error,
    )

    # Display bounds are estimated from the input-QC-retained population, NOT the balanced fit sample.
    # The reviewer's job is to judge the Step 1 separator and the Step 2 divisor, and a bound taken from
    # a handful of arm cells is set by the single brightest one -- which pushed the whole control region
    # into the bottom few percent of the axis. The divisor is then forced on-canvas so it can never be
    # clipped out of the panel it defines. Axes stay raw-linear (no log transform).
    display = np.column_stack([target_redsea, reference_redsea])
    display_bound_rows = np.flatnonzero(qc_retained)
    if not len(display_bound_rows):
        display_bound_rows = primary["fit_indices"]
    display_lower_q, display_upper_q, _adapted, _support = adaptive_bound_quantiles(
        len(display_bound_rows),
        lower_quantile=config.scale_lower_quantile,
        upper_quantile=config.scale_upper_quantile,
        min_support=config.min_quantile_support_cells,
        repair_support=config.degenerate_tail_support_cells,
    )
    _, display_lower, display_upper = robust_minmax_scale(
        display,
        display_bound_rows,
        lower_quantile=display_lower_q,
        upper_quantile=display_upper_q,
    )
    if target_full_stats is not None:
        display_upper[0] = max(
            float(display_upper[0]), 1.05 * float(target_full_stats["maximum"])
        )
    if controls.reference_separator is not None:
        display_upper[1] = max(
            float(display_upper[1]), 1.05 * float(controls.reference_separator)
        )
    full_labels = np.full(len(donor_df), -1, dtype=np.int16)
    full_labels[valid_idx] = primary["labels"]
    full_qc_retained = np.zeros(len(donor_df), dtype=bool)
    full_qc_retained[valid_idx] = qc_retained
    full_target_population = np.zeros(len(donor_df), dtype=bool)
    full_target_population[valid_idx] = target_population
    full_target_supported = np.zeros(len(donor_df), dtype=bool)
    full_target_supported[valid_idx] = target_supported
    full_reference_candidates = np.zeros(len(donor_df), dtype=bool)
    full_reference_candidates[valid_idx] = controls.reference_candidates
    full_reference_control = np.zeros(len(donor_df), dtype=bool)
    full_reference_control[valid_idx] = reference_control
    canonical_divisor = (
        float(target_full_stats["maximum"])
        if assessment.state is PairState.VALID and target_full_stats is not None
        else None
    )
    candidate_divisor = (
        float(target_full_stats["maximum"])
        if target_full_stats is not None
        else None
    )
    target_supported_additional = target_supported & ~target_population
    target_supported_below_divisor = np.zeros(len(target_redsea), dtype=bool)
    target_supported_additional_below_divisor = np.zeros(
        len(target_redsea), dtype=bool
    )
    potential_lineage_positive = np.zeros(len(target_redsea), dtype=bool)
    threshold_positive = np.zeros(len(target_redsea), dtype=bool)
    if candidate_divisor is not None:
        threshold_positive = target_redsea > candidate_divisor
        target_supported_below_divisor = target_supported & ~threshold_positive
        target_supported_additional_below_divisor = (
            target_supported_additional & ~threshold_positive
        )
        potential_lineage_positive = threshold_positive | target_supported
    full_target_supported_below_divisor = np.zeros(len(donor_df), dtype=bool)
    full_target_supported_below_divisor[valid_idx] = target_supported_below_divisor
    full_threshold_positive = np.zeros(len(donor_df), dtype=bool)
    full_threshold_positive[valid_idx] = threshold_positive
    full_target_supported_additional_below_divisor = np.zeros(
        len(donor_df), dtype=bool
    )
    full_target_supported_additional_below_divisor[valid_idx] = (
        target_supported_additional_below_divisor
    )
    full_normalized = np.full(len(donor_df), np.nan, dtype=float)
    full_lineage_positive = np.zeros(len(donor_df), dtype=bool)
    if assessment.state in {
        PairState.VALID,
        PairState.NO_TARGET_EXPRESSION_CONFIRMED,
    }:
        normalized, lineage_positive = normalize_for_assessment(
            target_redsea,
            assessment,
            canonical_divisor,
        )
        full_normalized[valid_idx] = normalized
        full_lineage_positive[valid_idx] = lineage_positive
    return {
        "method_version": METHOD_VERSION,
        "donor": str(donor),
        "target": target,
        "reference": reference,
        "expression_review": ExpressionReview(expression_review).value,
        "pair_review": PairReview(pair_review).value,
        "state": assessment.state.value,
        "state_reason": assessment.reason,
        "assessment": assessment,
        "config": config,
        "feature_columns": feature_columns,
        "valid_idx": valid_idx,
        "full_labels": full_labels,
        "full_qc_retained": full_qc_retained,
        "full_target_population": full_target_population,
        "full_target_supported": full_target_supported,
        "full_target_supported_below_divisor": full_target_supported_below_divisor,
        "full_threshold_positive": full_threshold_positive,
        "full_target_supported_additional_below_divisor": (
            full_target_supported_additional_below_divisor
        ),
        "full_reference_candidates": full_reference_candidates,
        "full_reference_control": full_reference_control,
        "full_normalized": full_normalized,
        "full_lineage_positive": full_lineage_positive,
        "target_raw": target_raw,
        "reference_raw": reference_raw,
        "target_redsea": target_redsea,
        "reference_redsea": reference_redsea,
        "qc_retained": qc_retained,
        "target_population": target_population,
        "target_supported": target_supported,
        "target_supported_additional": target_supported_additional,
        "threshold_positive": threshold_positive,
        "target_supported_below_divisor": target_supported_below_divisor,
        "target_supported_additional_below_divisor": (
            target_supported_additional_below_divisor
        ),
        "potential_lineage_positive": potential_lineage_positive,
        "reference_candidates": controls.reference_candidates,
        "reference_control": reference_control,
        "labels": primary["labels"],
        "reference_group": reference_group,
        "target_group": target_group,
        "distinct": reference_group != target_group,
        "directional_mean": direction.directional_mean,
        "directional_median": direction.directional_median,
        "converged": converged,
        "stable": stable,
        "min_control_jaccard": min_jaccard,
        "stability_rows": stability_rows,
        "target_fold": target_fold,
        "target_fold_by_seed": by_seed_target_folds,
        "reference_fold": reference_fold,
        "reference_n": int(reference_control.sum()),
        "target_n": int(target_population.sum()),
        # Reference-frequency rule (pure proportion, no absolute floor). Diagnostic only -- never gates
        # PairState, never alters a divisor.
        "reference_to_target_ratio": (
            float(int(reference_control.sum()) / int(target_population.sum()))
            if int(target_population.sum()) > 0
            else np.nan
        ),
        "reference_undersized": bool(
            int(target_population.sum()) > 0
            and int(reference_control.sum())
            < config.min_reference_to_target_ratio * int(target_population.sum())
        ),
        "target_anchor_n": int(target_population.sum()),
        "target_supported_n": int(target_supported.sum()),
        "threshold_positive_n": int(threshold_positive.sum()),
        "target_supported_additional_n": int(target_supported_additional.sum()),
        "target_supported_below_divisor_n": int(
            target_supported_below_divisor.sum()
        ),
        "target_supported_additional_below_divisor_n": int(
            target_supported_additional_below_divisor.sum()
        ),
        "target_supported_rescued_fraction": (
            float(
                target_supported_additional_below_divisor.sum()
                / target_supported_additional.sum()
            )
            if target_supported_additional.any() and candidate_divisor is not None
            else np.nan
        ),
        "potential_lineage_positive_n": int(potential_lineage_positive.sum()),
        "lineage_positive_n": int(full_lineage_positive.sum()),
        "reference_candidate_n": int(controls.reference_candidates.sum()),
        "reference_component_n": int(reference_component.sum()),
        "target_component_n": int(target_component.sum()),
        "reference_separator": controls.reference_separator,
        "reference_right_fraction": (
            float(reference_control.sum() / controls.reference_candidates.sum())
            if controls.reference_candidates.any()
            else np.nan
        ),
        "ordered_separation_verified": bool(
            controls.reference_separator is not None
            and np.all(
                reference_redsea[target_population]
                <= controls.reference_separator
            )
            and np.all(
                reference_redsea[target_supported]
                <= controls.reference_separator
            )
            and np.all(
                reference_redsea[reference_control]
                > controls.reference_separator
            )
        ),
        "control_fraction": float(reference_control.mean()),
        "n_input": int(len(donor_df)),
        "n_assigned": int(len(valid_idx)),
        "n_unavailable": unavailable_n,
        "unavailable_reason": unavailable_reason,
        "n_qc_retained": int(qc_retained.sum()),
        "qc_retained_fraction": float(qc_retained.mean()),
        "target_input_floor": target_floor,
        "target_input_floor_method": target_floor_result.method,
        "target_input_floor_linear_otsu": target_floor_result.linear_otsu,
        "target_input_floor_linear_triangle": target_floor_result.linear_triangle,
        "reference_input_floor": reference_floor,
        "reference_input_floor_method": reference_floor_result.method,
        "reference_input_floor_linear_otsu": reference_floor_result.linear_otsu,
        "reference_input_floor_linear_triangle": reference_floor_result.linear_triangle,
        **primary["arm_counts"],
        "target_arm_saturated": target_arm_saturated,
        "reference_arm_saturated": reference_arm_saturated,
        "scale_bound_adapted": primary["scale_bound_adapted"],
        "scale_bound_lower_quantile": primary["scale_bound_lower_quantile"],
        "scale_bound_upper_quantile": primary["scale_bound_upper_quantile"],
        "scale_bound_n": primary["scale_bound_n"],
        "scale_bound_support_cells": primary["scale_bound_support_cells"],
        "stability_informative": stability_informative,
        "probe_control_jaccard": probe_control_jaccard,
        # Fraction of the unambiguous raw target arm the model actually places in the target component.
        # ~0.98 on healthy fits; it collapses when a degenerate scaling bound lets the reference marker
        # dominate the factorization. Reported only -- no threshold gates on it yet.
        "anchor_recovery": (
            float(int(target_population.sum()) / int(primary["arm_counts"]["candidate_target_n"]))
            if int(primary["arm_counts"]["candidate_target_n"]) > 0
            else None
        ),
        "double_high_fraction": float(
            (
                (target_raw[qc_retained] > target_floor)
                & (reference_raw[qc_retained] > reference_floor)
            ).mean()
        ),
        "threshold_sample_n": threshold_sample_n,
        "target_full_stats": target_full_stats,
        "reference_full_stats": reference_full_stats,
        "target_background": target_background,
        "reference_background": reference_background,
        "target_threshold_grid": target_grid,
        "reference_threshold_grid": reference_grid,
        "candidate_full_maximum": candidate_divisor,
        "canonical_divisor": canonical_divisor,
        "display_target_bounds": (
            float(display_lower[0]),
            float(display_upper[0]),
        ),
        "display_reference_bounds": (
            float(display_lower[1]),
            float(display_upper[1]),
        ),
        "scale_lower": primary["scale_lower"],
        "scale_upper": primary["scale_upper"],
        "components": primary["components"],
    }


def locked_pair_metrics(evaluation: dict) -> dict:
    """Flatten one locked evaluation into the canonical donor/pair table row."""
    target_full = evaluation["target_full_stats"] or {}
    reference_full = evaluation["reference_full_stats"] or {}
    target_sample = evaluation["target_background"] or {}
    reference_sample = evaluation["reference_background"] or {}
    return {
        "method_version": evaluation["method_version"],
        "donor": evaluation["donor"],
        "target": evaluation["target"],
        "reference": evaluation["reference"],
        "expression_review": evaluation["expression_review"],
        "pair_review": evaluation["pair_review"],
        "state": evaluation["state"],
        "state_reason": evaluation["state_reason"],
        "directional_mean": evaluation["directional_mean"],
        "directional_median": evaluation["directional_median"],
        "converged": evaluation["converged"],
        "stable": evaluation["stable"],
        "min_control_jaccard": evaluation["min_control_jaccard"],
        "stability_informative": evaluation["stability_informative"],
        "probe_control_jaccard": evaluation["probe_control_jaccard"],
        "target_arm_saturated": evaluation["target_arm_saturated"],
        "reference_arm_saturated": evaluation["reference_arm_saturated"],
        "scale_bound_adapted": evaluation["scale_bound_adapted"],
        "scale_bound_lower_quantile": evaluation["scale_bound_lower_quantile"],
        "scale_bound_upper_quantile": evaluation["scale_bound_upper_quantile"],
        "scale_bound_n": evaluation["scale_bound_n"],
        "scale_bound_support_cells": evaluation["scale_bound_support_cells"],
        "anchor_recovery": evaluation["anchor_recovery"],
        "target_fold": evaluation["target_fold"],
        "reference_fold": evaluation["reference_fold"],
        "n_input": evaluation["n_input"],
        "n_assigned": evaluation["n_assigned"],
        "n_unavailable": evaluation["n_unavailable"],
        "unavailable_reason": evaluation["unavailable_reason"],
        "n_qc_retained": evaluation["n_qc_retained"],
        "qc_retained_fraction": evaluation["qc_retained_fraction"],
        "target_input_floor": evaluation["target_input_floor"],
        "target_input_floor_method": evaluation["target_input_floor_method"],
        "target_input_floor_linear_otsu": evaluation[
            "target_input_floor_linear_otsu"
        ],
        "target_input_floor_linear_triangle": evaluation[
            "target_input_floor_linear_triangle"
        ],
        "reference_input_floor": evaluation["reference_input_floor"],
        "reference_input_floor_method": evaluation["reference_input_floor_method"],
        "reference_input_floor_linear_otsu": evaluation[
            "reference_input_floor_linear_otsu"
        ],
        "reference_input_floor_linear_triangle": evaluation[
            "reference_input_floor_linear_triangle"
        ],
        "candidate_target_n": evaluation["candidate_target_n"],
        "candidate_reference_n": evaluation["candidate_reference_n"],
        "candidate_double_high_n": evaluation["candidate_double_high_n"],
        "fit_per_arm": evaluation["fit_per_arm"],
        "reference_n": evaluation["reference_n"],
        "target_n": evaluation["target_n"],
        "reference_to_target_ratio": evaluation["reference_to_target_ratio"],
        "reference_undersized": evaluation["reference_undersized"],
        "target_anchor_n": evaluation["target_anchor_n"],
        "target_supported_n": evaluation["target_supported_n"],
        "threshold_positive_n": evaluation["threshold_positive_n"],
        "target_supported_additional_n": evaluation[
            "target_supported_additional_n"
        ],
        "target_supported_below_divisor_n": evaluation[
            "target_supported_below_divisor_n"
        ],
        "target_supported_additional_below_divisor_n": evaluation[
            "target_supported_additional_below_divisor_n"
        ],
        "target_supported_rescued_fraction": evaluation[
            "target_supported_rescued_fraction"
        ],
        "potential_lineage_positive_n": evaluation[
            "potential_lineage_positive_n"
        ],
        "lineage_positive_n": evaluation["lineage_positive_n"],
        "reference_candidate_n": evaluation["reference_candidate_n"],
        "reference_component_n": evaluation["reference_component_n"],
        "target_component_n": evaluation["target_component_n"],
        "reference_separator": evaluation["reference_separator"],
        "reference_right_fraction": evaluation["reference_right_fraction"],
        "ordered_separation_verified": evaluation[
            "ordered_separation_verified"
        ],
        "threshold_sample_n": evaluation["threshold_sample_n"],
        "double_high_fraction": evaluation["double_high_fraction"],
        "candidate_full_maximum": evaluation["candidate_full_maximum"],
        "canonical_divisor": evaluation["canonical_divisor"],
        "target_full_mean_plus_3sd": target_full.get("mean_plus_3sd", np.nan),
        "target_sample_maximum": target_sample.get("maximum", np.nan),
        "target_sample_maximum_q05": target_sample.get("maximum_q05", np.nan),
        "target_sample_maximum_q95": target_sample.get("maximum_q95", np.nan),
        "reference_full_maximum": reference_full.get("maximum", np.nan),
        "reference_full_mean_plus_3sd": reference_full.get(
            "mean_plus_3sd", np.nan
        ),
        "reference_sample_maximum": reference_sample.get("maximum", np.nan),
        "reference_sample_maximum_q05": reference_sample.get(
            "maximum_q05", np.nan
        ),
        "reference_sample_maximum_q95": reference_sample.get(
            "maximum_q95", np.nan
        ),
        "feature_scale_lower": json.dumps(
            np.asarray(evaluation["scale_lower"]).tolist(), separators=(",", ":")
        ),
        "feature_scale_upper": json.dumps(
            np.asarray(evaluation["scale_upper"]).tolist(), separators=(",", ":")
        ),
    }


def evaluate_locked_cohort(
    sample_root: Path,
    raw_cells: Path,
    redsea_cells: Path,
    output_dir: Path,
    *,
    donors: Sequence[str],
    pairs: Sequence[tuple[str, str]] = BASELINE_PAIRS,
    config: PairValidationConfig | None = None,
    expression_reviews: Path | None = None,
    pair_reviews: Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Write one deterministic set of locked pair tables and method provenance."""
    config = config or PairValidationConfig()
    donors = ensure_eligible_donors(
        donors,
        context="RESTORE cohort evaluation",
    )
    pairs = tuple((str(target), str(reference)) for target, reference in pairs)
    if not pairs or len(set(pairs)) != len(pairs):
        raise ValueError("pairs must be a nonempty unique sequence")
    unknown_pairs = set(pairs).difference(CANDIDATE_PAIRS)
    if unknown_pairs:
        raise ValueError(f"pairs are outside the candidate specification: {unknown_pairs}")
    review_map, review_table = load_expression_reviews(expression_reviews)
    pair_review_map, pair_review_table = load_pair_reviews(pair_reviews)
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    temp_dir.mkdir(parents=True)

    audits = []
    pair_rows = []
    threshold_rows = []
    stability_rows = []
    try:
        for donor in donors:
            donor_df, audit = load_validation_sample(
                Path(sample_root), Path(raw_cells), Path(redsea_cells), donor
            )
            audits.append(audit)
            for target, reference in pairs:
                evaluation = evaluate_locked_pair(
                    donor_df,
                    donor,
                    target,
                    reference,
                    config=config,
                    expression_review=review_map.get(
                        (donor, target), ExpressionReview.UNREVIEWED
                    ),
                    pair_review=pair_review_map.get(
                        (target, reference), PairReview.UNREVIEWED
                    ),
                )
                pair_rows.append(locked_pair_metrics(evaluation))
                for axis, grid in (
                    ("target", evaluation["target_threshold_grid"]),
                    ("reference", evaluation["reference_threshold_grid"]),
                ):
                    for row in grid.to_dict("records"):
                        threshold_rows.append(
                            {
                                "donor": donor,
                                "target": target,
                                "reference": reference,
                                "axis": axis,
                                **row,
                            }
                        )
                for row in evaluation["stability_rows"]:
                    stability_rows.append(
                        {
                            "donor": donor,
                            "target": target,
                            "reference": reference,
                            **row,
                        }
                    )
            print(f"[locked-evaluate] donor {donor}: {len(donor_df):,} cells", flush=True)

        tables = {
            "data_audit": pd.DataFrame(audits),
            "pair_evaluations": pd.DataFrame(pair_rows),
            "threshold_sampling": pd.DataFrame(threshold_rows),
            "seed_stability": pd.DataFrame(stability_rows),
            "expression_reviews": review_table,
            "pair_reviews": pair_review_table,
        }
        for name, table in tables.items():
            table.to_csv(temp_dir / f"{name}.csv", index=False)
        (temp_dir / "method_spec.json").write_text(
            json.dumps(config.specification(), indent=2) + "\n"
        )
        (temp_dir / "run_manifest.json").write_text(
            json.dumps(
                {
                    "method_version": METHOD_VERSION,
                    "donors": list(donors),
                    "cohort_exclusions": DONOR_EXCLUSIONS,
                    "pairs": [list(pair) for pair in pairs],
                },
                indent=2,
            )
            + "\n"
        )
        temp_dir.rename(output_dir)
        return tables
    except Exception:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise


def evaluate_immune_joint(
    donor_df: pd.DataFrame,
    donor: str,
    *,
    fit_cap: int,
    seeds: Sequence[int],
    markers: Sequence[str] = IMMUNE_JOINT_MARKERS,
    compartments: Sequence[str] = LOCKED_FEATURE_COMPARTMENTS,
    n_components: Sequence[int] = (4, 5),
    scale_lower_quantile: float = 0.005,
    scale_upper_quantile: float = 0.995,
) -> list[dict]:
    """COMPARATOR ONLY -- fit one joint NNMF over the immune markers instead of one per pair.

    The pairwise fit takes ``per_arm = min(|target_only|, |reference_only|)`` cells, which for a sparse
    marker such as CD20 can be a few dozen; the arm is then too small for the robust scaling quantiles
    and the factorization is led by the abundant reference. Pooling defines the fitting set as the union
    of cells positive for ANY immune marker, which is two to three orders of magnitude larger, and lets
    one model describe all four targets consistently instead of refitting each target against every
    candidate reference.

    This writes comparator rows only. It does NOT emit a divisor or a call: RESTORE's divisor is defined
    per target against a mutually-exclusive negative population, and generalizing that (plus a
    directionality rule valid for more than two groups) is a method decision for the pair review, not a
    refactor. The locked production path is untouched and still requires exactly two components.
    """
    donor = ensure_eligible_donors(
        [donor], context="RESTORE immune-joint comparator"
    )[0]
    markers = list(markers)
    feature_cols = [
        compartment_column(compartment, marker)
        for compartment in compartments
        for marker in markers
    ]
    raw_cols = [f"raw__{marker}" for marker in markers]
    missing = set(feature_cols + raw_cols).difference(donor_df.columns)
    if missing:
        raise ValueError(f"donor {donor}: missing immune-joint columns {sorted(missing)}")

    matrix = donor_df[feature_cols].to_numpy(float)
    raw = donor_df[raw_cols].to_numpy(float)
    valid = np.isfinite(matrix).all(axis=1) & (matrix >= 0).all(axis=1) & np.isfinite(raw).all(axis=1)
    valid_idx = np.flatnonzero(valid)
    if not len(valid_idx):
        raise ValueError(f"donor {donor}: no projectable cells for the immune-joint comparator")
    matrix = matrix[valid_idx]
    raw = raw[valid_idx]

    floors = {
        marker: estimate_marker_input_floor(raw[:, j], marker).value
        for j, marker in enumerate(markers)
    }
    above = np.column_stack([raw[:, j] > floors[m] for j, m in enumerate(markers)])
    immune_any = above.any(axis=1)
    # "Exclusive arm" for a target: positive for it and negative for every other pooled marker. This is
    # the same unambiguous population the pairwise anchor uses, so anchor recovery is comparable.
    exclusive = {
        marker: above[:, j] & ~np.delete(above, j, axis=1).any(axis=1)
        for j, marker in enumerate(markers)
    }

    rows: list[dict] = []
    fit_pool = np.flatnonzero(immune_any)
    if len(fit_pool) < max(n_components) or not immune_any.any():
        return rows

    for k in n_components:
        for seed in seeds:
            fit_idx = _fit_indices(len(fit_pool), fit_cap, seed)
            fit_rows = np.sort(fit_pool[fit_idx])
            scaled, _lower, _upper = robust_minmax_scale(
                matrix,
                fit_rows,
                lower_quantile=scale_lower_quantile,
                upper_quantile=scale_upper_quantile,
            )
            try:
                labels, _components, metadata = fit_nnmf_predict(
                    scaled, n_components=k, seed=seed, fit_cap=None, fit_indices=fit_rows
                )
            except ValueError as exc:
                rows.append(
                    {"donor": donor, "method": "nnmf_immune_joint", "n_components": k,
                     "seed": seed, "target": None, "error": str(exc)}
                )
                continue
            # Assign each target the component with the highest mean raw intensity for that marker; a
            # component may win more than one marker, which is exactly the non-exclusivity this
            # comparator is meant to expose (CD11b spans myeloid subsets, CD163 overlaps CD68).
            comp_of = {}
            for j, marker in enumerate(markers):
                means = [
                    raw[labels == c, j].mean() if np.any(labels == c) else -np.inf
                    for c in range(k)
                ]
                comp_of[marker] = int(np.argmax(means))
            for j, marker in enumerate(markers):
                arm = exclusive[marker]
                recovery = (
                    float((labels[arm] == comp_of[marker]).mean()) if arm.any() else None
                )
                rows.append(
                    {
                        "donor": donor,
                        "method": "nnmf_immune_joint",
                        "n_components": k,
                        "seed": seed,
                        "target": marker,
                        "markers": "+".join(markers),
                        "n_fit": int(len(fit_rows)),
                        "n_immune_any": int(immune_any.sum()),
                        "exclusive_arm_n": int(arm.sum()),
                        "target_component": comp_of[marker],
                        "component_shared_with": ";".join(
                            other for other in markers
                            if other != marker and comp_of[other] == comp_of[marker]
                        ),
                        "anchor_recovery": recovery,
                        "target_input_floor": floors[marker],
                        "converged": bool(metadata["converged"]),
                        "error": None,
                    }
                )
    return rows


def evaluate_donor(
    donor_df: pd.DataFrame,
    donor: str,
    *,
    fit_cap: int,
    seeds: Sequence[int],
) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    """Evaluate NNMF feature sets and GMM comparators for one donor."""
    donor = ensure_eligible_donors(
        [donor],
        context="RESTORE exploratory donor evaluation",
    )[0]
    runs: list[dict] = []
    groups: list[dict] = []
    pair_rows: list[dict] = []
    stability: list[dict] = []
    assignments: dict[tuple[str, str, str, int, int], np.ndarray] = {}

    for target, reference in BASELINE_PAIRS:
        pair_rows.append(
            {
                "donor": donor,
                "target": target,
                "reference": reference,
                **_pair_metrics(donor_df, target, reference),
            }
        )
        nnmf_specs = {
            "nnmf_membrane_cell": ("Membrane", "Cell"),
            "nnmf_nonuclear": ("Cytoplasm", "Membrane", "Cell"),
            "nnmf_all_compartments": COMPARTMENTS,
        }
        for method, compartments in nnmf_specs.items():
            feature_cols = [
                compartment_column(compartment, marker)
                for compartment in compartments
                for marker in (target, reference)
            ]
            matrix = donor_df[feature_cols].to_numpy(float)
            valid = np.isfinite(matrix).all(axis=1) & ~(matrix == 0).all(axis=1)
            valid_idx = np.flatnonzero(valid)
            matrix = matrix[valid]
            method_seeds = seeds if method != "nnmf_all_compartments" else seeds[:1]
            for n_components in (2, 3):
                for seed in method_seeds:
                    labels, components, metadata = fit_nnmf_predict(
                        matrix,
                        n_components=n_components,
                        seed=seed,
                        fit_cap=fit_cap,
                    )
                    run, group_rows, control = _summarize_assignment(
                        donor_df,
                        valid_idx,
                        labels,
                        donor=donor,
                        target=target,
                        reference=reference,
                        method=method,
                        source="raw",
                        n_components=n_components,
                        seed=seed,
                        metadata=metadata,
                        components=components,
                    )
                    runs.append(run)
                    groups.extend(group_rows)
                    assignments[(target, reference, method, n_components, seed)] = control

        for method, source, method_seeds in (
            ("gmm_whole_raw", "raw", seeds),
            ("gmm_whole_redsea", "redsea", seeds[:1]),
        ):
            matrix = donor_df[
                [f"{source}__{target}", f"{source}__{reference}"]
            ].to_numpy(float)
            valid = np.isfinite(matrix).all(axis=1)
            valid_idx = np.flatnonzero(valid)
            matrix = matrix[valid]
            for n_components in (2, 3):
                for seed in method_seeds:
                    labels, components, metadata = fit_gmm_predict(
                        matrix,
                        n_components=n_components,
                        seed=seed,
                        fit_cap=fit_cap,
                    )
                    run, group_rows, control = _summarize_assignment(
                        donor_df,
                        valid_idx,
                        labels,
                        donor=donor,
                        target=target,
                        reference=reference,
                        method=method,
                        source=source,
                        n_components=n_components,
                        seed=seed,
                        metadata=metadata,
                        components=components,
                    )
                    runs.append(run)
                    groups.extend(group_rows)
                    assignments[(target, reference, method, n_components, seed)] = control

        comparisons = []
        for method in ("nnmf_membrane_cell", "nnmf_nonuclear", "gmm_whole_raw"):
            for n_components in (2, 3):
                for first_seed, second_seed in zip(seeds, seeds[1:]):
                    comparisons.append(
                        (
                            f"{method}_k{n_components}_seed",
                            (method, n_components, first_seed),
                            (method, n_components, second_seed),
                        )
                    )
            comparisons.append(
                (
                    f"{method}_component_count",
                    (method, 2, seeds[0]),
                    (method, 3, seeds[0]),
                )
            )
        comparisons.append(
            (
                "nnmf_vs_gmm_k2",
                ("nnmf_membrane_cell", 2, seeds[0]),
                ("gmm_whole_raw", 2, seeds[0]),
            )
        )
        for comparison, first_key, second_key in comparisons:
            first = assignments[(target, reference, *first_key)]
            second = assignments[(target, reference, *second_key)]
            stability.append(
                {
                    "donor": donor,
                    "target": target,
                    "reference": reference,
                    "comparison": comparison,
                    "first_method": first_key[0],
                    "first_components": first_key[1],
                    "first_seed": first_key[2],
                    "second_method": second_key[0],
                    "second_components": second_key[1],
                    "second_seed": second_key[2],
                    **_binary_agreement(first, second),
                }
            )
    return runs, groups, pair_rows, stability


def evaluate_cohort(
    sample_root: Path,
    raw_cells: Path,
    redsea_cells: Path,
    output_dir: Path,
    *,
    donors: Sequence[str],
    fit_cap: int = 50_000,
    seeds: Sequence[int] = (0, 1, 2),
) -> dict[str, pd.DataFrame]:
    """Run the raw/REDSEA method screen for all requested donors."""
    donors = ensure_eligible_donors(
        donors,
        context="RESTORE exploratory cohort evaluation",
    )
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"output already exists: {output_dir}")
    temp_dir = output_dir.with_name(output_dir.name + ".tmp")
    if temp_dir.exists():
        raise FileExistsError(f"temporary output already exists: {temp_dir}")
    temp_dir.mkdir(parents=True)

    audits: list[dict] = []
    runs: list[dict] = []
    groups: list[dict] = []
    pairs: list[dict] = []
    stability: list[dict] = []
    immune_joint: list[dict] = []
    try:
        for donor in donors:
            donor_df, audit = load_validation_sample(
                Path(sample_root), Path(raw_cells), Path(redsea_cells), donor
            )
            donor_runs, donor_groups, donor_pairs, donor_stability = evaluate_donor(
                donor_df, donor, fit_cap=fit_cap, seeds=seeds
            )
            audits.append(audit)
            runs.extend(donor_runs)
            groups.extend(donor_groups)
            pairs.extend(donor_pairs)
            stability.extend(donor_stability)
            immune_joint.extend(
                evaluate_immune_joint(donor_df, donor, fit_cap=fit_cap, seeds=seeds)
            )
            print(
                f"[evaluate] donor {donor}: {len(donor_df):,} cells, "
                f"{len(donor_runs)} model runs",
                flush=True,
            )

        tables = {
            "data_audit": pd.DataFrame(audits),
            "pair_metrics": pd.DataFrame(pairs),
            "model_runs": pd.DataFrame(runs),
            "group_metrics": pd.DataFrame(groups),
            "stability": pd.DataFrame(stability),
            "immune_joint_comparator": pd.DataFrame(immune_joint),
        }
        for name, table in tables.items():
            table.to_csv(temp_dir / f"{name}.csv", index=False)

        run_table = tables["model_runs"]
        summary = (
            run_table.groupby(["target", "reference", "method", "n_components"])
            .agg(
                donor_n=("donor", "nunique"),
                run_n=("donor", "size"),
                directional_mean_fraction=("directional_mean", "mean"),
                directional_median_fraction=("directional_median", "mean"),
                converged_fraction=("converged", "mean"),
                median_coverage=("coverage_fraction", "median"),
                median_target_fold=("target_mean_fold", "median"),
                median_reference_fold=("reference_mean_fold", "median"),
                median_spatial_excess=("spatial_neighbor_excess", "median"),
            )
            .reset_index()
        )
        summary.to_csv(temp_dir / "method_summary.csv", index=False)
        tables["method_summary"] = summary
        temp_dir.rename(output_dir)
        return tables
    except Exception:
        if temp_dir.exists():
            shutil.rmtree(temp_dir)
        raise


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract", help="extract raw compartment sample")
    extract.add_argument("--csv", type=Path, action="append", required=True)
    extract.add_argument("--canonical-cells", type=Path, required=True)
    extract.add_argument("--output", type=Path, required=True)
    extract.add_argument("--donor", action="append", required=True)
    extract.add_argument("--sample-modulus", type=int, default=20)
    extract.add_argument("--threads", type=int, default=8)
    evaluate = sub.add_parser("evaluate", help="evaluate candidate group separators")
    evaluate.add_argument("--sample", type=Path, required=True)
    evaluate.add_argument("--raw-cells", type=Path, required=True)
    evaluate.add_argument("--redsea-cells", type=Path, required=True)
    evaluate.add_argument("--output", type=Path, required=True)
    evaluate.add_argument("--donor", action="append", required=True)
    evaluate.add_argument("--fit-cap", type=int, default=50_000)
    evaluate.add_argument("--seed", type=int, action="append", default=None)
    locked = sub.add_parser(
        "evaluate-locked",
        help="run the frozen QC, balanced NNMF, and REDSEA divisor evaluation",
    )
    locked.add_argument("--sample", type=Path, required=True)
    locked.add_argument("--raw-cells", type=Path, required=True)
    locked.add_argument("--redsea-cells", type=Path, required=True)
    locked.add_argument("--output", type=Path, required=True)
    locked.add_argument("--donor", action="append", required=True)
    locked.add_argument("--fit-cap", type=int, default=50_000)
    locked.add_argument("--seed", type=int, action="append", default=None)
    locked.add_argument("--expression-reviews", type=Path)
    locked.add_argument("--pair-reviews", type=Path)
    locked.add_argument(
        "--pair-set",
        choices=tuple(PAIR_SETS),
        default="baseline",
        help="predeclared candidate pair set to evaluate",
    )

    joint = sub.add_parser(
        "immune-joint",
        help="COMPARATOR: one pooled NNMF per donor over the immune markers (emits no divisor or call)",
    )
    joint.add_argument("--sample", type=Path, required=True)
    joint.add_argument("--raw-cells", type=Path, required=True)
    joint.add_argument("--redsea-cells", type=Path, required=True)
    joint.add_argument("--output", type=Path, required=True)
    joint.add_argument("--donor", action="append", required=True)
    joint.add_argument("--fit-cap", type=int, default=50_000)
    joint.add_argument("--seed", type=int, action="append", default=None)
    joint.add_argument("--n-components", type=int, action="append", default=None)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = _parse_args(argv)
    if args.command == "extract":
        audit = extract_compartment_sample(
            args.csv,
            args.canonical_cells,
            args.output,
            donors=args.donor,
            sample_modulus=args.sample_modulus,
            threads=args.threads,
        )
        print(audit.to_string(index=False))
        return 0
    if args.command == "evaluate":
        tables = evaluate_cohort(
            args.sample,
            args.raw_cells,
            args.redsea_cells,
            args.output,
            donors=args.donor,
            fit_cap=args.fit_cap,
            seeds=tuple(args.seed or [0]),
        )
        print(tables["method_summary"].to_string(index=False))
        return 0
    if args.command == "immune-joint":
        donors = ensure_eligible_donors(args.donor, context="immune-joint comparator")
        rows: list[dict] = []
        for donor in donors:
            donor_df, _audit = load_validation_sample(
                args.sample, args.raw_cells, args.redsea_cells, donor
            )
            donor_rows = evaluate_immune_joint(
                donor_df,
                donor,
                fit_cap=args.fit_cap,
                seeds=tuple(args.seed or [0, 1, 2]),
                n_components=tuple(args.n_components or (4, 5)),
            )
            rows.extend(donor_rows)
            print(f"[immune-joint] donor {donor}: {len(donor_rows)} comparator rows", flush=True)
        table = pd.DataFrame(rows)
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(out, index=False)
        print(f"[immune-joint] wrote {out} ({len(table)} rows)")
        ok = table[table["anchor_recovery"].notna()]
        if len(ok):
            print(
                ok.groupby(["n_components", "target"])["anchor_recovery"]
                .agg(["count", "median", "min"])
                .round(3)
                .to_string()
            )
        return 0
    if args.command == "evaluate-locked":
        tables = evaluate_locked_cohort(
            args.sample,
            args.raw_cells,
            args.redsea_cells,
            args.output,
            donors=args.donor,
            pairs=PAIR_SETS[args.pair_set],
            config=PairValidationConfig(
                fit_cap=args.fit_cap,
                seeds=tuple(args.seed or [0, 1, 2]),
            ),
            expression_reviews=args.expression_reviews,
            pair_reviews=args.pair_reviews,
        )
        print(tables["pair_evaluations"].to_string(index=False))
        return 0
    raise AssertionError(f"unhandled command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
