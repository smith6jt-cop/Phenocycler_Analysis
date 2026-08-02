"""Donor-disjoint fitting and calibration for probabilistic typing bundles.

This module is a development API, not a pipeline stage.  It can read only the
new label namespace authorized by :mod:`phenocycler.refinement_contracts` and
uses development donors for coefficient/bootstrap fitting and held-out
calibration donors for temperature and operating-point diagnostics.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp

from .artifacts import RunManifest
from .hierarchical_typing import (
    FEATURE_SCHEMA_VERSION,
    AdditiveMultinomialModel,
    AssignmentThresholds,
    TypingRegistry,
    balanced_seed_weights,
    classifier_from_registry,
    fit_constrained_classifier,
)
from .marker_registry import MarkerRegistry, load_registry
from .refinement_contracts import (
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
    validate_development_labels,
    validate_evidence_source,
    validate_training_evidence,
)
from .threshold_selection import ThresholdSelectionProvenance
from .typing_model_bundle import (
    MINIMUM_VALID_BOOTSTRAP_FRACTION,
    PRODUCTION_MAX_ABS_WEIGHT,
    TEMPERATURE_MAX,
    TEMPERATURE_MIN,
    FrozenModelEnsemble,
    TypingModelBundle,
    load_typing_model_bundle,
)

TEMPERATURE_BOUNDS = (TEMPERATURE_MIN, TEMPERATURE_MAX)


def _labels_for(
    labels: pd.DataFrame,
    *,
    level: str,
    parent: str | None,
) -> pd.DataFrame:
    parent_value = "" if parent is None else str(parent)
    selected = labels.loc[
        labels["level"].eq(level) & labels["parent"].eq(parent_value),
        ["donor_id", "object_id", "label"],
    ].copy()
    return selected.reset_index(drop=True)


def _require_class_coverage(
    labels: pd.DataFrame,
    classes: Sequence[str],
    *,
    donors: Sequence[str],
    split: str,
    level: str,
    parent: str | None,
) -> pd.DataFrame:
    selected = labels.loc[labels["donor_id"].isin(set(map(str, donors)))].copy()
    observed = set(selected["label"])
    missing = sorted(set(classes) - observed)
    unknown = sorted(observed - set(classes))
    context = level if parent is None else f"{level}/{parent}"
    if missing or unknown:
        raise ValueError(
            f"{context} {split} labels do not cover the model vocabulary: "
            f"missing={missing}, unknown={unknown}"
        )
    return selected.reset_index(drop=True)


def _aligned_probabilities(
    model: AdditiveMultinomialModel,
    evidence: pd.DataFrame,
    labels: pd.DataFrame,
) -> pd.DataFrame:
    probability = model.predict_proba(evidence)
    aligned = labels.merge(
        probability,
        on=["donor_id", "object_id"],
        how="left",
        validate="1:1",
        indicator=True,
    )
    absent = aligned["_merge"].ne("both")
    if absent.any():
        examples = aligned.loc[absent, ["donor_id", "object_id"]].head(5)
        raise ValueError(
            "calibration labels have no marker evidence: "
            f"{examples.to_dict(orient='records')}"
        )
    return aligned.drop(columns="_merge")


def fit_temperature_scaling(
    model: AdditiveMultinomialModel,
    evidence: pd.DataFrame,
    calibration_labels: pd.DataFrame,
) -> float:
    """Fit one positive softmax temperature on donor-balanced held-out labels."""

    aligned = _aligned_probabilities(model, evidence, calibration_labels)
    classes = list(model.classes)
    unknown = sorted(set(aligned["label"]) - set(classes))
    if unknown:
        raise ValueError(f"calibration labels are outside model classes: {unknown}")
    probabilities = aligned.loc[:, classes].to_numpy(float)
    if not np.isfinite(probabilities).all() or (probabilities <= 0).any():
        # Softmax probabilities should be finite and strictly positive. Refuse
        # underflowed inputs instead of silently clipping a broken model.
        raise ValueError("base model produced invalid calibration probabilities")
    if len(classes) == 1:
        # A one-class parent (currently the configured mesenchymal subtype)
        # has an identically-one softmax. Temperature is not identifiable and
        # cannot change either the probability or ranking, so record the
        # neutral transform explicitly rather than pretending to calibrate it.
        return 1.0
    log_probability = np.log(probabilities)
    centered_log_probability = log_probability - log_probability.mean(
        axis=1, keepdims=True
    )
    if np.max(np.abs(centered_log_probability)) <= 1e-8:
        raise ValueError(
            "temperature scaling is unidentified for uniform base probabilities"
        )
    y = aligned["label"].map({name: i for i, name in enumerate(classes)}).to_numpy(int)
    weights = balanced_seed_weights(
        aligned,
        class_col="label",
        donor_col="donor_id",
    ).to_numpy(float)

    def objective(log_temperature: float) -> float:
        inverse = np.exp(-float(log_temperature))
        logits = log_probability * inverse
        log_calibrated = logits - logsumexp(logits, axis=1, keepdims=True)
        return -float(np.sum(weights * log_calibrated[np.arange(len(y)), y]))

    lower_log_temperature, upper_log_temperature = map(
        np.log, TEMPERATURE_BOUNDS
    )
    result = minimize_scalar(
        objective,
        method="bounded",
        bounds=(lower_log_temperature, upper_log_temperature),
        options={"xatol": 1e-8, "maxiter": 500},
    )
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError(f"temperature scaling failed: {result.message}")
    boundary_tolerance = 1e-4
    if (
        result.x <= lower_log_temperature + boundary_tolerance
        or result.x >= upper_log_temperature - boundary_tolerance
    ):
        raise ValueError(
            "temperature calibration optimum lies on the fitted domain boundary"
        )
    curvature_step = 1e-3
    curvature = (
        objective(float(result.x) - curvature_step)
        + objective(float(result.x) + curvature_step)
        - 2.0 * objective(float(result.x))
    )
    if not np.isfinite(curvature) or curvature <= 1e-10:
        raise ValueError("temperature calibration objective is locally unidentified")
    temperature = float(np.exp(result.x))
    if not np.isfinite(temperature) or temperature <= 0:
        raise RuntimeError("temperature scaling returned an invalid value")
    return temperature


def _fit_ensemble(
    evidence: pd.DataFrame,
    registry: TypingRegistry,
    labels: pd.DataFrame,
    splits: DonorSplitManifest,
    *,
    level: str,
    parent: str | None,
    n_bootstrap: int,
    rng: np.random.Generator,
    l2: float,
    max_abs_weight: float,
    maxiter: int,
) -> FrozenModelEnsemble:
    template = classifier_from_registry(registry, level=level, parent=parent)
    development = _require_class_coverage(
        labels,
        template.classes,
        donors=splits.development,
        split="development",
        level=level,
        parent=parent,
    )
    calibration = _require_class_coverage(
        labels,
        template.classes,
        donors=splits.calibration,
        split="calibration",
        level=level,
        parent=parent,
    )
    training_donors = tuple(
        donor for donor in splits.development if donor in set(development["donor_id"])
    )
    calibration_donors = tuple(
        donor for donor in splits.calibration if donor in set(calibration["donor_id"])
    )
    if len(training_donors) < 2:
        raise ValueError(f"{level} fitting requires labels from at least two donors")
    class_donor_support = development.groupby("label", observed=True)[
        "donor_id"
    ].nunique()
    under_supported = {
        class_name: int(class_donor_support.get(class_name, 0))
        for class_name in template.classes
        if int(class_donor_support.get(class_name, 0)) < 2
    }
    if under_supported:
        raise ValueError(
            f"{level} classes require development labels from at least two donors: "
            f"{under_supported}"
        )

    base = fit_constrained_classifier(
        evidence,
        registry,
        development,
        level=level,
        parent=parent,
        l2=l2,
        max_abs_weight=max_abs_weight,
        maxiter=maxiter,
    )
    if not base.converged:
        raise RuntimeError(f"{level} base fit did not converge: {base.message}")
    temperature = fit_temperature_scaling(base, evidence, calibration)

    base_coefficients = tuple(
        tuple(float(item) for item in row)
        for row in np.asarray(base.coefficients, dtype=float)
    )
    bootstrap_coefficients: list[tuple[tuple[float, ...], ...]] = []
    bootstrap_draws: list[tuple[str, ...]] = []
    bootstrap_valid: list[bool] = []
    bootstrap_failures: list[str] = []
    class_set = set(template.classes)
    for _replicate in range(int(n_bootstrap)):
        draw = tuple(
            str(value)
            for value in rng.choice(
                np.asarray(training_donors, dtype=object),
                size=len(training_donors),
                replace=True,
            )
        )
        parts: list[pd.DataFrame] = []
        for draw_index, donor in enumerate(draw):
            part = development.loc[development["donor_id"].eq(donor)].copy()
            part["_bootstrap_group"] = f"{draw_index}:{donor}"
            parts.append(part)
        boot_labels = pd.concat(parts, ignore_index=True)
        missing_classes = sorted(class_set - set(boot_labels["label"]))
        bootstrap_draws.append(draw)
        if missing_classes:
            bootstrap_coefficients.append(base_coefficients)
            bootstrap_valid.append(False)
            bootstrap_failures.append(
                "missing_classes:" + "|".join(missing_classes)
            )
            continue
        fitted = fit_constrained_classifier(
            evidence,
            registry,
            boot_labels,
            level=level,
            parent=parent,
            l2=l2,
            max_abs_weight=max_abs_weight,
            maxiter=maxiter,
            balance_group_col="_bootstrap_group",
        )
        if not fitted.converged:
            bootstrap_coefficients.append(base_coefficients)
            bootstrap_valid.append(False)
            bootstrap_failures.append("nonconverged:" + str(fitted.message))
            continue
        bootstrap_coefficients.append(
            tuple(
                tuple(float(item) for item in row)
                for row in np.asarray(fitted.coefficients, dtype=float)
            )
        )
        bootstrap_valid.append(True)
        bootstrap_failures.append("")

    valid_fraction = float(np.mean(bootstrap_valid))
    if valid_fraction < MINIMUM_VALID_BOOTSTRAP_FRACTION:
        raise RuntimeError(
            "too few valid whole-donor bootstrap fits: "
            f"{sum(bootstrap_valid)}/{n_bootstrap} < "
            f"{MINIMUM_VALID_BOOTSTRAP_FRACTION:.0%}"
        )

    return FrozenModelEnsemble(
        level=level,
        parent=parent,
        classes=base.classes,
        markers=base.markers,
        coefficients=base_coefficients,
        bootstrap_coefficients=tuple(bootstrap_coefficients),
        bootstrap_valid=tuple(bootstrap_valid),
        bootstrap_failures=tuple(bootstrap_failures),
        training_donors=training_donors,
        calibration_donors=calibration_donors,
        bootstrap_draws=tuple(bootstrap_draws),
        temperature=temperature,
        n_development_labels=len(development),
        n_calibration_labels=len(calibration),
        fit_converged=True,
    )


def fit_probabilistic_typing_bundle(
    evidence: pd.DataFrame,
    labels: pd.DataFrame,
    *,
    marker_registry: MarkerRegistry,
    typing_registry: TypingRegistry,
    splits: DonorSplitManifest,
    label_provenance: DevelopmentLabelProvenance,
    evidence_provenance: DevelopmentEvidenceProvenance,
    bundle_version: str,
    thresholds: AssignmentThresholds | None = None,
    n_bootstrap: int = 100,
    random_state: int | np.random.Generator = 0,
    l2: float = 0.02,
    max_abs_weight: float = 8.0,
    maxiter: int = 500,
    require_all_subtype_parents: bool = False,
    created_at_utc: str | None = None,
) -> TypingModelBundle:
    """Fit a complete bundle without reading validation or locked-test labels."""

    if (
        isinstance(n_bootstrap, bool)
        or int(n_bootstrap) != n_bootstrap
        or n_bootstrap < 1
    ):
        raise ValueError("n_bootstrap must be an integer >= 1")
    if not np.isfinite(l2) or float(l2) < 0:
        raise ValueError("l2 must be finite and >= 0")
    if (
        not np.isfinite(max_abs_weight)
        or float(max_abs_weight) <= 0
        or float(max_abs_weight) > PRODUCTION_MAX_ABS_WEIGHT
    ):
        raise ValueError(
            "max_abs_weight must be finite, > 0, and <= "
            f"{PRODUCTION_MAX_ABS_WEIGHT:g}"
        )
    if isinstance(maxiter, bool) or int(maxiter) != maxiter or int(maxiter) < 1:
        raise ValueError("maxiter must be an integer >= 1")
    if thresholds is None:
        raise ValueError(
            "an unfrozen candidate requires explicit assignment thresholds"
        )
    work = validate_development_labels(
        labels,
        provenance=label_provenance,
        splits=splits,
    )
    if label_provenance.source_run_id != splits.source_run_id:
        raise ValueError("label and split source runs differ")
    validate_training_evidence(
        evidence,
        provenance=evidence_provenance,
        splits=splits,
        feature_schema_version=FEATURE_SCHEMA_VERSION,
    )
    missing_evidence_keys = sorted(
        {"donor_id", "object_id"} - set(evidence.columns)
    )
    if missing_evidence_keys:
        raise KeyError(f"training evidence is missing keys: {missing_evidence_keys}")
    label_keys = {
        (str(donor), str(object_id))
        for donor, object_id in work[["donor_id", "object_id"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    }
    evidence_keys = {
        (str(donor), str(object_id))
        for donor, object_id in evidence[["donor_id", "object_id"]]
        .drop_duplicates()
        .itertuples(index=False, name=None)
    }
    missing_keys = sorted(label_keys - evidence_keys)[:5]
    extra_keys = sorted(evidence_keys - label_keys)[:5]
    if missing_keys or extra_keys:
        raise ValueError(
            "training evidence must be the exact labeled-cell projection: "
            f"missing_examples={missing_keys}, extra_examples={extra_keys}"
        )
    if isinstance(random_state, np.random.Generator):
        rng = random_state
        recorded_seed: int | None = None
    else:
        recorded_seed = int(random_state)
        rng = np.random.default_rng(recorded_seed)

    broad_labels = _labels_for(work, level="broad", parent=None)
    if broad_labels.empty:
        raise ValueError("development label bundle has no broad labels")
    if len(classifier_from_registry(typing_registry).classes) < 2:
        raise ValueError("probabilistic broad typing requires at least two classes")
    broad = _fit_ensemble(
        evidence,
        typing_registry,
        broad_labels,
        splits,
        level="broad",
        parent=None,
        n_bootstrap=int(n_bootstrap),
        rng=rng,
        l2=l2,
        max_abs_weight=max_abs_weight,
        maxiter=maxiter,
    )

    subtype_models: list[FrozenModelEnsemble] = []
    configured_parents = tuple(
        dict.fromkeys(str(rule.parent) for rule in typing_registry.subtype_rules)
    )
    final_specific_classes = sorted(
        {
            leaf
            for broad_class in broad.classes
            for leaf in (
                tuple(
                    rule.name
                    for rule in typing_registry.rules_for(
                        level="subtype", parent=broad_class
                    )
                )
                or (broad_class,)
            )
        }
    )
    missing_parents: list[str] = []
    structural_rules_only_parents: list[str] = []
    for parent in configured_parents:
        subtype_template = classifier_from_registry(
            typing_registry, level="subtype", parent=parent
        )
        if len(subtype_template.classes) < 2:
            structural_rules_only_parents.append(parent)
            continue
        parent_labels = _labels_for(work, level="subtype", parent=parent)
        if parent_labels.empty:
            missing_parents.append(parent)
            continue
        subtype_models.append(
            _fit_ensemble(
                evidence,
                typing_registry,
                parent_labels,
                splits,
                level="subtype",
                parent=parent,
                n_bootstrap=int(n_bootstrap),
                rng=rng,
                l2=l2,
                max_abs_weight=max_abs_weight,
                maxiter=maxiter,
            )
        )
    if require_all_subtype_parents and missing_parents:
        raise ValueError(
            "development label bundle has no subtype labels for parents: "
            f"{sorted(missing_parents)}"
        )

    bundle = TypingModelBundle.create(
        bundle_version=bundle_version,
        marker_registry_fingerprint=marker_registry.fingerprint,
        typing_rules_fingerprint=typing_registry.rules_fingerprint,
        split_manifest=splits,
        label_provenance=label_provenance,
        evidence_provenance=evidence_provenance,
        thresholds=thresholds,
        broad=broad,
        subtypes=subtype_models,
        fit_parameters={
            "n_bootstrap": int(n_bootstrap),
            "random_state": recorded_seed,
            "l2": float(l2),
            "max_abs_weight": float(max_abs_weight),
            "maxiter": int(maxiter),
            "require_all_subtype_parents": bool(require_all_subtype_parents),
            "class_priors": "uniform_no_intercept",
            "donor_weighting": "class_then_donor_balanced",
            "minimum_valid_bootstrap_fraction": MINIMUM_VALID_BOOTSTRAP_FRACTION,
            "structural_rules_only_parents": sorted(
                structural_rules_only_parents
            ),
            "missing_probabilistic_subtype_parents": sorted(missing_parents),
            "final_specific_classes": final_specific_classes,
        },
        created_at_utc=created_at_utc,
    )
    bundle.validate_current(
        marker_registry=marker_registry,
        typing_registry=typing_registry,
    )
    return bundle


def acceptance_operating_points(
    reference: Sequence[object],
    probabilities: pd.DataFrame | np.ndarray,
    classes: Sequence[str],
    stability: Sequence[float],
    *,
    probability_thresholds: Sequence[float] = (0.70, 0.75, 0.80, 0.85, 0.90),
    margin_thresholds: Sequence[float] = (0.10, 0.20, 0.25, 0.30),
    stability_thresholds: Sequence[float] = (0.80, 0.85, 0.90, 0.95),
    weights: Sequence[float] | None = None,
    false_positive_cost: float = 1.0,
    false_negative_cost: float = 1.0,
    current_prediction: Sequence[object] | None = None,
    current_assignment_status: Sequence[object] | None = None,
    current_reason: Sequence[object] | None = None,
    candidate_prediction: Sequence[object] | None = None,
    threshold_eligible: Sequence[bool] | None = None,
    valid_evidence: Sequence[bool] | None = None,
) -> pd.DataFrame:
    """Evaluate a production-aware probability/margin/stability grid.

    Wrong accepted calls count as both a false positive for the called class and
    a false negative for a named reference class. Abstentions count only as
    false negatives. In production-aware mode, fixed anchors and ``Other`` keep
    their disposition, and only structurally eligible rows with valid positive
    evidence are rethresholded. No class-prevalence target enters this
    calculation.
    """

    classes = tuple(str(value) for value in classes)
    truth = pd.Series(reference, dtype="string").astype(str)
    values = (
        probabilities.loc[:, list(classes)].to_numpy(float)
        if isinstance(probabilities, pd.DataFrame)
        else np.asarray(probabilities, dtype=float)
    )
    if values.shape != (len(truth), len(classes)):
        raise ValueError("probability matrix shape does not match labels/classes")
    if not np.isfinite(values).all() or (values < 0).any() or (values > 1).any():
        raise ValueError("probabilities must be finite and lie in [0, 1]")
    if not np.allclose(values.sum(axis=1), 1.0, atol=1e-6):
        raise ValueError("each probability row must sum to one")
    production_arguments = (
        current_prediction,
        current_assignment_status,
        current_reason,
        candidate_prediction,
        threshold_eligible,
        valid_evidence,
    )
    production_mode = any(value is not None for value in production_arguments)
    if production_mode and not all(
        value is not None for value in production_arguments
    ):
        raise ValueError(
            "current_prediction, current_assignment_status, current_reason, "
            "candidate_prediction, threshold_eligible, and valid_evidence must "
            "be provided together"
        )
    allowed_truth = set(classes) | ({"Other"} if production_mode else set())
    unknown = sorted(set(truth) - allowed_truth)
    if unknown:
        raise ValueError(f"reference contains unknown classes: {unknown}")
    stability_values = np.asarray(stability, dtype=float)
    if stability_values.shape != (len(truth),):
        raise ValueError("stability must be a vector aligned to reference")
    finite_stability = np.isfinite(stability_values)
    if np.isinf(stability_values).any() or (
        (stability_values[finite_stability] < 0)
        | (stability_values[finite_stability] > 1)
    ).any():
        raise ValueError("finite stability values must lie in [0, 1]")
    sample_weight = (
        np.ones(len(truth), dtype=float)
        if weights is None
        else np.asarray(weights, dtype=float)
    )
    if sample_weight.shape != (len(truth),) or not np.isfinite(sample_weight).all():
        raise ValueError("weights must be a finite vector aligned to reference")
    if (sample_weight < 0).any() or sample_weight.sum() <= 0:
        raise ValueError("weights must be non-negative with positive total")
    for name, thresholds in (
        ("probability_thresholds", probability_thresholds),
        ("margin_thresholds", margin_thresholds),
        ("stability_thresholds", stability_thresholds),
    ):
        threshold_values = np.asarray(thresholds, dtype=float)
        if (
            threshold_values.ndim != 1
            or len(threshold_values) == 0
            or not np.isfinite(threshold_values).all()
            or ((threshold_values < 0) | (threshold_values > 1)).any()
        ):
            raise ValueError(f"{name} must be a nonempty finite sequence in [0, 1]")
    if (
        not np.isfinite(false_positive_cost)
        or not np.isfinite(false_negative_cost)
        or false_positive_cost < 0
        or false_negative_cost < 0
    ):
        raise ValueError("false-positive and false-negative costs must be finite and >= 0")
    order = np.argsort(-values, axis=1, kind="stable")
    best_index = order[:, 0]
    best_probability = values[np.arange(len(values)), best_index]
    second = values[np.arange(len(values)), order[:, 1]] if len(classes) > 1 else 0.0
    margin = best_probability - second
    model_candidate = np.asarray(classes, dtype=object)[best_index]
    if production_mode:
        current = pd.Series(current_prediction, dtype="string").astype(str).to_numpy(
            object
        )
        status = (
            pd.Series(current_assignment_status, dtype="string")
            .fillna("")
            .str.strip()
            .str.lower()
            .to_numpy(object)
        )
        reason = (
            pd.Series(current_reason, dtype="string")
            .fillna("")
            .str.strip()
            .to_numpy(object)
        )
        candidate = pd.Series(candidate_prediction, dtype="string").astype(
            str
        ).to_numpy(object)
        if any(
            value.shape != (len(truth),)
            for value in (current, status, reason, candidate)
        ):
            raise ValueError("production predictions must align to reference")
        eligible_series = pd.Series(threshold_eligible)
        evidence_series = pd.Series(valid_evidence)
        for name, series in (
            ("threshold_eligible", eligible_series),
            ("valid_evidence", evidence_series),
        ):
            if len(series) != len(truth) or series.isna().any() or not series.isin(
                [True, False, 0, 1]
            ).all():
                raise ValueError(f"{name} must contain aligned explicit booleans")
        eligible = eligible_series.astype(bool).to_numpy()
        evidence_pass = evidence_series.astype(bool).to_numpy()
        allowed_status = {"anchor", "inferred", "other", "ambiguous", "unavailable"}
        unknown_status = sorted(set(status) - allowed_status)
        if unknown_status:
            raise ValueError(
                f"current_assignment_status contains unknown values: {unknown_status}"
            )
        expected_disposition = {
            "other": "Other",
            "ambiguous": "Ambiguous",
            "unavailable": "Unavailable",
        }
        for disposition_status, expected_label in expected_disposition.items():
            mismatched = (status == disposition_status) & (current != expected_label)
            if mismatched.any():
                raise ValueError(
                    f"{disposition_status} status requires {expected_label!r} prediction"
                )
        named_current = np.isin(current, classes)
        invalid_named_status = named_current & ~np.isin(status, ["anchor", "inferred"])
        if invalid_named_status.any():
            raise ValueError("fixed named predictions must be anchors or inferred calls")
        if ((status == "anchor") & ~named_current).any():
            raise ValueError("anchor status requires a modeled named prediction")
        if ((status == "inferred") & ~named_current).any():
            raise ValueError("inferred status requires a modeled named prediction")
        rescue_reasons = {
            "nonunique_model_winner",
            "valid_positive_support_required",
            "posterior_margin_or_stability_failed",
        }
        should_be_eligible = (status == "inferred") | (
            (status == "ambiguous") & np.isin(reason, list(rescue_reasons))
        )
        if not np.array_equal(eligible, should_be_eligible):
            raise ValueError(
                "threshold_eligible differs from the assignment status/reason contract"
            )
        invalid_candidate = eligible & ~np.isin(candidate, classes)
        if invalid_candidate.any():
            raise ValueError("every threshold-eligible row needs a modeled candidate")
        if not np.array_equal(candidate[eligible], model_candidate[eligible]):
            raise ValueError(
                "candidate_prediction differs from the probability argmax"
            )
    else:
        current = np.full(len(truth), "Ambiguous", dtype=object)
        candidate = model_candidate
        eligible = np.ones(len(truth), dtype=bool)
        evidence_pass = np.ones(len(truth), dtype=bool)
    required_stability = eligible if production_mode else np.ones(len(truth), dtype=bool)
    if (~finite_stability & required_stability).any():
        raise ValueError("threshold-eligible rows require finite stability")
    total = float(sample_weight.sum())
    rows: list[dict[str, float]] = []
    for probability_threshold in probability_thresholds:
        for margin_threshold in margin_thresholds:
            for stability_threshold in stability_thresholds:
                accepted_lane = (
                    eligible
                    & evidence_pass
                    & (best_probability >= float(probability_threshold))
                    & (margin >= float(margin_threshold))
                    & (margin > 1e-12)
                    & (stability_values >= float(stability_threshold))
                )
                final_prediction = current.copy()
                final_prediction[eligible] = "Ambiguous"
                final_prediction[accepted_lane] = candidate[accepted_lane]
                named_call = np.isin(final_prediction, classes)
                truth_named = np.isin(truth.to_numpy(object), classes)
                correct_named = named_call & (
                    final_prediction == truth.to_numpy(object)
                )
                wrong_named = named_call & ~correct_named
                false_negative = truth_named & ~correct_named
                correct_weight = float(sample_weight[correct_named].sum())
                wrong_weight = float(sample_weight[wrong_named].sum())
                false_negative_weight = float(sample_weight[false_negative].sum())
                called_weight = float(sample_weight[named_call].sum())
                eligible_weight = float(sample_weight[eligible].sum())
                rescue_called_weight = float(sample_weight[accepted_lane].sum())
                rescue_correct_weight = float(
                    sample_weight[
                        accepted_lane
                        & (candidate == truth.to_numpy(object))
                    ].sum()
                )
                rows.append(
                    {
                        "probability_threshold": float(probability_threshold),
                        "margin_threshold": float(margin_threshold),
                        "stability_threshold": float(stability_threshold),
                        "coverage": called_weight / total,
                        "selective_error": (
                            wrong_weight / called_weight if called_weight else np.nan
                        ),
                        "correct_call_rate": correct_weight / total,
                        "false_positive_weight": wrong_weight,
                        "false_negative_weight": false_negative_weight,
                        "rescue_eligible_weight": eligible_weight,
                        "rescue_called_weight": rescue_called_weight,
                        "rescue_coverage": (
                            rescue_called_weight / eligible_weight
                            if eligible_weight
                            else np.nan
                        ),
                        "rescue_selective_error": (
                            1.0 - rescue_correct_weight / rescue_called_weight
                            if rescue_called_weight
                            else np.nan
                        ),
                        "error_cost": (
                            float(false_positive_cost) * wrong_weight
                            + float(false_negative_cost) * false_negative_weight
                        )
                        / total,
                    }
                )
    return pd.DataFrame(rows)


def choose_operating_point(
    table: pd.DataFrame,
    *,
    maximum_selective_error: float,
    minimum_coverage: float = 0.0,
    maximum_rescue_selective_error: float | None = None,
    minimum_rescue_coverage: float | None = None,
) -> pd.Series:
    """Select the lowest-cost feasible frozen operating point deterministically.

    Rescue-specific constraints are optional for legacy stand-alone grids but
    should be supplied for production-aware grids. They prevent fixed anchors
    from satisfying total-coverage/error gates while the probabilistic rescue
    lane makes no useful calls.
    """

    required = {
        "probability_threshold",
        "margin_threshold",
        "stability_threshold",
        "coverage",
        "selective_error",
        "error_cost",
    }
    missing = sorted(required - set(table.columns))
    if missing:
        raise KeyError(f"operating-point table is missing columns: {missing}")
    for value, name in (
        (maximum_selective_error, "maximum_selective_error"),
        (minimum_coverage, "minimum_coverage"),
        (maximum_rescue_selective_error, "maximum_rescue_selective_error"),
        (minimum_rescue_coverage, "minimum_rescue_coverage"),
    ):
        if value is not None and (
            not np.isfinite(value) or float(value) < 0 or float(value) > 1
        ):
            raise ValueError(f"{name} must lie in [0, 1]")
    feasible_mask = (
        table["selective_error"].le(float(maximum_selective_error))
        & table["coverage"].ge(float(minimum_coverage))
    )
    rescue_constraints = {
        "rescue_selective_error": maximum_rescue_selective_error,
        "rescue_coverage": minimum_rescue_coverage,
    }
    requested_rescue_columns = {
        column for column, value in rescue_constraints.items() if value is not None
    }
    missing_rescue = sorted(requested_rescue_columns - set(table.columns))
    if missing_rescue:
        raise KeyError(
            "operating-point table is missing rescue columns: "
            f"{missing_rescue}"
        )
    if maximum_rescue_selective_error is not None:
        feasible_mask &= table["rescue_selective_error"].le(
            float(maximum_rescue_selective_error)
        )
    if minimum_rescue_coverage is not None:
        feasible_mask &= table["rescue_coverage"].ge(
            float(minimum_rescue_coverage)
        )
    feasible = table.loc[feasible_mask].copy()
    if feasible.empty:
        raise ValueError("no operating point meets the frozen error/coverage constraints")
    tie_break_columns = ["error_cost", "coverage"]
    ascending = [True, False]
    if "rescue_coverage" in feasible:
        tie_break_columns.append("rescue_coverage")
        ascending.append(False)
    if "rescue_selective_error" in feasible:
        tie_break_columns.append("rescue_selective_error")
        ascending.append(True)
    tie_break_columns.extend(
        [
            "probability_threshold",
            "margin_threshold",
            "stability_threshold",
        ]
    )
    ascending.extend([False, False, False])
    feasible = feasible.sort_values(
        tie_break_columns,
        ascending=ascending,
        kind="mergesort",
    )
    return feasible.iloc[0]


def _read_frame(path: Path | str) -> pd.DataFrame:
    source = Path(path).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    if source.suffix.lower() in {".parquet", ".pq"}:
        return pd.read_parquet(source)
    if source.suffix.lower() in {".csv", ".tsv"}:
        return pd.read_csv(source, sep="\t" if source.suffix.lower() == ".tsv" else ",")
    raise ValueError(f"unsupported table format: {source.suffix}")


def load_labeled_evidence(
    evidence_path: Path | str,
    labels: pd.DataFrame,
) -> pd.DataFrame:
    """Project only labeled cells from a Parquet evidence file/dataset.

    This avoids materializing a many-million-cell production donor set during
    model development. CSV/TSV is supported for small pilots; production should
    use the partitioned Parquet marker-evidence artifact.
    """

    source = Path(evidence_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(source)
    required = {"donor_id", "object_id"}
    missing = sorted(required - set(labels.columns))
    if missing:
        raise KeyError(f"labels are missing evidence keys: {missing}")
    keys = labels.loc[:, ["donor_id", "object_id"]].copy()
    keys["donor_id"] = keys["donor_id"].astype(str)
    keys["object_id"] = keys["object_id"].astype(str)
    keys = keys.drop_duplicates(ignore_index=True)
    if source.is_file() and source.suffix.lower() in {".csv", ".tsv"}:
        evidence = _read_frame(source)
        for column in ("donor_id", "object_id"):
            if column not in evidence:
                raise KeyError(f"evidence is missing {column!r}")
            evidence[column] = evidence[column].astype(str)
        selected = evidence.merge(
            keys,
            on=["donor_id", "object_id"],
            how="inner",
            validate="m:1",
        )
    else:
        if source.is_dir():
            parquet_files = sorted(source.rglob("*.parquet"))
            if not parquet_files:
                raise FileNotFoundError(f"no Parquet evidence below {source}")
            parquet_source = (source / "**" / "*.parquet").as_posix()
        elif source.suffix.lower() in {".parquet", ".pq"}:
            parquet_source = source.as_posix()
        else:
            raise ValueError(f"unsupported evidence source: {source}")
        connection = duckdb.connect()
        try:
            connection.register("development_label_keys", keys)
            selected = connection.execute(
                """
                SELECT e.*
                FROM read_parquet(?, hive_partitioning=true) e
                JOIN development_label_keys k
                  ON CAST(e.donor_id AS VARCHAR) = k.donor_id
                 AND CAST(e.object_id AS VARCHAR) = k.object_id
                """,
                [parquet_source],
            ).fetchdf()
        finally:
            connection.close()
        selected["donor_id"] = selected["donor_id"].astype(str)
        selected["object_id"] = selected["object_id"].astype(str)
    observed = selected.loc[:, ["donor_id", "object_id"]].drop_duplicates()
    aligned = keys.merge(
        observed,
        on=["donor_id", "object_id"],
        how="left",
        indicator=True,
        validate="1:1",
    )
    absent = aligned["_merge"].ne("both")
    if absent.any():
        examples = aligned.loc[absent, ["donor_id", "object_id"]].head(5)
        raise ValueError(
            "development labels have no marker evidence: "
            f"{examples.to_dict(orient='records')}"
        )
    return selected


def _freeze_main(argv: Sequence[str]) -> int:
    """Attach reviewed threshold provenance to a preliminary fit without refitting."""

    parser = argparse.ArgumentParser(
        description="Freeze calibration-selected thresholds onto a preliminary fit."
    )
    parser.add_argument("--model-bundle", type=Path, required=True)
    parser.add_argument("--threshold-provenance", type=Path, required=True)
    parser.add_argument("--marker-registry", type=Path)
    parser.add_argument("--typing-rules", type=Path)
    parser.add_argument("--bundle-version", required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    markers = load_registry(args.marker_registry)
    rules = TypingRegistry.from_marker_registry(
        markers,
        typing_rules=args.typing_rules,
    )
    preliminary = load_typing_model_bundle(
        args.model_bundle,
        marker_registry=markers,
        typing_registry=rules,
    )
    provenance = ThresholdSelectionProvenance.read_json(
        args.threshold_provenance
    )
    frozen = preliminary.freeze_threshold_selection(
        provenance,
        bundle_version=args.bundle_version,
    )
    frozen.validate_current(marker_registry=markers, typing_registry=rules)
    output = frozen.write_json(args.out)
    print(
        f"wrote {output} content_id={frozen.content_id} "
        f"model_fit_content_id={frozen.model_fit_content_id} thresholds=frozen"
    )
    return 0


def main(argv=None) -> int:
    """Build one immutable model bundle from manifest-authorized labels."""

    raw = list(sys.argv[1:] if argv is None else argv)
    if raw and raw[0] == "freeze":
        return _freeze_main(raw[1:])

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fit", nargs="?", choices=("fit",), default="fit")
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--source-run-manifest", type=Path, required=True)
    parser.add_argument("--marker-registry", type=Path)
    parser.add_argument("--typing-rules", type=Path)
    parser.add_argument("--label-provenance", type=Path, required=True)
    parser.add_argument("--evidence-provenance", type=Path, required=True)
    parser.add_argument("--bundle-version", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--l2", type=float, default=0.02)
    parser.add_argument("--max-abs-weight", type=float, default=8.0)
    parser.add_argument("--maxiter", type=int, default=500)
    parser.add_argument(
        "--unfrozen-candidate",
        action="store_true",
        required=True,
        help="build a preliminary non-promotable candidate with explicit thresholds",
    )
    parser.add_argument("--probability-threshold", type=float)
    parser.add_argument("--margin-threshold", type=float)
    parser.add_argument("--stability-threshold", type=float)
    parser.add_argument("--require-all-subtypes", action="store_true")
    args = parser.parse_args(raw)

    raw_thresholds = (
        args.probability_threshold,
        args.margin_threshold,
        args.stability_threshold,
    )
    if any(value is None for value in raw_thresholds):
        parser.error(
            "--unfrozen-candidate requires --probability-threshold, "
            "--margin-threshold, and --stability-threshold"
        )
    assignment_thresholds = AssignmentThresholds(
        inferred_probability=args.probability_threshold,
        inferred_margin=args.margin_threshold,
        inferred_stability=args.stability_threshold,
    )

    splits = DonorSplitManifest.read_json(args.splits)
    labels = _read_frame(args.labels)
    provenance = DevelopmentLabelProvenance.read_json(args.label_provenance)
    labels = validate_development_labels(
        labels,
        provenance=provenance,
        splits=splits,
    )
    evidence_provenance = DevelopmentEvidenceProvenance.read_json(
        args.evidence_provenance
    )
    source_run_manifest = RunManifest.read_json(args.source_run_manifest)
    validate_evidence_source(
        args.evidence,
        provenance=evidence_provenance,
        splits=splits,
        source_run_manifest=source_run_manifest,
    )
    evidence = load_labeled_evidence(args.evidence, labels)
    validate_training_evidence(
        evidence,
        provenance=evidence_provenance,
        splits=splits,
        feature_schema_version=FEATURE_SCHEMA_VERSION,
    )
    markers = load_registry(args.marker_registry)
    rules = TypingRegistry.from_marker_registry(
        markers,
        typing_rules=args.typing_rules,
    )
    bundle = fit_probabilistic_typing_bundle(
        evidence,
        labels,
        marker_registry=markers,
        typing_registry=rules,
        splits=splits,
        label_provenance=provenance,
        evidence_provenance=evidence_provenance,
        bundle_version=args.bundle_version,
        thresholds=assignment_thresholds,
        n_bootstrap=args.bootstrap,
        random_state=args.seed,
        l2=args.l2,
        max_abs_weight=args.max_abs_weight,
        maxiter=args.maxiter,
        require_all_subtype_parents=args.require_all_subtypes,
    )
    output = bundle.write_json(args.out)
    print(
        f"wrote {output} content_id={bundle.content_id} "
        f"model_fit_content_id={bundle.model_fit_content_id} "
        "thresholds=unfrozen "
        f"bootstrap={bundle.broad.bootstrap_replicates}"
    )
    return 0


__all__ = [
    "acceptance_operating_points",
    "choose_operating_point",
    "fit_probabilistic_typing_bundle",
    "fit_temperature_scaling",
    "load_labeled_evidence",
    "main",
]


if __name__ == "__main__":
    raise SystemExit(main())
