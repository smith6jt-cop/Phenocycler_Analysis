from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    InputArtifact,
    RunManifest,
    StageManifest,
    StaleArtifactError,
    snapshot_parquet_dataset,
)
from phenocycler.candidate_evaluation import (
    CandidateEvaluationManifest,
    candidate_assignment_run_id,
)
from phenocycler.evaluation_reporting import (
    EvaluationTargetReview,
    build_split_evaluation_report,
)
from phenocycler.hierarchical_typing import (
    FEATURE_SCHEMA_VERSION,
    AssignmentThresholds,
    TypingRegistry,
    classifier_from_registry,
)
from phenocycler.marker_registry import MarkerRegistry, load_registry
from phenocycler.refinement_contracts import (
    DEVELOPMENT_LABEL_PURPOSE,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
    ReviewLabelProvenance,
    blinded_review_sample_frame,
    build_review_sampling_frame,
)
from phenocycler.threshold_selection import (
    CalibrationDecisionRow,
    CalibrationTargetDecisionInput,
    FrozenOperatingPoint,
    ThresholdSelectionPolicy,
    ThresholdSelectionProvenance,
    _decision_inputs_fingerprint,
    _derive_operating_grid,
    operating_grid_fingerprint,
    threshold_grid_content_id,
)
from phenocycler.typing_model_bundle import (
    EvaluatedDecisionRow,
    FrozenModelEnsemble,
    PromotionGatePolicy,
    SplitEvaluationReport,
    TypingModelBundle,
    _target_simultaneous_intervals,
    _two_stage_survey_bootstrap,
)
from phenocycler.typing_stage import type_candidate_cells

_REVIEW_BLINDING_KEY = "evaluation-test-review-key"


def _specific_reference(object_id: str) -> str:
    if object_id.startswith("s-"):
        return object_id.removeprefix("s-")
    if object_id.startswith("x-"):
        return object_id.removeprefix("x-")
    if object_id == "b-Neural":
        return "Neural"
    return "Other"


def _ingest_output_content_id(source_run: RunManifest) -> str:
    references = [reference for reference in source_run.stages if reference.stage == "ingest"]
    assert len(references) == 1
    return references[0].output_content_sha256


@dataclass(frozen=True)
class _TargetSpec:
    level: str
    parent: str | None
    classes: tuple[str, ...]


@dataclass(frozen=True)
class EvaluationFixture:
    root: Path
    evidence: Path
    splits: DonorSplitManifest
    markers: MarkerRegistry
    rules: TypingRegistry
    bundle: TypingModelBundle
    source_run: RunManifest
    references: dict[tuple[str, str, str], str]

    def candidate(
        self,
        name: str,
        *,
        split_name: str = "validation",
        transform=None,
    ) -> tuple[Path, pd.DataFrame]:
        donors = tuple(getattr(self.splits, split_name))
        run_id = candidate_assignment_run_id(
            model_bundle_content_id=self.bundle.content_id,
            split_manifest_content_id=self.splits.content_id,
            split_name=split_name,
            source_run_id=self.splits.source_run_id,
            source_run_manifest_content_id=(
                self.bundle.evidence_provenance.source_run_manifest_content_id
            ),
            source_stage_manifest_content_id=(
                self.bundle.evidence_provenance.source_stage_manifest_content_id
            ),
            source_artifact_content_id=(self.bundle.evidence_provenance.source_artifact_content_id),
            marker_registry_fingerprint=self.markers.fingerprint,
            typing_rules_fingerprint=self.rules.rules_fingerprint,
            donors=donors,
        )
        frames: list[pd.DataFrame] = []
        for donor in donors:
            evidence = pd.read_parquet(self.evidence / f"donor_id={donor}" / "evidence.parquet")
            assignments = type_candidate_cells(
                evidence,
                model_bundle=self.bundle,
                marker_registry=self.markers,
                typing_registry=self.rules,
            )
            assignments["candidate_assignment_run_id"] = run_id
            assignments["candidate_evaluation_split_name"] = split_name
            frames.append(assignments)
        assignments = pd.concat(frames, ignore_index=True)
        if transform is not None:
            assignments = transform(assignments.copy())
        output_root = self.root / f"{name}-{split_name}-assignments"
        for donor in donors:
            partition = output_root / f"donor_id={donor}"
            partition.mkdir(parents=True)
            assignments.loc[assignments["donor_id"].eq(donor)].to_parquet(
                partition / "assignments.parquet", index=False
            )
        snapshot = snapshot_parquet_dataset(output_root, expected_donors=donors)
        manifest = CandidateEvaluationManifest.create(
            model_bundle=self.bundle,
            split_name=split_name,
            source_run_manifest=self.source_run,
            marker_registry=self.markers,
            typing_registry=self.rules,
            output=snapshot,
        )
        manifest_path = self.root / f"{name}-{split_name}.manifest.json"
        manifest.write_json(manifest_path)
        return manifest_path, assignments

    def reviews(
        self,
        assignments: pd.DataFrame,
        *,
        split_name: str,
        candidate_manifest: CandidateEvaluationManifest,
    ) -> tuple[EvaluationTargetReview, ...]:
        reviews = []
        targets = (
            self.bundle.broad,
            _TargetSpec(
                level="specific",
                parent=None,
                classes=tuple(self.bundle.fit_parameters["final_specific_classes"]),
            ),
            *self.bundle.subtypes,
        )
        for ensemble in targets:
            sample = build_review_sampling_frame(
                assignments,
                splits=self.splits,
                split_name=split_name,
                candidate_assignment_run_id=candidate_manifest.assignment_run_id,
                level=ensemble.level,
                parent=ensemble.parent,
                blinding_key=_REVIEW_BLINDING_KEY,
            )
            references = []
            for donor, object_id in sample[["donor_id", "object_id"]].itertuples(
                index=False, name=None
            ):
                if ensemble.level == "broad":
                    references.append(
                        self.references.get(
                            (split_name, donor, object_id),
                            "Immune",
                        )
                        if object_id.startswith("b-")
                        else "Immune"
                    )
                elif ensemble.level == "subtype":
                    references.append(
                        self.references[(split_name, donor, object_id)]
                        if object_id.startswith("s-")
                        else "Other"
                    )
                else:
                    references.append(_specific_reference(object_id))
            reviewer_frame = blinded_review_sample_frame(sample)
            reviews.append(
                _review_for_sample(
                    sample,
                    reviewer_frame,
                    references,
                    ensemble=ensemble,
                    splits=self.splits,
                    split_name=split_name,
                    candidate_manifest=candidate_manifest,
                    bundle=self.bundle,
                    review_context_artifact_content_id=(
                        _ingest_output_content_id(self.source_run)
                    ),
                )
            )
        return tuple(reviews)


def _fit_parameters(rules: TypingRegistry) -> dict[str, object]:
    structural = sorted(
        {
            str(rule.parent)
            for rule in rules.subtype_rules
            if len(rules.rules_for(level="subtype", parent=str(rule.parent))) == 1
        }
    )
    return {
        "n_bootstrap": 1,
        "random_state": 7,
        "l2": 0.02,
        "max_abs_weight": 8.0,
        "maxiter": 500,
        "require_all_subtype_parents": False,
        "structural_rules_only_parents": structural,
        "missing_probabilistic_subtype_parents": ["Endothelial", "Epithelial"],
        "final_specific_classes": sorted(
            {
                leaf
                for broad_class in classifier_from_registry(rules).classes
                for leaf in (
                    tuple(
                        rule.name
                        for rule in rules.rules_for(
                            level="subtype", parent=broad_class
                        )
                    )
                    or (broad_class,)
                )
            }
        ),
        "minimum_valid_bootstrap_fraction": 0.8,
        "class_priors": "uniform_no_intercept",
        "donor_weighting": "class_then_donor_balanced",
    }


def _ensemble(
    rules: TypingRegistry,
    *,
    level: str,
    parent: str | None = None,
) -> FrozenModelEnsemble:
    model = classifier_from_registry(rules, level=level, parent=parent)
    coefficients = tuple(
        tuple(float(item) for item in row) for row in np.asarray(model.coefficients, dtype=float)
    )
    return FrozenModelEnsemble(
        level=level,
        parent=parent,
        classes=model.classes,
        markers=model.markers,
        coefficients=coefficients,
        bootstrap_coefficients=(coefficients,),
        training_donors=("D1", "D2"),
        calibration_donors=("D3",),
        bootstrap_draws=(("D1", "D2"),),
        bootstrap_valid=(True,),
        bootstrap_failures=("",),
        temperature=1.0,
        n_development_labels=20,
        n_calibration_labels=10,
    )


def _fixture_threshold_policy() -> ThresholdSelectionPolicy:
    return ThresholdSelectionPolicy(
        policy_version="evaluation-fixture-v1",
        false_positive_cost=1.0,
        false_negative_cost=1.0,
        maximum_selective_error=1.0,
        minimum_coverage=0.0,
        maximum_rescue_selective_error=1.0,
        minimum_rescue_coverage=0.0,
    )


def _freeze_bundle(bundle: TypingModelBundle) -> TypingModelBundle:
    assignment_run_id = "9" * 64
    targets = (bundle.broad, *bundle.subtypes)
    inputs = []
    reviews = []
    for target_index, ensemble in enumerate(targets):
        rows = []
        for class_index, class_name in enumerate(ensemble.classes):
            probabilities = np.full(len(ensemble.classes), 0.01 / (len(ensemble.classes) - 1))
            probabilities[class_index] = 0.99
            rows.append(
                CalibrationDecisionRow(
                    donor_id=bundle.split_manifest.calibration[0],
                        object_id=f"calibration-{target_index}-{class_index}",
                        reference_label=class_name,
                        reviewer_1="reviewer-a",
                        reviewer_2="reviewer-b",
                        adjudicated=True,
                        root_cause="independent evidence review",
                        evidence_sufficient=True,
                        sampling_stratum=f"S{target_index}-{class_index}",
                        sampling_method="hmac_sha256_candidate_bound_v2",
                        sampling_salt="calibration-fixture-salt",
                        stratum_population=1,
                        stratum_sampled=1,
                        sampling_weight=1.0,
                        selection_probability=1.0,
                    current_prediction=class_name,
                    current_assignment_status="inferred",
                    current_reason="posterior_margin_stability_pass",
                    candidate_prediction=class_name,
                    threshold_eligible=True,
                    valid_evidence=True,
                    stability=1.0,
                    probabilities=tuple(float(value) for value in probabilities),
                )
            )
        target_input = CalibrationTargetDecisionInput(
            level=ensemble.level,
            parent="" if ensemble.parent is None else str(ensemble.parent),
            classes=ensemble.classes,
            rows=tuple(rows),
        )
        inputs.append(target_input)
        sample = pd.DataFrame(
            [
                {
                    "donor_id": row.donor_id,
                    "object_id": row.object_id,
                    "sampling_stratum": row.sampling_stratum,
                    "sampling_method": row.sampling_method,
                    "sampling_salt": row.sampling_salt,
                    "stratum_population": row.stratum_population,
                    "stratum_sampled": row.stratum_sampled,
                    "sampling_weight": row.sampling_weight,
                    "selection_probability": row.selection_probability,
                }
                for row in rows
            ]
        )
        labels = pd.DataFrame(
            [
                {
                    "donor_id": row.donor_id,
                    "object_id": row.object_id,
                    "reference_label": row.reference_label,
                    "reviewer_1": row.reviewer_1,
                    "reviewer_2": row.reviewer_2,
                    "adjudicated": row.adjudicated,
                    "root_cause": row.root_cause,
                    "evidence_sufficient": row.evidence_sufficient,
                }
                for row in rows
            ]
        )
        calibration_reviewer_frame = blinded_review_sample_frame(sample)
        calibration_reviewer_frame.attrs[
            "review_context_artifact_content_id"
        ] = bundle.evidence_provenance.source_artifact_content_id
        reviews.append(
            ReviewLabelProvenance.from_labels(
                labels,
                bundle_id=f"calibration-{ensemble.level}-{ensemble.parent or 'broad'}",
                splits=bundle.split_manifest,
                split_name="calibration",
                sample_frame=sample,
                reviewer_frame=calibration_reviewer_frame,
                assignment_run_id=assignment_run_id,
                typing_model_bundle_content_id=bundle.content_id,
                candidate_evaluation_manifest_content_id="a" * 64,
                candidate_assignment_output_content_id="b" * 64,
                level=ensemble.level,
                parent=ensemble.parent,
            )
        )
    policy = _fixture_threshold_policy()
    probability_thresholds = (bundle.thresholds.inferred_probability,)
    margin_thresholds = (bundle.thresholds.inferred_margin,)
    stability_thresholds = (bundle.thresholds.inferred_stability,)
    target_inputs = tuple(inputs)
    grid = _derive_operating_grid(
        target_inputs,
        probability_thresholds=probability_thresholds,
        margin_thresholds=margin_thresholds,
        stability_thresholds=stability_thresholds,
        policy=policy,
    )
    points = tuple(FrozenOperatingPoint.from_series(row) for _, row in grid.iterrows())
    selected = grid.iloc[0]
    provenance = ThresholdSelectionProvenance(
        selection_version="evaluation-fixture-v1",
        preliminary_bundle_content_id=bundle.content_id,
        model_fit_content_id=bundle.model_fit_content_id,
        preliminary_assignment_thresholds={
            name: float(getattr(bundle.thresholds, name))
            for name in (
                "inferred_probability",
                "inferred_margin",
                "inferred_stability",
                "anchor_probability",
                "negative_probability",
            )
        },
        split_manifest=bundle.split_manifest,
        calibration_candidate_manifest_content_id="a" * 64,
        calibration_candidate_output_content_id="b" * 64,
        calibration_review_provenance=tuple(reviews),
        assignment_run_id=assignment_run_id,
        calibration_donors=bundle.split_manifest.calibration,
        calibration_inputs=target_inputs,
        calibration_row_count=sum(len(value.rows) for value in target_inputs),
        calibration_frame_fingerprint=_decision_inputs_fingerprint(target_inputs),
        probability_thresholds=probability_thresholds,
        margin_thresholds=margin_thresholds,
        stability_thresholds=stability_thresholds,
        operating_grid=points,
        operating_grid_row_count=len(points),
        operating_grid_fingerprint=operating_grid_fingerprint(grid),
        policy=policy,
        selected_probability_threshold=float(selected["probability_threshold"]),
        selected_margin_threshold=float(selected["margin_threshold"]),
        selected_stability_threshold=float(selected["stability_threshold"]),
        selected_metrics={
            name: float(selected[name])
            for name in (
                "coverage",
                "selective_error",
                "correct_call_rate",
                "false_positive_weight",
                "false_negative_weight",
                "rescue_eligible_weight",
                "rescue_called_weight",
                "rescue_coverage",
                "rescue_selective_error",
                "error_cost",
            )
        },
    )
    return bundle.freeze_threshold_selection(provenance)


@pytest.fixture
def evaluation_fixture(tmp_path: Path) -> EvaluationFixture:
    threshold_policy = _fixture_threshold_policy()
    promotion_policy = _report_policy()
    splits = DonorSplitManifest(
        split_version="heldout-evaluation-v1",
        source_run_id="source-run",
        cohort_content_id="source-cohort",
        development=("D1", "D2"),
        calibration=("D3",),
        validation=("V1", "V2", "V3", "V4", "V5"),
        locked_test=("L1", "L2", "L3", "L4", "L5"),
        review_n_per_stratum={
            "calibration": 1,
            "validation": 2,
            "locked_test": 2,
        },
        review_blinding_key_sha256=hashlib.sha256(
            _REVIEW_BLINDING_KEY.encode()
        ).hexdigest(),
        threshold_selection_policy_content_id=threshold_policy.content_id,
        threshold_grid_content_id=threshold_grid_content_id(
            probability_thresholds=(0.70,),
            margin_thresholds=(0.20,),
            stability_thresholds=(0.90,),
        ),
        promotion_policy_content_id=promotion_policy.content_id,
    )
    markers = load_registry()
    rules = TypingRegistry.from_marker_registry(markers)
    broad_classes = classifier_from_registry(rules).classes
    subtype_classes = classifier_from_registry(rules, level="subtype", parent="Immune").classes
    evidence = tmp_path / "evidence"
    for donor in splits.donors:
        partition = evidence / f"donor_id={donor}"
        partition.mkdir(parents=True)
        object_ids = [
            *(f"b-{label}" for label in broad_classes),
            *(f"s-{label}" for label in subtype_classes),
            "x-Lymphatic_endothelial",
            "x-CD34_endothelial",
            "x-Beta",
            "x-Alpha",
            "x-Delta",
            "x-SMA_positive_mesenchymal",
        ]
        frame = pd.DataFrame(
            {
                "donor_id": donor,
                "object_id": object_ids,
                "qc_analysis_eligible": True,
            }
        )
        for marker in markers.active_markers:
            frame[f"{marker.name}__state"] = "negative"
            frame[f"{marker.name}__log2_threshold_ratio"] = -1.0
            frame[f"{marker.name}__model_valid"] = True
        positive_markers = {
            "b-Immune": ("CD4", "CD8", "CD11b", "CD206"),
            "b-Endothelial": ("CD34", "Podoplanin"),
            "b-Epithelial": ("Ker8_18", "Keratin_5", "EpCAM"),
            "b-Neural": ("B3TUBB",),
            "b-Mesenchymal": ("Vimentin",),
            "s-T_cell": ("CD3e", "CD4", "CD8"),
            "s-B_lineage": ("CD79a",),
            "s-Macrophage": ("CD68",),
            "s-Myeloid_APC": ("CD11c", "CD11b", "CD206"),
            "s-Neutrophil": ("MPO",),
            "x-Lymphatic_endothelial": ("CD31", "Podoplanin"),
            "x-CD34_endothelial": ("CD31", "CD34"),
            "x-Beta": ("INS",),
            "x-Alpha": ("GCG",),
            "x-Delta": ("E_cadherin", "SST"),
            "x-SMA_positive_mesenchymal": ("SMA",),
        }
        for row_index, object_id in enumerate(frame["object_id"]):
            for marker_name in positive_markers[object_id]:
                frame.loc[row_index, f"{marker_name}__state"] = "positive"
                frame.loc[row_index, f"{marker_name}__log2_threshold_ratio"] = 2.0
            if object_id == "s-T_cell":
                # Numeric CD3e support below the authoritative-anchor threshold
                # exercises the probabilistic rescue lane without creating an
                # anchor conflict, while CD4/CD8 make T_cell the clear winner.
                frame.loc[row_index, "CD3e__state"] = ""
        frame.to_parquet(partition / "evidence.parquet", index=False)
    producer = tmp_path / "producer.py"
    producer.write_text("METHOD_VERSION = 'evaluation-test-v1'\n", encoding="utf-8")
    raw = tmp_path / "raw.txt"
    raw.write_text("immutable\n", encoding="utf-8")
    code = CodeMetadata.capture(tmp_path, source_paths=(producer,))
    ingest_stage = StageManifest.create(
        stage="ingest",
        method_version="evaluation-test-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=evidence,
        config={"test": "evaluation-review-context"},
        code=code,
    )
    ingest_path = ingest_stage.write_json(tmp_path / "ingest.manifest.json")
    stage = StageManifest.create(
        stage="calibrate",
        method_version="evaluation-test-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=evidence,
        config={"test": "evaluation-reporting"},
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "calibrate.manifest.json")
    source_run = RunManifest.create(
        run_name=splits.source_run_id,
        expected_donors=splits.donors,
        stage_manifest_paths=(ingest_path, stage_path),
        config={"test": "evaluation-reporting"},
        code=code,
    )
    source_run.write_json(tmp_path / "source-run.manifest.json")

    label_provenance = DevelopmentLabelProvenance(
        bundle_id="independent-labels-v1",
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
        fingerprint="a" * 64,
    )
    evidence_provenance = DevelopmentEvidenceProvenance(
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        source_run_manifest_content_id=source_run.content_id,
        source_stage_manifest_content_id=stage.content_id,
        source_artifact_content_id=stage.output.content_sha256,
        labeled_frame_fingerprint="b" * 64,
        feature_schema_version=FEATURE_SCHEMA_VERSION,
        row_count=6,
    )
    bundle = TypingModelBundle.create(
        bundle_version="heldout-candidate-v1",
        marker_registry_fingerprint=markers.fingerprint,
        typing_rules_fingerprint=rules.rules_fingerprint,
        split_manifest=splits,
        label_provenance=label_provenance,
        evidence_provenance=evidence_provenance,
        thresholds=AssignmentThresholds(
            inferred_probability=0.70,
            inferred_margin=0.20,
        ),
        broad=_ensemble(rules, level="broad"),
        subtypes=(_ensemble(rules, level="subtype", parent="Immune"),),
        fit_parameters=_fit_parameters(rules),
    )
    bundle = _freeze_bundle(bundle)
    references: dict[tuple[str, str, str], str] = {}
    for split_name in ("calibration", "validation", "locked_test"):
        for donor in getattr(splits, split_name):
            for label in bundle.broad.classes:
                references[(split_name, donor, f"b-{label}")] = label
            for label in bundle.subtypes[0].classes:
                references[(split_name, donor, f"s-{label}")] = label
            references[(split_name, donor, "x-Lymphatic_endothelial")] = "Endothelial"
            references[(split_name, donor, "x-CD34_endothelial")] = "Endothelial"
            references[(split_name, donor, "x-Beta")] = "Epithelial"
            references[(split_name, donor, "x-Alpha")] = "Epithelial"
            references[(split_name, donor, "x-Delta")] = "Epithelial"
            references[(split_name, donor, "x-SMA_positive_mesenchymal")] = "Mesenchymal"
    return EvaluationFixture(
        root=tmp_path,
        evidence=evidence,
        splits=splits,
        markers=markers,
        rules=rules,
        bundle=bundle,
        source_run=source_run,
        references=references,
    )


def _review_for_sample(
    sample: pd.DataFrame,
    reviewer_frame: pd.DataFrame,
    references: list[str],
    *,
    ensemble: FrozenModelEnsemble,
    splits: DonorSplitManifest,
    split_name: str,
    candidate_manifest: CandidateEvaluationManifest,
    bundle: TypingModelBundle,
    review_context_artifact_content_id: str,
) -> EvaluationTargetReview:
    sample = sample.reset_index(drop=True).copy()
    reviewer_frame = reviewer_frame.reset_index(drop=True).copy()
    reviewer_frame.attrs["review_context_artifact_content_id"] = (
        review_context_artifact_content_id
    )
    labels = sample.loc[:, ["donor_id", "object_id"]].copy()
    labels["reference_label"] = references
    labels["reviewer_1"] = "reviewer-a"
    labels["reviewer_2"] = "reviewer-b"
    labels["adjudicated"] = True
    labels["root_cause"] = "independent evidence review"
    labels["evidence_sufficient"] = True
    provenance = ReviewLabelProvenance.from_labels(
        labels,
        bundle_id=f"{split_name}-{ensemble.level}-{ensemble.parent or 'broad'}",
        splits=splits,
        split_name=split_name,
        sample_frame=sample,
        reviewer_frame=reviewer_frame,
        assignment_run_id=candidate_manifest.assignment_run_id,
        typing_model_bundle_content_id=bundle.content_id,
        candidate_evaluation_manifest_content_id=candidate_manifest.content_id,
        candidate_assignment_output_content_id=candidate_manifest.output_content_sha256,
        level=ensemble.level,
        parent=ensemble.parent,
    )
    return EvaluationTargetReview(
        sample_frame=sample,
        reviewer_frame=reviewer_frame,
        review_labels=labels,
        provenance=provenance,
    )


def _report_policy(n_bootstrap: int = 1000) -> PromotionGatePolicy:
    return PromotionGatePolicy(
        minimum_donors=5,
        minimum_rows=20,
        minimum_coverage_lower=1e-6,
        maximum_selective_error_upper=0.999999,
        minimum_rescue_coverage_lower=1e-6,
        maximum_rescue_selective_error_upper=0.999999,
        maximum_class_false_positive_rate_upper=0.999999,
        maximum_class_false_negative_rate_upper=0.999999,
        minimum_class_positive_rows=3,
        minimum_class_negative_rows=3,
        minimum_class_positive_donors=3,
        minimum_class_negative_donors=3,
        confidence_level=0.95,
        bootstrap_replicates=n_bootstrap,
    )


def _build(
    fixture: EvaluationFixture,
    manifest_path: Path,
    reviews: tuple[EvaluationTargetReview, ...],
    *,
    allow_locked_test: bool = False,
    n_bootstrap: int = 1000,
):
    return build_split_evaluation_report(
        candidate_evaluation=manifest_path,
        model_bundle=fixture.bundle,
        marker_registry=fixture.markers,
        typing_registry=fixture.rules,
        source_run_manifest=fixture.source_run,
        target_reviews=reviews,
        promotion_policy=_report_policy(n_bootstrap),
        review_blinding_key=_REVIEW_BLINDING_KEY,
        allow_locked_test=allow_locked_test,
    )


def _run_and_reviews(
    fixture: EvaluationFixture,
    name: str,
    *,
    split_name: str = "validation",
    transform=None,
):
    manifest_path, assignments = fixture.candidate(name, split_name=split_name, transform=transform)
    manifest = CandidateEvaluationManifest.read_json(manifest_path)
    reviews = fixture.reviews(
        assignments,
        split_name=split_name,
        candidate_manifest=manifest,
    )
    return manifest_path, assignments, reviews


def test_builds_weighted_report_from_exact_artifacts_and_is_deterministic(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "happy")
    report = _build(evaluation_fixture, manifest_path, tuple(reversed(reviews)))
    repeated = _build(evaluation_fixture, manifest_path, reviews)

    manifest = CandidateEvaluationManifest.read_json(manifest_path)
    assert report.content_id == repeated.content_id
    assert report.candidate_evaluation_manifest_content_id == manifest.content_id
    assert report.assignment_artifact_content_id == manifest.output_content_sha256
    assert report.evaluation_donors == evaluation_fixture.splits.validation
    assert [target.target_key for target in report.targets] == [
        ("broad", None),
        ("specific", None),
        ("subtype", "Immune"),
    ]
    broad = report.targets[0]
    expected_population = reviews[0].sample_frame["sampling_weight"].sum()
    assert broad.population_weight == pytest.approx(expected_population)
    assert broad.coverage.estimate == pytest.approx(1.0)
    assert broad.rescue_coverage.estimate == pytest.approx(1.0)
    assert sum(error.false_positive_weight for error in broad.class_errors) == pytest.approx(
        broad.wrong_call_weight
    )
    assert all(
        error.false_positive_rate.lower <= error.false_positive_rate.upper
        for error in broad.class_errors
    )
    policy_path = evaluation_fixture.root / "promotion-policy.json"
    report_path = evaluation_fixture.root / "validation-report.json"
    report.promotion_policy.write_json(policy_path)
    report.write_json(report_path)
    assert PromotionGatePolicy.read_json(policy_path).content_id == (
        report.promotion_policy.content_id
    )
    assert SplitEvaluationReport.read_json(report_path).content_id == report.content_id
    report.validate_candidate_source_replay(
        model_bundle=evaluation_fixture.bundle,
        marker_registry=evaluation_fixture.markers,
        typing_registry=evaluation_fixture.rules,
        review_blinding_key=_REVIEW_BLINDING_KEY,
    )


def test_preliminary_bundle_cannot_be_used_for_heldout_evaluation(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "preliminary-denied")
    preliminary = replace(
        evaluation_fixture.bundle,
        threshold_selection_provenance=None,
        model_fit_content_id="",
        content_id="",
    )
    with pytest.raises(ValueError, match="calibration-frozen"):
        build_split_evaluation_report(
            candidate_evaluation=manifest_path,
            model_bundle=preliminary,
            marker_registry=evaluation_fixture.markers,
            typing_registry=evaluation_fixture.rules,
            source_run_manifest=evaluation_fixture.source_run,
            target_reviews=reviews,
            promotion_policy=_report_policy(),
            review_blinding_key=_REVIEW_BLINDING_KEY,
        )


def test_locked_test_requires_explicit_authorization(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(
        evaluation_fixture, "locked", split_name="locked_test"
    )
    with pytest.raises(PermissionError, match="explicit"):
        _build(evaluation_fixture, manifest_path, reviews)
    report = _build(
        evaluation_fixture,
        manifest_path,
        reviews,
        allow_locked_test=True,
    )
    assert report.split_name == "locked_test"


def test_calibration_assignments_cannot_become_heldout_reports(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(
        evaluation_fixture, "calibration", split_name="calibration"
    )
    with pytest.raises(ValueError, match="reserved for threshold selection"):
        _build(evaluation_fixture, manifest_path, reviews)


def test_requires_exact_unique_model_target_reviews(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "targets")
    with pytest.raises(ValueError, match="incomplete"):
        _build(evaluation_fixture, manifest_path, reviews[:1])
    with pytest.raises(ValueError, match="duplicate"):
        _build(evaluation_fixture, manifest_path, (reviews[0], reviews[0], reviews[1]))
    extra_provenance = replace(
        reviews[2].provenance,
        parent="Epithelial",
        content_id="",
    )
    extra = replace(reviews[2], provenance=extra_provenance)
    with pytest.raises(ValueError, match="extra"):
        _build(evaluation_fixture, manifest_path, (*reviews, extra))


def test_rejects_donor_leakage_and_incomplete_review(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "governance")
    leaked_sample = reviews[0].sample_frame.copy()
    leaked_labels = reviews[0].review_labels.copy()
    leaked_reviewer = reviews[0].reviewer_frame.copy()
    leaked_sample.loc[0, "donor_id"] = "D1"
    leaked_labels.loc[0, "donor_id"] = "D1"
    leaked_reviewer.loc[0, "donor_id"] = "D1"
    leaked = replace(
        reviews[0],
        sample_frame=leaked_sample,
        reviewer_frame=leaked_reviewer,
        review_labels=leaked_labels,
    )
    with pytest.raises(ValueError, match="sample differs|escape|mis-stratified"):
        _build(evaluation_fixture, manifest_path, (leaked, *reviews[1:]))

    incomplete = replace(
        reviews[0],
        provenance=replace(
            reviews[0].provenance,
            complete_evidence_sufficient=False,
            content_id="",
        ),
    )
    with pytest.raises(ValueError, match="completeness|evidence-sufficient"):
        _build(evaluation_fixture, manifest_path, (incomplete, *reviews[1:]))


def test_rejects_missing_reviewed_cells_and_assignment_columns(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "missing-cell")
    broad = reviews[0]
    malformed_sample = broad.sample_frame.copy()
    malformed_sample.loc[0, "object_id"] = "not-assigned"
    manifest = CandidateEvaluationManifest.read_json(manifest_path)
    missing_cell_review = _review_for_sample(
        malformed_sample,
        blinded_review_sample_frame(malformed_sample),
        broad.review_labels["reference_label"].tolist(),
        ensemble=evaluation_fixture.bundle.broad,
        splits=evaluation_fixture.splits,
        split_name="validation",
        candidate_manifest=manifest,
        bundle=evaluation_fixture.bundle,
        review_context_artifact_content_id=(
            _ingest_output_content_id(evaluation_fixture.source_run)
        ),
    )
    with pytest.raises(ValueError, match="absent|mis-stratified"):
        _build(
            evaluation_fixture,
            manifest_path,
            (missing_cell_review, *reviews[1:]),
        )

    def drop_reason(frame):
        return frame.drop(columns="broad_reason")

    missing_column_path, _, missing_column_reviews = _run_and_reviews(
        evaluation_fixture,
        "missing-column",
        transform=drop_reason,
    )
    with pytest.raises(KeyError, match="missing .*columns"):
        _build(
            evaluation_fixture,
            missing_column_path,
            missing_column_reviews,
        )


def test_rejects_unsupported_labels_and_missing_class_support(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "labels")
    broad = reviews[0]
    unsupported_labels = broad.review_labels.copy()
    unsupported_labels.loc[0, "reference_label"] = "Invented_cell"
    unsupported_provenance = ReviewLabelProvenance.from_labels(
        unsupported_labels,
        bundle_id="unsupported",
        splits=evaluation_fixture.splits,
        split_name="validation",
        sample_frame=broad.sample_frame,
        reviewer_frame=broad.reviewer_frame,
        assignment_run_id=broad.provenance.assignment_run_id,
        typing_model_bundle_content_id=evaluation_fixture.bundle.content_id,
        candidate_evaluation_manifest_content_id=(
            broad.provenance.candidate_evaluation_manifest_content_id
        ),
        candidate_assignment_output_content_id=(
            broad.provenance.candidate_assignment_output_content_id
        ),
    )
    unsupported = replace(
        broad,
        review_labels=unsupported_labels,
        provenance=unsupported_provenance,
    )
    with pytest.raises(ValueError, match="unsupported labels"):
        _build(evaluation_fixture, manifest_path, (unsupported, *reviews[1:]))

    one_class_labels = broad.review_labels.copy()
    one_class_labels["reference_label"] = evaluation_fixture.bundle.broad.classes[0]
    one_class_provenance = ReviewLabelProvenance.from_labels(
        one_class_labels,
        bundle_id="one-class",
        splits=evaluation_fixture.splits,
        split_name="validation",
        sample_frame=broad.sample_frame,
        reviewer_frame=broad.reviewer_frame,
        assignment_run_id=broad.provenance.assignment_run_id,
        typing_model_bundle_content_id=evaluation_fixture.bundle.content_id,
        candidate_evaluation_manifest_content_id=(
            broad.provenance.candidate_evaluation_manifest_content_id
        ),
        candidate_assignment_output_content_id=(
            broad.provenance.candidate_assignment_output_content_id
        ),
    )
    one_class = replace(
        broad,
        review_labels=one_class_labels,
        provenance=one_class_provenance,
    )
    with pytest.raises(ValueError, match="no (positive|negative) support"):
        _build(evaluation_fixture, manifest_path, (one_class, *reviews[1:]))


def test_reference_other_is_retained_as_a_named_call_false_positive(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "other-negative")
    baseline = _build(evaluation_fixture, manifest_path, reviews)
    broad = reviews[0]
    labels = broad.review_labels.copy()
    first_class = evaluation_fixture.bundle.broad.classes[0]
    row = labels.index[labels["donor_id"].eq("V2") & labels["reference_label"].eq(first_class)][0]
    added_weight = float(broad.sample_frame.loc[row, "sampling_weight"])
    labels.loc[row, "reference_label"] = "Other"
    provenance = ReviewLabelProvenance.from_labels(
        labels,
        bundle_id="broad-with-other-negative",
        splits=evaluation_fixture.splits,
        split_name="validation",
        sample_frame=broad.sample_frame,
        reviewer_frame=broad.reviewer_frame,
        assignment_run_id=broad.provenance.assignment_run_id,
        typing_model_bundle_content_id=evaluation_fixture.bundle.content_id,
        candidate_evaluation_manifest_content_id=(
            broad.provenance.candidate_evaluation_manifest_content_id
        ),
        candidate_assignment_output_content_id=(
            broad.provenance.candidate_assignment_output_content_id
        ),
    )
    modified = replace(broad, review_labels=labels, provenance=provenance)
    report = _build(evaluation_fixture, manifest_path, (modified, *reviews[1:]))

    baseline_broad = baseline.targets[0]
    broad_result = report.targets[0]
    first_error = next(error for error in broad_result.class_errors if error.label == first_class)
    baseline_first_error = next(
        error for error in baseline_broad.class_errors if error.label == first_class
    )
    assert broad_result.wrong_call_weight == pytest.approx(
        baseline_broad.wrong_call_weight + added_weight
    )
    assert first_error.false_positive_weight == pytest.approx(
        baseline_first_error.false_positive_weight + added_weight
    )


def test_rejects_targets_without_actual_rescue_calls(
    evaluation_fixture: EvaluationFixture,
):
    def anchors_only(frame):
        broad = frame["object_id"].str.startswith("b-")
        frame.loc[broad, "broad_assignment_status"] = "anchor"
        frame.loc[broad, "broad_reason"] = "authoritative_anchor"
        frame.loc[broad, "broad_threshold_eligible"] = False
        return frame

    manifest_path, _, reviews = _run_and_reviews(
        evaluation_fixture, "no-rescue", transform=anchors_only
    )
    with pytest.raises(ValueError, match="deterministic replay"):
        _build(evaluation_fixture, manifest_path, reviews)


def test_zero_denominator_bootstrap_draws_are_retained_as_fail_closed_bounds(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, assignments, reviews = _run_and_reviews(evaluation_fixture, "sparse-support")
    first_class = evaluation_fixture.bundle.broad.classes[0]
    broad = reviews[0]
    donor = broad.review_labels["donor_id"]
    sparse_labels = broad.review_labels.copy()
    sparse_labels.loc[donor == "V1", "reference_label"] = first_class
    other_classes = evaluation_fixture.bundle.broad.classes[1:]
    v2_rows = sparse_labels.index[donor == "V2"]
    sparse_labels.loc[v2_rows, "reference_label"] = [
        other_classes[index % len(other_classes)] for index in range(len(v2_rows))
    ]
    sparse_provenance = ReviewLabelProvenance.from_labels(
        sparse_labels,
        bundle_id="sparse-broad-review",
        splits=evaluation_fixture.splits,
        split_name="validation",
        sample_frame=broad.sample_frame,
        reviewer_frame=broad.reviewer_frame,
        assignment_run_id=broad.provenance.assignment_run_id,
        typing_model_bundle_content_id=evaluation_fixture.bundle.content_id,
        candidate_evaluation_manifest_content_id=(
            broad.provenance.candidate_evaluation_manifest_content_id
        ),
        candidate_assignment_output_content_id=(
            broad.provenance.candidate_assignment_output_content_id
        ),
    )
    sparse = replace(broad, review_labels=sparse_labels, provenance=sparse_provenance)
    report = _build(
        evaluation_fixture,
        manifest_path,
        (sparse, *reviews[1:]),
        n_bootstrap=1000,
    )
    first_error = report.targets[0].class_errors[0]
    assert first_error.false_negative_rate.upper == 1.0
    assert report.n_bootstrap == 1000
    assert set(assignments["donor_id"]) == {"V1", "V2", "V3", "V4", "V5"}


def test_assignment_content_tampering_fails_before_metrics(
    evaluation_fixture: EvaluationFixture,
):
    manifest_path, _, reviews = _run_and_reviews(evaluation_fixture, "tamper")
    manifest = CandidateEvaluationManifest.read_json(manifest_path)
    assignment_path = Path(manifest.output.donors[0].files[0].path)
    assignments = pd.read_parquet(assignment_path)
    assignments.loc[0, "broad_label"] = "tampered"
    assignments.to_parquet(assignment_path, index=False)
    with pytest.raises(StaleArtifactError, match="changed"):
        _build(evaluation_fixture, manifest_path, reviews)


def _bootstrap_row(
    object_id: str,
    *,
    population: int,
    sampled: int,
) -> EvaluatedDecisionRow:
    sample = {
        "donor_id": "D1",
        "object_id": object_id,
        "sampling_stratum": "S1",
        "sampling_method": "hmac_sha256_candidate_bound_v2",
        "sampling_salt": "unit-test-salt",
        "stratum_population": population,
        "stratum_sampled": sampled,
        "sampling_weight": population / sampled,
        "selection_probability": sampled / population,
    }
    return EvaluatedDecisionRow(
        donor_id="D1",
        object_id=object_id,
        sample_record=sample,
        reviewer_record={
            "donor_id": "D1",
            "object_id": object_id,
            "sampling_stratum": "S1",
        },
        reference_label="A",
        reviewer_1="reviewer-a",
        reviewer_2="reviewer-b",
        adjudicated=True,
        root_cause="unit test",
        evidence_sufficient=True,
        prediction="A",
        assignment_status="inferred",
        assignment_reason="posterior_margin_stability_pass",
        target_scope=True,
        called=True,
        rescue_eligible=True,
        rescue_called=True,
    )


def test_two_stage_bootstrap_keeps_census_exact_and_adds_noncensus_variance():
    values = np.asarray([[1.0], [3.0]])
    census_rows = tuple(
        _bootstrap_row(f"census-{index}", population=2, sampled=2)
        for index in range(2)
    )
    noncensus_rows = tuple(
        _bootstrap_row(f"sample-{index}", population=4, sampled=2)
        for index in range(2)
    )
    census = _two_stage_survey_bootstrap(
        rows=census_rows,
        values=values,
        donors=("D1",),
        n_bootstrap=200,
        random_state=17,
        target_key=("broad", None),
    )
    noncensus = _two_stage_survey_bootstrap(
        rows=noncensus_rows,
        values=values,
        donors=("D1",),
        n_bootstrap=200,
        random_state=17,
        target_key=("broad", None),
    )
    assert np.all(census[:, 0] == 4.0)
    assert np.std(noncensus[:, 0]) > 0


def test_simultaneous_intervals_use_finite_b_order_and_boundary_guard():
    unresolved = _target_simultaneous_intervals(
        {"coverage": 0.8},
        {"coverage": np.full(100, 0.8)},
        {"coverage": np.ones(5)},
        donor_ids=np.asarray(["D1", "D2", "D3", "D4", "D5"]),
        confidence_level=0.95,
        multiplicity_target_count=3,
    )
    assert unresolved["coverage"].lower == 0.0
    assert unresolved["coverage"].upper == 1.0

    replicates = 4000
    boundary = _target_simultaneous_intervals(
        {"coverage": 1.0, "selective_error": 0.0},
        {
            "coverage": np.ones(replicates),
            "selective_error": np.zeros(replicates),
        },
        {
            "coverage": np.ones(5),
            "selective_error": np.ones(5),
        },
        donor_ids=np.asarray(["D1", "D2", "D3", "D4", "D5"]),
        confidence_level=0.95,
        multiplicity_target_count=1,
    )
    alpha_boundary = 0.05 / (2 * 1 * 2)
    expected_upper = 1.0 - alpha_boundary ** (1.0 / 5)
    assert boundary["selective_error"].upper == pytest.approx(expected_upper)
    assert boundary["coverage"].lower == pytest.approx(1.0 - expected_upper)
