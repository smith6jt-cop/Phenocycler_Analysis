from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, replace
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    InputArtifact,
    RunManifest,
    StageManifest,
    sha256_json,
    snapshot_parquet_dataset,
)
from phenocycler.candidate_evaluation import (
    CandidateEvaluationManifest,
    candidate_assignment_run_id,
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
    OPERATING_POINT_COLUMNS,
    ThresholdSelectionPolicy,
    ThresholdSelectionProvenance,
    ThresholdSelectionReplayInputs,
    ThresholdSelectionSourceReplayProof,
    build_threshold_selection_provenance,
    operating_grid_fingerprint,
    threshold_grid_content_id,
    validate_threshold_selection_source_replay,
)
from phenocycler.typing_model_bundle import FrozenModelEnsemble, TypingModelBundle
from phenocycler.typing_stage import RULE_SCORE_SEMANTICS, type_candidate_cells

Target = tuple[str, str]
REVIEW_BLINDING_KEY = "threshold-tests-secret-key-v1"


@dataclass(frozen=True)
class ThresholdFixture:
    bundle: TypingModelBundle
    manifest: CandidateEvaluationManifest
    samples: dict[Target, pd.DataFrame]
    ledgers: dict[Target, pd.DataFrame]
    reviews: tuple[ReviewLabelProvenance, ...]
    assignment_path: Path
    markers: MarkerRegistry
    rules: TypingRegistry
    source_run: RunManifest
    subtype_parent: str


def _splits() -> DonorSplitManifest:
    return DonorSplitManifest(
        split_version="threshold-split-v2",
        source_run_id="source-run",
        cohort_content_id="cohort-content",
        development=("D0", "D1"),
        calibration=("D2",),
        validation=("D3",),
        review_blinding_key_sha256=hashlib.sha256(
            REVIEW_BLINDING_KEY.encode()
        ).hexdigest(),
        threshold_selection_policy_content_id=_policy().content_id,
        threshold_grid_content_id=threshold_grid_content_id(
            probability_thresholds=(0.50, 0.80),
            margin_thresholds=(0.05,),
            stability_thresholds=(0.80, 0.90),
        ),
        locked_test=("D4",),
    )


def _fit_parameters(
    rules: TypingRegistry, *, modeled_parents: tuple[str, ...] = ()
) -> dict[str, object]:
    structural = sorted(
        {
            str(rule.parent)
            for rule in rules.subtype_rules
            if len(rules.rules_for(level="subtype", parent=str(rule.parent))) == 1
        }
    )
    broad_classes = classifier_from_registry(rules, level="broad").classes
    final_specific_classes = sorted(
        {
            leaf
            for broad_class in broad_classes
            for leaf in (
                tuple(
                    rule.name
                    for rule in rules.rules_for(level="subtype", parent=broad_class)
                )
                or (broad_class,)
            )
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
        "missing_probabilistic_subtype_parents": sorted(
            {str(rule.parent) for rule in rules.subtype_rules}
            - set(structural)
            - set(modeled_parents)
        ),
        "final_specific_classes": final_specific_classes,
        "minimum_valid_bootstrap_fraction": 0.8,
        "class_priors": "uniform_no_intercept",
        "donor_weighting": "class_then_donor_balanced",
    }


def _ensemble(model, *, splits: DonorSplitManifest) -> FrozenModelEnsemble:
    coefficients = tuple(
        tuple(float(value) for value in row) for row in np.asarray(model.coefficients, dtype=float)
    )
    return FrozenModelEnsemble(
        level=model.level,
        parent=model.parent,
        classes=model.classes,
        markers=model.markers,
        coefficients=coefficients,
        bootstrap_coefficients=(coefficients,),
        training_donors=splits.development,
        calibration_donors=splits.calibration,
        bootstrap_draws=(splits.development,),
        bootstrap_valid=(True,),
        bootstrap_failures=("",),
        temperature=1.0,
        n_development_labels=20,
        n_calibration_labels=10,
    )


def _ranking(classes: tuple[str, ...], winner: str, probability: float) -> str:
    rest = (1.0 - probability) / (len(classes) - 1)
    values = [(name, probability if name == winner else rest) for name in classes]
    values.sort(key=lambda item: (-item[1], classes.index(item[0])))
    return json.dumps(values, separators=(",", ":"))


def _probability_vector(
    classes: tuple[str, ...], winner: str, probability: float
) -> dict[str, float]:
    rest = (1.0 - probability) / (len(classes) - 1)
    return {name: probability if name == winner else rest for name in classes}


def _evidence_frame(donor: str, object_ids: list[str], *, markers: MarkerRegistry) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "donor_id": [donor] * len(object_ids),
            "object_id": object_ids,
            "qc_analysis_eligible": [True] * len(object_ids),
        }
    )
    for marker in markers.active_markers:
        frame[f"{marker.name}__state"] = "negative"
        frame[f"{marker.name}__log2_threshold_ratio"] = -1.0
        frame[f"{marker.name}__model_valid"] = True
    # The default registry has an Immune/CD3e anchor.  It gives the template
    # output a real probabilistic Immune subtype evaluation on every row.
    frame["CD3e__state"] = "positive"
    frame["CD3e__log2_threshold_ratio"] = 3.0
    return frame


def _ledger(sample: pd.DataFrame, references: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": sample["donor_id"].to_numpy(),
            "object_id": sample["object_id"].to_numpy(),
            "reference_label": references,
            "reviewer_1": ["reviewer-A"] * len(sample),
            "reviewer_2": ["reviewer-B"] * len(sample),
            "adjudicated": [True] * len(sample),
            "root_cause": ["independent_blinded_review"] * len(sample),
            "evidence_sufficient": [True] * len(sample),
        }
    )


def _reviewer_frame(sample: pd.DataFrame, *, context_content_id: str) -> pd.DataFrame:
    frame = blinded_review_sample_frame(sample)
    frame.attrs["review_context_artifact_content_id"] = context_content_id
    return frame


def _review_provenances(
    *,
    bundle: TypingModelBundle,
    manifest: CandidateEvaluationManifest,
    samples: dict[Target, pd.DataFrame],
    ledgers: dict[Target, pd.DataFrame],
) -> tuple[ReviewLabelProvenance, ...]:
    return tuple(
        ReviewLabelProvenance.from_labels(
            ledgers[target],
            bundle_id=f"calibration-{target[0]}-{target[1] or 'none'}-v1",
            splits=bundle.split_manifest,
            split_name="calibration",
            sample_frame=samples[target],
            reviewer_frame=_reviewer_frame(
                samples[target],
                context_content_id=bundle.evidence_provenance.source_artifact_content_id,
            ),
            assignment_run_id=manifest.assignment_run_id,
            typing_model_bundle_content_id=bundle.content_id,
            candidate_evaluation_manifest_content_id=manifest.content_id,
            candidate_assignment_output_content_id=manifest.output_content_sha256,
            level=target[0],
            parent=target[1],
        )
        for target in samples
    )


def _stamp_level(
    assignments: pd.DataFrame,
    index: int,
    *,
    prefix: str,
    classes: tuple[str, ...],
    status: str,
    winner: str,
    probability: float,
) -> None:
    if status in {"anchor", "inferred"}:
        label = winner
        reason = (
            "single_authoritative_anchor"
            if status == "anchor"
            else "posterior_margin_stability_pass"
        )
    elif status == "ambiguous":
        label = "Ambiguous" if prefix == "broad" else "Immune-unclassified"
        reason = "posterior_margin_or_stability_failed"
    elif status == "other":
        label = "Other"
        reason = "all_required_models_valid_and_negative"
    else:
        label = "Unavailable"
        reason = "required_models_unavailable"
    vector = _probability_vector(classes, winner, probability)
    second = max(value for name, value in vector.items() if name != winner)
    assignments.loc[index, f"{prefix}_label"] = label
    assignments.loc[index, f"{prefix}_assignment_status"] = status
    assignments.loc[index, f"{prefix}_reason"] = reason
    assignments.loc[index, f"{prefix}_best_candidate"] = winner
    assignments.loc[index, f"{prefix}_best_probability"] = probability
    assignments.loc[index, f"{prefix}_margin"] = probability - second
    assignments.loc[index, f"{prefix}_stability"] = 1.0 if status != "other" else np.nan
    assignments.loc[index, f"{prefix}_threshold_eligible"] = status in {
        "inferred",
        "ambiguous",
    }
    assignments.loc[index, f"{prefix}_valid_evidence_pass"] = status not in {
        "other",
        "unavailable",
    }
    assignments.loc[index, f"ranked_{prefix}_probabilities"] = _ranking(
        classes, winner, probability
    )
    if prefix == "broad":
        for class_name, value in vector.items():
            assignments.loc[index, f"broad_probability__{class_name}"] = value


def _fixture(tmp_path: Path) -> ThresholdFixture:
    splits = _splits()
    markers = load_registry()
    rules = TypingRegistry.from_marker_registry(markers)
    subtype_parent = "Immune"
    broad_model = classifier_from_registry(rules, level="broad")
    subtype_model = classifier_from_registry(rules, level="subtype", parent=subtype_parent)
    broad_classes = broad_model.classes
    subtype_classes = subtype_model.classes

    subtype_anchor_keys = [f"subtype-{index}" for index in range(len(subtype_classes))]
    broad_anchor_keys = [
        f"broad-{index}"
        for index, class_name in enumerate(broad_classes)
        if class_name != subtype_parent
    ]
    special_keys = [
        "broad-rescue-correct",
        "broad-rescue-wrong",
        "subtype-rescue-correct",
        "subtype-rescue-wrong",
        "other",
    ]
    calibration_keys = [*subtype_anchor_keys, *broad_anchor_keys, *special_keys]

    evidence_root = tmp_path / "marker-evidence"
    for donor in splits.donors:
        keys = calibration_keys if donor == "D2" else [f"{donor}-cell"]
        partition = evidence_root / f"donor_id={donor}"
        partition.mkdir(parents=True)
        _evidence_frame(donor, keys, markers=markers).to_parquet(
            partition / "evidence.parquet", index=False
        )
    producer = tmp_path / "calibration_producer.py"
    producer.write_text("METHOD_VERSION = 'calibrate-threshold-test-v1'\n", encoding="utf-8")
    raw = tmp_path / "raw-input.txt"
    raw.write_text("immutable input\n", encoding="utf-8")
    code = CodeMetadata.capture(tmp_path, source_paths=(producer,))
    ingest_stage = StageManifest.create(
        stage="ingest",
        method_version="ingest-threshold-test-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=evidence_root,
        config={"test": "threshold-selection-review-context"},
        code=code,
    )
    ingest_stage_path = ingest_stage.write_json(tmp_path / "ingest.manifest.json")
    stage = StageManifest.create(
        stage="calibrate",
        method_version="calibrate-threshold-test-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=evidence_root,
        config={"test": "threshold-selection"},
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "calibrate.manifest.json")
    source_run = RunManifest.create(
        run_name=splits.source_run_id,
        expected_donors=splits.donors,
        stage_manifest_paths=(ingest_stage_path, stage_path),
        config={"test": "threshold-selection"},
        code=code,
    )

    bundle = TypingModelBundle.create(
        bundle_version="unpromoted-calibration-candidate-v1",
        marker_registry_fingerprint=markers.fingerprint,
        typing_rules_fingerprint=rules.rules_fingerprint,
        split_manifest=splits,
        label_provenance=DevelopmentLabelProvenance(
            bundle_id="development-labels-v1",
            source_run_id=splits.source_run_id,
            split_manifest_content_id=splits.content_id,
            purpose=DEVELOPMENT_LABEL_PURPOSE,
            fingerprint="4" * 64,
        ),
        evidence_provenance=DevelopmentEvidenceProvenance(
            source_run_id=splits.source_run_id,
            split_manifest_content_id=splits.content_id,
            source_run_manifest_content_id=source_run.content_id,
            source_stage_manifest_content_id=stage.content_id,
            source_artifact_content_id=stage.output.content_sha256,
            labeled_frame_fingerprint="5" * 64,
            feature_schema_version=FEATURE_SCHEMA_VERSION,
            row_count=20,
        ),
        thresholds=AssignmentThresholds(
            inferred_probability=0.50,
            inferred_margin=0.05,
            inferred_stability=0.80,
        ),
        broad=_ensemble(broad_model, splits=splits),
        subtypes=(_ensemble(subtype_model, splits=splits),),
        fit_parameters=_fit_parameters(rules, modeled_parents=(subtype_parent,)),
    )
    run_id = candidate_assignment_run_id(
        model_bundle_content_id=bundle.content_id,
        split_manifest_content_id=splits.content_id,
        split_name="calibration",
        source_run_id=splits.source_run_id,
        source_run_manifest_content_id=source_run.content_id,
        source_stage_manifest_content_id=stage.content_id,
        source_artifact_content_id=stage.output.content_sha256,
        marker_registry_fingerprint=bundle.marker_registry_fingerprint,
        typing_rules_fingerprint=bundle.typing_rules_fingerprint,
        donors=splits.calibration,
    )
    calibration_evidence = pd.read_parquet(evidence_root / "donor_id=D2" / "evidence.parquet")
    assignments = type_candidate_cells(
        calibration_evidence,
        model_bundle=bundle,
        marker_registry=markers,
        typing_registry=rules,
    )
    assignments["candidate_assignment_run_id"] = run_id
    assignments["candidate_evaluation_split_name"] = "calibration"
    row_index = {value: index for index, value in enumerate(assignments["object_id"])}

    broad_references: dict[str, str] = {}
    subtype_references: dict[str, str] = {}
    for index, class_name in enumerate(subtype_classes):
        key = subtype_anchor_keys[index]
        row = row_index[key]
        _stamp_level(
            assignments,
            row,
            prefix="broad",
            classes=broad_classes,
            status="anchor",
            winner=subtype_parent,
            probability=0.95,
        )
        _stamp_level(
            assignments,
            row,
            prefix="subtype",
            classes=subtype_classes,
            status="anchor",
            winner=class_name,
            probability=0.95,
        )
        broad_references[key] = subtype_parent
        subtype_references[key] = class_name
    for key, class_name in zip(
        broad_anchor_keys,
        (name for name in broad_classes if name != subtype_parent),
        strict=True,
    ):
        row = row_index[key]
        _stamp_level(
            assignments,
            row,
            prefix="broad",
            classes=broad_classes,
            status="anchor",
            winner=class_name,
            probability=0.95,
        )
        broad_references[key] = class_name

    broad_other = next(name for name in broad_classes if name != subtype_parent)
    for key, winner, probability, reference in (
        ("broad-rescue-correct", subtype_parent, 0.90, subtype_parent),
        ("broad-rescue-wrong", broad_other, 0.55, subtype_parent),
    ):
        _stamp_level(
            assignments,
            row_index[key],
            prefix="broad",
            classes=broad_classes,
            status="ambiguous",
            winner=winner,
            probability=probability,
        )
        broad_references[key] = reference
    subtype_first, subtype_second = subtype_classes[:2]
    for key, winner, probability, reference in (
        ("subtype-rescue-correct", subtype_first, 0.90, subtype_first),
        ("subtype-rescue-wrong", subtype_first, 0.55, subtype_second),
    ):
        row = row_index[key]
        _stamp_level(
            assignments,
            row,
            prefix="broad",
            classes=broad_classes,
            status="anchor",
            winner=subtype_parent,
            probability=0.95,
        )
        _stamp_level(
            assignments,
            row,
            prefix="subtype",
            classes=subtype_classes,
            status="ambiguous",
            winner=winner,
            probability=probability,
        )
        broad_references[key] = subtype_parent
        subtype_references[key] = reference
    _stamp_level(
        assignments,
        row_index["other"],
        prefix="broad",
        classes=broad_classes,
        status="other",
        winner=subtype_parent,
        probability=0.30,
    )
    broad_references["other"] = "Other"

    accepted = assignments["broad_assignment_status"].isin({"anchor", "inferred"})
    configured_parents = {str(rule.parent) for rule in rules.subtype_rules}
    applicable = accepted & assignments["broad_label"].isin(configured_parents)
    bundled = applicable & assignments["broad_label"].eq(subtype_parent)
    fallback = applicable & ~bundled
    assignments["typing_subtype_model_source"] = "not_applicable"
    assignments.loc[bundled, "typing_subtype_model_source"] = "probabilistic_bundle"
    assignments.loc[fallback, "typing_subtype_model_source"] = "rules_fallback"
    assignments["typing_subtype_score_semantics"] = ""
    assignments.loc[bundled, "typing_subtype_score_semantics"] = bundle.score_semantics
    assignments.loc[fallback, "typing_subtype_score_semantics"] = RULE_SCORE_SEMANTICS
    assignments["typing_subtype_bootstrap_replicates"] = np.where(
        bundled, bundle.subtypes[0].bootstrap_replicates, 0
    ).astype(int)
    assignments["typing_subtype_valid_bootstrap_replicates"] = np.where(
        bundled, bundle.subtypes[0].valid_bootstrap_replicates, 0
    ).astype(int)

    output_root = tmp_path / "calibration-assignments"
    partition = output_root / "donor_id=D2"
    partition.mkdir(parents=True)
    assignment_path = partition / "assignments.parquet"
    assignments.to_parquet(assignment_path, index=False)
    snapshot = snapshot_parquet_dataset(output_root, expected_donors=splits.calibration)
    manifest = CandidateEvaluationManifest(
        assignment_run_id=run_id,
        model_bundle_content_id=bundle.content_id,
        split_manifest_content_id=splits.content_id,
        split_name="calibration",
        source_run_id=splits.source_run_id,
        source_run_manifest_content_id=source_run.content_id,
        source_stage_manifest_content_id=stage.content_id,
        source_artifact_content_id=stage.output.content_sha256,
        marker_registry_fingerprint=bundle.marker_registry_fingerprint,
        typing_rules_fingerprint=bundle.typing_rules_fingerprint,
        donors=splits.calibration,
        output=snapshot,
        output_content_sha256=snapshot.content_sha256,
    )
    broad_keys = list(assignments["object_id"].astype(str))
    subtype_keys = [*subtype_anchor_keys, "subtype-rescue-correct", "subtype-rescue-wrong"]
    broad_sample = build_review_sampling_frame(
        assignments,
        splits=splits,
        split_name="calibration",
        candidate_assignment_run_id=manifest.assignment_run_id,
        level="broad",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )
    subtype_sample = build_review_sampling_frame(
        assignments,
        splits=splits,
        split_name="calibration",
        candidate_assignment_run_id=manifest.assignment_run_id,
        level="subtype",
        parent=subtype_parent,
        blinding_key=REVIEW_BLINDING_KEY,
    )
    assert set(broad_sample["object_id"]) == set(broad_keys)
    assert set(subtype_sample["object_id"]) == set(subtype_keys)
    samples = {
        ("broad", ""): broad_sample,
        ("subtype", subtype_parent): subtype_sample,
    }
    ledgers = {
        ("broad", ""): _ledger(
            samples[("broad", "")],
            [
                broad_references[key]
                for key in samples[("broad", "")]["object_id"].astype(str)
            ],
        ),
        ("subtype", subtype_parent): _ledger(
            samples[("subtype", subtype_parent)],
            [
                subtype_references[key]
                for key in samples[("subtype", subtype_parent)]["object_id"].astype(str)
            ],
        ),
    }
    reviews = _review_provenances(
        bundle=bundle, manifest=manifest, samples=samples, ledgers=ledgers
    )
    return ThresholdFixture(
        bundle=bundle,
        manifest=manifest,
        samples=samples,
        ledgers=ledgers,
        reviews=reviews,
        assignment_path=assignment_path,
        markers=markers,
        rules=rules,
        source_run=source_run,
        subtype_parent=subtype_parent,
    )


def _policy(**changes) -> ThresholdSelectionPolicy:
    values = {
        "policy_version": "calibration-policy-v2",
        "false_positive_cost": 2.0,
        "false_negative_cost": 1.0,
        "maximum_selective_error": 0.0,
        "minimum_coverage": 0.5,
        "maximum_rescue_selective_error": 0.0,
        "minimum_rescue_coverage": 0.5,
    }
    values.update(changes)
    return ThresholdSelectionPolicy(**values)


def _build(fixture: ThresholdFixture, **changes) -> ThresholdSelectionProvenance:
    selected_samples = changes.get("calibration_sample_frames", fixture.samples)
    arguments = {
        "selection_version": "threshold-selection-v2",
        "model_bundle": fixture.bundle,
        "calibration_candidate_manifest": fixture.manifest,
        "marker_registry": fixture.markers,
        "typing_registry": fixture.rules,
        "source_run_manifest": fixture.source_run,
        "calibration_sample_frames": selected_samples,
        "calibration_reviewer_frames": {
            target: _reviewer_frame(
                sample,
                context_content_id=fixture.bundle.evidence_provenance.source_artifact_content_id,
            )
            for target, sample in selected_samples.items()
        },
        "calibration_review_ledgers": fixture.ledgers,
        "calibration_review_provenance": fixture.reviews,
        "review_blinding_key": REVIEW_BLINDING_KEY,
        "policy": _policy(),
        "probability_thresholds": (0.50, 0.80),
        "margin_thresholds": (0.05,),
        "stability_thresholds": (0.80, 0.90),
    }
    arguments.update(changes)
    # This fixture deliberately fabricates reviewed operating points that its
    # trivial classifier cannot emit.  Deterministic candidate replay has its
    # own adversarial coverage in test_candidate_evaluation.py; keep these
    # tests focused on the threshold-selection contract.
    with patch(
        "phenocycler.candidate_evaluation._validate_deterministic_assignment_replay"
    ):
        return build_threshold_selection_provenance(**arguments)


def _replay_inputs(
    fixture: ThresholdFixture,
    **changes,
) -> ThresholdSelectionReplayInputs:
    selected_samples = changes.get("calibration_sample_frames", fixture.samples)
    arguments = {
        "preliminary_bundle": fixture.bundle,
        "calibration_candidate_manifest": fixture.manifest,
        "marker_registry": fixture.markers,
        "typing_registry": fixture.rules,
        "source_run_manifest": fixture.source_run,
        "calibration_sample_frames": selected_samples,
        "calibration_reviewer_frames": {
            target: _reviewer_frame(
                sample,
                context_content_id=(
                    fixture.bundle.evidence_provenance.source_artifact_content_id
                ),
            )
            for target, sample in selected_samples.items()
        },
        "calibration_review_ledgers": fixture.ledgers,
        "review_blinding_key": REVIEW_BLINDING_KEY,
    }
    arguments.update(changes)
    return ThresholdSelectionReplayInputs(**arguments)


def _source_replay(
    frozen: TypingModelBundle,
    replay_inputs: ThresholdSelectionReplayInputs,
) -> ThresholdSelectionSourceReplayProof:
    # The fixture stamps reviewed operating points that its trivial classifier
    # cannot emit; deterministic candidate replay is independently covered by
    # candidate-evaluation tests.  All source/snapshot/lineage checks remain on.
    with patch(
        "phenocycler.candidate_evaluation._validate_deterministic_assignment_replay"
    ):
        return validate_threshold_selection_source_replay(frozen, replay_inputs)


def _readdress(payload: dict) -> dict:
    payload = json.loads(json.dumps(payload))
    payload["content_id"] = sha256_json(
        {key: value for key, value in payload.items() if key != "content_id"}
    )
    return payload


def _manifest_for_output(
    fixture: ThresholdFixture, assignments: pd.DataFrame, root: Path
) -> tuple[CandidateEvaluationManifest, Path]:
    partition = root / "donor_id=D2"
    partition.mkdir(parents=True)
    path = partition / "assignments.parquet"
    assignments.to_parquet(path, index=False)
    snapshot = snapshot_parquet_dataset(root, expected_donors=("D2",))
    manifest = replace(
        fixture.manifest,
        output=snapshot,
        output_content_sha256=snapshot.content_sha256,
        content_id="",
    )
    return manifest, path


def test_builder_derives_target_gated_grid_and_roundtrips(tmp_path):
    fixture = _fixture(tmp_path)
    result = _build(fixture)

    assert result.preliminary_bundle_content_id == fixture.bundle.content_id
    assert result.model_fit_content_id == fixture.bundle.model_fit_content_id
    assert result.calibration_candidate_manifest_content_id == fixture.manifest.content_id
    assert result.calibration_candidate_output_content_id == fixture.manifest.output_content_sha256
    assert result.assignment_run_id == fixture.manifest.assignment_run_id
    assert result.calibration_donors == ("D2",)
    assert result.evaluated_targets == (
        ("broad", ""),
        ("subtype", fixture.subtype_parent),
    )
    assert result.operating_grid_row_count == 4
    assert result.selected_probability_threshold == 0.8
    assert result.selected_metrics["selective_error"] == pytest.approx(0.0)

    frozen = fixture.bundle.freeze_threshold_selection(
        result, bundle_version="frozen-after-calibration-v1"
    )
    assert frozen.model_fit_content_id == fixture.bundle.model_fit_content_id
    assert frozen.content_id != fixture.bundle.content_id
    assert frozen.thresholds.inferred_probability == result.selected_probability_threshold

    path = tmp_path / "threshold-selection.json"
    result.write_json(path)
    assert ThresholdSelectionProvenance.read_json(path) == result
    with pytest.raises(FileExistsError, match="immutable"):
        result.write_json(path)


def test_frozen_threshold_selection_replays_from_exact_protected_sources(tmp_path):
    fixture = _fixture(tmp_path)
    provenance = _build(fixture)
    frozen = fixture.bundle.freeze_threshold_selection(
        provenance, bundle_version="source-replayed-frozen-v1"
    )

    proof = _source_replay(frozen, _replay_inputs(fixture))

    assert proof.threshold_selection_provenance_content_id == provenance.content_id
    assert proof.preliminary_bundle_content_id == fixture.bundle.content_id
    assert proof.model_fit_content_id == fixture.bundle.model_fit_content_id
    assert proof.calibration_candidate_manifest_content_id == fixture.manifest.content_id
    assert (
        proof.calibration_candidate_output_content_id
        == fixture.manifest.output_content_sha256
    )
    assert proof.assignment_run_id == fixture.manifest.assignment_run_id
    assert (
        proof.source_run_manifest_content_id
        == fixture.bundle.evidence_provenance.source_run_manifest_content_id
    )
    assert (
        proof.source_stage_manifest_content_id
        == fixture.bundle.evidence_provenance.source_stage_manifest_content_id
    )
    assert (
        proof.source_artifact_content_id
        == fixture.bundle.evidence_provenance.source_artifact_content_id
    )
    assert not hasattr(proof, "to_dict")


def test_threshold_source_replay_rejects_candidate_output_tampering(tmp_path):
    fixture = _fixture(tmp_path)
    frozen = fixture.bundle.freeze_threshold_selection(_build(fixture))
    assignments = pd.read_parquet(fixture.assignment_path)
    assignments.loc[0, "broad_best_candidate"] = "tampered-after-freeze"
    assignments.to_parquet(fixture.assignment_path, index=False)

    with pytest.raises(Exception, match="changed|differs|snapshot"):
        _source_replay(frozen, _replay_inputs(fixture))


def test_threshold_source_replay_rejects_source_dataset_tampering(tmp_path):
    fixture = _fixture(tmp_path)
    frozen = fixture.bundle.freeze_threshold_selection(_build(fixture))
    source_path = tmp_path / "marker-evidence" / "donor_id=D2" / "evidence.parquet"
    evidence = pd.read_parquet(source_path)
    evidence.loc[0, "CD3e__log2_threshold_ratio"] += 0.125
    evidence.to_parquet(source_path, index=False)

    with pytest.raises(Exception, match="changed|differs|snapshot|hash"):
        _source_replay(frozen, _replay_inputs(fixture))


def test_threshold_source_replay_rejects_sampling_audit_tampering(tmp_path):
    fixture = _fixture(tmp_path)
    frozen = fixture.bundle.freeze_threshold_selection(_build(fixture))
    samples = {target: frame.copy() for target, frame in fixture.samples.items()}
    samples[("broad", "")]["sampling_salt"] = "f" * 64

    with pytest.raises(ValueError, match="salt differs|sample differs"):
        _source_replay(
            frozen,
            _replay_inputs(fixture, calibration_sample_frames=samples),
        )


def test_threshold_source_replay_proof_constructor_is_guarded(tmp_path):
    fixture = _fixture(tmp_path)
    provenance = _build(fixture)
    digest = "a" * 64
    with pytest.raises(TypeError, match="created only by successful"):
        ThresholdSelectionSourceReplayProof(
            threshold_selection_provenance_content_id=provenance.content_id,
            preliminary_bundle_content_id=fixture.bundle.content_id,
            model_fit_content_id=fixture.bundle.model_fit_content_id,
            calibration_candidate_manifest_content_id=fixture.manifest.content_id,
            calibration_candidate_output_content_id=(
                fixture.manifest.output_content_sha256
            ),
            assignment_run_id=fixture.manifest.assignment_run_id,
            source_run_manifest_content_id=digest,
            source_stage_manifest_content_id=digest,
            source_artifact_content_id=digest,
            review_context_artifact_content_id=digest,
            _verification_token=object(),
        )


def test_public_builder_has_no_caller_authored_grid_or_metrics_api(tmp_path):
    fixture = _fixture(tmp_path)
    with pytest.raises(TypeError, match="unexpected keyword argument"):
        _build(
            fixture,
            operating_points=pd.DataFrame(),
            selected_metrics={"coverage": 1.0},
        )


def test_grid_cannot_admit_broad_parents_absent_from_preliminary_subtype_scoring(
    tmp_path,
):
    fixture = _fixture(tmp_path)
    with pytest.raises(ValueError, match="grid differs from donor-split precommitment"):
        _build(fixture, probability_thresholds=(0.49, 0.80))


def test_policy_must_match_donor_split_precommitment(tmp_path):
    fixture = _fixture(tmp_path)
    with pytest.raises(ValueError, match="policy differs from donor-split precommitment"):
        _build(fixture, policy=_policy(false_positive_cost=3.0))


def test_resigned_sampling_population_cannot_reweight_calibration_errors(tmp_path):
    fixture = _fixture(tmp_path)
    samples = {key: value.copy() for key, value in fixture.samples.items()}
    broad = samples[("broad", "")]
    stratum = broad.loc[0, "sampling_stratum"]
    selected = broad["sampling_stratum"].eq(stratum)
    declared_sampled = int(selected.sum())
    broad.loc[selected, "stratum_population"] = declared_sampled + 100
    broad.loc[selected, "sampling_weight"] = (declared_sampled + 100) / declared_sampled
    broad.loc[selected, "selection_probability"] = declared_sampled / (
        declared_sampled + 100
    )
    reviews = _review_provenances(
        bundle=fixture.bundle,
        manifest=fixture.manifest,
        samples=samples,
        ledgers=fixture.ledgers,
    )
    with pytest.raises(ValueError, match="populations differ"):
        _build(
            fixture,
            calibration_sample_frames=samples,
            calibration_review_provenance=reviews,
        )


def test_exact_target_set_rejects_missing_and_unrelated_reviews(tmp_path):
    fixture = _fixture(tmp_path)
    with pytest.raises(ValueError, match="every probabilistic subtype"):
        _build(
            fixture,
            calibration_sample_frames={("broad", ""): fixture.samples[("broad", "")]},
        )
    unrelated = dict(fixture.ledgers)
    unrelated[("subtype", "Unrelated")] = fixture.ledgers[("subtype", fixture.subtype_parent)]
    with pytest.raises(ValueError, match="every probabilistic subtype"):
        _build(fixture, calibration_review_ledgers=unrelated)


def test_unrelated_cell_key_fails_after_exact_ledger_validation(tmp_path):
    fixture = _fixture(tmp_path)
    samples = {key: value.copy() for key, value in fixture.samples.items()}
    ledgers = {key: value.copy() for key, value in fixture.ledgers.items()}
    samples[("broad", "")].loc[0, "object_id"] = "unrelated"
    ledgers[("broad", "")].loc[0, "object_id"] = "unrelated"
    reviews = _review_provenances(
        bundle=fixture.bundle,
        manifest=fixture.manifest,
        samples=samples,
        ledgers=ledgers,
    )
    with pytest.raises(ValueError, match="absent from or mis-stratified"):
        _build(
            fixture,
            calibration_sample_frames=samples,
            calibration_review_ledgers=ledgers,
            calibration_review_provenance=reviews,
        )


def test_review_binds_preliminary_full_bundle_not_fit_identity(tmp_path):
    fixture = _fixture(tmp_path)
    bad_reviews = tuple(
        replace(
            review,
            typing_model_bundle_content_id=fixture.bundle.model_fit_content_id,
            content_id="",
        )
        for review in fixture.reviews
    )
    with pytest.raises(ValueError, match="different typing model bundle|preliminary full bundle"):
        _build(fixture, calibration_review_provenance=bad_reviews)


def test_assignment_tampering_is_detected_before_metrics(tmp_path):
    fixture = _fixture(tmp_path)
    assignments = pd.read_parquet(fixture.assignment_path)
    assignments.loc[0, "broad_best_candidate"] = "tampered"
    assignments.to_parquet(fixture.assignment_path, index=False)
    with pytest.raises(Exception, match="changed|differs|snapshot"):
        _build(fixture)


def test_resigned_subset_manifest_fails_exact_source_universe_check(tmp_path):
    fixture = _fixture(tmp_path)
    assignments = pd.read_parquet(fixture.assignment_path).iloc[:-1].copy()
    subset_manifest, _ = _manifest_for_output(fixture, assignments, tmp_path / "resigned-subset")
    with pytest.raises(ValueError, match="different row count|object universe"):
        _build(fixture, calibration_candidate_manifest=subset_manifest)


def test_targetwise_gate_rejects_subtype_with_zero_rescue_despite_aggregate(tmp_path):
    fixture = _fixture(tmp_path)
    assignments = pd.read_parquet(fixture.assignment_path)
    subtype_rows = assignments["typing_subtype_model_source"].eq("probabilistic_bundle")
    assignments.loc[subtype_rows, "subtype_assignment_status"] = "anchor"
    assignments.loc[subtype_rows, "subtype_label"] = assignments.loc[
        subtype_rows, "subtype_best_candidate"
    ]
    assignments.loc[subtype_rows, "subtype_reason"] = "single_authoritative_anchor"
    assignments.loc[subtype_rows, "subtype_threshold_eligible"] = False
    manifest, _ = _manifest_for_output(fixture, assignments, tmp_path / "zero-subtype-rescue")
    samples = {
        ("broad", ""): build_review_sampling_frame(
            assignments,
            splits=fixture.bundle.split_manifest,
            split_name="calibration",
            candidate_assignment_run_id=manifest.assignment_run_id,
            level="broad",
            parent=None,
            blinding_key=REVIEW_BLINDING_KEY,
        ),
        ("subtype", fixture.subtype_parent): build_review_sampling_frame(
            assignments,
            splits=fixture.bundle.split_manifest,
            split_name="calibration",
            candidate_assignment_run_id=manifest.assignment_run_id,
            level="subtype",
            parent=fixture.subtype_parent,
            blinding_key=REVIEW_BLINDING_KEY,
        ),
    }
    old_references = {
        target: dict(
            zip(
                fixture.ledgers[target]["object_id"].astype(str),
                fixture.ledgers[target]["reference_label"].astype(str),
                strict=True,
            )
        )
        for target in fixture.ledgers
    }
    ledgers = {
        target: _ledger(
            sample,
            [old_references[target][key] for key in sample["object_id"].astype(str)],
        )
        for target, sample in samples.items()
    }
    reviews = _review_provenances(
        bundle=fixture.bundle,
        manifest=manifest,
        samples=samples,
        ledgers=ledgers,
    )
    with pytest.raises(ValueError, match="every probabilistic target"):
        _build(
            fixture,
            calibration_candidate_manifest=manifest,
            calibration_sample_frames=samples,
            calibration_review_ledgers=ledgers,
            calibration_review_provenance=reviews,
        )


def test_every_modeled_class_requires_positive_and_negative_review_support(tmp_path):
    fixture = _fixture(tmp_path)
    ledgers = {key: value.copy() for key, value in fixture.ledgers.items()}
    subtype_target = ("subtype", fixture.subtype_parent)
    ledgers[subtype_target]["reference_label"] = fixture.bundle.subtypes[0].classes[0]
    reviews = _review_provenances(
        bundle=fixture.bundle,
        manifest=fixture.manifest,
        samples=fixture.samples,
        ledgers=ledgers,
    )
    with pytest.raises(ValueError, match="positive or negative support"):
        _build(
            fixture,
            calibration_review_ledgers=ledgers,
            calibration_review_provenance=reviews,
        )


def test_already_frozen_bundle_cannot_be_retuned(tmp_path):
    fixture = _fixture(tmp_path)
    provenance = _build(fixture)
    frozen = fixture.bundle.freeze_threshold_selection(provenance, bundle_version="frozen-v1")
    with pytest.raises(ValueError, match="preliminary/unfrozen"):
        _build(fixture, model_bundle=frozen)


def test_absent_threshold_and_readdressed_false_metric_fail(tmp_path):
    result = _build(_fixture(tmp_path))
    absent = result.to_dict()
    absent.pop("selected_probability_threshold")
    absent = _readdress(absent)
    with pytest.raises(ValueError, match="missing required fields"):
        ThresholdSelectionProvenance.from_dict(absent)

    false_metric = result.to_dict()
    false_metric["selected_metrics"]["coverage"] -= 0.1
    false_metric = _readdress(false_metric)
    with pytest.raises(ValueError, match="selected_metrics differs"):
        ThresholdSelectionProvenance.from_dict(false_metric)


def test_readdressed_false_grid_metric_cannot_override_scientific_rows(tmp_path):
    result = _build(_fixture(tmp_path))
    payload = result.to_dict()
    payload["operating_grid"][0]["error_cost"] += 0.1
    grid = pd.DataFrame(
        [{name: row[name] for name in OPERATING_POINT_COLUMNS} for row in payload["operating_grid"]]
    )
    payload["operating_grid_fingerprint"] = operating_grid_fingerprint(grid)
    payload = _readdress(payload)
    with pytest.raises(ValueError, match="embedded operating grid differs"):
        ThresholdSelectionProvenance.from_dict(payload)


def test_resigned_embedded_reviewer_row_cannot_bypass_review_fingerprint(tmp_path):
    fixture = _fixture(tmp_path)
    payload = _build(fixture).to_dict()
    payload["calibration_inputs"][0]["rows"][0]["reviewer_1"] = "forged-reviewer"
    payload = _readdress(payload)
    with pytest.raises(
        ValueError,
        match="calibration_frame_fingerprint|embedded calibration ledger differs",
    ):
        ThresholdSelectionProvenance.from_dict(payload)


def test_resigned_provenance_cannot_make_grid_more_permissive_than_seed(tmp_path):
    payload = _build(_fixture(tmp_path)).to_dict()
    payload["preliminary_assignment_thresholds"]["inferred_probability"] = 0.51
    payload = _readdress(payload)
    with pytest.raises(ValueError, match="cannot be more permissive"):
        ThresholdSelectionProvenance.from_dict(payload)


def test_strict_json_rejects_unknown_duplicate_and_nonfinite(tmp_path):
    result = _build(_fixture(tmp_path))
    unknown = result.to_dict()
    unknown["passed"] = True
    with pytest.raises(ValueError, match="unknown fields"):
        ThresholdSelectionProvenance.from_dict(unknown)

    duplicate_path = tmp_path / "duplicate.json"
    duplicate_path.write_text('{"schema_version": 2, "schema_version": 2}\n')
    with pytest.raises(ValueError, match="duplicate key"):
        ThresholdSelectionProvenance.read_json(duplicate_path)

    nonfinite_path = tmp_path / "nonfinite.json"
    nonfinite_path.write_text('{"value": NaN}\n')
    with pytest.raises(ValueError, match="non-finite"):
        ThresholdSelectionProvenance.read_json(nonfinite_path)


def test_policy_is_strict_content_addressed_and_immutable(tmp_path):
    policy = _policy()
    assert ThresholdSelectionPolicy.from_dict(policy.to_dict()) == policy
    assert replace(policy, minimum_coverage=0.25, content_id="").content_id != policy.content_id
    with pytest.raises(ValueError, match="at least one"):
        _policy(false_positive_cost=0, false_negative_cost=0)

    path = tmp_path / "policy.json"
    policy.write_json(path)
    assert ThresholdSelectionPolicy.read_json(path) == policy
    with pytest.raises(FileExistsError, match="immutable"):
        policy.write_json(path)
