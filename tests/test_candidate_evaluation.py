from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    DuplicateObjectError,
    InputArtifact,
    RunManifest,
    StageManifest,
    StaleArtifactError,
    snapshot_parquet_dataset,
)
from phenocycler.candidate_evaluation import (
    CandidateEvaluationManifest,
    load_candidate_evaluation_manifest,
    run_candidate_evaluation,
)
from phenocycler.hierarchical_typing import (
    FEATURE_SCHEMA_VERSION,
    AssignmentThresholds,
    TypingRegistry,
    classifier_from_registry,
)
from phenocycler.marker_registry import load_registry
from phenocycler.refinement_contracts import (
    DEVELOPMENT_LABEL_PURPOSE,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
)
from phenocycler.typing_model_bundle import (
    FrozenModelEnsemble,
    TypingModelBundle,
    load_typing_model_bundle,
)


@dataclass(frozen=True)
class CandidateFixture:
    evidence: Path
    source_run_path: Path
    bundle_path: Path
    splits: DonorSplitManifest


def _evidence_frame(donor: str) -> pd.DataFrame:
    registry = load_registry()
    frame = pd.DataFrame(
        {
            "donor_id": [donor, donor],
            "object_id": [f"{donor}-immune", f"{donor}-other"],
            "qc_analysis_eligible": [True, True],
        }
    )
    for marker in registry.active_markers:
        frame[f"{marker.name}__state"] = "negative"
        frame[f"{marker.name}__log2_threshold_ratio"] = -1.0
        frame[f"{marker.name}__model_valid"] = True
    frame.loc[0, "CD3e__state"] = "positive"
    frame.loc[0, "CD3e__log2_threshold_ratio"] = 3.0
    return frame


def _fit_parameters(rules: TypingRegistry) -> dict[str, object]:
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
            {str(rule.parent) for rule in rules.subtype_rules} - set(structural)
        ),
        "final_specific_classes": final_specific_classes,
        "minimum_valid_bootstrap_fraction": 0.8,
        "class_priors": "uniform_no_intercept",
        "donor_weighting": "class_then_donor_balanced",
    }


def _candidate_fixture(tmp_path: Path) -> CandidateFixture:
    splits = DonorSplitManifest(
        split_version="candidate-split-v1",
        source_run_id="source-run",
        cohort_content_id="source-cohort",
        development=("D1", "D2"),
        calibration=("D3",),
        validation=("D4",),
        locked_test=("D5",),
    )
    evidence = tmp_path / "marker_evidence"
    for donor in splits.donors:
        partition = evidence / f"donor_id={donor}"
        partition.mkdir(parents=True)
        _evidence_frame(donor).to_parquet(partition / "evidence.parquet", index=False)

    producer = tmp_path / "calibration_producer.py"
    producer.write_text("METHOD_VERSION = 'calibrate-test-v1'\n", encoding="utf-8")
    raw = tmp_path / "raw-input.txt"
    raw.write_text("immutable source\n", encoding="utf-8")
    code = CodeMetadata.capture(tmp_path, source_paths=(producer,))
    stage = StageManifest.create(
        stage="calibrate",
        method_version="calibrate-test-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=evidence,
        config={"test": "candidate-evaluation"},
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "calibrate.manifest.json")
    source_run = RunManifest.create(
        run_name=splits.source_run_id,
        expected_donors=splits.donors,
        stage_manifest_paths=(stage_path,),
        config={"test": "candidate-evaluation"},
        code=code,
    )
    source_run_path = source_run.write_json(tmp_path / "source-run.manifest.json")

    markers = load_registry()
    rules = TypingRegistry.from_marker_registry(markers)
    model = classifier_from_registry(rules)
    coefficients = tuple(
        tuple(float(item) for item in row) for row in np.asarray(model.coefficients, dtype=float)
    )
    ensemble = FrozenModelEnsemble(
        level="broad",
        parent=None,
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
    label_provenance = DevelopmentLabelProvenance(
        bundle_id="separate-development-labels-v1",
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
        bundle_version="unpromoted-candidate-v1",
        marker_registry_fingerprint=markers.fingerprint,
        typing_rules_fingerprint=rules.rules_fingerprint,
        split_manifest=splits,
        label_provenance=label_provenance,
        evidence_provenance=evidence_provenance,
        thresholds=AssignmentThresholds(),
        broad=ensemble,
        fit_parameters=_fit_parameters(rules),
    )
    bundle_path = bundle.write_json(tmp_path / "candidate-model.json")
    return CandidateFixture(
        evidence=evidence,
        source_run_path=source_run_path,
        bundle_path=bundle_path,
        splits=splits,
    )


def _loaded_dependencies(fixture: CandidateFixture):
    markers = load_registry()
    rules = TypingRegistry.from_marker_registry(markers)
    bundle = load_typing_model_bundle(
        fixture.bundle_path,
        marker_registry=markers,
        typing_registry=rules,
    )
    source_run = RunManifest.read_json(fixture.source_run_path)
    return markers, rules, bundle, source_run


def _run(
    fixture: CandidateFixture,
    tmp_path: Path,
    *,
    name: str,
    split_name: str = "validation",
    allow_locked_test: bool = False,
):
    return run_candidate_evaluation(
        model_bundle_path=fixture.bundle_path,
        evidence_path=fixture.evidence,
        source_run_manifest_path=fixture.source_run_path,
        output_root=tmp_path / f"{name}-assignments",
        manifest_path=tmp_path / f"{name}.manifest.json",
        split_name=split_name,
        allow_locked_test=allow_locked_test,
    )


def _resnapshot_candidate_manifest(manifest_path: Path, output_root: Path) -> None:
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    snapshot = snapshot_parquet_dataset(
        output_root,
        expected_donors=payload["donors"],
    )
    payload["output"] = snapshot.to_dict()
    payload["output_content_sha256"] = snapshot.content_sha256
    payload["content_id"] = ""
    resigned = CandidateEvaluationManifest.from_dict(payload)
    manifest_path.write_text(
        json.dumps(resigned.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def test_candidate_evaluation_reads_and_writes_validation_only(tmp_path, monkeypatch):
    fixture = _candidate_fixture(tmp_path)
    from phenocycler import candidate_evaluation

    observed: list[str] = []
    original = candidate_evaluation._read_evidence_partition

    def recording_reader(root, donor_snapshot):
        observed.append(donor_snapshot.donor_id)
        return original(root, donor_snapshot)

    monkeypatch.setattr(candidate_evaluation, "_read_evidence_partition", recording_reader)
    manifest = _run(fixture, tmp_path, name="validation")

    # Generation and the independent persisted-output replay both remain
    # confined to the selected donor.
    assert observed == ["D4", "D4"]
    assert manifest.split_name == "validation"
    assert manifest.donors == fixture.splits.validation
    assert manifest.output.donor_ids == ("D4",)
    output = pd.read_parquet(
        tmp_path / "validation-assignments" / "donor_id=D4" / "assignments.parquet"
    )
    assert set(output["donor_id"].astype(str)) == {"D4"}
    assert output["candidate_assignment_run_id"].eq(manifest.assignment_run_id).all()
    assert output["typing_model_bundle_content_id"].eq(manifest.model_bundle_content_id).all()
    assert output["typing_mode"].eq("probabilistic_candidate").all()
    assert not (tmp_path / "validation-assignments" / "donor_id=D5").exists()

    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    loaded = load_candidate_evaluation_manifest(
        tmp_path / "validation.manifest.json",
        model_bundle=bundle,
        marker_registry=markers,
        typing_registry=rules,
        source_run_manifest=source_run,
    )
    assert loaded.content_id == manifest.content_id


def test_candidate_evaluation_allows_calibration_without_locked_flag(tmp_path, monkeypatch):
    fixture = _candidate_fixture(tmp_path)
    from phenocycler import candidate_evaluation

    observed: list[str] = []
    original = candidate_evaluation._read_evidence_partition

    def recording_reader(root, donor_snapshot):
        observed.append(donor_snapshot.donor_id)
        return original(root, donor_snapshot)

    monkeypatch.setattr(candidate_evaluation, "_read_evidence_partition", recording_reader)
    manifest = _run(fixture, tmp_path, name="calibration", split_name="calibration")

    assert observed == ["D3", "D3"]
    assert manifest.split_name == "calibration"
    assert manifest.donors == fixture.splits.calibration
    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    loaded = load_candidate_evaluation_manifest(
        tmp_path / "calibration.manifest.json",
        model_bundle=bundle,
        marker_registry=markers,
        typing_registry=rules,
        source_run_manifest=source_run,
    )
    assert loaded.assignment_run_id == manifest.assignment_run_id


def test_locked_test_requires_explicit_flag_for_generation_and_loading(tmp_path):
    fixture = _candidate_fixture(tmp_path)
    with pytest.raises(PermissionError, match="explicit"):
        _run(fixture, tmp_path, name="locked-denied", split_name="locked_test")
    assert not (tmp_path / "locked-denied-assignments").exists()

    manifest = _run(
        fixture,
        tmp_path,
        name="locked-allowed",
        split_name="locked_test",
        allow_locked_test=True,
    )
    assert manifest.donors == ("D5",)
    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    with pytest.raises(PermissionError, match="explicit"):
        load_candidate_evaluation_manifest(
            tmp_path / "locked-allowed.manifest.json",
            model_bundle=bundle,
            marker_registry=markers,
            typing_registry=rules,
            source_run_manifest=source_run,
        )
    loaded = load_candidate_evaluation_manifest(
        tmp_path / "locked-allowed.manifest.json",
        model_bundle=bundle,
        marker_registry=markers,
        typing_registry=rules,
        source_run_manifest=source_run,
        allow_locked_test=True,
    )
    assert loaded.assignment_run_id == manifest.assignment_run_id


def test_manifest_and_assignment_tampering_fail_closed(tmp_path):
    fixture = _candidate_fixture(tmp_path)
    _run(fixture, tmp_path, name="tamper")
    manifest_path = tmp_path / "tamper.manifest.json"
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["split_name"] = "locked_test"
    with pytest.raises(ValueError, match="assignment_run_id|content_id"):
        CandidateEvaluationManifest.from_dict(payload)

    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload["unexpected"] = True
    with pytest.raises(ValueError, match="unknown fields"):
        CandidateEvaluationManifest.from_dict(payload)

    assignment_path = tmp_path / "tamper-assignments" / "donor_id=D4" / "assignments.parquet"
    changed = pd.read_parquet(assignment_path)
    changed.loc[0, "broad_label"] = "tampered"
    changed.to_parquet(assignment_path, index=False)
    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    with pytest.raises(StaleArtifactError, match="changed"):
        load_candidate_evaluation_manifest(
            manifest_path,
            model_bundle=bundle,
            marker_registry=markers,
            typing_registry=rules,
            source_run_manifest=source_run,
        )


def test_outputs_and_manifest_are_immutable(tmp_path):
    fixture = _candidate_fixture(tmp_path)
    first = _run(fixture, tmp_path, name="immutable")
    with pytest.raises(FileExistsError, match="immutable"):
        _run(fixture, tmp_path, name="immutable")
    with pytest.raises(FileExistsError, match="immutable"):
        first.write_json(tmp_path / "immutable.manifest.json")


def test_assignment_run_id_is_deterministic_and_promotion_independent(tmp_path):
    fixture = _candidate_fixture(tmp_path)
    first = _run(fixture, tmp_path, name="first")
    second = _run(fixture, tmp_path, name="second")

    assert first.assignment_run_id == second.assignment_run_id
    assert first.model_bundle_content_id == second.model_bundle_content_id
    assert first.output.root != second.output.root
    assert "promotion" not in first.to_dict()
    # Output paths are validation caches, not scientific identity.  Relocating
    # byte-identical assignments must not change downstream review draws or
    # bootstrap seeds.
    assert first.content_id == second.content_id


def test_evidence_must_be_the_exact_source_calibrate_dataset(tmp_path):
    fixture = _candidate_fixture(tmp_path)
    copied = tmp_path / "copied-evidence"
    copied.mkdir()
    with pytest.raises(ValueError, match="not the calibration output"):
        run_candidate_evaluation(
            model_bundle_path=fixture.bundle_path,
            evidence_path=copied,
            source_run_manifest_path=fixture.source_run_path,
            output_root=tmp_path / "wrong-source-assignments",
            manifest_path=tmp_path / "wrong-source.manifest.json",
        )
    assert not (tmp_path / "wrong-source-assignments").exists()


def test_generation_rejects_monkeypatched_row_dropping(tmp_path, monkeypatch):
    fixture = _candidate_fixture(tmp_path)
    from phenocycler import candidate_evaluation

    original = candidate_evaluation.type_candidate_cells

    def dropping_typer(*args, **kwargs):
        return original(*args, **kwargs).iloc[:-1].copy()

    monkeypatch.setattr(candidate_evaluation, "type_candidate_cells", dropping_typer)
    with pytest.raises(ValueError, match="source evidence has"):
        _run(fixture, tmp_path, name="dropped")
    assert not (tmp_path / "dropped.manifest.json").exists()


def test_generation_rejects_duplicate_assignment_objects(tmp_path, monkeypatch):
    fixture = _candidate_fixture(tmp_path)
    from phenocycler import candidate_evaluation

    original = candidate_evaluation.type_candidate_cells

    def duplicate_typer(*args, **kwargs):
        output = original(*args, **kwargs)
        output.loc[1, "object_id"] = output.loc[0, "object_id"]
        return output

    monkeypatch.setattr(candidate_evaluation, "type_candidate_cells", duplicate_typer)
    with pytest.raises(DuplicateObjectError, match="duplicate object_id"):
        _run(fixture, tmp_path, name="duplicate")
    assert not (tmp_path / "duplicate.manifest.json").exists()


@pytest.mark.parametrize(
    ("mutation", "expected_exception", "message"),
    [
        ("subset", ValueError, "row count"),
        ("missing_schema", KeyError, "missing required columns"),
        ("wrong_boolean_schema", TypeError, "non-null booleans"),
        ("candidate_assignment_run_id", ValueError, "candidate_assignment_run_id"),
        ("candidate_evaluation_split_name", ValueError, "candidate_evaluation_split_name"),
        ("typing_model_bundle_content_id", ValueError, "typing_model_bundle_content_id"),
        ("typing_model_bundle_version", ValueError, "typing_model_bundle_version"),
        ("typing_mode", ValueError, "typing_mode"),
        ("typing_feature_schema_version", ValueError, "typing_feature_schema_version"),
        ("typing_score_semantics", ValueError, "typing_score_semantics"),
    ],
)
def test_resnapshotted_invalid_output_fails_current_validation(
    tmp_path,
    mutation,
    expected_exception,
    message,
):
    fixture = _candidate_fixture(tmp_path)
    _run(fixture, tmp_path, name=f"resnapshot-{mutation}")
    output_root = tmp_path / f"resnapshot-{mutation}-assignments"
    assignment_path = output_root / "donor_id=D4" / "assignments.parquet"
    frame = pd.read_parquet(assignment_path)
    if mutation == "subset":
        frame = frame.iloc[:-1].copy()
    elif mutation == "missing_schema":
        frame = frame.drop(columns=["broad_reason"])
    elif mutation == "wrong_boolean_schema":
        frame["broad_threshold_eligible"] = frame["broad_threshold_eligible"].astype(str)
    else:
        frame[mutation] = "wrong-metadata"
    frame.to_parquet(assignment_path, index=False)
    manifest_path = tmp_path / f"resnapshot-{mutation}.manifest.json"
    _resnapshot_candidate_manifest(manifest_path, output_root)

    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    with pytest.raises(expected_exception, match=message):
        load_candidate_evaluation_manifest(
            manifest_path,
            model_bundle=bundle,
            marker_registry=markers,
            typing_registry=rules,
            source_run_manifest=source_run,
        )


@pytest.mark.parametrize(
    ("column", "replacement"),
    [
        ("broad_label", "fabricated_label"),
        ("broad_assignment_status", "ambiguous"),
        ("broad_reason", "fabricated_reason"),
        ("broad_best_probability", 0.123456789),
        ("broad_probability_pass", True),
    ],
)
def test_resnapshotted_scientific_edits_fail_deterministic_replay(
    tmp_path,
    column,
    replacement,
):
    fixture = _candidate_fixture(tmp_path)
    name = f"replay-{column}"
    _run(fixture, tmp_path, name=name)
    output_root = tmp_path / f"{name}-assignments"
    assignment_path = output_root / "donor_id=D4" / "assignments.parquet"
    frame = pd.read_parquet(assignment_path)
    target = frame["object_id"].eq("D4-other")
    assert target.sum() == 1
    if column == "broad_probability_pass":
        replacement = not bool(frame.loc[target, column].iloc[0])
    frame.loc[target, column] = replacement
    frame.to_parquet(assignment_path, index=False)
    manifest_path = tmp_path / f"{name}.manifest.json"
    _resnapshot_candidate_manifest(manifest_path, output_root)

    markers, rules, bundle, source_run = _loaded_dependencies(fixture)
    with pytest.raises(StaleArtifactError, match="deterministic replay"):
        load_candidate_evaluation_manifest(
            manifest_path,
            model_bundle=bundle,
            marker_registry=markers,
            typing_registry=rules,
            source_run_manifest=source_run,
        )


@pytest.mark.parametrize(
    "invalid_json",
    [
        '{"schema_version": 1, "schema_version": 1}',
        '{"schema_version": NaN}',
        '{"schema_version": Infinity}',
        '{"schema_version": -Infinity}',
        '{"schema_version": 1e999}',
    ],
)
def test_manifest_json_rejects_duplicate_keys_and_nonfinite_numbers(
    tmp_path,
    invalid_json,
):
    path = tmp_path / "invalid-candidate.manifest.json"
    path.write_text(invalid_json, encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate JSON key|non-finite JSON number"):
        CandidateEvaluationManifest.read_json(path)
