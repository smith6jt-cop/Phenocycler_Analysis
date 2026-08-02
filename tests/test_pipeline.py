from __future__ import annotations

from types import SimpleNamespace

import pytest

from phenocycler.artifacts import ContractError, StaleArtifactError
from phenocycler.config import PipelineConfig
from phenocycler.pipeline import (
    STAGE_ORDER,
    _identity_config,
    _selected_stages,
    _validate_export,
)


def test_stage_selection_is_explicit_and_ordered():
    assert _selected_stages(None, None) == STAGE_ORDER
    assert _selected_stages(None, "expression") == (
        "ingest",
        "geometry",
        "redsea",
        "expression",
    )
    assert _selected_stages(["type"], None) == ("type",)


def test_probabilistic_evaluate_dispatches_to_candidate_evaluation(monkeypatch):
    from phenocycler import candidate_evaluation, pipeline

    observed = []
    monkeypatch.setattr(
        candidate_evaluation,
        "main",
        lambda arguments: observed.append(arguments) or 17,
    )

    assert pipeline.main(["probabilistic", "evaluate", "--model-bundle", "x"]) == 17
    assert observed == [["--model-bundle", "x"]]


def test_typing_bundle_semantics_are_part_of_run_identity_config():
    cfg = PipelineConfig()
    rules_only = _identity_config(cfg)
    bundle = SimpleNamespace(
        content_id="a" * 64,
        method_version="probabilistic-v1",
        feature_schema_version="features-v1",
    )
    promotion = SimpleNamespace(content_id="b" * 64)
    probabilistic = _identity_config(
        cfg,
        typing_model_bundle=bundle,
        typing_model_promotion_manifest=promotion,
    )

    assert rules_only["typing_model"] == {"mode": "rules_only"}
    assert probabilistic["typing_model"] == {
        "mode": "probabilistic_bundle",
        "content_id": "a" * 64,
        "method_version": "probabilistic-v1",
        "feature_schema_version": "features-v1",
        "promotion_manifest_content_id": "b" * 64,
    }
    assert probabilistic != rules_only


def test_type_stage_uses_captured_bundle_input_artifact(monkeypatch, tmp_path):
    from phenocycler import pipeline

    dependency = tmp_path / "calibrate.json"
    dependency.write_text("{}\n")
    captured_bundle = object()
    captured_promotion = object()
    generated = []
    monkeypatch.setattr(
        pipeline.InputArtifact,
        "from_stage_manifest",
        lambda name, path: ("stage", name, path),
    )
    monkeypatch.setattr(
        pipeline.InputArtifact,
        "from_path",
        lambda name, path: generated.append((name, path)) or (name, path),
    )
    monkeypatch.setattr(
        pipeline,
        "_validate_typing_calibration_compatibility",
        lambda _context: None,
    )
    context = SimpleNamespace(
        base_config=SimpleNamespace(
            marker_registry=tmp_path / "markers.json",
            typing_rules=tmp_path / "rules.json",
        ),
        typing_model_bundle=SimpleNamespace(input_artifact=captured_bundle),
        typing_model_promotion_manifest=SimpleNamespace(
            input_artifact=captured_promotion,
            validate_for_bundle=lambda _bundle: None,
        ),
        stage_manifest_path=lambda stage: tmp_path / f"{stage}.json",
    )

    inputs = pipeline._stage_inputs(context, "type")
    assert inputs[-2:] == (captured_bundle, captured_promotion)
    assert generated == [
        ("marker_registry", tmp_path / "markers.json"),
        ("typing_rules", tmp_path / "rules.json"),
    ]


def test_type_stage_passes_verified_promotion_capability(monkeypatch, tmp_path):
    from phenocycler import pipeline

    evidence = object()
    assignments = object()
    bundle = object()
    promotion = object()
    observed: dict[str, object] = {}
    monkeypatch.setattr(
        pipeline,
        "read_single_partition",
        lambda *_args, **_kwargs: evidence,
    )

    def fake_type(frame, **kwargs):
        observed.update(frame=frame, **kwargs)
        return assignments

    monkeypatch.setattr(pipeline, "type_calibrated_cells", fake_type)
    monkeypatch.setattr(
        pipeline,
        "write_donor_partition",
        lambda frame, *_args, **_kwargs: observed.update(written=frame),
    )
    context = SimpleNamespace(
        donors=("D1",),
        config=SimpleNamespace(
            marker_evidence_dir=tmp_path / "evidence",
            assignments_dir=tmp_path / "assignments",
        ),
        registry="markers",
        typing_registry="rules",
        typing_model_bundle=bundle,
        typing_model_promotion_manifest=promotion,
    )

    pipeline._run_type(context)

    assert observed["frame"] is evidence
    assert observed["model_bundle"] is bundle
    assert observed["model_promotion_manifest"] is promotion
    assert observed["written"] is assignments


def test_resolve_rejects_typing_bundle_from_another_cohort(monkeypatch, tmp_path):
    from phenocycler import pipeline

    cohort = SimpleNamespace(
        content_id="current-cohort",
        expected_donors=("D1",),
        validate_current=lambda *, mode: None,
    )
    marker_registry = SimpleNamespace(fingerprint="m" * 64)
    typing_registry = SimpleNamespace(rules_fingerprint="r" * 64)
    validations = []
    promotion_validations = []
    bundle = SimpleNamespace(
        input_artifact=SimpleNamespace(
            validate_current=lambda *, mode: validations.append(mode)
        ),
        split_manifest=SimpleNamespace(cohort_content_id="different-cohort"),
    )
    promotion = SimpleNamespace(
        input_artifact=SimpleNamespace(
            validate_current=lambda *, mode: promotion_validations.append(mode)
        ),
        validate_for_bundle=lambda _bundle, **_kwargs: None,
    )
    monkeypatch.setattr(
        pipeline.QuPathCohortManifest,
        "read_json",
        lambda _path: cohort,
    )
    monkeypatch.setattr(pipeline, "load_registry", lambda _path: marker_registry)
    monkeypatch.setattr(
        pipeline.TypingRegistry,
        "from_marker_registry",
        lambda *_args, **_kwargs: typing_registry,
    )
    monkeypatch.setattr(
        pipeline,
        "load_typing_model_bundle",
        lambda *_args, **_kwargs: bundle,
    )
    monkeypatch.setattr(
        pipeline,
        "load_typing_model_promotion_manifest",
        lambda *_args, **_kwargs: promotion,
    )
    cfg = PipelineConfig(
        data_dir=tmp_path / "data",
        qupath_manifest=tmp_path / "cohort.json",
        marker_registry=tmp_path / "markers.json",
        typing_rules=tmp_path / "rules.json",
        typing_model_bundle=tmp_path / "typing-model.json",
        typing_model_promotion_manifest=tmp_path / "typing-promotion.json",
    )
    release_key_path = tmp_path / "release-signing.key"
    release_key_path.write_bytes(b"pipeline-test-release-signing-key")
    monkeypatch.setenv(
        "PHENOCYCLER_RELEASE_SIGNING_KEY_FILE", str(release_key_path)
    )

    with pytest.raises(ContractError, match="different cohort"):
        pipeline.resolve_run_context(cfg)
    assert validations == ["content"]
    assert promotion_validations == ["content"]
    assert not cfg.runs_dir.exists()


def test_malformed_export_manifest_fails_with_a_contract_error(tmp_path):
    manifests = tmp_path / "manifests"
    manifests.mkdir()
    (manifests / "qupath_export.json").write_text("{}\n")
    context = SimpleNamespace(
        config=SimpleNamespace(manifests_dir=manifests),
        run_id="run",
    )

    with pytest.raises(ContractError, match="missing content_id"):
        _validate_export(context)


def test_stage_inputs_are_validated_before_the_runner(monkeypatch, tmp_path):
    from phenocycler import pipeline

    class StaleInput:
        def validate_current(self, *, mode):
            raise StaleArtifactError("upstream changed")

    context = SimpleNamespace(
        stage_manifest_path=lambda stage: tmp_path / f"{stage}.json",
        output_root=lambda stage: tmp_path / stage,
    )
    ran = False

    def runner(_context):
        nonlocal ran
        ran = True

    monkeypatch.setattr(pipeline, "_stage_inputs", lambda *_: (StaleInput(),))
    monkeypatch.setitem(pipeline.STAGE_RUNNERS, "ingest", runner)

    with pytest.raises(StaleArtifactError, match="upstream changed"):
        pipeline.run_stage(context, "ingest")
    assert not ran
    assert not (tmp_path / "ingest").exists()


def test_probabilistic_type_requires_exact_current_calibration_content(
    monkeypatch, tmp_path
):
    from phenocycler import pipeline

    manifest_path = tmp_path / "calibrate.json"
    manifest_path.write_text("{}\n")
    validations = []
    manifest = SimpleNamespace(
        stage="calibrate",
        output=SimpleNamespace(content_sha256="current-content"),
        validate_current=lambda **kwargs: validations.append(kwargs),
    )
    monkeypatch.setattr(
        pipeline.StageManifest,
        "read_json",
        lambda _path: manifest,
    )
    monkeypatch.setattr(
        pipeline,
        "_stage_config",
        lambda _context, stage: {"stage": stage},
    )
    context = SimpleNamespace(
        typing_model_bundle=SimpleNamespace(
            evidence_provenance=SimpleNamespace(
                source_artifact_content_id="training-content"
            )
        ),
        stage_manifest_path=lambda _stage: manifest_path,
    )

    with pytest.raises(ContractError, match="differs from model training evidence"):
        pipeline._validate_typing_calibration_compatibility(context)
    assert validations == [
        {
            "config": {"stage": "calibrate"},
            "validation_mode": "content",
        }
    ]

    manifest.output.content_sha256 = "training-content"
    assert (
        pipeline._validate_typing_calibration_compatibility(context)
        is manifest
    )


def test_pipelined_cli_dispatches_dependency_prefix_and_worker_overrides(
    monkeypatch,
):
    from phenocycler import donor_pipeline, pipeline

    context = object()
    loaded = {}
    dispatched = {}

    def fake_load(path, **overrides):
        loaded.update(overrides)
        return object()

    monkeypatch.setattr(pipeline, "load_config", fake_load)
    monkeypatch.setattr(
        pipeline,
        "resolve_run_context",
        lambda _config: context,
    )
    monkeypatch.setattr(
        donor_pipeline,
        "run_pipelined_pipeline",
        lambda passed, *, stages, export: dispatched.update(
            context=passed,
            stages=stages,
            export=export,
        ),
    )

    assert (
        pipeline.main(
            [
                "run",
                "--pipelined",
                "--through",
                "geometry",
                "--jobs",
                "4",
                "--geometry-workers",
                "2",
                "--redsea-workers",
                "1",
                "--downstream-workers",
                "3",
                "--no-export",
            ]
        )
        == 0
    )
    assert loaded == {
        "n_jobs": 4,
        "pipeline_geometry_workers": 2,
        "pipeline_redsea_workers": 1,
        "pipeline_downstream_workers": 3,
    }
    assert dispatched == {
        "context": context,
        "stages": ("ingest", "geometry"),
        "export": False,
    }


def test_probabilistic_pipelined_run_seals_calibration_before_type(
    monkeypatch,
):
    from phenocycler import donor_pipeline, pipeline

    context = SimpleNamespace(typing_model_bundle=object())
    calls = []
    monkeypatch.setattr(pipeline, "load_config", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(pipeline, "resolve_run_context", lambda _cfg: context)
    monkeypatch.setattr(
        pipeline,
        "_validate_typing_calibration_compatibility",
        lambda passed: calls.append(("validated", passed)),
    )
    monkeypatch.setattr(
        donor_pipeline,
        "run_pipelined_pipeline",
        lambda passed, *, stages, export: calls.append(
            ("run", passed, stages, export)
        ),
    )

    assert (
        pipeline.main(
            ["run", "--pipelined", "--through", "type", "--no-export"]
        )
        == 0
    )
    assert calls == [
        ("run", context, pipeline.STAGE_ORDER[:6], False),
        ("validated", context),
        ("run", context, pipeline.STAGE_ORDER[:7], False),
    ]


def test_status_reports_donor_receipt_progress(capsys, tmp_path):
    from phenocycler.pipeline import status

    manifests = tmp_path / "manifests"
    receipt = (
        manifests
        / "donors"
        / "donor_id=1"
        / "geometry.json"
    )
    receipt.parent.mkdir(parents=True)
    receipt.write_text("{}\n")
    context = SimpleNamespace(
        donors=("1", "2"),
        run_id="run",
        config=SimpleNamespace(manifests_dir=manifests),
        stage_manifest_path=lambda stage: manifests / f"{stage}.json",
    )

    assert status(context) == 1
    assert "geometry    IN PROGRESS   1/2" in capsys.readouterr().out
