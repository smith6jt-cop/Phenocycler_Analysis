from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    ContractError,
    InputArtifact,
    PartialArtifactError,
    StageManifest,
    StaleArtifactError,
)
from phenocycler.donor_pipeline import (
    DONOR_STAGES,
    DonorStageReceipt,
    DonorTask,
    ResourceLimits,
    ResourceState,
    _aggregate_donor_audit,
    _initial_task_state,
    _input_from_receipt,
    _run_donor_stage_body,
    _seal_ready_stages,
    _terminate_executor_workers,
    _validate_task_receipt,
    donor_receipt_path,
    donor_stage_inputs,
    run_donor_queue,
    select_ready_task,
    snapshot_donor_partition,
    task_ready,
    task_resources,
)
from phenocycler.expression import write_donor_partition

REPOSITORY = Path(__file__).resolve().parents[1]


def _frame(donor: str, objects=("a", "b")) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": [donor] * len(objects),
            "object_id": list(objects),
            "value": range(len(objects)),
        }
    )


def test_donor_receipt_validates_output_input_config_and_code(tmp_path):
    output_root = tmp_path / "states"
    write_donor_partition(_frame("1"), output_root, donor_id="1")
    source = tmp_path / "source.json"
    source.write_text('{"version": 1}\n')
    inputs = (InputArtifact.from_path("source", source),)
    config = {"stage": "states", "donor_id": "1"}
    code = CodeMetadata.capture(
        REPOSITORY,
        source_paths=("phenocycler/state_markers.py",),
    )

    receipt = DonorStageReceipt.create(
        stage="states",
        method_version="test-v1",
        donor_id="1",
        inputs=inputs,
        output_root=output_root,
        sidecar_paths=(),
        config=config,
        code=code,
    )
    path = tmp_path / "manifests" / "states.json"
    receipt.write_json(path)
    current = DonorStageReceipt.read_json(path)
    current.validate_current(config=config, validation_mode="content")
    assert current.output.row_count == 2
    assert current.output.object_id_sha256 == snapshot_donor_partition(
        output_root, "1"
    ).object_id_sha256

    _frame("1", objects=("a", "changed")).to_parquet(
        output_root / "donor_id=1" / "data.parquet",
        index=False,
    )
    with pytest.raises(StaleArtifactError):
        current.validate_current(
            config=config,
            validate_code=False,
            validate_inputs=False,
        )


def test_atomic_donor_write_is_immutable_and_leaves_no_temporary_file(tmp_path):
    root = tmp_path / "output"
    path = write_donor_partition(_frame("1"), root, donor_id="1")
    assert path.is_file()
    assert list(path.parent.glob("*.tmp")) == []
    with pytest.raises(FileExistsError, match="immutable donor artifact"):
        write_donor_partition(_frame("1"), root, donor_id="1")


def test_dependency_receipt_must_match_the_requested_donor(tmp_path):
    output_root = tmp_path / "geometry"
    write_donor_partition(_frame("1"), output_root, donor_id="1")
    receipt = DonorStageReceipt.create(
        stage="geometry",
        method_version="test-v1",
        donor_id="1",
        inputs=(),
        output_root=output_root,
        sidecar_paths=(),
        config={"donor_id": "1"},
        code=CodeMetadata.capture(
            REPOSITORY,
            source_paths=("phenocycler/geometry_qc.py",),
        ),
    )
    path = tmp_path / "geometry.json"
    receipt.write_json(path)

    with pytest.raises(ContractError, match="expected geometry/2"):
        _input_from_receipt("geometry", path, donor_id="2")


def test_task_receipt_must_own_the_canonical_stage_output(tmp_path):
    wrong_root = tmp_path / "wrong"
    write_donor_partition(_frame("1"), wrong_root, donor_id="1")
    receipt = DonorStageReceipt.create(
        stage="geometry",
        method_version="geometry-raster-qc-v1",
        donor_id="1",
        inputs=(),
        output_root=wrong_root,
        sidecar_paths=(),
        config={},
        code=CodeMetadata.capture(
            REPOSITORY,
            source_paths=("phenocycler/geometry_qc.py",),
        ),
    )
    context = SimpleNamespace(
        output_root=lambda _stage: tmp_path / "geometry"
    )

    with pytest.raises(ContractError, match="output root differs"):
        _validate_task_receipt(context, "geometry", "1", receipt)


def test_task_readiness_follows_the_donor_local_dag():
    completed = {DonorTask("geometry", "1")}
    assert task_ready(DonorTask("redsea", "1"), completed)
    assert not task_ready(DonorTask("expression", "1"), completed)
    completed.add(DonorTask("redsea", "1"))
    assert task_ready(DonorTask("expression", "1"), completed)
    completed.add(DonorTask("expression", "1"))
    assert task_ready(DonorTask("controls", "1"), completed)
    assert task_ready(DonorTask("states", "1"), completed)


def test_resource_selection_keeps_one_gpu_redsea_job_and_feeds_geometry():
    context = SimpleNamespace(
        donors=("1", "2"),
        config=SimpleNamespace(
            n_jobs=4,
            use_gpu=True,
            pipeline_geometry_workers=2,
            pipeline_redsea_workers=4,
            pipeline_downstream_workers=4,
        ),
    )
    limits = ResourceLimits.from_context(context, task_count=10)
    assert limits.gpu_slots == 1
    assert limits.redsea_workers == 1
    resources = ResourceState()
    completed = {DonorTask("geometry", "1")}
    pending = {
        DonorTask("redsea", "1"),
        DonorTask("geometry", "2"),
    }

    first = select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    )
    assert first == DonorTask("redsea", "1")
    resources.reserve(task_resources(context, first))
    pending.remove(first)
    second = select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    )
    assert second == DonorTask("geometry", "2")
    assert not resources.can_start(
        task_resources(context, DonorTask("redsea", "2")),
        limits,
    )


def test_queue_abort_terminates_owned_workers_without_waiting():
    class Process:
        def __init__(self):
            self.alive = True
            self.terminated = False
            self.killed = False
            self.joins = []

        def is_alive(self):
            return self.alive

        def terminate(self):
            self.terminated = True

        def kill(self):
            self.killed = True
            self.alive = False

        def join(self, timeout):
            self.joins.append(timeout)

    class Executor:
        def __init__(self, process):
            self._processes = {1: process}
            self.shutdown_call = None

        def shutdown(self, *, wait, cancel_futures):
            self.shutdown_call = (wait, cancel_futures)

    process = Process()
    executor = Executor(process)
    _terminate_executor_workers({"geometry": executor})
    assert process.terminated
    assert process.killed
    assert process.joins == [2.0, 2.0]
    assert executor.shutdown_call == (False, True)


def test_serial_queue_obeys_dependencies_without_stage_barriers(
    monkeypatch,
):
    from phenocycler import donor_pipeline

    context = SimpleNamespace(
        donors=("1", "2"),
        config=SimpleNamespace(
            n_jobs=1,
            use_gpu=False,
            pipeline_geometry_workers=1,
            pipeline_redsea_workers=1,
            pipeline_downstream_workers=1,
        ),
    )
    events = []
    monkeypatch.setattr(
        donor_pipeline,
        "_initial_task_state",
        lambda *_: (set(), set()),
    )
    monkeypatch.setattr(
        donor_pipeline,
        "_seal_ready_stages",
        lambda *_: None,
    )
    monkeypatch.setattr(
        donor_pipeline,
        "execute_donor_stage",
        lambda _context, stage, donor: events.append((stage, donor)),
    )

    run_donor_queue(
        context,
        ("geometry", "redsea", "expression"),
    )
    positions = {
        event: index for index, event in enumerate(events)
    }
    for donor in context.donors:
        assert positions[("geometry", donor)] < positions[("redsea", donor)]
        assert positions[("redsea", donor)] < positions[("expression", donor)]
    assert positions[("redsea", "1")] < positions[("geometry", "2")]

    with pytest.raises(ContractError, match="dependency prefix"):
        run_donor_queue(context, ("geometry", "expression"))


def test_donor_audits_aggregate_in_cohort_order(tmp_path):
    context = SimpleNamespace(
        donors=("2", "1"),
        run_root=tmp_path,
        config=SimpleNamespace(audit_dir=tmp_path / "audit"),
    )
    for donor in context.donors:
        root = (
            context.config.audit_dir
            / "donors"
            / "expression_availability"
        )
        write_donor_partition(
            pd.DataFrame(
                {"donor_id": [donor], "marker": ["CD3e"]}
            ),
            root,
            donor_id=donor,
        )

    path = _aggregate_donor_audit(context, "expression")
    assert path == tmp_path / "audit" / "expression_availability.parquet"
    assert pd.read_parquet(path)["donor_id"].tolist() == ["2", "1"]
    assert _aggregate_donor_audit(context, "expression") == path


def test_all_declared_donor_stages_have_resource_requests():
    context = SimpleNamespace(
        config=SimpleNamespace(use_gpu=False)
    )
    assert {
        task_resources(context, DonorTask(stage, "1")).group
        for stage in DONOR_STAGES
    } == {"geometry", "redsea", "downstream"}


def test_donor_type_inputs_reuse_captured_bundle_artifact(monkeypatch, tmp_path):
    from phenocycler import donor_pipeline, pipeline

    captured_bundle = object()
    captured_promotion = object()
    generated = []
    validated = []
    monkeypatch.setattr(
        donor_pipeline,
        "_dependency_input",
        lambda _context, stage, donor: (stage, donor),
    )
    monkeypatch.setattr(
        donor_pipeline.InputArtifact,
        "from_path",
        lambda name, path: generated.append((name, path)) or (name, path),
    )
    monkeypatch.setattr(
        pipeline,
        "_validate_typing_calibration_compatibility",
        lambda context: validated.append(context),
    )
    promotion = SimpleNamespace(
        input_artifact=captured_promotion,
        validate_for_bundle=lambda bundle: validated.append(bundle),
    )
    bundle = SimpleNamespace(input_artifact=captured_bundle)
    context = SimpleNamespace(
        base_config=SimpleNamespace(
            marker_registry=tmp_path / "markers.json",
            typing_rules=tmp_path / "rules.json",
        ),
        typing_model_bundle=bundle,
        typing_model_promotion_manifest=promotion,
    )

    inputs = donor_stage_inputs(context, "type", "1")
    assert inputs[-2:] == (captured_bundle, captured_promotion)
    assert validated == [context, bundle]
    assert generated == [
        ("marker_registry", tmp_path / "markers.json"),
        ("typing_rules", tmp_path / "rules.json"),
    ]


def test_probabilistic_donor_typing_waits_for_sealed_calibration(tmp_path):
    context = SimpleNamespace(
        donors=("1",),
        typing_model_bundle=object(),
        stage_manifest_path=lambda stage: tmp_path / f"{stage}.json",
        config=SimpleNamespace(use_gpu=False),
    )
    task = DonorTask("type", "1")
    pending = {task}
    completed = {DonorTask("calibrate", "1")}
    resources = ResourceState()
    limits = ResourceLimits(
        total_cpu_slots=1,
        gpu_slots=0,
        geometry_workers=1,
        redsea_workers=1,
        downstream_workers=1,
    )

    assert (
        select_ready_task(
            pending,
            completed,
            context=context,
            resources=resources,
            limits=limits,
        )
        is None
    )
    context.stage_manifest_path("calibrate").write_text("{}\n")
    assert select_ready_task(
        pending,
        completed,
        context=context,
        resources=resources,
        limits=limits,
    ) == task


def test_donor_type_body_passes_configured_bundle(monkeypatch, tmp_path):
    from phenocycler import donor_pipeline

    evidence = pd.DataFrame({"donor_id": ["1"], "object_id": ["a"]})
    assignments = evidence.assign(cell_type="Immune")
    bundle = object()
    promotion = object()
    observed = {}
    monkeypatch.setattr(
        donor_pipeline,
        "read_single_partition",
        lambda *_args, **_kwargs: evidence,
    )

    def fake_type(frame, **kwargs):
        observed.update(frame=frame, **kwargs)
        return assignments

    monkeypatch.setattr(donor_pipeline, "type_calibrated_cells", fake_type)
    monkeypatch.setattr(
        donor_pipeline,
        "write_donor_partition",
        lambda *_args, **_kwargs: tmp_path / "data.parquet",
    )
    context = SimpleNamespace(
        config=SimpleNamespace(
            marker_evidence_dir=tmp_path / "evidence",
            assignments_dir=tmp_path / "assignments",
        ),
        cohort=SimpleNamespace(images=()),
        registry="markers",
        typing_registry="rules",
        typing_model_bundle=bundle,
        typing_model_promotion_manifest=promotion,
    )

    assert _run_donor_stage_body(context, "type", "1") == ()
    assert observed["frame"] is evidence
    assert observed["marker_registry"] == "markers"
    assert observed["typing_registry"] == "rules"
    assert observed["model_bundle"] is bundle
    assert observed["model_promotion_manifest"] is promotion


def test_state_donor_body_writes_output_and_donor_audit(tmp_path):
    config = SimpleNamespace(
        selected_expression_dir=tmp_path / "expression",
        state_dir=tmp_path / "states",
        audit_dir=tmp_path / "audit",
    )
    context = SimpleNamespace(config=config, cohort=SimpleNamespace(images=()))
    expression = pd.DataFrame(
        {
            "donor_id": ["1", "1"],
            "object_id": ["a", "b"],
            "qc_analysis_eligible": [True, False],
            "qc_estimation_eligible": [True, False],
            "Ki67": [1.0, 2.0],
            "PCNA": [3.0, 4.0],
        }
    )
    write_donor_partition(
        expression,
        config.selected_expression_dir,
        donor_id="1",
    )

    sidecars = _run_donor_stage_body(context, "states", "1")
    assert len(sidecars) == 1
    assert sidecars[0].is_file()
    assert (
        config.state_dir / "donor_id=1" / "data.parquet"
    ).is_file()


def test_completed_receipts_seal_the_cohort_stage(monkeypatch, tmp_path):
    from phenocycler import pipeline

    code = CodeMetadata.capture(
        REPOSITORY,
        source_paths=("phenocycler/state_markers.py",),
    )
    stage_config = {"run_id": "run", "stage": "geometry"}
    monkeypatch.setattr(
        pipeline,
        "_stage_config",
        lambda _context, _stage: stage_config,
    )
    monkeypatch.setattr(
        pipeline,
        "_stage_inputs",
        lambda _context, _stage: (),
    )
    monkeypatch.setattr(
        pipeline,
        "_stage_code",
        lambda _context, _stage: code,
    )
    context = SimpleNamespace(
        donors=("1", "2"),
        run_root=tmp_path,
        config=SimpleNamespace(
            manifests_dir=tmp_path / "manifests",
        ),
        output_root=lambda stage: tmp_path / stage,
        stage_manifest_path=lambda stage: (
            tmp_path / "manifests" / f"{stage}.json"
        ),
    )
    for donor in context.donors:
        write_donor_partition(
            _frame(donor),
            context.output_root("geometry"),
            donor_id=donor,
        )
        receipt = DonorStageReceipt.create(
            stage="geometry",
            method_version="geometry-raster-qc-v1",
            donor_id=donor,
            inputs=(),
            output_root=context.output_root("geometry"),
            sidecar_paths=(),
            config={
                **stage_config,
                "donor_id": donor,
                "donor_receipt_version": 1,
            },
            code=code,
        )
        receipt.write_json(
            donor_receipt_path(context, "geometry", donor)
        )

    sealed = set()
    _seal_ready_stages(context, ("geometry",), sealed)
    assert sealed == {"geometry"}
    manifest = StageManifest.read_json(
        context.stage_manifest_path("geometry")
    )
    assert manifest.completed_donors == context.donors
    assert {
        Path(item.path)
        for item in manifest.sidecars
    } == {
        donor_receipt_path(context, "geometry", donor).resolve()
        for donor in context.donors
    }


def test_orphan_donor_output_is_not_silently_adopted(tmp_path):
    context = SimpleNamespace(
        donors=("1",),
        config=SimpleNamespace(
            manifests_dir=tmp_path / "manifests",
            audit_dir=tmp_path / "audit",
        ),
        output_root=lambda stage: tmp_path / stage,
        stage_manifest_path=lambda stage: (
            tmp_path / "manifests" / f"{stage}.json"
        ),
    )
    write_donor_partition(
        _frame("1"),
        context.output_root("geometry"),
        donor_id="1",
    )

    with pytest.raises(
        PartialArtifactError, match="output exists without a donor receipt"
    ):
        _initial_task_state(context, ("geometry",))
