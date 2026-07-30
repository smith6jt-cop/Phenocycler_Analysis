from __future__ import annotations

from types import SimpleNamespace

import pytest

from phenocycler.artifacts import ContractError
from phenocycler.artifacts import StaleArtifactError
from phenocycler.pipeline import (
    STAGE_ORDER,
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
