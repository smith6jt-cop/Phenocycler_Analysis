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
