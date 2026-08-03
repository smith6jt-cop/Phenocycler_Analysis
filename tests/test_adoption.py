from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from phenocycler import adoption
from phenocycler.artifacts import (
    CodeMetadata,
    ContractError,
    StageManifest,
    StaleArtifactError,
    sha256_file,
)
from phenocycler.pipeline import STAGE_DEPENDENCIES, STAGE_ORDER, STAGE_SOURCE_PATHS
from phenocycler.resume import _content_payload


@dataclass
class _TinyContext:
    base_config: object
    config: object
    donors: tuple[str, ...]
    repository: Path
    run_id: str
    run_root: Path

    def output_root(self, stage: str) -> Path:
        assert stage == "ingest"
        return self.run_root / "cells"

    def stage_manifest_path(self, stage: str) -> Path:
        return self.config.manifests_dir / f"{stage}.json"


def _write_partition(root: Path, donor: str) -> None:
    path = root / f"donor_id={donor}"
    path.mkdir(parents=True)
    pd.DataFrame(
        {"object_id": [f"{donor}-cell"], "donor_id": [donor]}
    ).to_parquet(path / "data.parquet", index=False)


def _tiny_adoption(tmp_path: Path) -> tuple[_TinyContext, dict, StageManifest]:
    repository = Path(__file__).resolve().parents[1]
    runs = tmp_path / "runs"
    source_id = "a" * 20
    target_id = "b" * 20
    source_root = runs / source_id
    target_root = runs / target_id
    source_output = source_root / "cells"
    for donor in ("1", "2"):
        _write_partition(source_output, donor)
    rejects = source_root / "ingest_rejects.parquet"
    pd.DataFrame({"reason": []}).to_parquet(rejects, index=False)
    panel = source_root / "panel_availability.json"
    panel.write_text("{}\n", encoding="utf-8")
    source = StageManifest.create(
        stage="ingest",
        method_version="qupath-ingest-union-v2",
        expected_donors=("1", "2"),
        inputs=(),
        output_root=source_output,
        sidecar_paths=(rejects, panel),
        config={"source": True},
        code=CodeMetadata.capture(
            repository, source_paths=STAGE_SOURCE_PATHS["ingest"]
        ),
    )
    source_manifest_path = source_root / "manifests" / "ingest.json"
    source_manifest_path.parent.mkdir(parents=True)
    source.write_json(source_manifest_path)

    context = _TinyContext(
        base_config=SimpleNamespace(runs_dir=runs),
        config=SimpleNamespace(manifests_dir=target_root / "manifests"),
        donors=("1", "2"),
        repository=repository,
        run_id=target_id,
        run_root=target_root,
    )
    stages = []
    for stage in STAGE_ORDER:
        if stage == "ingest":
            stages.append(
                {
                    "stage": stage,
                    "action": "adopt_exact",
                    "reason_codes": [],
                    "source_manifest_content_id": source.content_id,
                    "source_manifest_file_sha256": sha256_file(
                        source_manifest_path
                    ),
                    "source_output_content_sha256": source.output.content_sha256,
                    "source_schema_sha256": source.output.schema_sha256,
                    "source_object_id_sha256": source.output.object_id_sha256,
                    "configuration_keys": [
                        "run_schema_version",
                        "cells_min_cell_area",
                    ],
                    "source_configuration": {"cells_min_cell_area": 10},
                    "target_configuration": {"cells_min_cell_area": 10},
                    "compatibility_policy_id": None,
                    "compatibility_policy_content_id": None,
                }
            )
        else:
            stages.append(
                {
                    "stage": stage,
                    "action": "recompute",
                    "reason_codes": ["test"],
                }
            )
    plan = _content_payload(
        {
            "adoption_contract_version": 1,
            "source_run_id": source_id,
            "target_run_id": target_id,
            "validation_mode": "fast",
            "stage_order": list(STAGE_ORDER),
            "stages": stages,
            "summary": {"recompute_stages": list(STAGE_ORDER[1:])},
        }
    )
    return context, plan, source


def test_configuration_projection_is_stage_scoped_and_legacy_rules_only():
    configuration = {
        "run_schema_version": 1,
        "cells_min_cell_area": 10,
        "geometry_qc": {"duplicate_policy": "deterministic_keep"},
        "typing_model": {"mode": "probabilistic_bundle"},
    }
    ingest = adoption.project_stage_configuration(configuration, "ingest")
    geometry = adoption.project_stage_configuration(configuration, "geometry")
    typing = adoption.project_stage_configuration(
        {"run_schema_version": 1}, "type"
    )

    assert "typing_model" not in ingest
    assert "typing_model" not in geometry
    assert typing["typing_model"] == {"mode": "rules_only"}


def test_stage_graph_keeps_states_independent_of_type():
    assert STAGE_DEPENDENCIES["type"] == ("calibrate",)
    assert STAGE_DEPENDENCIES["states"] == ("expression",)
    assert "type" not in STAGE_DEPENDENCIES["states"]
    decisions = {
        "expression": {"action": "adopt_exact"},
        "calibrate": {"action": "recompute"},
    }
    assert adoption._upstream_incompatibilities("states", decisions) == []
    assert adoption._upstream_incompatibilities("type", decisions) == [
        "upstream_not_adoptable:calibrate"
    ]


def test_probabilistic_bundle_is_checked_against_planned_calibration():
    context = SimpleNamespace(
        typing_model_bundle=SimpleNamespace(
            evidence_provenance=SimpleNamespace(
                source_artifact_content_id="expected-calibration"
            )
        )
    )
    detail, reason = adoption._typing_calibration_plan(
        context,
        {
            "calibrate": {
                "action": "adopt_reviewed",
                "source_output_content_sha256": "different-calibration",
            }
        },
    )
    assert detail == {
        "bundle_evidence_content_sha256": "expected-calibration",
        "planned_calibrate_content_sha256": "different-calibration",
    }
    assert reason == "typing_bundle_calibration_mismatch"


def test_plan_identity_rejects_tampering():
    plan = _content_payload(
        {
            "adoption_contract_version": 1,
            "stage_order": list(STAGE_ORDER),
            "stages": [],
        }
    )
    adoption.verify_adoption_plan(plan)
    plan["stages"] = [{"stage": "ingest"}]
    with pytest.raises(ContractError, match="content ID is invalid"):
        adoption.verify_adoption_plan(plan)


def test_checked_in_policy_pins_the_reviewed_donor_pipeline_transition():
    policy, content_id = adoption._load_policy()
    assert len(content_id) == 64
    assert {item["stage"] for item in policy["transitions"]} == {
        "geometry",
        "redsea",
        "expression",
        "controls",
        "calibrate",
        "states",
    }
    assert {
        tuple(item["changed_paths"]) for item in policy["transitions"]
    } == {("phenocycler/donor_pipeline.py",)}


def test_changed_source_files_supports_multi_file_and_path_set_drift(tmp_path):
    repository = tmp_path / "repository"
    repository.mkdir()
    subprocess.run(["git", "init", "-q", str(repository)], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.email", "test@example.org"],
        check=True,
    )
    subprocess.run(
        ["git", "-C", str(repository), "config", "user.name", "Test"],
        check=True,
    )
    (repository / "shared.py").write_text("old shared\n", encoding="utf-8")
    (repository / "removed.py").write_text("removed\n", encoding="utf-8")
    subprocess.run(["git", "-C", str(repository), "add", "."], check=True)
    subprocess.run(
        ["git", "-C", str(repository), "commit", "-q", "-m", "source"],
        check=True,
    )
    commit = subprocess.run(
        ["git", "-C", str(repository), "rev-parse", "HEAD"],
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    ).stdout.strip()
    (repository / "shared.py").write_text("new shared\n", encoding="utf-8")
    (repository / "removed.py").unlink()
    (repository / "added.py").write_text("added\n", encoding="utf-8")

    changed = adoption._changed_source_files(
        repository,
        commit,
        ("shared.py", "removed.py"),
        ("shared.py", "added.py"),
    )
    assert set(changed) == {"added.py", "removed.py", "shared.py"}
    assert changed["added.py"][0] is None
    assert changed["removed.py"][1] is None


def test_apply_hardlinks_reseals_and_is_idempotent(tmp_path, monkeypatch):
    context, plan, source = _tiny_adoption(tmp_path)
    monkeypatch.setattr(adoption, "build_adoption_plan", lambda *args, **kwargs: plan)
    monkeypatch.setattr(adoption, "_stage_inputs", lambda *args, **kwargs: ())
    monkeypatch.setattr(
        adoption,
        "_stage_config",
        lambda context, stage: {"target": True, "stage": stage},
    )
    monkeypatch.setattr(
        adoption.StageManifest,
        "create",
        classmethod(
            lambda cls, **kwargs: pytest.fail(
                "adoption must relocate the verified snapshot, not rescan Parquet"
            )
        ),
    )

    adopted = adoption.apply_adoption_plan(context, plan)
    assert adopted == ("ingest",)
    target = StageManifest.read_json(context.stage_manifest_path("ingest"))
    assert target.output.content_sha256 == source.output.content_sha256
    source_file = Path(source.output.donors[0].files[0].path)
    target_file = context.output_root("ingest") / "donor_id=1" / "data.parquet"
    assert source_file.samefile(target_file)
    assert (context.run_root / "ingest_rejects.parquet").samefile(
        Path(source.sidecars[0].path)
    )
    receipt = context.config.manifests_dir / "adoption" / "ingest.json"
    assert receipt.resolve() in {Path(item.path).resolve() for item in target.sidecars}

    manifest_id = target.content_id
    adopted_again = adoption.apply_adoption_plan(context, plan)
    assert adopted_again == ("ingest",)
    assert StageManifest.read_json(
        context.stage_manifest_path("ingest")
    ).content_id == manifest_id


def test_apply_rejects_source_manifest_swap_before_writing(tmp_path, monkeypatch):
    context, plan, _ = _tiny_adoption(tmp_path)
    monkeypatch.setattr(adoption, "build_adoption_plan", lambda *args, **kwargs: plan)
    source_path = (
        context.base_config.runs_dir
        / plan["source_run_id"]
        / "manifests"
        / "ingest.json"
    )
    source_path.write_text(
        source_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )

    with pytest.raises(StaleArtifactError, match="differs from the reviewed plan"):
        adoption.apply_adoption_plan(context, plan)
    assert not context.run_root.exists()
