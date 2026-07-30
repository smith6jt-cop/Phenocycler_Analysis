from __future__ import annotations

import pytest

from phenocycler import load_config
from phenocycler.cohort import DONOR_EXCLUSIONS, ensure_eligible_donors


def test_defaults_expose_only_the_current_stage_roots():
    cfg = load_config()
    assert cfg.data_dir.is_absolute()
    assert cfg.cells_dir == cfg.data_dir / "cells"
    assert cfg.geometry_qc_dir == cfg.data_dir / "geometry_qc"
    assert cfg.redsea_dir == cfg.data_dir / "redsea"
    assert cfg.selected_expression_dir == cfg.data_dir / "expression"
    assert cfg.marker_evidence_dir == cfg.data_dir / "marker_evidence"
    assert cfg.assignments_dir == cfg.data_dir / "assignments"
    for legacy in (
        "cells_redsea_dir",
        "restore_dir",
        "restore_gated_dir",
        "broad_dir",
        "restore_pair_reviews_csv",
    ):
        assert not hasattr(cfg, legacy)


def test_fixed_scientific_defaults():
    cfg = load_config()
    assert cfg.cells_min_cell_area == 5.0
    assert cfg.cell_qc_min_cell_area == 20.0
    assert cfg.redsea_comp_mode == 1
    assert cfg.redsea_norm_form == "donor"
    assert cfg.redsea_alpha == 1.0


def test_keyword_and_unknown_overrides():
    assert load_config(n_jobs=8, redsea_alpha=0.5).n_jobs == 8
    with pytest.raises(TypeError):
        load_config(not_a_field=1)


def test_relative_paths_resolve_against_config_file(tmp_path):
    ini = tmp_path / "config.ini"
    ini.write_text(
        "[paths]\n"
        "data_dir = output\n"
        "qupath_manifest = inputs/qupath.json\n"
        "marker_registry = policy/markers.json\n"
        "typing_rules = policy/types.json\n"
    )
    cfg = load_config(ini)
    assert cfg.data_dir == (tmp_path / "output").resolve()
    assert cfg.qupath_manifest == (tmp_path / "inputs/qupath.json").resolve()
    assert cfg.marker_registry == (tmp_path / "policy/markers.json").resolve()
    assert cfg.typing_rules == (tmp_path / "policy/types.json").resolve()


def test_environment_overrides(monkeypatch, tmp_path):
    monkeypatch.setenv("PHENOCYCLER_JOBS", "4")
    monkeypatch.setenv("PHENOCYCLER_DATA_DIR", str(tmp_path / "data"))
    cfg = load_config()
    assert cfg.n_jobs == 4
    assert cfg.data_dir == tmp_path / "data"


def test_discover_donors_keeps_central_fail_closed_exclusions(tmp_path):
    cfg = load_config(data_dir=tmp_path)
    for donor in ("6579", "6539", "6457", "6450"):
        (cfg.cells_dir / f"donor_id={donor}").mkdir(parents=True)
    assert cfg.discover_donors() == ["6450", "6539"]
    for donor in DONOR_EXCLUSIONS:
        with pytest.raises(ValueError, match=rf"excluded donor.*{donor}"):
            ensure_eligible_donors(["6539", donor], context="test")


def test_stale_redsea_modes_are_rejected():
    with pytest.raises(ValueError, match="mass-conserving"):
        load_config(redsea_norm_form="recipient")
    with pytest.raises(TypeError, match="unknown"):
        load_config(redsea_emit_compartments=False)
