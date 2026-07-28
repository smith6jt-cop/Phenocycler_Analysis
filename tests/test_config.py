"""Unit tests for config loading, overrides, and derived paths (8-class pipeline)."""

from __future__ import annotations


from phenocycler import load_config
from phenocycler.config import (LINEAGES, STRUCT_LINEAGES,
                                DEFAULT_MARKER_PAIRS, EXTRA_MARKER_PAIRS, EXTRA_MARKERS)


def test_defaults_and_derived_paths():
    cfg = load_config()
    assert cfg.data_dir.is_absolute()
    assert cfg.cells_dir == cfg.data_dir / "cells"
    assert cfg.cells_redsea_dir == cfg.data_dir / "cells_redsea"
    assert cfg.restore_gated_dir == cfg.data_dir / "restore_gated_redsea"
    assert cfg.restore_thresholds_csv == cfg.data_dir / "restore_thresholds_redsea.csv"
    assert cfg.broad_dir == cfg.data_dir / "phenotype" / "broad"
    assert cfg.qupath_class_dir == cfg.data_dir / "phenotype" / "qupath_class"
    # 8-class additions
    assert cfg.restore_gated_extra_dir == cfg.data_dir / "restore_gated_redsea_extra"
    assert cfg.restore_redsea_extra_dir == cfg.data_dir / "restore_redsea_extra"
    assert cfg.restore_thresholds_extra_csv == cfg.data_dir / "restore_thresholds_extra.csv"
    assert cfg.restore_gated_prefloor_dir == cfg.data_dir / "restore_gated_redsea.pre_hormonefloor"
    assert cfg.redsea_reassess_dir == cfg.data_dir / "redsea_reassess"


def test_scientific_defaults():
    cfg = load_config()
    assert cfg.redsea_comp_mode == 0
    assert cfg.redsea_alpha == 1.0
    assert cfg.redsea_edge_radius == 0
    assert cfg.restore_model == "SSC"
    assert cfg.restore_robust is True
    assert cfg.restore_robust_factor == 3.0
    assert cfg.restore_min_cell_area == 5.0
    # lineage knobs (the false-endocrine floor + CD99 bright gate)
    assert cfg.hormone_min_norm == 5.0
    assert cfg.cd99_bright == 3.0


def test_keyword_override():
    cfg = load_config(n_jobs=8, redsea_alpha=0.5)
    assert cfg.n_jobs == 8
    assert cfg.redsea_alpha == 0.5


def test_unknown_override_raises():
    import pytest
    with pytest.raises(TypeError):
        load_config(not_a_field=1)


def test_relative_paths_resolve_against_config_dir(tmp_path):
    ini = tmp_path / "config.ini"
    ini.write_text("[paths]\ndata_dir = mydata\nimages_dir = imgs\n")
    cfg = load_config(ini)
    assert cfg.data_dir == (tmp_path / "mydata").resolve()
    assert cfg.images_dir == (tmp_path / "imgs").resolve()


def test_absolute_paths_preserved(tmp_path):
    ini = tmp_path / "config.ini"
    abs_imgs = tmp_path / "abs" / "images"
    ini.write_text(f"[paths]\nimages_dir = {abs_imgs}\n")
    cfg = load_config(ini)
    assert cfg.images_dir == abs_imgs


def test_lineage_ini_section(tmp_path):
    ini = tmp_path / "config.ini"
    ini.write_text("[lineage]\nhormone_min_norm = 4\ncd99_bright = 2.5\n")
    cfg = load_config(ini)
    assert cfg.hormone_min_norm == 4.0
    assert cfg.cd99_bright == 2.5


def test_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("PHENOCYCLER_JOBS", "4")
    monkeypatch.setenv("PHENOCYCLER_DATA_DIR", str(tmp_path / "envdata"))
    monkeypatch.setenv("PHENOCYCLER_HORMONE_MIN_NORM", "6")
    cfg = load_config()
    assert cfg.n_jobs == 4
    assert cfg.data_dir == tmp_path / "envdata"
    assert cfg.hormone_min_norm == 6.0


def test_discover_donors(tmp_path):
    cfg = load_config(data_dir=tmp_path)
    for d in ("6539", "6450", "6414"):
        (tmp_path / "cells" / f"donor_id={d}").mkdir(parents=True)
    assert cfg.discover_donors() == ["6414", "6450", "6539"]   # sorted


def test_constants_shape():
    assert len(LINEAGES) == 8
    assert {"Neural", "Neutrophil"}.issubset(LINEAGES)
    assert "CD99" in LINEAGES["Endocrine"]
    assert LINEAGES["Neural"] == ["B3TUBB"] and LINEAGES["Neutrophil"] == ["MPO"]
    assert set(STRUCT_LINEAGES) == {"Epithelial", "Fibroblast", "Muscle"}
    assert len(DEFAULT_MARKER_PAIRS) == 10          # the validated 10-marker gates stay unchanged
    assert len(EXTRA_MARKER_PAIRS) == 3
    assert set(EXTRA_MARKERS) == {"CD99", "B3TUBB", "MPO"}


def test_config_matches_science_modules():
    """Guard the two independent definitions against 6-vs-8 / CD99 drift."""
    from phenocycler import lineage, marker_taxonomy
    assert list(lineage.LNAMES) == list(LINEAGES)                 # lineage uses config's LINEAGES
    assert load_config().cd99_bright == marker_taxonomy.CD99_BRIGHT
