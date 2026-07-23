"""Unit tests for config loading, overrides, and derived paths (5-compartment pipeline)."""

from __future__ import annotations

from phenocycler import load_config
from phenocycler.cohort import DONOR_EXCLUSIONS, ensure_eligible_donors
from phenocycler.config import (COMPARTMENT_ORDER, COMPARTMENT_GATES, OTHER_LABEL,
                                ENDOCRINE_MARKERS)


def test_defaults_and_derived_paths():
    cfg = load_config()
    assert cfg.data_dir.is_absolute()
    assert cfg.cells_dir == cfg.data_dir / "cells"
    assert cfg.cells_redsea_dir == cfg.data_dir / "cells_redsea"
    assert cfg.restore_gated_dir == cfg.data_dir / "restore_gated_redsea"
    assert cfg.restore_thresholds_csv == cfg.data_dir / "restore_thresholds_redsea.csv"
    assert cfg.broad_dir == cfg.data_dir / "phenotype" / "broad"
    assert cfg.qupath_class_dir == cfg.data_dir / "phenotype" / "qupath_class"


def test_scientific_defaults():
    cfg = load_config()
    assert cfg.redsea_comp_mode == 0
    assert cfg.redsea_alpha == 1.0
    assert cfg.redsea_edge_radius == 0
    assert cfg.restore_model == "SSC"
    assert cfg.restore_robust is True
    assert cfg.restore_robust_factor == 3.0
    assert cfg.restore_min_cell_area == 5.0
    # the former _norm-floor + CD99-bright knobs are removed in the curated redesign
    assert not hasattr(cfg, "hormone_min_norm")
    assert not hasattr(cfg, "immune_min_norm")
    assert not hasattr(cfg, "cd99_bright")


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


def test_env_override(monkeypatch, tmp_path):
    monkeypatch.setenv("PHENOCYCLER_JOBS", "4")
    monkeypatch.setenv("PHENOCYCLER_DATA_DIR", str(tmp_path / "envdata"))
    cfg = load_config()
    assert cfg.n_jobs == 4
    assert cfg.data_dir == tmp_path / "envdata"


def test_discover_donors(tmp_path):
    cfg = load_config(data_dir=tmp_path)
    for d in ("6579", "6539", "6457", "6450", "6414"):
        (tmp_path / "cells" / f"donor_id={d}").mkdir(parents=True)
    assert cfg.discover_donors() == ["6414", "6450", "6539"]   # sorted


def test_excluded_donors_fail_closed_when_requested_explicitly():
    import pytest

    assert "channel registration" in DONOR_EXCLUSIONS["6457"]
    assert "pancreatitis outlier" in DONOR_EXCLUSIONS["6579"]
    for donor in DONOR_EXCLUSIONS:
        with pytest.raises(ValueError, match=rf"excluded donor.*{donor}"):
            ensure_eligible_donors(["6539", donor], context="test analysis")


def test_compartment_constants_shape():
    assert COMPARTMENT_ORDER == ["Epithelial", "Endothelial", "Neural", "Immune", "Mesenchymal"]
    assert OTHER_LABEL == "Other"
    assert COMPARTMENT_GATES["Epithelial"] == ["E_cadherin"]      # E-cadherin anchor (stays + on endocrine)
    assert COMPARTMENT_GATES["Endothelial"] == ["CD31"]
    assert COMPARTMENT_GATES["Neural"] == ["B3TUBB"]
    assert COMPARTMENT_GATES["Mesenchymal"] == ["SMA", "Vimentin"]   # ordered: SMA then Vimentin
    assert {"CD3e", "CD68", "MPO"}.issubset(set(COMPARTMENT_GATES["Immune"]))
    # CD20 is deliberately NOT a compartment anchor: B cells are sparse in pancreas, so a lone CD20
    # call would more often be stray signal pulling an abundant cell into Immune than a real B cell.
    assert "CD20" not in COMPARTMENT_GATES["Immune"]
    assert "CD79a" in COMPARTMENT_GATES["Immune"]       # B cells still reach Immune via CD79a
    assert ENDOCRINE_MARKERS == ["INS", "GCG", "SST"]      # IAPP removed 2026-07-10 (failed marker)


def test_legacy_pair_webs_are_retired():
    """MARKER_PAIRS / FUNCTIONAL_PAIRS were deleted with the vendored-SSC driver: ten of their pairs
    referenced CD20, which is too sparse in pancreas to define a negative control, and the method is
    being rebuilt in restore_validation.py. Nothing may quietly resurrect them."""
    import phenocycler
    from phenocycler import config as cfg_mod
    for name in ("MARKER_PAIRS", "FUNCTIONAL_PAIRS"):
        assert not hasattr(cfg_mod, name), f"{name} was retired; do not reintroduce it"
        assert not hasattr(phenocycler, name)
    from phenocycler import restore
    for name in ("run_restore", "run_restore_functional"):
        assert not hasattr(restore, name), f"restore.{name} was retired with the pair web"
    assert hasattr(restore, "run_restore_proliferation")   # the standalone pass survives


def test_config_matches_science_modules():
    """Guard the module-level constants against drift."""
    from phenocycler import lineage, marker_taxonomy
    assert list(lineage.CLASSES) == COMPARTMENT_ORDER + [OTHER_LABEL]
    assert not hasattr(marker_taxonomy, "CD99_BRIGHT")            # CD99 demoted to a state marker
    for gates in COMPARTMENT_GATES.values():                     # gate markers must be TYPE (heatmap)
        for m in gates:
            assert m in marker_taxonomy.TYPE, f"gate marker {m} missing from marker_taxonomy.TYPE"
