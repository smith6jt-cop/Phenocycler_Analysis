"""Unit tests for config loading, overrides, and derived paths (5-compartment pipeline)."""

from __future__ import annotations

from phenocycler import load_config
from phenocycler.cohort import DONOR_EXCLUSIONS, ensure_eligible_donors
from phenocycler.config import (COMPARTMENT_ORDER, COMPARTMENT_GATES, OTHER_LABEL,
                                ENDOCRINE_MARKERS, MARKER_PAIRS, FUNCTIONAL_PAIRS)


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
    assert {"CD3e", "CD20", "CD68", "MPO"}.issubset(set(COMPARTMENT_GATES["Immune"]))
    assert ENDOCRINE_MARKERS == ["INS", "GCG", "SST"]      # IAPP removed 2026-07-10 (failed marker)


def test_marker_pairs_directional_no_target_collision():
    """RESTORE keys threshs on the TARGET only -> two pairs sharing a target silently overwrite."""
    for pairs in (MARKER_PAIRS, FUNCTIONAL_PAIRS):
        targets = [p[0] for p in pairs]
        assert len(targets) == len(set(targets)), f"duplicate target in {pairs}"
        assert all(len(p) == 2 for p in pairs)
    pm = {t: r for t, r in MARKER_PAIRS}
    assert pm["E_cadherin"] == "CD31"          # epithelial <- endothelial (spatially separated)
    assert pm["B3TUBB"] == "CD3e"              # neural <- immune (TUBB3-negative)
    assert ["INS", "GCG"] in MARKER_PAIRS and ["GCG", "INS"] in MARKER_PAIRS   # intra-islet reciprocal


def test_config_matches_science_modules():
    """Guard the module-level constants against drift."""
    from phenocycler import lineage, marker_taxonomy
    assert list(lineage.CLASSES) == COMPARTMENT_ORDER + [OTHER_LABEL]
    assert not hasattr(marker_taxonomy, "CD99_BRIGHT")            # CD99 demoted to a state marker
    for gates in COMPARTMENT_GATES.values():                     # gate markers must be TYPE (heatmap)
        for m in gates:
            assert m in marker_taxonomy.TYPE, f"gate marker {m} missing from marker_taxonomy.TYPE"
