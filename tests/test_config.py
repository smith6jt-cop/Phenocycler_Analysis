"""Unit tests for config loading, overrides, and derived paths (5-compartment pipeline)."""

from __future__ import annotations

from phenocycler import load_config
from phenocycler.cohort import DONOR_EXCLUSIONS, ensure_eligible_donors
from phenocycler.config import (COMPARTMENT_ORDER, COMPARTMENT_GATES, OTHER_LABEL,
                                UNRESOLVED_LABEL, ENDOCRINE_MARKERS)


def test_defaults_and_derived_paths():
    cfg = load_config()
    assert cfg.data_dir.is_absolute()
    assert cfg.cells_dir == cfg.data_dir / "cells"
    assert cfg.cells_redsea_dir == cfg.data_dir / "cells_redsea"
    assert cfg.restore_gated_dir == cfg.data_dir / "restore_gated_redsea"
    assert cfg.restore_thresholds_csv == cfg.data_dir / "restore_thresholds_redsea.csv"
    assert cfg.broad_dir == cfg.data_dir / "phenotype" / "broad"
    assert cfg.qupath_class_dir == cfg.data_dir / "phenotype" / "qupath_class"


def test_restore_apply_input_paths():
    """Gate-4 restore_apply inputs are derived under data_dir/restore_pair_validation."""
    cfg = load_config()
    pv = cfg.data_dir / "restore_pair_validation"
    assert cfg.restore_pair_validation_dir == pv
    assert cfg.restore_allcells_sample_dir == pv / "qupath_compartment_full_20donor"
    assert cfg.restore_pair_reviews_csv == pv / "pair_reviews_accepted.csv"
    assert cfg.restore_expression_reviews_csv == pv / "expression_reviews.csv"
    assert cfg.restore_apply_manifest_csv == cfg.restore_gated_dir / "apply_manifest.csv"


def test_scientific_defaults():
    cfg = load_config()
    assert cfg.redsea_comp_mode == 1          # subtract + reinforce (METHOD v10, 2026-07-25)
    assert cfg.redsea_norm_form == "donor"    # published mass-conserving contact normalisation
    assert cfg.redsea_exclude_no_spillover is True
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


def test_restore_pair_method_policy_is_configurable():
    """The METHOD v10 knobs must be reachable from config.ini, or the arm under review and the arm
    production runs cannot be made to agree (they could not, before 2026-07-27)."""
    from phenocycler.restore_validation import (CONTROL_DEFINITIONS, DIVISOR_STATISTICS,
                                                PairValidationConfig)
    cfg = load_config()
    # UNSET by default, deliberately. These used to mirror the frozen dataclass "so exposing them
    # changes nothing on its own" -- but PairValidationConfig defaults to divisor_statistic="max",
    # which restore_apply cannot run, while the frozen production p99 lives only in the PARENT repo's
    # config.ini. Mirroring therefore made "forgot --config" silently select the unrunnable policy.
    # See test_restore_apply.test_unset_policy_refuses_rather_than_defaulting.
    assert cfg.restore_control_definition == ""
    assert cfg.restore_divisor_statistic == ""
    assert PairValidationConfig.control_definition in CONTROL_DEFINITIONS
    assert PairValidationConfig.divisor_statistic in DIVISOR_STATISTICS
    over = load_config(restore_control_definition="reference_and_target", restore_divisor_statistic="p99")
    assert (over.restore_control_definition, over.restore_divisor_statistic) == (
        "reference_and_target", "p99")


def test_restore_policy_reaches_config_ini_schema():
    """The keys must be in _INI_SCHEMA['restore'] — a dataclass field alone would be silently ignored
    by config.ini."""
    from phenocycler.config import _INI_SCHEMA
    assert _INI_SCHEMA["restore"]["control_definition"] == ("restore_control_definition", str)
    assert _INI_SCHEMA["restore"]["divisor_statistic"] == ("restore_divisor_statistic", str)


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
    # Gate 4 (2026-07-24): Immune is FIRST, and every gate marker is a divisor-backed ACCEPTED_PAIRS target.
    # Endothelial precedes Epithelial as of 2026-07-27 (see config.py for the measurement). Immune
    # stays first (CD31 is on leukocytes); Neural still follows Epithelial (endocrine is TUBB3+).
    assert COMPARTMENT_ORDER == ["Immune", "Endothelial", "Epithelial", "Neural", "Mesenchymal"]
    assert COMPARTMENT_ORDER.index("Immune") < COMPARTMENT_ORDER.index("Endothelial")
    assert COMPARTMENT_ORDER.index("Epithelial") < COMPARTMENT_ORDER.index("Neural")
    assert OTHER_LABEL == "Other"
    assert UNRESOLVED_LABEL == "Unresolved"
    # Expanded 2026-07-27 from 7 gate markers to 24 (maintainer sign-off). `Other` is the ordered
    # residual's terminal bucket, so a cell whose only lineage marker was ungated could not be typed.
    assert COMPARTMENT_GATES["Immune"] == ["CD3e", "CD68", "CD11b", "CD79a", "CD11c", "Iba1",
                                           "CD206", "MPO", "CD4", "CD8"]
    assert COMPARTMENT_GATES["Epithelial"] == ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM",
                                               "Keratin_5", "INS", "GCG", "SST"]
    assert COMPARTMENT_GATES["Endothelial"] == ["CD31", "CD34", "Podoplanin"]
    assert COMPARTMENT_GATES["Neural"] == ["B3TUBB"]
    assert COMPARTMENT_GATES["Mesenchymal"] == ["Vimentin", "SMA"]
    # Still NOT gates, each for a stated reason rather than by omission:
    #   CD20   — takes no part in RESTORE at all (RESTORE_EXCLUDED_MARKERS)
    #   CD163  — screened 2026-07-27 and REJECTED: arm separation 1.73-2.25x, at or below the 2x floor
    #   CD56   — marks NK, neural AND endocrine cells, so it cannot discriminate three compartments
    #   CD66   — epithelial and granulocyte
    for rejected in ("CD20", "CD163", "CD56", "CD66"):
        assert rejected not in COMPARTMENT_GATES["Immune"]
        assert rejected not in COMPARTMENT_GATES["Epithelial"]
    assert ENDOCRINE_MARKERS == ["INS", "GCG", "SST"]      # retained for the deferred level-2 pass


def test_gate_markers_are_the_accepted_restore_targets():
    """Every compartment gate marker must be a divisor-backed ACCEPTED_PAIRS target (Gate 4), so a
    tri-state 'unavailable' always means a real per-donor QC gap, never 'never screened'."""
    from phenocycler.restore_validation import ACCEPTED_PAIRS
    gate_markers = {m for gates in COMPARTMENT_GATES.values() for m in gates}
    accepted_targets = {t for t, _ in ACCEPTED_PAIRS}
    assert gate_markers == accepted_targets


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
    assert list(lineage.CLASSES) == COMPARTMENT_ORDER + [OTHER_LABEL, UNRESOLVED_LABEL]
    assert not hasattr(marker_taxonomy, "CD99_BRIGHT")            # CD99 demoted to a state marker
    for gates in COMPARTMENT_GATES.values():                     # gate markers must be TYPE (heatmap)
        for m in gates:
            assert m in marker_taxonomy.TYPE, f"gate marker {m} missing from marker_taxonomy.TYPE"
