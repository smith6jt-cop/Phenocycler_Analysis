"""Tissue as a first-class dataset feature: derivation, scoping, and per-tissue divergence.

The cohort is 20 donors x {PhenoCycler, Xenium} x {pancreas, pancreatic lymph node}. Every
test here pins one of the two ways those tissues differ — the structural unit and the
identity panels — or one of the places a tissue-blind assumption would silently lose half
the dataset.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration import tissues as T
from phenocycler.integration.contract import (CellTable, ContractError, coerce_cell_table,
                                              validate_cell_table)
from phenocycler.integration.donor import COMPARABLE, comparable_for, composition, disease_trend
from phenocycler.integration.pipeline import NEEDS_STRUCTURES, applicable_tissues
from phenocycler.integration.qc import FAIL, PASS, WARN, verdict_for
from phenocycler.integration.structures import extract_donor


@pytest.fixture(scope="module")
def cfg():
    return load_config()


# --------------------------------------------------------------------------- #
# Derivation
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("roi,expected", [
    ("panc", T.PANCREAS), ("Panc", T.PANCREAS), ("pancreas", T.PANCREAS),
    ("pln_1", T.LYMPH_NODE), ("pln_3", T.LYMPH_NODE), ("pLN", T.LYMPH_NODE),
    ("ln_2", T.LYMPH_NODE),
    ("region_1", ""), ("", ""),
])
def test_tissue_derived_from_roi(roi, expected):
    """ROI is what the data carries; tissue is what the analysis needs. One derives the other.

    An unrecognised ROI returns '' rather than guessing a tissue — a wrong guess would route
    a section through the wrong identity panels and produce confident nonsense.
    """
    assert T.for_roi(roi) == expected


@pytest.mark.parametrize("given", ["pancreatic_lymph_node", "pLN", "lymph node", "LYMPH-NODE"])
def test_normalize_is_lenient_on_input_strict_on_output(given):
    assert T.normalize(given) == T.LYMPH_NODE


def test_normalize_rejects_an_unknown_tissue():
    """A typo must not become a third tissue that silently matches nothing downstream."""
    with pytest.raises(T.TissueError, match="unknown tissue"):
        T.normalize("thymus")
    with pytest.raises(T.TissueError, match="empty"):
        T.normalize("  ")


def test_cell_table_carries_tissue_from_creation(cfg):
    """`tissue` is set when the dataset is imported, not inferred later by each consumer."""
    raw = pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(5)],
        "x_um": np.linspace(0, 100, 5), "y_um": np.linspace(0, 100, 5),
        "lineage_native": ["T"] * 5,
    })
    panc = coerce_cell_table(raw, modality="xenium", donor_id="1", roi="panc")
    pln = coerce_cell_table(raw, modality="xenium", donor_id="1", roi="pln_2")
    assert (panc["tissue"] == T.PANCREAS).all()
    assert (pln["tissue"] == T.LYMPH_NODE).all()

    # An explicit tissue wins over the ROI, for data whose ROI naming does not follow suit.
    forced = coerce_cell_table(raw, modality="xenium", donor_id="1", roi="region_1",
                               tissue="pLN")
    assert (forced["tissue"] == T.LYMPH_NODE).all()


def test_untissued_table_is_rejected():
    """An ROI that names no tissue must fail validation rather than default to pancreas.

    Defaulting is the dangerous option: a mislabelled lymph node would be scored against
    pancreatic endocrine and exocrine identity cores, which is not a degraded answer but a
    wrong one. `coerce` fills and `validate` rejects — the same split the pixel guard uses —
    so the failure surfaces before anything is written, since `write_cell_table` validates.
    """
    raw = pd.DataFrame({
        "cell_id": ["a", "b"], "x_um": [0.0, 1.0], "y_um": [0.0, 1.0],
        "lineage_native": ["T", "B"],
    })
    df = coerce_cell_table(raw, modality="xenium", donor_id="1", roi="region_1")
    assert (df["tissue"] == "").all(), "an unmappable ROI must not be guessed at"
    with pytest.raises(ContractError, match="no tissue"):
        validate_cell_table(df)


def test_an_unknown_explicit_tissue_is_refused():
    """`tissue='thymus'` is a typo or a scope error, never a silent third tissue."""
    raw = pd.DataFrame({
        "cell_id": ["a", "b"], "x_um": [0.0, 1.0], "y_um": [0.0, 1.0],
        "lineage_native": ["T", "B"],
    })
    with pytest.raises(T.TissueError, match="unknown tissue"):
        coerce_cell_table(raw, modality="xenium", donor_id="1", roi="panc", tissue="thymus")


# --------------------------------------------------------------------------- #
# Scoping: separately, and combined
# --------------------------------------------------------------------------- #

def test_resolve_rois_covers_the_three_run_shapes(cfg):
    """The user's requirement verbatim: usable separately and when combined."""
    assert cfg.resolve_rois(None, ["pancreas"]) == ["panc"]
    assert cfg.resolve_rois(None, ["pancreatic_lymph_node"]) == ["pln_1", "pln_2", "pln_3"]
    # Combined is the default, so a bare run does not quietly leave the lymph nodes out.
    assert cfg.resolve_rois(None, None) == ["panc", "pln_1", "pln_2", "pln_3"]
    # A single ROI is the escape hatch and beats the tissue selection.
    assert cfg.resolve_rois("pln_2", ["pancreas"]) == ["pln_2"]


def test_tissue_list_and_roi_mapping_round_trip(cfg):
    for tissue in cfg.tissue_list:
        for roi in cfg.rois_for_tissue(tissue):
            assert cfg.tissue_for_roi(roi) == tissue


# --------------------------------------------------------------------------- #
# Divergence 1: the structural unit
# --------------------------------------------------------------------------- #

def test_only_pancreas_has_a_structural_unit():
    assert T.structure_kinds(T.PANCREAS) == ("islet", "duct", "vessel")
    assert T.structure_kinds(T.LYMPH_NODE) == ()
    assert T.has_structures(T.PANCREAS)
    assert not T.has_structures(T.LYMPH_NODE)


def test_structure_stages_are_skipped_per_tissue_not_per_run():
    """A combined run must still match islets in the pancreas half.

    Gating the whole stage off because *one* selected tissue lacks islets would throw away
    the pancreas result, which is the opposite of what skipping pLN is for.
    """
    both = [T.PANCREAS, T.LYMPH_NODE]
    for stage in NEEDS_STRUCTURES:
        assert applicable_tissues(stage, both) == [T.PANCREAS]
        assert applicable_tissues(stage, [T.LYMPH_NODE]) == []
    # Every other stage applies to both tissues untouched.
    assert applicable_tissues("register", both) == both
    assert applicable_tissues("grid", both) == both
    assert applicable_tissues("donor", both) == both


def test_extract_donor_returns_empty_for_a_tissue_without_structures(cfg):
    """No islets is not an error and not zero islets found — the level does not apply."""
    n = 200
    rng = np.random.default_rng(0)
    raw = pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(n)],
        "x_um": rng.uniform(0, 2000, n), "y_um": rng.uniform(0, 2000, n),
        "lineage_native": ["T_NK"] * n, "lineage_common": ["T_NK"] * n,
    })
    df = coerce_cell_table(raw, modality="xenium", donor_id="1", roi="pln_1")
    table = CellTable(df=df, modality="xenium", donor_id="1", roi="pln_1",
                      tissue=T.LYMPH_NODE)
    frames, labels = extract_donor(table, cfg)
    assert frames == {}
    assert len(labels) == len(df)
    assert (labels == "").all(), "no cell may be assigned to a structure that cannot exist"


def test_match_skips_a_tissue_with_no_structural_unit(cfg, tmp_path):
    """`match_donor` reports a skip, not an error: nothing is wrong with a lymph node."""
    from phenocycler.integration.match import match_donor

    pairs, stats = match_donor(cfg, "6539", "pln_1")
    assert pairs.empty
    assert stats.get("skipped") is True
    assert "no structural unit" in stats["reason"]
    # The pancreas path is untouched and still reports the ordinary missing-input error.
    _pairs, panc_stats = match_donor(cfg, "6539", "panc")
    assert not panc_stats.get("skipped")


# --------------------------------------------------------------------------- #
# Divergence 2: identity panels and comparable lineages
# --------------------------------------------------------------------------- #

def test_pancreas_specific_panels_do_not_apply_to_lymph_node():
    """Scoring a lymph-node cell on an endocrine axis is worse than uninformative.

    The broad call is an argmax over identity-core scores, so an axis that cannot exist in
    the tissue is not a harmless zero — it is a competitor that can win.
    """
    for panel in T.PANCREAS_ONLY_PANELS:
        assert T.applies(T.PANCREAS, panel)
        assert not T.applies(T.LYMPH_NODE, panel)
    for shared in ("03_t_cell", "05_myeloid", "09_endothelial", "12_neural"):
        assert T.applies(T.PANCREAS, shared)
        assert T.applies(T.LYMPH_NODE, shared)


def test_lymph_node_borrows_the_shared_core_pool():
    """XeniumPanelExplorer ships pancreas and thymus, no lymph node — shared cores are the
    correct answer for pLN rather than a degraded fallback."""
    assert T.panel_tissue_dir(T.PANCREAS) == "pancreas"
    assert T.panel_tissue_dir(T.LYMPH_NODE) is None


def test_comparable_lineages_exclude_compartments_the_tissue_lacks():
    """A near-zero endocrine fraction on both sides is two instruments correctly finding
    nothing — not agreement about biology. Letting it into the Spearman would inflate the
    concordance with a free correct answer."""
    panc = comparable_for(T.PANCREAS)
    pln = comparable_for(T.LYMPH_NODE)
    assert "Endocrine" in panc and "Endocrine" not in pln
    assert set(pln) < set(panc)
    assert {"T_NK", "B_Plasma", "Myeloid", "Vascular", "Stromal"} <= set(pln)
    # An unknown or absent tissue falls back to the instrument-only baseline, never to empty.
    assert comparable_for("") == COMPARABLE


# --------------------------------------------------------------------------- #
# QC gating
# --------------------------------------------------------------------------- #

def test_lymph_node_registration_is_graded_on_dice_alone():
    """The islet-RMSE gate is undefined without islets, so it is dropped, not failed.

    Counting it as a missed gate would cap every lymph-node donor at WARN forever and
    permanently withhold `crossmodal` from half the cohort — an artefact of the gate design,
    not a property of the registration.
    """
    good = dict(dice_min=0.80, rmse_max=50.0)
    v, reasons = verdict_for(0.91, np.nan, **good, require_structures=False)
    assert v == PASS
    assert any("n/a" in r for r in reasons), "the dropped gate must be stated, not hidden"

    v, _ = verdict_for(0.42, np.nan, **good, require_structures=False)
    assert v == FAIL

    # Same inputs in a tissue that does have islets stay at WARN — the gate is real there.
    assert verdict_for(0.91, np.nan, **good)[0] == WARN


def test_dice_gate_still_binds_for_both_tissues():
    good = dict(dice_min=0.80, rmse_max=50.0)
    assert verdict_for(0.91, 20.0, **good)[0] == PASS
    assert verdict_for(0.50, 20.0, **good)[0] == WARN        # rmse rescues it to WARN
    assert verdict_for(0.50, 900.0, **good)[0] == FAIL


# --------------------------------------------------------------------------- #
# Donor-level: two tissues per donor must not collapse into one
# --------------------------------------------------------------------------- #

def _table(roi, modality="xenium", lineage="T_NK", n=50):
    rng = np.random.default_rng(1)
    raw = pd.DataFrame({
        "cell_id": [f"{roi}_{i}" for i in range(n)],
        "x_um": rng.uniform(0, 500, n), "y_um": rng.uniform(0, 500, n),
        "lineage_native": [lineage] * n, "lineage_common": [lineage] * n,
    })
    df = coerce_cell_table(raw, modality=modality, donor_id="6539", roi=roi)
    return CellTable(df=df, modality=modality, donor_id="6539", roi=roi,
                     tissue=T.for_roi(roi))


def test_composition_is_reported_per_section_with_its_tissue():
    panc = composition(_table("panc"))
    pln = composition(_table("pln_1"))
    assert (panc["tissue"] == T.PANCREAS).all()
    assert (pln["tissue"] == T.LYMPH_NODE).all()
    assert (panc["donor_id"] == pln["donor_id"]).all(), "same donor, two sections"


def test_disease_trend_groups_by_tissue():
    """Pooling a pancreas and a lymph node into one mean would halve the pancreatic signal
    with structural zeros and turn the ND -> AAB -> T1D gradient into noise."""
    comp = pd.concat([composition(_table("panc")), composition(_table("pln_1"))],
                     ignore_index=True)
    trend = disease_trend(comp, {"6539": "ND"})
    assert "tissue" in trend.columns
    assert set(trend["tissue"]) == {T.PANCREAS, T.LYMPH_NODE}
    # Each (tissue, lineage) keeps its own row rather than being averaged together.
    t_nk = trend[trend["lineage_common"] == "T_NK"]
    assert len(t_nk) == 2


def test_run_donor_keeps_both_sections_of_a_donor(cfg, tmp_path):
    """The regression this guards: `{donor: roi for donor, roi in partitions}` keeps only the
    last ROI, so a donor with a pancreas and a lymph node would silently lose one tissue."""
    from phenocycler.integration.contract import write_cell_table
    from phenocycler.integration.donor import run_donor

    data = tmp_path / "data"
    local = load_config()
    local.data_dir = data
    for roi in ("panc", "pln_1"):
        for modality, base in (("phenocycler", local.cells_pheno_dir),
                               ("xenium", local.cells_xen_dir)):
            t = _table(roi, modality=modality)
            write_cell_table(t.df, base / "donor_id=6539" / f"roi={roi}" / "cells.parquet")

    stats = run_donor(local, ["6539"])
    assert set(stats["roi"]) == {"panc", "pln_1"}, "both sections must be compared"
    assert set(stats["tissue"]) == {T.PANCREAS, T.LYMPH_NODE}


# --------------------------------------------------------------------------- #
# PhenoCycler source resolution
# --------------------------------------------------------------------------- #

def test_pheno_source_override_is_still_available(cfg, tmp_path):
    """`pheno_tissue_dirs` survives as an override for a tissue processed through a separate
    core-pipeline run (a different QuPath export). It is no longer the normal path: section
    keys separate the tissues within one run."""
    local = load_config()
    local.pheno_tissue_dirs = (f"pancreas={tmp_path/'panc'},"
                               f"pancreatic_lymph_node={tmp_path/'pln'}")
    assert local.pheno_dir_for_tissue(T.PANCREAS) == tmp_path / "panc"
    assert local.pheno_dir_for_tissue(T.LYMPH_NODE) == tmp_path / "pln"
    # Unmapped tissues fall back to the shared data dir, which is the default and is correct.
    local.pheno_tissue_dirs = f"pancreas={tmp_path/'panc'}"
    assert local.pheno_dir_for_tissue(T.LYMPH_NODE) == local.data_dir


def test_one_source_dir_yields_both_tissues(cfg, tmp_path):
    """The whole point of the section key: one core-pipeline run covers both tissues.

    Previously `6539` and `6539pLN` collided in `donor_id=6539`, and the workaround was a
    separate data dir per tissue. Now the partitions are distinct, so a single directory
    exports a pancreas section *and* a lymph-node section for the same donor — each read
    from its own partition and written to its own `(donor, roi)`.
    """
    from phenocycler.integration.export_pheno import _sections_in

    local = load_config()
    local.data_dir = tmp_path
    for key in ("6539", "6539pln", "6414"):
        (local.broad_dir / f"donor_id={key}").mkdir(parents=True, exist_ok=True)

    both = _sections_in(local, wanted_rois={"panc", "pln_1"}, keep_donors=None)
    assert sorted(both) == ["6414", "6539", "6539pln"]

    # And each tissue can still be selected on its own, from the same directory.
    assert sorted(_sections_in(local, {"panc"}, None)) == ["6414", "6539"]
    assert sorted(_sections_in(local, {"pln_1"}, None)) == ["6539pln"]
    assert sorted(_sections_in(local, {"panc", "pln_1"}, {"6539"})) == ["6539", "6539pln"]


def test_unparseable_partition_is_skipped_not_guessed(cfg, tmp_path, capsys):
    """An unknown region token must not be exported as a pancreas.

    Folding `6539spleen` into the pancreas would put a third tissue's cells into pancreatic
    composition and islet statistics with no error — the failure the section key exists to
    prevent, reintroduced one level up.
    """
    from phenocycler.integration.export_pheno import _sections_in

    local = load_config()
    local.data_dir = tmp_path
    for key in ("6539", "6539spleen"):
        (local.broad_dir / f"donor_id={key}").mkdir(parents=True, exist_ok=True)

    found = _sections_in(local, wanted_rois={"panc"}, keep_donors=None)
    assert found == ["6539"], "the unknown region must not be exported at all"
    assert "6539spleen" in capsys.readouterr().out, "and the skip must be announced"


# --------------------------------------------------------------------------- #
# Stage scoping
# --------------------------------------------------------------------------- #

def test_stages_span_the_selection_rather_than_being_called_per_roi():
    """Stage outputs are stage-level paths, so a per-ROI call would clobber the last ROI's.

    Each per-ROI runner therefore takes `roi`/`tissue` and iterates the manifest itself. This
    pins the signature so a future refactor cannot reintroduce the orchestrator-side loop.
    """
    import inspect

    from phenocycler.integration.crossmodal import run_crossmodal
    from phenocycler.integration.figures import run_figures
    from phenocycler.integration.grid import run_grid
    from phenocycler.integration.match import run_match
    from phenocycler.integration.qc import run_qc
    from phenocycler.integration.register import run_register
    from phenocycler.integration.sameslide import run_sameslide

    for fn in (run_grid, run_match, run_register, run_crossmodal, run_sameslide, run_figures,
               run_qc):
        params = inspect.signature(fn).parameters
        assert "tissue" in params, f"{fn.__name__} cannot be scoped to a tissue"
        assert params["roi"].default is None, (
            f"{fn.__name__} defaults to a single ROI; it must default to the whole selection")


def test_paired_rows_is_the_single_scoping_gate():
    """Every per-ROI stage must scope through `paired_rows`, not a hand-rolled roi equality.

    The old `man["roi"] == roi` pattern is what made the lymph node invisible; a new stage
    copying it would silently reintroduce that.
    """
    import pathlib
    import re

    base = pathlib.Path("phenocycler/integration")
    offenders = []
    for path in base.glob("*.py"):
        text = path.read_text()
        if re.search(r'man\["roi"\]\s*==\s*roi', text):
            offenders.append(path.name)
    assert not offenders, f"hand-rolled roi filter in {offenders}; use paired_rows(...)"


# --------------------------------------------------------------------------- #
# Registry integrity
# --------------------------------------------------------------------------- #

def test_every_configured_tissue_is_fully_described(cfg):
    """A tissue in config.ini with no registry entry would fail late and obscurely."""
    for tissue in cfg.tissue_list:
        assert tissue in T.TISSUES
        assert T.DISPLAY_NAME[tissue]
        assert cfg.rois_for_tissue(tissue), f"{tissue} has no ROIs"
        assert T.comparable_lineages(tissue), f"{tissue} has nothing to compare"
        assert T.describe(tissue)


def test_structure_seeds_cover_every_declared_kind():
    """Adding a structure kind without its seed definition would fail at extraction time."""
    declared = {k for kinds in T.STRUCTURE_KINDS.values() for k in kinds}
    assert declared <= set(T.STRUCTURE_SEEDS), f"undefined seeds: {declared - set(T.STRUCTURE_SEEDS)}"
