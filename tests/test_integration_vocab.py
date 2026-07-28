"""Vocabulary harmonisation: totality of the crosswalk, and the off-panel marker trap."""

from __future__ import annotations

import pytest

from phenocycler import load_config
from phenocycler.config import LINEAGES
from phenocycler.integration.vocab import (
    COMMON_LINEAGES, GATE_MARKERS, PHENO_TO_COMMON,
    PROTEIN_GENE_CANDIDATES, RESOLVABLE, SURROGATE_PANEL, XENIUM_TO_COMMON, check_panel_roles_agreement, coverage_report, crosswalk_frame,
    load_panel_taxonomy, pheno_to_common, protein_gene_pairs, usable_anchors,
    xenium_to_common,
)

#: The Xenium pancreas broad vocabulary (scripts/tissue_markers/pancreas.py BROAD_TO_FINE).
XENIUM_BROAD = ("Endocrine", "Exocrine", "Stromal", "Vascular",
                "Immune_T", "Immune_B", "Myeloid", "Neural")

#: Markers verified absent from the deployed 5,081-gene panel.
KNOWN_OFF_PANEL = ("INS", "GCG", "SST", "Vimentin", "HLA_DR", "HLA_A",
                   "Iba1", "Caveolin", "M2Gal3")


@pytest.fixture(scope="module")
def cfg():
    return load_config()


@pytest.fixture(scope="module")
def tax(cfg):
    try:
        return load_panel_taxonomy(cfg)
    except Exception as exc:  # pragma: no cover - submodule not checked out
        pytest.skip(f"XeniumPanelExplorer submodule unavailable: {exc}")


# --------------------------------------------------------------------------- #
# Totality — nothing may fall through silently
# --------------------------------------------------------------------------- #

def test_every_phenocycler_lineage_maps():
    """All eight PhenoCycler classes reach the common vocabulary.

    `Immune` is intentionally absent from PHENO_TO_COMMON because it resolves through
    `immune_subclass`; every other class must map directly.
    """
    for lineage in LINEAGES:
        if lineage == "Immune":
            continue
        assert lineage in PHENO_TO_COMMON, lineage
        assert PHENO_TO_COMMON[lineage] in COMMON_LINEAGES


def test_every_xenium_lineage_maps():
    for lineage in XENIUM_BROAD:
        assert lineage in XENIUM_TO_COMMON, lineage
        assert XENIUM_TO_COMMON[lineage] in COMMON_LINEAGES


def test_immune_splits_three_ways():
    """PhenoCycler collapses T/B/myeloid; the split is recovered from the per-marker gates."""
    assert pheno_to_common("Immune", "T") == "T_NK"
    assert pheno_to_common("Immune", "B") == "B_Plasma"
    assert pheno_to_common("Immune", "Myeloid") == "Myeloid"
    # Multi-gate-positive cells go to the least-specific class rather than being resolved
    # by column order.
    assert pheno_to_common("Immune", "Mixed") == "Myeloid"
    assert pheno_to_common("Immune", "") == "Myeloid"


def test_fibroblast_and_muscle_collapse_to_stromal():
    """Xenium cannot split them: ACTA2 is not in the 5K reference and 12_pericyte
    contributes one SMC-hinted gene."""
    assert pheno_to_common("Fibroblast") == pheno_to_common("Muscle") == "Stromal"
    assert XENIUM_TO_COMMON["Stromal"] == "Stromal"


def test_xenium_reaches_granulocyte_only_through_a_fine_label():
    """Xenium folds `07_immune_granulocyte` into `Myeloid` at broad level, whatever the
    protein side does — so a fine label is the only route."""
    assert RESOLVABLE["Granulocyte"]["xenium"] == "fine_only"
    assert xenium_to_common("Myeloid", "Mast") == "Granulocyte"
    assert xenium_to_common("Myeloid", "Myeloid") == "Myeloid"


def test_phenocycler_granulocyte_route_follows_the_core_lineages():
    """With a `Neutrophil` class it maps to Granulocyte; folded into `Immune` it does not
    exist, and MPO-only cells reach Myeloid via the empty immune subclass."""
    from phenocycler.config import LINEAGES

    if "Neutrophil" in LINEAGES:
        assert pheno_to_common("Neutrophil") == "Granulocyte"
    else:
        assert pheno_to_common("Immune", "") == "Myeloid"


def test_epithelial_is_a_default_sink_and_is_flagged():
    assert pheno_to_common("Epithelial") == "Exocrine"
    assert RESOLVABLE["Exocrine"]["phenocycler"] == "default_sink"


def test_unknown_labels_do_not_crash():
    assert pheno_to_common("SomethingNew") == "Unassigned"
    assert xenium_to_common("SomethingNew") == "Unassigned"
    assert pheno_to_common("") == "Unassigned"


def test_every_common_lineage_has_a_resolvability_entry():
    for lineage in COMMON_LINEAGES:
        assert lineage in RESOLVABLE, lineage
        assert set(RESOLVABLE[lineage]) == {"phenocycler", "xenium"}


# --------------------------------------------------------------------------- #
# Panel taxonomy (needs the pinned submodule)
# --------------------------------------------------------------------------- #

def test_panel_roles_drift_guard(tax):
    """Upstream renaming or adding a lineage must be caught, not silently Unassigned."""
    assert check_panel_roles_agreement(tax) == []


def test_identity_cores_loaded(tax):
    assert tax.panel_to_lineage, "no broad_mask panels found"
    for stem in tax.panel_to_lineage:
        assert tax.core_genes.get(stem), f"{stem} has no identity_core genes"
    # The endocrine core is TF/secretory machinery — the hormones are not on the panel.
    endo = set(tax.lineage_genes("Endocrine"))
    assert {"PDX1", "ISL1", "NEUROD1", "ABCC8"} <= endo
    assert not ({"INS", "GCG", "SST"} & endo)


def test_granulocyte_surrogate_genes_available(tax):
    """Granulocytes have no Xenium broad class; the hinted genes are the only route."""
    genes = tax.genes_by_hint("07_immune_granulocyte", "neutrophil")
    assert {"CSF3R", "CXCR1", "CXCR2", "FCGR3B", "MPO"} & set(genes) or genes
    assert len(genes) >= 5


# --------------------------------------------------------------------------- #
# Protein <-> gene crosswalk
# --------------------------------------------------------------------------- #

def test_pairs_are_validated_against_the_actual_panel(tax):
    """A gene absent from the run must not become a phantom anchor."""
    pairs = protein_gene_pairs(["PECAM1", "CD3E"], tax=tax)
    on = {p.protein for p in pairs if p.on_panel}
    assert on == {"CD31", "CD3e"}
    assert all(not p.on_panel for p in pairs if p.protein not in on)


def test_off_panel_markers_route_to_surrogate_panels(tax):
    """INS/GCG/SST/Vimentin define PhenoCycler's Endocrine and Fibroblast classes and are
    absent from the Xenium panel — so those comparisons must go through surrogate panels."""
    panel = sorted({g for gs in tax.core_genes.values() for g in gs})
    pairs = {p.protein: p for p in protein_gene_pairs(panel, tax=tax)}
    for marker in ("INS", "GCG", "SST", "Vimentin"):
        p = pairs[marker]
        assert not p.on_panel
        assert not p.usable_anchor, "an off-panel marker must never count as an anchor"
        assert p.surrogate_panel, f"{marker} has no surrogate panel"
        assert p.surrogate_genes, f"{marker} surrogate panel resolved to no genes"


def test_surrogate_genes_are_filtered_to_the_run_panel(tax):
    """A surrogate is only useful if its genes are actually present in the run.

    Filtering here rather than at use-time is what stops a surrogate panel from looking
    available when none of its genes were measured.
    """
    pairs = {p.protein: p for p in protein_gene_pairs(["PECAM1"], tax=tax)}
    assert pairs["INS"].surrogate_panel == "01_pancreas_endocrine"
    assert pairs["INS"].surrogate_genes == []


def test_known_off_panel_set(tax):
    """Guards the empirical finding against silent drift in the candidate table."""
    real_panel = sorted({g for gs in tax.core_genes.values() for g in gs})
    # Add the markers we know ARE on the deployed panel via the custom-100 add-on.
    real_panel += ["KRT19", "ACTA2", "MS4A1", "CD163", "CD68", "CD99", "TUBB3", "MPO",
                   "PECAM1", "CD3E", "CD3G"]
    pairs = {p.protein: p for p in protein_gene_pairs(real_panel, tax=tax)}
    for marker in KNOWN_OFF_PANEL:
        assert not pairs[marker].on_panel, f"{marker} unexpectedly resolved on-panel"


def test_gate_markers_are_flagged():
    for marker in GATE_MARKERS:
        assert marker in PROTEIN_GENE_CANDIDATES, marker
    pairs = protein_gene_pairs(["PECAM1"])
    assert {p.protein for p in pairs if p.is_gate} == set(GATE_MARKERS)


def test_usable_anchors_requires_a_real_gene(tax):
    pairs = protein_gene_pairs(["PECAM1", "MS4A1", "CD163"], tax=tax)
    anchors = usable_anchors(pairs)
    assert {a.protein for a in anchors} == {"CD31", "CD20", "CD163"}
    assert all(a.genes for a in anchors)


def test_type_only_anchors_exclude_state_markers(tax):
    pairs = protein_gene_pairs(["PECAM1", "MKI67"], tax=tax)
    assert {a.protein for a in usable_anchors(pairs)} == {"CD31", "Ki67"}
    assert {a.protein for a in usable_anchors(pairs, type_only=True)} == {"CD31"}


def test_surrogate_panels_reference_real_cores(tax):
    for marker, stem in SURROGATE_PANEL.items():
        if stem:
            assert stem in tax.core_genes, f"{marker} -> unknown identity core {stem}"


def test_crosswalk_frame_and_report(tax):
    pairs = protein_gene_pairs(["PECAM1", "CD3E", "MS4A1"], tax=tax)
    df = crosswalk_frame(pairs, tax=tax)
    assert set(df["kind"]) == {"lineage", "marker"}
    assert len(df[df["kind"] == "lineage"]) == len(COMMON_LINEAGES)
    report = coverage_report(pairs, tax)
    assert "OFF-PANEL LINEAGE GATES" in report
    assert "default_sink" in report


def test_granulocyte_resolvability_tracks_the_core_pipeline():
    """The crosswalk must not advertise a capability the pipeline has stopped having.

    While a dedicated `Neutrophil` lineage exists, PhenoCycler resolves granulocytes at broad
    level. Fold MPO into `Immune` and it does not — `immune_subclass` gates CD3e/CD20/CD163
    only, so an MPO-only cell lands in Myeloid. Deriving from `config.LINEAGES` makes the
    crosswalk follow that instead of drifting behind it.
    """
    from phenocycler.config import LINEAGES
    from phenocycler.integration.vocab import PHENO_TO_COMMON, RESOLVABLE

    has_neutrophil = "Neutrophil" in LINEAGES
    assert ("Neutrophil" in PHENO_TO_COMMON) == has_neutrophil
    assert RESOLVABLE["Granulocyte"]["phenocycler"] == ("yes" if has_neutrophil else "fine_only")
    # Xenium never calls it at broad level, whatever the protein side does.
    assert RESOLVABLE["Granulocyte"]["xenium"] == "fine_only"


def test_crosswalk_covers_exactly_the_core_lineages():
    """The drift guard: a class added to or removed from `LINEAGES` must not leave the
    crosswalk stale in either direction."""
    from phenocycler.config import LINEAGES
    from phenocycler.integration.vocab import PHENO_TO_COMMON

    # `Immune` is deliberately absent — it resolves via immune_subclass, not the broad label.
    expected = set(LINEAGES) - {"Immune"}
    assert set(PHENO_TO_COMMON) == expected, (
        f"crosswalk drifted from phenocycler.config.LINEAGES: "
        f"missing={expected - set(PHENO_TO_COMMON)}, stale={set(PHENO_TO_COMMON) - expected}")


def test_mpo_only_cells_reach_myeloid_when_neutrophil_is_folded_in():
    """Broad types only: with MPO inside `Immune`, a neutrophil has no subclass and maps to
    Myeloid — which is exactly where Xenium puts a granulocyte at broad level."""
    from phenocycler.integration.vocab import pheno_to_common

    assert pheno_to_common("Immune", "") == "Myeloid"
    assert pheno_to_common("Immune", "Myeloid") == "Myeloid"


def test_granulocyte_never_enters_a_comparable_set():
    """Xenium cannot call it at broad level, so it is not evidence either way."""
    from phenocycler.integration.donor import COMPARABLE, comparable_for
    from phenocycler.integration.tissues import TISSUES

    assert "Granulocyte" not in COMPARABLE
    for tissue in TISSUES:
        assert "Granulocyte" not in comparable_for(tissue), tissue
