"""
What each tissue supports — the single place that knows pancreas from lymph node.

The cohort is 20 donors × {PhenoCycler, Xenium} × {pancreas, pancreatic lymph node}. The two
tissues share most of the pipeline and diverge in two specific ways, and both divergences are
biological rather than technical:

**Structures.** Islets exist in pancreas. Lymph nodes have no islets, so hormone-seeded DBSCAN
has nothing to seed on and islet matching (S6) has nothing to match. Rather than invent a
substitute landmark, the structure level is simply **not run** for lymph node — the
registration-free levels (donor, niche) and the registered grid still apply, and those are
the levels that were always going to carry the lymph-node analysis anyway.

That is a deliberate choice, not an omission. B-cell follicles are the obvious candidate
landmark and the code would support them (``call_landmarks`` is already generic over a seed
lineage), but a follicle called from a 5K RNA panel on one side and a 55-plex protein panel on
the other is a much weaker correspondence than an islet, and shipping it would invite
structure-level claims the data cannot support. ``STRUCTURE_SEEDS`` below is where it would
go if that changes.

**Identity panels.** XeniumPanelExplorer ships curated identity cores per tissue, and has
``pancreas`` and ``thymus`` — no lymph node. But the cores a lymph node needs are exactly the
*shared* ones: T cell, B/plasma, NK/ILC, myeloid, granulocyte, endothelial, pericyte, neural,
epithelial. The three pancreas-specific cores (endocrine, exocrine, fibroblast/stellate) do
not apply. So a tissue with no directory of its own resolves to the shared pool, which is the
correct answer rather than a fallback.
"""

from __future__ import annotations

from typing import Optional

#: Canonical tissue identifiers. Anything else is rejected at import rather than silently
#: creating a third tissue nobody configured.
PANCREAS = "pancreas"
LYMPH_NODE = "pancreatic_lymph_node"
TISSUES: tuple[str, ...] = (PANCREAS, LYMPH_NODE)

#: Human-readable, for figures and reports.
DISPLAY_NAME = {PANCREAS: "pancreas", LYMPH_NODE: "pancreatic lymph node"}

#: ROI prefix -> tissue. One donor has one pancreas ROI and up to three lymph-node ROIs.
ROI_PREFIX_TO_TISSUE = {"pan": PANCREAS, "pln": LYMPH_NODE, "ln": LYMPH_NODE}

#: The XeniumPanelExplorer tissue directory backing each tissue's identity cores.
#: ``None`` means "no tissue-specific cores — use the shared pool", which is the right answer
#: for lymph node rather than a degraded one.
PANEL_TISSUE_DIR: dict[str, Optional[str]] = {PANCREAS: "pancreas", LYMPH_NODE: None}

#: Identity cores that are pancreas-specific and must NOT be applied to lymph node.
PANCREAS_ONLY_PANELS = frozenset({
    "01_pancreas_endocrine", "02_pancreas_exocrine", "13_fibroblast_stellate",
})

#: Structure kinds each tissue supports. Empty tuple = the structure level does not run, and
#: everything gated on it (S6 islet matching, the islet QC metric) is skipped for that tissue.
STRUCTURE_KINDS: dict[str, tuple[str, ...]] = {
    PANCREAS: ("islet", "duct", "vessel"),
    LYMPH_NODE: (),
}

#: Seed definition per structure kind, as (common lineage, description). The registry exists
#: so adding follicles later is a one-line change here plus a config flag, not a redesign.
STRUCTURE_SEEDS: dict[str, tuple[str, str]] = {
    "islet": ("Endocrine", "hormone-positive endocrine cells"),
    "duct": ("Exocrine", "ductal/acinar epithelium"),
    "vessel": ("Vascular", "endothelium"),
    # "follicle": ("B_Plasma", "B-cell follicles") — deliberately not enabled; see module docstring.
}

#: Common lineages worth comparing per tissue. A lymph node has no endocrine or exocrine
#: compartment, so reporting a near-zero fraction for those would be noise dressed as a
#: finding; a pancreas has all of them.
COMPARABLE_LINEAGES: dict[str, tuple[str, ...]] = {
    PANCREAS: ("Endocrine", "Exocrine", "Stromal", "Vascular",
               "T_NK", "B_Plasma", "Myeloid", "Neural"),
    LYMPH_NODE: ("Stromal", "Vascular", "T_NK", "B_Plasma", "Myeloid", "Granulocyte"),
}


class TissueError(ValueError):
    """An unknown or unsupported tissue."""


def normalize(tissue: str) -> str:
    """Canonicalise a tissue name, accepting the obvious spellings.

    Lenient on input (``pLN``, ``lymph node``, ``pancreatic_lymph_node`` all work) and strict
    on output, so config files and CLI flags stay forgiving without letting a typo become a
    third tissue that quietly matches nothing.
    """
    t = str(tissue).strip().lower().replace(" ", "_").replace("-", "_")
    if not t:
        raise TissueError("tissue is empty; it must be declared when the dataset is imported")
    if t in TISSUES:
        return t
    if t in {"pln", "lymph_node", "ln", "pancreatic_ln", "panc_lymph_node"}:
        return LYMPH_NODE
    if t in {"panc", "pan", "pancreatic"}:
        return PANCREAS
    raise TissueError(f"unknown tissue {tissue!r}; expected one of {TISSUES}")


def for_roi(roi: str) -> str:
    """ROI name -> tissue. Returns '' for an unrecognised ROI rather than guessing."""
    r = str(roi).strip().lower()
    for prefix, tissue in ROI_PREFIX_TO_TISSUE.items():
        if r.startswith(prefix):
            return tissue
    return ""


def rois(tissue: str) -> tuple[str, ...]:
    """Known ROI names for a tissue."""
    t = normalize(tissue)
    return ("panc",) if t == PANCREAS else ("pln_1", "pln_2", "pln_3")


def structure_kinds(tissue: str) -> tuple[str, ...]:
    """Structure kinds to extract. Empty means the structure level does not apply."""
    return STRUCTURE_KINDS.get(normalize(tissue), ())


def has_structures(tissue: str) -> bool:
    """Whether structure extraction and islet matching are meaningful for this tissue."""
    return bool(structure_kinds(tissue))


def comparable_lineages(tissue: str) -> tuple[str, ...]:
    """Lineages worth reporting a cross-modal comparison for, in this tissue."""
    return COMPARABLE_LINEAGES.get(normalize(tissue), ())


def panel_tissue_dir(tissue: str) -> Optional[str]:
    """XeniumPanelExplorer tissue directory, or ``None`` for shared-cores-only."""
    return PANEL_TISSUE_DIR.get(normalize(tissue))


def applies(tissue: str, panel_stem: str) -> bool:
    """Whether an identity-core panel applies to this tissue.

    Lymph node gets the shared immune / vascular / stromal / neural cores and none of the
    pancreas-specific ones. Scoring a lymph-node cell against ``01_pancreas_endocrine`` would
    not merely be uninformative — it would put weight on an endocrine axis that cannot exist
    there, and the broad call is an argmax over those axes.
    """
    if normalize(tissue) == PANCREAS:
        return True
    return panel_stem not in PANCREAS_ONLY_PANELS


def describe(tissue: str) -> str:
    """One-line summary of what this tissue supports — used in stage logs."""
    t = normalize(tissue)
    kinds = structure_kinds(t)
    struct = ", ".join(kinds) if kinds else "none (registration-free + grid levels only)"
    return (f"{DISPLAY_NAME[t]}: rois={'/'.join(rois(t))} | structures={struct} | "
            f"{len(comparable_lineages(t))} comparable lineages")
