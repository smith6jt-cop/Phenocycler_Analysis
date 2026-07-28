"""
S2 — harmonising two label vocabularies and two molecular alphabets.

Both pipelines call broad classes, which makes them look comparable. They are not the same
partition of the tissue:

    PhenoCycler   Epithelial  Fibroblast  Muscle  Immune              Endocrine  Endothelial  Neural
    Xenium        Exocrine    Stromal ----------  Immune_T/B/Myeloid  Endocrine  Vascular     Neural

(PhenoCycler additionally carries a ``Neutrophil`` class while MPO is gated separately; see
the granulocyte note below for what happens when it is folded into ``Immune``.)

Four specific mismatches, each handled explicitly rather than averaged away:

*Immune.* PhenoCycler gates CD3e/CD20/CD163 but collapses all three into one ``Immune``
class. Xenium resolves T, B and myeloid as separate broad lineages. The split is recoverable
on the PhenoCycler side without re-running anything, because the per-marker ``_pos`` columns
survive in ``restore_gated_redsea/`` — ``export_pheno`` derives ``immune_subclass`` from
them, and this module maps that, not the broad label, onto the common vocabulary.

*Fibroblast vs Muscle.* PhenoCycler separates them by Vimentin/SMA argmax. Xenium folds
stellate cells and pericyte/smooth-muscle into one ``Stromal`` lineage. Worse, ``ACTA2`` is
not in the Xenium 5K reference at all — it arrives only via the custom-100 add-on, and
``12_pericyte_smooth_muscle`` contributes exactly one SMC-hinted gene (``MYLK``). So the
common vocabulary keeps ``Stromal`` as the comparable unit and marks the finer
Fibroblast/Muscle split PhenoCycler-only.

*Granulocytes.* This one tracks the core pipeline rather than asserting a fixed answer.
While a dedicated ``Neutrophil`` lineage exists (MPO gate), PhenoCycler resolves granulocytes
at broad level and Xenium does not — Xenium folds ``07_immune_granulocyte`` into ``Myeloid``,
reaching ``Granulocyte`` only through a fine label. When MPO is instead folded into ``Immune``
(neutrophils are immune cells), PhenoCycler stops resolving them too: ``immune_subclass``
gates CD3e/CD20/CD163 only, so an MPO-only cell has no subclass and lands in ``Myeloid`` —
exactly where Xenium puts it. The two sides then agree, at the cost of the distinction.

So ``PHENO_TO_COMMON`` and ``RESOLVABLE`` are *derived* from ``config.LINEAGES`` instead of
hardcoded. A mapping for a class the pipeline no longer emits is not harmless: it survives as
a dead entry while the crosswalk CSV keeps advertising a capability that is gone. Either way
``Granulocyte`` never enters the comparable set, because Xenium cannot call it at broad level
— a lineage only one side can see is not evidence of a difference between modalities.

Restoring the distinction later is additive on both sides: an ``MPO`` entry in
``export_pheno.IMMUNE_GATES`` for protein, the granulocyte identity-core score for RNA.

*Epithelial.* Xenium deliberately has no epithelial mask — ``16_epithelial_general`` is
``embedding_extra``, not ``broad_mask``, because KRT/EPCAM are shared by ductal, acinar and
endocrine cells and scoring them as one lineage biases the argmax. PhenoCycler's
``Epithelial``, meanwhile, is a *default sink*. Mapping it to Exocrine is defensible but the
asymmetry is real, which is why the contract carries ``evidence_positive``.

The mapping is *derived* from ``external/XeniumPanelExplorer/data/panel_roles.csv`` and the
``identity_core/*.csv`` files rather than hardcoded — the upstream repo ships those two
artifacts specifically as a machine-readable contract for Python consumers, and warns that
unpinned masks drift silently. We pin the submodule and read them.

The molecular crosswalk has its own trap. 46 of 55 PhenoCycler markers have a gene on the
deployed 5,081-gene panel, which sounds like plenty — but **INS, GCG, SST and Vimentin are
off-panel**, and those are precisely the markers that define PhenoCycler's Endocrine and
Fibroblast classes. So for those lineages there is no 1:1 protein<->gene pair and the
comparison has to route through Xenium's curated *surrogate* identity panels (PDX1, ISL1,
NEUROD1, ABCC8... for endocrine). That is orthogonal evidence rather than circular
confirmation, which is scientifically a feature — but it must be stated, not assumed away,
and anything that consumes ``protein_gene_pairs()`` must handle the reduced anchor set.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import pandas as pd

from ..config import LINEAGES as _CORE_LINEAGES
from ..config import PipelineConfig
from ..marker_taxonomy import PROCESS, TYPE

# --------------------------------------------------------------------------- #
# The common vocabulary
# --------------------------------------------------------------------------- #

#: Harmonised broad classes. Deliberately finer than PhenoCycler on the immune axis and
#: coarser on the stromal axis — each level is one that *both* modalities can actually
#: support, given `Stromal` cannot be split by Xenium and `Immune` cannot be left unsplit
#: without discarding Xenium's resolution.
COMMON_LINEAGES: tuple[str, ...] = (
    "Endocrine",
    "Exocrine",       # PhenoCycler Epithelial (Pan-CK) <-> Xenium Acinar/Ductal
    "Stromal",        # PhenoCycler Fibroblast + Muscle <-> Xenium Stellate + pericyte/SMC
    "Vascular",
    "T_NK",
    "B_Plasma",
    "Myeloid",
    "Granulocyte",
    "Neural",
    "Unassigned",
)

#: Every PhenoCycler broad class this repo has ever emitted. `PHENO_TO_COMMON` below is this
#: filtered to the classes the *current* `config.LINEAGES` actually produces, so the crosswalk
#: tracks the core pipeline instead of drifting from it. `Neutrophil` is the live example:
#: folding MPO into `Immune` deletes the class, and a hardcoded mapping would have been left
#: behind pointing at a lineage nothing emits.
_PHENO_BROAD_TO_COMMON: dict[str, str] = {
    "Epithelial": "Exocrine",
    "Fibroblast": "Stromal",
    "Muscle": "Stromal",
    "Endocrine": "Endocrine",
    "Endothelial": "Vascular",
    "Neural": "Neural",
    "Neutrophil": "Granulocyte",
}

#: PhenoCycler broad_lineage -> common. `Immune` is intentionally absent: it is resolved
#: through `immune_subclass` by :func:`pheno_to_common`, because mapping it to a single
#: common class would throw away information both modalities have.
PHENO_TO_COMMON: dict[str, str] = {
    k: v for k, v in _PHENO_BROAD_TO_COMMON.items() if k in _CORE_LINEAGES
}

#: PhenoCycler `immune_subclass` (from CD3e/CD20/CD163 `_pos`) -> common.
PHENO_IMMUNE_TO_COMMON: dict[str, str] = {
    "T": "T_NK",
    "B": "B_Plasma",
    "Myeloid": "Myeloid",
    "Mixed": "Myeloid",     # multi-gate-positive; myeloid is the least-specific bucket
    "": "Myeloid",           # includes MPO-only neutrophils — see the module docstring
}

#: Xenium celltype_broad -> common. Mirrors `panel_roles.csv` `lineage_or_group`, which is
#: verified at load time by :func:`check_panel_roles_agreement`.
XENIUM_TO_COMMON: dict[str, str] = {
    "Endocrine": "Endocrine",
    "Exocrine": "Exocrine",
    "Stromal": "Stromal",
    "Vascular": "Vascular",
    "Immune_T": "T_NK",
    "Immune_B": "B_Plasma",
    "Myeloid": "Myeloid",
    "Neural": "Neural",
}

#: Xenium fine labels that should override the broad mapping. Granulocytes have no Xenium
#: broad class (folded into Myeloid), so the only route to a Granulocyte call is a fine
#: label or the granulocyte identity-core score.
XENIUM_FINE_TO_COMMON: dict[str, str] = {
    "Mast": "Granulocyte",
    "Neutrophil": "Granulocyte",
    "Eosinophil": "Granulocyte",
    "Basophil": "Granulocyte",
}

#: Which modality can resolve which common lineage, and why not when it cannot. Consumers
#: check this before reporting a lineage — a class only one side can see is not evidence of
#: a difference between modalities.
RESOLVABLE: dict[str, dict[str, str]] = {
    # Panel-only default: INS/GCG/SST are absent from the 5K panel, so Xenium infers
    # endocrine identity from surrogate cores (PDX1, ISL1, NEUROD1, ABCC8...). Post-Xenium IF
    # staining for INS/GCG makes it a direct protein measurement on the same cells — see
    # `resolvable_for`, which upgrades this when [integration] xenium_hormone_if is set.
    "Endocrine":   {"phenocycler": "yes", "xenium": "surrogate"},
    "Exocrine":    {"phenocycler": "default_sink", "xenium": "yes"},
    "Stromal":     {"phenocycler": "yes", "xenium": "yes"},
    "Vascular":    {"phenocycler": "yes", "xenium": "yes"},
    "T_NK":        {"phenocycler": "subclass", "xenium": "yes"},
    "B_Plasma":    {"phenocycler": "subclass", "xenium": "yes"},
    "Myeloid":     {"phenocycler": "subclass", "xenium": "yes"},
    # Tracks the core pipeline. While a dedicated `Neutrophil` lineage exists, PhenoCycler
    # resolves granulocytes at broad level; once MPO is folded into `Immune` it does not, and
    # an MPO-only cell lands in Myeloid — exactly where Xenium puts a granulocyte. Hardcoding
    # "yes" would publish a capability the pipeline had stopped having.
    "Granulocyte": {"phenocycler": "yes" if "Neutrophil" in _CORE_LINEAGES else "fine_only",
                    "xenium": "fine_only"},
    "Neural":      {"phenocycler": "yes", "xenium": "yes"},
    "Unassigned":  {"phenocycler": "never", "xenium": "no"},
}

_RESOLVABLE_NOTES = {
    "yes": "resolved by a positive broad-level call",
    "surrogate": "no 1:1 marker; resolved via curated surrogate identity panel",
    "subclass": "recovered from per-marker gates, not the broad label",
    "default_sink": "assigned by fallback when all positive gates fail (epi_default)",
    "fine_only": "no broad class; only reachable as a fine subtype",
    "never": "this modality has no Unassigned bucket by construction",
    "no": "not resolvable",
}

#: PhenoCycler-only distinctions that the common vocabulary deliberately collapses. Kept so
#: the loss is documented in the crosswalk CSV rather than being invisible.
PHENO_ONLY_SPLITS = {
    "Stromal": ["Fibroblast (Vimentin)", "Muscle (SMA)"],
}

# --------------------------------------------------------------------------- #
# Protein <-> gene crosswalk
# --------------------------------------------------------------------------- #

#: PhenoCycler marker -> candidate gene symbol(s), best first. Candidates are *proposals*;
#: :func:`protein_gene_pairs` keeps only those actually present in the run's panel, so a
#: panel change silently shrinks the anchor set instead of producing phantom features.
PROTEIN_GENE_CANDIDATES: dict[str, list[str]] = {
    # epithelial / exocrine
    "Pan_Cytokeratin": ["KRT19", "KRT8", "KRT18", "KRT7"],
    "EpCAM": ["EPCAM"],
    "E_cadherin": ["CDH1"],
    "TP63": ["TP63"],
    "CD66": ["CEACAM1", "CEACAM5", "CEACAM6"],
    # mesenchyme / muscle
    "Vimentin": ["VIM"],                    # off-panel -> surrogate
    "Collagen_IV": ["COL4A1", "COL4A2"],
    "SMA": ["ACTA2"],                       # custom-100 only
    # neural
    "B3TUBB": ["TUBB3"],
    "CD56": ["NCAM1"],
    # vascular
    "CD31": ["PECAM1"],
    "CD34": ["CD34"],
    "Podoplanin": ["PDPN"],
    "Caveolin": ["CAV1"],                   # off-panel
    # endocrine
    "INS": ["INS"],                         # off-panel -> surrogate
    "GCG": ["GCG"],                         # off-panel -> surrogate
    "SST": ["SST"],                         # off-panel -> surrogate
    "CD99": ["CD99"],
    # immune
    "CD3e": ["CD3E", "CD3D", "CD3G"],
    "CD8": ["CD8A", "CD8B"],
    "CD4": ["CD4"],
    "FOXP3": ["FOXP3"],
    "CD20": ["MS4A1"],
    "CD79a": ["CD79A"],
    "CD68": ["CD68"],
    "CD163": ["CD163"],
    "Iba1": ["AIF1"],                       # off-panel
    "CD206": ["MRC1"],
    "M2Gal3": ["LGALS3"],                   # off-panel
    "CD11c": ["ITGAX"],
    "CD209": ["CD209"],
    "CD11b": ["ITGAM"],
    "MPO": ["MPO"],
    "b_Catenin1": ["CTNNB1"],
    # process / state
    "Ki67": ["MKI67"],
    "PCNA": ["PCNA"],
    "Granzyme_B": ["GZMB"],
    "CD107a": ["LAMP1"],
    "PD_1": ["PDCD1"],
    "PD_L1": ["CD274"],
    "LAG3": ["LAG3"],
    "TOX": ["TOX"],
    "VISTA": ["VSIR"],
    "ICOS": ["ICOS"],
    "TCF_1": ["TCF7"],
    "CD39": ["ENTPD1"],
    "HLA_DR": ["HLA-DRA", "HLA-DRB1"],      # off-panel
    "HLA_A": ["HLA-A"],                     # off-panel
    "iNOS": ["NOS2"],
    "IDO1": ["IDO1"],
    "SOX2": ["SOX2"],
    "Bcl_2": ["BCL2"],
    "CD44": ["CD44"],
    "CD57": ["B3GAT1"],
    "CD38": ["CD38"],
}

#: Where a protein has no gene on the panel, the comparison routes through the Xenium
#: identity-core panel for the lineage that protein defines. `identity_core/<stem>.csv`
#: supplies the actual gene list; this only names which panel to use.
SURROGATE_PANEL: dict[str, str] = {
    "INS": "01_pancreas_endocrine",
    "GCG": "01_pancreas_endocrine",
    "SST": "01_pancreas_endocrine",
    "Vimentin": "13_fibroblast_stellate",
    "Caveolin": "11_endothelial_vascular",
    "Iba1": "06_immune_myeloid",
    "M2Gal3": "06_immune_myeloid",
    "HLA_DR": "06_immune_myeloid",
    "HLA_A": "",                            # ubiquitous; no lineage surrogate is meaningful
}

#: Markers that gate PhenoCycler's broad lineage call. Anchors drawn from this set carry the
#: most weight for cross-modal identity comparison.
GATE_MARKERS = ("Pan_Cytokeratin", "Vimentin", "SMA", "CD31", "INS", "GCG", "SST",
                "CD3e", "CD20", "CD163", "CD99", "B3TUBB", "MPO")


class VocabError(RuntimeError):
    """The panel taxonomy could not be loaded or is internally inconsistent."""


@dataclass
class MarkerPair:
    """One resolved protein <-> gene(s) correspondence."""

    protein: str
    genes: list[str]                  # present on the panel; empty when off-panel
    on_panel: bool
    is_gate: bool
    role: str                         # 'TYPE' | 'PROCESS' | 'OTHER'
    surrogate_panel: str = ""         # identity_core stem used when off-panel
    surrogate_genes: list[str] = field(default_factory=list)

    @property
    def usable_anchor(self) -> bool:
        """Whether this pair can anchor a shared feature space.

        Requires a real gene on the panel — a surrogate panel supports *lineage-level*
        comparison but is not a per-cell feature correspondence, so it must not be counted
        toward the pseudo-cell linking anchor budget.
        """
        return self.on_panel and bool(self.genes)


@dataclass
class PanelTaxonomy:
    """The XeniumPanelExplorer identity layer, loaded for one tissue."""

    tissue: str
    root: Path
    #: identity_core stem -> Xenium broad lineage (from panel_roles.csv broad_mask rows)
    panel_to_lineage: dict[str, str] = field(default_factory=dict)
    #: identity_core stem -> [gene]
    core_genes: dict[str, list[str]] = field(default_factory=dict)
    #: identity_core stem -> {gene: subtype_hint}
    core_hints: dict[str, dict[str, str]] = field(default_factory=dict)
    #: stems flagged embedding_extra rather than broad_mask (Epithelial)
    embedding_extra: list[str] = field(default_factory=list)

    def lineage_genes(self, lineage: str) -> list[str]:
        """Union of identity-core genes mapped to a Xenium broad lineage."""
        out: list[str] = []
        for stem, lin in self.panel_to_lineage.items():
            if lin == lineage:
                out.extend(self.core_genes.get(stem, []))
        return sorted(set(out))

    def genes_by_hint(self, stem: str, hint: str) -> list[str]:
        """Genes in one identity core carrying a given ``subtype_hint``.

        This is how a Granulocyte score is recovered from ``07_immune_granulocyte``, whose
        neutrophil-hinted genes are the only route to a class Xenium has no broad label for.
        """
        return sorted(g for g, h in self.core_hints.get(stem, {}).items() if h == hint)


def _read_csv_rows(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def load_panel_taxonomy(cfg: PipelineConfig, tissue: Optional[str] = None) -> PanelTaxonomy:
    """Load ``panel_roles.csv`` + the ``identity_core`` CSVs from the pinned submodule.

    ``tissue`` selects which cores apply. XeniumPanelExplorer ships ``pancreas`` and
    ``thymus``; a tissue with no directory of its own (pancreatic lymph node) resolves to the
    **shared** cores — T cell, B/plasma, NK/ILC, myeloid, granulocyte, endothelial, pericyte,
    neural, epithelial — and drops the three pancreas-specific ones. That is the correct
    answer for a lymph node rather than a degraded fallback: scoring lymph-node cells against
    ``01_pancreas_endocrine`` would put weight on an endocrine axis that cannot exist there,
    and the broad call is an argmax over those axes.

    Raises with an actionable message when the submodule is not checked out — this is the
    single most likely first-run failure, and 'FileNotFoundError: panel_roles.csv' would not
    tell anyone to run ``git submodule update``.
    """
    from .tissues import applies as tissue_applies
    from .tissues import panel_tissue_dir

    tissue = tissue or cfg.tissue
    root = Path(cfg.panel_explorer)
    roles_csv = root / "data" / "panel_roles.csv"
    if not roles_csv.exists():
        raise VocabError(
            f"XeniumPanelExplorer panel taxonomy not found at {roles_csv}.\n"
            f"  git submodule update --init --recursive\n"
            f"or point [integration] panel_explorer at a checkout."
        )

    tax = PanelTaxonomy(tissue=tissue, root=root)
    # Which row set in panel_roles.csv describes this tissue's cores. A tissue with no
    # directory borrows the pancreas rows and then filters out the pancreas-specific stems,
    # which is how the shared pool is reached without duplicating the role table upstream.
    panel_dir = panel_tissue_dir(tissue)
    roles_tissue = panel_dir or "pancreas"

    for row in _read_csv_rows(roles_csv):
        if row.get("tissue") != roles_tissue:
            continue
        stem = (row.get("subpanel") or "").strip()
        use = (row.get("phenotyping_use") or "").strip()
        lineage = (row.get("lineage_or_group") or "").strip()
        if not stem or not tissue_applies(tissue, stem):
            continue
        if use == "broad_mask":
            tax.panel_to_lineage[stem] = lineage
        elif use == "embedding_extra":
            tax.embedding_extra.append(stem)

    # identity_core lives per-tissue for tissue-specific panels and in a shared pool for the
    # rest; panel_roles.phenotyping_source records which, but checking both is simpler and
    # equivalent (the tissue copy wins).
    search_dirs = [
        *( [root / "data" / "tissues" / panel_dir / "identity_core"] if panel_dir else [] ),
        root / "data" / "identity_core",
    ]
    for stem in list(tax.panel_to_lineage) + tax.embedding_extra:
        for d in search_dirs:
            path = d / f"{stem}.csv"
            if not path.exists():
                continue
            rows = _read_csv_rows(path)
            genes, hints = [], {}
            for r in rows:
                g = (r.get("gene") or "").strip()
                if not g:
                    continue
                genes.append(g)
                hints[g] = (r.get("subtype_hint") or "").strip()
            tax.core_genes[stem] = genes
            tax.core_hints[stem] = hints
            break

    missing = [s for s in tax.panel_to_lineage if s not in tax.core_genes]
    if missing:
        raise VocabError(
            f"identity_core CSVs missing for {missing} under {root}/data — the pinned "
            f"submodule looks incomplete"
        )
    return tax


def check_panel_roles_agreement(tax: PanelTaxonomy) -> list[str]:
    """Verify every upstream lineage has a place in the common vocabulary.

    A drift guard, not a formality: XeniumPanelExplorer is a separate repo on its own release
    cycle, and a renamed or added lineage would otherwise land silently in `Unassigned`.
    """
    problems = []
    for stem, lineage in sorted(tax.panel_to_lineage.items()):
        if lineage not in XENIUM_TO_COMMON:
            problems.append(
                f"panel_roles.csv maps {stem} -> lineage {lineage!r}, which XENIUM_TO_COMMON "
                f"does not cover (known: {sorted(XENIUM_TO_COMMON)})"
            )
    return problems


def pheno_to_common(lineage_native: str, immune_subclass: str = "") -> str:
    """PhenoCycler broad lineage (+ immune subclass) -> common lineage."""
    lineage_native = (lineage_native or "").strip()
    if lineage_native == "Immune":
        return PHENO_IMMUNE_TO_COMMON.get((immune_subclass or "").strip(), "Myeloid")
    return PHENO_TO_COMMON.get(lineage_native, "Unassigned")


def xenium_to_common(celltype_broad: str, celltype_fine: str = "") -> str:
    """Xenium broad (+ fine) -> common lineage. Fine wins where it is more specific."""
    fine = (celltype_fine or "").strip()
    if fine in XENIUM_FINE_TO_COMMON:
        return XENIUM_FINE_TO_COMMON[fine]
    return XENIUM_TO_COMMON.get((celltype_broad or "").strip(), "Unassigned")


def map_pheno_series(lineage_native: pd.Series, immune_subclass: Optional[pd.Series] = None
                     ) -> pd.Series:
    """Vectorised :func:`pheno_to_common` over a column pair."""
    if immune_subclass is None:
        immune_subclass = pd.Series([""] * len(lineage_native), index=lineage_native.index)
    return pd.Series(
        [pheno_to_common(a, b) for a, b in zip(lineage_native.fillna(""),
                                               immune_subclass.fillna(""))],
        index=lineage_native.index, dtype="object",
    )


def map_xenium_series(celltype_broad: pd.Series, celltype_fine: Optional[pd.Series] = None
                      ) -> pd.Series:
    """Vectorised :func:`xenium_to_common` over a column pair."""
    if celltype_fine is None:
        celltype_fine = pd.Series([""] * len(celltype_broad), index=celltype_broad.index)
    return pd.Series(
        [xenium_to_common(a, b) for a, b in zip(celltype_broad.fillna(""),
                                                celltype_fine.fillna(""))],
        index=celltype_broad.index, dtype="object",
    )


def protein_gene_pairs(
    panel_genes: Iterable[str],
    tax: Optional[PanelTaxonomy] = None,
    markers: Optional[Iterable[str]] = None,
) -> list[MarkerPair]:
    """Resolve the protein <-> gene crosswalk against a run's actual panel.

    Parameters
    ----------
    panel_genes
        The gene symbols actually present in the Xenium run (``adata.var_names``). Passing
        the *nominal* panel here rather than the run's own genes is the mistake this
        signature exists to prevent: the deployed panel is the 5K base plus a custom add-on,
        and several textbook markers live only in the add-on.
    tax
        Panel taxonomy, used to attach surrogate gene lists for off-panel proteins.
    markers
        Restrict to these PhenoCycler markers; default is every marker with a candidate.
    """
    present = {str(g).strip() for g in panel_genes}
    wanted = list(markers) if markers is not None else list(PROTEIN_GENE_CANDIDATES)
    type_set, process_set = set(TYPE), set(PROCESS)

    out: list[MarkerPair] = []
    for protein in wanted:
        candidates = PROTEIN_GENE_CANDIDATES.get(protein, [])
        hits = [g for g in candidates if g in present]
        role = "TYPE" if protein in type_set else ("PROCESS" if protein in process_set
                                                   else "OTHER")
        stem = SURROGATE_PANEL.get(protein, "") if not hits else ""
        surrogate_genes: list[str] = []
        if stem and tax is not None:
            surrogate_genes = [g for g in tax.core_genes.get(stem, []) if g in present]
        out.append(MarkerPair(
            protein=protein,
            genes=hits,
            on_panel=bool(hits),
            is_gate=protein in GATE_MARKERS,
            role=role,
            surrogate_panel=stem,
            surrogate_genes=surrogate_genes,
        ))
    return out


def usable_anchors(pairs: Iterable[MarkerPair], *, type_only: bool = False) -> list[MarkerPair]:
    """Pairs that can anchor a shared feature space.

    ``type_only`` keeps identity markers and drops state/process ones — appropriate when the
    anchors must discriminate cell *type* rather than activation state.
    """
    out = [p for p in pairs if p.usable_anchor]
    if type_only:
        out = [p for p in out if p.role == "TYPE"]
    return out


def resolvable_for(cfg: Optional[PipelineConfig] = None) -> dict[str, dict[str, str]]:
    """`RESOLVABLE`, adjusted for what this run's Xenium sections actually carry.

    Post-Xenium immunofluorescence for INS/GCG turns Xenium's endocrine call from a surrogate
    inference into a direct protein measurement *on the same cells as the RNA* — the single
    biggest gap in the panel crosswalk, since the hormones are the markers that define the
    compartment a T1D study is about. Reported rather than assumed: with the flag off, the
    crosswalk keeps saying `surrogate`, which is the truthful answer for panel-only data.
    """
    out = {k: dict(v) for k, v in RESOLVABLE.items()}
    if cfg is not None and getattr(cfg, "xenium_hormone_if", False):
        out["Endocrine"]["xenium"] = "yes"
    return out


def crosswalk_frame(
    pairs: Iterable[MarkerPair],
    tax: Optional[PanelTaxonomy] = None,
) -> pd.DataFrame:
    """The full crosswalk as a table — this is what gets written to disk and reviewed.

    Two blocks: one row per common lineage (how each modality resolves it), then one row per
    PhenoCycler marker (its gene, or the surrogate panel standing in for it).
    """
    rows: list[dict] = []

    for lineage in COMMON_LINEAGES:
        res = RESOLVABLE.get(lineage, {})
        p, x = res.get("phenocycler", "no"), res.get("xenium", "no")
        xen_native = sorted({k for k, v in XENIUM_TO_COMMON.items() if v == lineage})
        pheno_native = sorted({k for k, v in PHENO_TO_COMMON.items() if v == lineage})
        if lineage in PHENO_IMMUNE_TO_COMMON.values():
            pheno_native += [f"Immune/{k}" for k, v in PHENO_IMMUNE_TO_COMMON.items()
                             if v == lineage and k]
        rows.append({
            "kind": "lineage",
            "common": lineage,
            "phenocycler": "|".join(pheno_native),
            "xenium": "|".join(xen_native),
            "pheno_resolution": p,
            "xenium_resolution": x,
            "pheno_note": _RESOLVABLE_NOTES.get(p, ""),
            "xenium_note": _RESOLVABLE_NOTES.get(x, ""),
            "collapsed": "|".join(PHENO_ONLY_SPLITS.get(lineage, [])),
            "identity_core_genes": (len(tax.lineage_genes(
                next((k for k, v in XENIUM_TO_COMMON.items() if v == lineage), "")))
                if tax is not None else ""),
        })

    for mp in pairs:
        rows.append({
            "kind": "marker",
            "common": "",
            "phenocycler": mp.protein,
            "xenium": "|".join(mp.genes),
            "pheno_resolution": mp.role,
            "xenium_resolution": "on_panel" if mp.on_panel else "off_panel",
            "pheno_note": "lineage gate" if mp.is_gate else "",
            "xenium_note": (f"surrogate: {mp.surrogate_panel} "
                            f"({len(mp.surrogate_genes)} genes)" if mp.surrogate_panel else ""),
            "collapsed": "",
            "identity_core_genes": len(mp.surrogate_genes) if mp.surrogate_panel else "",
        })

    return pd.DataFrame(rows)


def coverage_report(pairs: Iterable[MarkerPair], tax: Optional[PanelTaxonomy] = None) -> str:
    """Human-readable summary — printed by the stage and pasted into notebooks."""
    pairs = list(pairs)
    on = [p for p in pairs if p.on_panel]
    off = [p for p in pairs if not p.on_panel]
    gates_on = [p for p in on if p.is_gate]
    gates_off = [p for p in off if p.is_gate]
    anchors = usable_anchors(pairs, type_only=True)

    lines = [
        "PhenoCycler <-> Xenium vocabulary crosswalk",
        "=" * 60,
        f"common lineages      : {len(COMMON_LINEAGES) - 1} (+Unassigned)",
        f"markers with a gene  : {len(on)}/{len(pairs)}",
        f"  lineage gates      : {len(gates_on)}/{len(GATE_MARKERS)}",
        f"identity anchors     : {len(anchors)} (TYPE markers with a panel gene)",
        "",
    ]
    if gates_off:
        lines.append("OFF-PANEL LINEAGE GATES — these route through surrogate panels, so the")
        lines.append("corresponding lineage comparison is orthogonal evidence, not a 1:1 check:")
        for p in gates_off:
            sg = f"{len(p.surrogate_genes)} surrogate genes" if p.surrogate_genes else "NO surrogate"
            lines.append(f"    {p.protein:<18} -> {p.surrogate_panel or '(none)'}  [{sg}]")
        lines.append("")
    other_off = [p for p in off if not p.is_gate]
    if other_off:
        lines.append(f"other off-panel markers ({len(other_off)}): "
                     f"{', '.join(p.protein for p in other_off)}")
        lines.append("")

    lines.append("lineage resolution by modality:")
    lines.append(f"  {'common':<13} {'phenocycler':<15} {'xenium':<15} note")
    for lin in COMMON_LINEAGES:
        r = RESOLVABLE.get(lin, {})
        p, x = r.get("phenocycler", "no"), r.get("xenium", "no")
        note = ""
        if lin in PHENO_ONLY_SPLITS:
            note = "collapses " + " + ".join(PHENO_ONLY_SPLITS[lin])
        elif x == "fine_only":
            note = "no Xenium broad class"
        elif p == "default_sink":
            note = "PhenoCycler fallback bucket — compare evidence-positive only"
        lines.append(f"  {lin:<13} {p:<15} {x:<15} {note}")

    if tax is not None:
        probs = check_panel_roles_agreement(tax)
        if probs:
            lines.append("")
            lines.append("DRIFT WARNINGS (XeniumPanelExplorer changed under us):")
            lines.extend(f"    {p}" for p in probs)
    return "\n".join(lines)


def run_vocab(cfg: PipelineConfig, panel_genes: Optional[Iterable[str]] = None,
              *, tissue: Optional[str] = None, write: bool = True) -> pd.DataFrame:
    """Stage entry point: build the crosswalk, print coverage, write the CSV.

    ``panel_genes`` defaults to the union of the tissue's identity cores, which is enough to
    build the lineage half of the crosswalk without any Xenium run mounted. Pass the run's
    real ``var_names`` for an accurate marker-level report.

    The crosswalk is tissue-specific because the identity cores are: a lymph node is scored
    against the shared immune/vascular/stromal cores only, so its lineage rows carry
    different core-gene counts than the pancreas rows. Each tissue therefore gets its own
    ``vocab_crosswalk_{tissue}.csv``, and the canonical ``vocab_crosswalk.csv`` is the
    concatenation of them, tagged by tissue.
    """
    tissue = tissue or cfg.tissue
    tax = load_panel_taxonomy(cfg, tissue)
    if panel_genes is None:
        panel_genes = sorted({g for genes in tax.core_genes.values() for g in genes})
        print(f"[vocab] {tissue}: no run panel supplied — marker coverage is computed against "
              "the identity-core union only and will understate on-panel markers", flush=True)

    pairs = protein_gene_pairs(panel_genes, tax=tax)
    print(coverage_report(pairs, tax), flush=True)

    df = crosswalk_frame(pairs, tax=tax)
    df.insert(0, "tissue", tissue)
    if write:
        out = cfg.vocab_crosswalk_csv
        out.parent.mkdir(parents=True, exist_ok=True)
        per_tissue = out.with_name(f"{out.stem}_{tissue}{out.suffix}")
        df.to_csv(per_tissue, index=False)
        # Merge into the canonical file rather than overwrite it, so a per-tissue re-run does
        # not silently delete the other tissue's rows.
        combined = df
        if out.exists():
            try:
                prev = pd.read_csv(out)
                if "tissue" in prev.columns:
                    prev = prev[prev["tissue"] != tissue]
                    combined = pd.concat([prev, df], ignore_index=True)
            except Exception:  # noqa: BLE001 - a corrupt/legacy file is replaced, not fatal
                combined = df
        combined.to_csv(out, index=False)
        print(f"\n[vocab] wrote {per_tissue} ({len(df)} rows) and {out} "
              f"({len(combined)} rows)", flush=True)
    return df


def main(argv: Optional[list[str]] = None) -> int:
    import argparse

    from ..config import load_config

    ap = argparse.ArgumentParser(description="Build the PhenoCycler<->Xenium vocab crosswalk")
    ap.add_argument("--config", default=None)
    ap.add_argument("--panel", default=None,
                    help="text file of gene symbols (one per line) from the actual run")
    ap.add_argument("--tissue", action="append", dest="tissues",
                    help="build for this tissue (repeatable); default is every configured one")
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    genes = None
    if args.panel:
        genes = [ln.strip() for ln in Path(args.panel).read_text().splitlines() if ln.strip()]
    for tissue in (args.tissues or cfg.tissue_list):
        run_vocab(cfg, genes, tissue=tissue, write=not args.no_write)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
