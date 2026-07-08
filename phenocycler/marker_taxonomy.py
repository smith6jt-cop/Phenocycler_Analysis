"""
Marker taxonomy — cell-TYPE vs cell-PROCESS markers (single source of truth).

Faithful port of ``scripts/senior/marker_taxonomy.py`` (Islet-Explorer-Senior).

Lineage *identity* is defined only by **TYPE markers** (INS/GCG/SST→endocrine, CD3e→T, CD20→B,
CD163/CD68→myeloid, CD31→endothelial, Pan_CK→epithelial, Vimentin→fibroblast, SMA→muscle,
B3TUBB→neural, MPO→neutrophil). **PROCESS/state markers are type-agnostic** — they describe a
function on *any* lineage (Ki67/PCNA=proliferation, Granzyme_B=cytotoxicity, PD-1/LAG3/TOX=exhaustion,
HLA-DR=activation, iNOS/IDO1=inflammation, SOX2=stress/dediff): a Ki67+ cell is a proliferating cell
of whatever type, never a cell type. So TYPE markers drive phenotyping and the identity heatmap;
PROCESS markers drive within-lineage state axes (Phase 5) and are excluded from identity views.

This module is the ONE place the TYPE/PROCESS split is defined so phenotyping and every plot agree.
`DAPI` (nuclear stain — not a marker), `IAPP` (failed marker), and `Ker8_18`/`Keratin_5` (redundant
epithelial keratins — `Pan_Cytokeratin` is the epithelial identity marker) are EXCLUDED from everything.

NB: the *lineage-gate* set (the RESTORE-gated markers that actually drive ``lineage.py`` /
``assign_broad_lineage.py``) is a SUBSET of TYPE — e.g. FOXP3, Caveolin, CD66 are TYPE (shown in the
identity heatmap) but are NOT used to call a broad lineage.
"""

# --- PROCESS / state markers: excluded from identity calls and the identity heatmap ---
PROCESS = [
    "Ki67", "PCNA",                                             # proliferation
    "Granzyme_B", "CD107a",                                     # cytotoxicity / degranulation
    "PD_1", "PD_L1", "LAG3", "TOX", "VISTA", "ICOS", "TCF_1", "CD39",  # checkpoint / exhaustion
    "HLA_DR", "HLA_A",                                          # activation / antigen-presentation state
    "iNOS", "IDO1",                                             # inflammation
    "SOX2",                                                     # stress / de-differentiation
    "Bcl_2",                                                    # survival / apoptosis
    "Beta_actin",                                               # housekeeping
    "CD44", "CD57",                                             # activation / senescence (weak identity)
    "CD38",                                                     # activation / plasma-cell state (batch-1 only)
]

# --- Never used for phenotyping OR any plot ---
EXCLUDED = ["DAPI", "IAPP", "Ker8_18", "Keratin_5"]

# --- TYPE / lineage-identity markers: drive phenotyping and the identity heatmap ---
# Ordered in CLASS BLOCKS following the heatmap row order (LIN_ORDER: Epithelial, Fibroblast, Muscle,
# Neural, Endothelial, Endocrine, Immune, Neutrophil) so the heatmap reads as a rough diagonal and
# markers of the same class sit together (e.g. CD56 next to B3TUBB in the neural block).
TYPE = [
    # epithelial / exocrine (acinar + ductal)
    "Pan_Cytokeratin", "EpCAM", "E_cadherin", "TP63", "CD66",
    # fibroblast / mesenchyme
    "Vimentin", "Collagen_IV",
    # muscle / pericyte
    "SMA",
    # neural (CD56/NCAM groups with B3TUBB — both light up the Neural class)
    "B3TUBB", "CD56",
    # endothelial / vascular
    "CD31", "CD34", "Podoplanin", "Caveolin",
    # endocrine
    "INS", "GCG", "SST", "CD99",
    # immune lineage / subtype (lymphoid then myeloid)
    "CD3e", "CD8", "CD4", "FOXP3", "CD20", "CD79a", "CD68", "CD163", "Iba1", "CD206",
    "M2Gal3", "CD11c", "CD209", "CD11b",
    # neutrophil / myeloid granulocyte
    "MPO",
    # batch-2-only structural (present only in batch-2 donors; harmless if absent)
    "b_Catenin1",
]

# CD99 is a lineage marker for Endocrine, but it is broadly expressed (~96% detectable), so the
# standard RESTORE mean+3σ gate (CD99_norm ≥ 1) over-calls a spatially-diffuse, non-islet population.
# CD99 marks the "other" endocrine cell types we lack specific markers for (PP/epsilon/EC), so we take
# only BRIGHT CD99 = ≥ this multiple of the per-image RESTORE threshold. Validated (cohort-wide): at
# ≥3× the hormone-negative CD99-bright cells are 36% islet-core / 42% core+peri — as islet-coherent as
# INS/GCG/SST (37% core) and 11× enriched over baseline; the permissive ≥1 gate is 0.7× (diffuse).
# NB: must equal PipelineConfig.cd99_bright (guarded by tests/test_config.py).
CD99_BRIGHT = 3.0

_TYPE = set(TYPE)
_PROCESS = set(PROCESS)
_EXCLUDED = set(EXCLUDED)


def heatmap_markers(available):
    """Return the TYPE markers present in `available`, preserving TYPE order.

    Drops PROCESS and EXCLUDED (DAPI/IAPP/Ker8_18/Keratin_5) markers. Warns about any available marker
    that is not classified in TYPE/PROCESS/EXCLUDED, so a newly added panel channel is never silently
    shown.
    """
    avail = set(available)
    unknown = sorted(avail - _TYPE - _PROCESS - _EXCLUDED)
    if unknown:
        print(f"[marker_taxonomy] WARNING: {len(unknown)} unclassified marker(s) hidden from heatmap "
              f"(add them to TYPE or PROCESS): {unknown}", flush=True)
    return [m for m in TYPE if m in avail]
