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
`DAPI` (nuclear stain — not a marker) is EXCLUDED from everything. `Ker8_18`/`Keratin_5` (basal-ductal
keratins) are now USED; `IAPP` FAILED (removed 2026-07-10 → EXCLUDED); `CD99` is a broad state marker (PROCESS), NOT a gate.
`Collagen_IV` is a **MASK** (basement-membrane ECM, not a per-cell signal) — out of the identity heatmap.
`Ki67`/`PCNA` (**PROLIFERATION**) are lineage-agnostic with no mutually-exclusive counterpart, so they are
binarized by a standalone per-image bimodal fit (`restore.bimodal_thresh`) rather than a RESTORE pair.

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
    "CD99",                                                     # broad (~96% detectable) — state marker, NOT a lineage gate
]

# --- Never used for phenotyping OR any plot ---
EXCLUDED = ["DAPI", "IAPP"]   # DAPI = nuclear stain; IAPP = failed marker (removed 2026-07-10). Ker8_18/Keratin_5 used; CD99 -> PROCESS

# --- MASK markers: basement-membrane / ECM, NOT a per-cell signal -> used as a spatial mask; kept OUT of
#     the identity heatmap and phenotyping (per PhenotypingPlan.txt: "use as a mask"). ---
MASK = ["Collagen_IV"]

# --- PROLIFERATION markers: lineage-agnostic, no mutually-exclusive counterpart, so RESTORE pairing does
#     not apply -> binarized by a STANDALONE per-image bimodal fit (restore.bimodal_thresh). They remain
#     PROCESS/state markers (excluded from identity), but DO get a positivity call. ---
PROLIFERATION = ["Ki67", "PCNA"]

# --- TYPE / lineage-identity markers: drive phenotyping and the identity heatmap ---
# Ordered in CLASS BLOCKS following the heatmap row order (LIN_ORDER: Epithelial, Fibroblast, Muscle,
# Neural, Endothelial, Endocrine, Immune) so the heatmap reads as a rough diagonal and
# markers of the same class sit together (e.g. CD56 next to B3TUBB in the neural block).
# Ordered in the 5 compartment blocks (Epithelial[+endocrine] / Endothelial / Neural / Immune / Mesenchymal)
# so the identity heatmap reads as a rough diagonal.
TYPE = [
    # Epithelial — exocrine (acinar/ductal) incl. the basal-ductal keratins now in use
    "E_cadherin", "Pan_Cytokeratin", "Ker8_18", "Keratin_5", "EpCAM", "TP63", "CD66",
    #   endocrine sub-branch of Epithelial (hormones: beta INS, alpha GCG, delta SST; IAPP removed — failed)
    "INS", "GCG", "SST",
    # Endothelial / vascular (+ lymphatic Podoplanin, CD34)
    "CD31", "CD34", "Podoplanin", "Caveolin",
    # Neural (CD56/NCAM groups with B3TUBB, but CD56 is promiscuous NK/neural/endocrine)
    "B3TUBB", "CD56",
    # Immune — lymphoid then myeloid then granulocyte
    "CD3e", "CD8", "CD4", "FOXP3", "CD20", "CD79a", "CD68", "CD163", "Iba1", "CD206",
    "M2Gal3", "CD11c", "CD209", "CD11b", "MPO",
    # Mesenchymal / stromal  (Collagen_IV -> MASK, out of the identity heatmap)
    "SMA", "Vimentin",
    # batch-2-only structural (present only in batch-2 donors; harmless if absent)
    "b_Catenin1",
]

# CD99 is broadly expressed (~96% detectable) and NOT lineage-clean -> demoted to a state marker
# (PROCESS above), no longer an Endocrine gate. Endocrine = E-cadherin+ epithelial that is hormone+
# (INS/GCG/SST); the old bright-CD99 gate + CD99_BRIGHT constant are removed.

# --- NO_SPILLOVER markers: REDSEA must NOT correct these -------------------------------------------
# REDSEA compensates *lateral* spillover across a shared cell boundary: it subtracts a contact-weighted
# share of each neighbour's 1-px boundary-band signal. That model only makes sense for a marker whose
# signal actually reaches the boundary. A nuclear marker's signal sits at the centre of the cell, so the
# neighbour's band cannot contain its spillover — subtracting there removes real signal and nothing else.
# Measured on the comp_mode=0 run: DAPI lost 23–30% of its cohort median and 4–10% of cells were clamped
# to exactly zero, purely as an artifact of correcting a marker the method was never meant to touch.
# Bai et al. 2021 expose ``markers_of_interest`` for exactly this reason; this list is our equivalent.
#
# Collagen IV is here for a different reason: it is basement-membrane ECM, already declared a ``MASK``
# above, and is not a per-cell signal at all — its "spillover" is the whole measurement.
#
# Scope decision (maintainer, 2026-07-25): exclude nuclear + ECM ONLY. Cytoskeletal, filament, granule
# and cytoplasmic markers KEEP the correction, because their signal is sub-membranous or fills the cell
# body and therefore does reach the shared boundary. That also leaves all seven RESTORE gate markers
# corrected, so a REDSEA change cannot confound a RESTORE threshold sweep.
NO_SPILLOVER = [
    # nuclear — signal is centred, never at the shared membrane
    "DAPI", "Ki67", "PCNA", "FOXP3", "TCF_1", "TOX", "SOX2", "TP63",
    # ECM — not a per-cell signal
    "Collagen_IV",
]

_TYPE = set(TYPE)
_PROCESS = set(PROCESS)
_EXCLUDED = set(EXCLUDED)
_MASK = set(MASK)
_NO_SPILLOVER = set(NO_SPILLOVER)


def spillover_corrected_mask(cols):
    """Boolean mask over `cols`: True where REDSEA compensation should be applied.

    `cols` are sanitized panel channel names in acquisition order (``redsea.SANITIZE`` output), so the
    returned mask lines up with the columns of REDSEA's ``data``/``edge`` matrices. Unknown channels
    default to True (corrected) — the conservative choice, matching the pre-2026-07-25 behaviour where
    every channel was corrected.
    """
    import numpy as np

    return np.array([c not in _NO_SPILLOVER for c in cols], dtype=bool)


def heatmap_markers(available):
    """Return the TYPE markers present in `available`, preserving TYPE order.

    Drops PROCESS and EXCLUDED (DAPI, IAPP) markers. Warns about any available marker
    that is not classified in TYPE/PROCESS/EXCLUDED, so a newly added panel channel is never silently
    shown.
    """
    avail = set(available)
    unknown = sorted(avail - _TYPE - _PROCESS - _EXCLUDED - _MASK)
    if unknown:
        print(f"[marker_taxonomy] WARNING: {len(unknown)} unclassified marker(s) hidden from heatmap "
              f"(add them to TYPE or PROCESS): {unknown}", flush=True)
    return [m for m in TYPE if m in avail]
