# Level-2 marker hierarchy — the deferred second phenotyping pass (59-plex)

> **STATUS: DEFERRED SPEC — NOT WIRED.** This is the maintainer's design (2026-07-24) for the *second*
> phenotyping pass, recorded here verbatim so it is not lost and so the Gate-4 architecture stays
> extension-ready. **None of the level-2 pairs below has been RESTORE-screened.** Wiring any of them into
> production requires putting it through the same Gate-1→3 validation the level-1 pairs went through
> (see `docs/restore_faithful_rebuild_plan.md`). The **level-1 broad deliverable** (Gate 4) uses ONLY the
> 7 frozen `restore_validation.ACCEPTED_PAIRS`; `phenocycler/config.py::COMPARTMENT_GATES` gates on those
> 7 targets alone. Do **not** add anything here to `COMPARTMENT_GATES` or `ACCEPTED_PAIRS` without a
> validated screen.

## How the two levels compose

RESTORE thresholds are **global per scene** (one divisor per image per pair), so level 2 is applied the
same way as level 1: each level-2 marker gets its own `[target, counterpart]` pair, is thresholded
independently, and the **hierarchy is applied afterward to interpret the calls** within the broad
compartment a cell already landed in. Level 2 does not re-threshold or move a cell between broad
compartments — it *subtypes* within one. A level-2 call is only read inside its parent node (e.g. CD4 is
read only within CD3e⁺).

**Level-1 anchors (unchanged, not repeated below):** `E-cadherin`, `CD31`, `CD3e`, `CD20`, `CD68`,
`CD11b`, `B3TUBB`, `Vimentin`. (In the *current* production run CD20 is deferred out of RESTORE entirely
— `RESTORE_EXCLUDED_MARKERS` — so the wired level-1 Immune union is `CD3e/CD68/CD11b`; CD20/CD79a are the
first immune additions when this pass runs.)

**Counterpart rules (carried over from level 1):**
- epithelial-family target → **CD31**
- endocrine-confounded target → **EpCAM**
- hormone target → **the other hormone**
- everything else → **E-cadherin**

## Level-2 pair list (verbatim)

```python
level2_pairs = [
    # === EPITHELIAL (E-cadherin+) ===
    ["EpCAM",           "CD31"],   # exocrine-high; must be reliable -- it counterparts CD56
    ["Pan-Cytokeratin", "CD31"],
    ["Ker8-18",         "CD31"],
    ["Keratin 5",       "CD31"],   # basal/metaplastic -- near-zero prevalence, see (2)
    ["TP63",            "CD31"],   # basal, nuclear -- near-zero prevalence, see (2)
    ["SOX2",            "CD31"],   # near-zero prevalence, see (2)
    #   endocrine (intra-islet, microenvironment-matched)
    ["CD56",            "EpCAM"],  # pan-endocrine proxy inside ECAD+; NK inside CD45-union
    ["INS",             "GCG"],    # beta
    ["GCG",             "INS"],    # alpha
    ["SST",             "INS"],    # delta -- sparse, see (2)
    ["IAPP",            "GCG"],    # beta-associated
    # === IMMUNE (CD3e / CD20 / CD68 / CD11b union) ===
    #   T
    ["CD4",             "E-cadherin"],  # interpret ONLY within CD3e+ (monocytes are CD4-dim)
    ["CD8",             "E-cadherin"],
    ["FOXP3",           "E-cadherin"],  # Treg -- sparse, see (2)
    #   B / plasma
    ["CD79a",           "E-cadherin"],  # recovers CD20- plasma cells
    ["CD38",            "E-cadherin"],
    #   myeloid
    ["Iba1",            "E-cadherin"],  # pan-macrophage, confirms CD68
    ["CD163",           "E-cadherin"],  # polarization, nested under CD68/Iba1
    ["CD206",           "E-cadherin"],
    ["CD11c",           "E-cadherin"],  # DC
    ["CD209",           "E-cadherin"],  # DC -- sparse, see (2)
    ["MPO",             "E-cadherin"],  # neutrophil -- sparse, see (2)
    #   NK
    ["CD57",            "E-cadherin"],
    #   functional / checkpoint (all sparse without insulitis, see (2))
    ["PD-1",            "E-cadherin"],
    ["LAG3",            "E-cadherin"],
    ["TOX",             "E-cadherin"],
    ["TCF-1",           "CD68"],        # myeloid counterpart: epithelium can be TCF7+ via Wnt
    ["ICOS",            "E-cadherin"],
    ["VISTA",           "E-cadherin"],
    ["CD39",            "E-cadherin"],
    ["Granzyme B",      "E-cadherin"],
    ["CD107a",          "E-cadherin"],  # degranulation
    ["HLA-DR",          "E-cadherin"],  # weakest pair in the set, see (1)
    ["PD-L1",           "E-cadherin"],  # inducible, see (1)
    ["IDO1",            "E-cadherin"],  # inducible, see (1)
    ["iNOS",            "E-cadherin"],  # inducible, see (1)
    # === ENDOTHELIAL (CD31+) ===
    ["Podoplanin",      "E-cadherin"],  # NOT CD31 -- lymphatic EC are CD31+
    ["CD34",            "E-cadherin"],
    # === MESENCHYMAL (Vimentin+) ===
    ["SMA",             "E-cadherin"],  # mural / myofibroblast
    #   Podoplanin and CD34 above are reused here -- see dual-use note
    # === STATE, not lineage ===
    ["CD44",            "E-cadherin"],
    ["CD99",            "E-cadherin"],
    ["Bcl-2",           "E-cadherin"],
    ["Caveolin",        "CD31"],        # corrected from E-cadherin, see (1)
    ["M2Gal3",          "CD31"],        # corrected from E-cadherin, see (1)
    ["CD66",            "CD31"],        # assumes CEACAM5/6; if CD66b, use "E-cadherin"
]
```

## Dual-use markers (one threshold, interpreted by parent node)

Three markers do double duty, disambiguated by **where the cell already landed**, not by a second
threshold:
- **CD56** — endocrine inside `E-cadherin⁺`; NK inside the immune union. (NK cells were removed by the
  immune gate and nerve is `E-cadherin⁻`, so `ECAD⁺/CD56⁺` is endocrine.)
- **Podoplanin** — lymphatic/capillary endothelium inside `CD31⁺`; fibroblast subset inside `Vimentin⁺`.
- **CD34** — capillary endothelium inside `CD31⁺`; resting-fibroblast subset inside `Vimentin⁺`.

## What resolves, and what stays residual

- **Epithelial** — materially weaker than a 35-plex here: no chromogranin/synaptophysin, so endocrine has
  no dedicated pan-marker, but **CD56 substitutes cleanly inside `E-cadherin⁺`** (its confounders are
  already removed). This matters: hormone-gating alone would misfile γ/PP and ε cells as exocrine (no
  PPY/ghrelin in the panel); CD56 recovers them as "endocrine, hormone-negative." **Exocrine stays
  unsplittable** (no KRT19/CFTR/SOX9/MUC1/acinar enzyme) — ductal and acinar collapse to one bin, with
  Keratin 5/TP63 marking only basal/metaplastic epithelium.
- **Immune** — much stronger: full T subsetting (CD4/CD8/FOXP3) with an 8-marker exhaustion/checkpoint
  layer, B+plasma via CD79a/CD38, 4-marker macrophage phenotyping, DC via CD11c/CD209, neutrophil via
  MPO, NK via CD56/CD57. **CD4 must be read only inside `CD3e⁺`** (monocytes are CD4-dim).
- **Endothelial** — splits on Podoplanin (lymphatic vs blood vascular), CD34 at capillary level.
- **Neural** — **terminal at level 1.** B3TUBB is the only neural marker (no GAP43/PGP9.5/S100/SOX10) →
  no subtyping, no sprouting-vs-mature axis, no Schwann call. Largest structural gap vs a 35-plex.
- **Mesenchymal** — splits on SMA (mural/pericyte/myofibroblast/activated stellate) against a
  `CD34⁺/SMA⁻` resting-fibroblast subset and a `Podoplanin⁺` inflammatory-CAF-like subset.

## Runtime issues specific to level 2

**(1) Four pairs are structurally compromised (two corrected from an earlier 59-plex list).** **M2Gal3**
and **Caveolin** are counterpart **CD31**, not E-cadherin — both are expressed by ductal epithelium, so an
E-cadherin-high negative would contain true positives and bias those thresholds high (same reasoning as
LGALS3 in the 35-plex). **HLA-DR** is contaminated the same way (ductal HLA-II induction in T1D) but has
no better abundant counterpart, so **flag its call as conservative**. **PD-L1, IDO1, iNOS** are inducible
on epithelium and myeloid alike — treat all of these as **state annotations, never class calls**.

**(2) Prevalence collapse is the dominant failure mode at level 2** and is worse than in a 35-plex,
because the checkpoint layer is only populated in inflamed tissue. RESTORE's 2-component fit will **split
the negative distribution rather than fail** when no positive mode exists. Affected: Keratin 5, TP63, SOX2
(absent from normal adult pancreas — basal/squamous/progenitor programs appear only in ADM/metaplasia);
FOXP3, Granzyme B, CD57, TOX, LAG3, PD-1, ICOS, VISTA, CD39, MPO, CD209 (absent from any scene without
insulitis); SST (δ cells are a small fraction of an already-small compartment). Compounding it, the
`counterpart > median` prefilter selects roughly half the tissue while these positive modes can sit under
0.1% of cells. **Estimate these on scenes known to contain positives and apply globally, or check
positive-cluster size before trusting any per-scene value.**

**(3) Acinar autofluorescence** biases every E-cadherin-paired threshold conservative (unchanged from
level 1), and it lands hardest on exactly the sparse immune targets in (2), where a conservatively-set cut
on a thin positive population is where real cells are lost. The islet-ROI `scene` mitigation applies if
running level 2 on islet-centered acquisitions rather than whole sections.

## Excluded — no valid counterpart (6)

- **DAPI**, **Beta-actin** — segmentation / reference, not lineage.
- **HLA-A** — MHC-I on all nucleated cells, hyperexpressed in T1D islets → no stable negative.
- **Collagen IV** — ECM mask, not a per-cell signal.
- **Ki67**, **PCNA** — lineage-agnostic; binarize by a standalone bimodal fit and apply as an **orthogonal
  proliferation annotation** across all five compartments.

## Relationship to the code

This maps onto the existing `marker_taxonomy` split: **level-2 lineage markers ⊂ `TYPE`** (minus the 7
level-1 gate markers), and the **state markers** here (`CD44/CD99/Bcl-2/Caveolin/M2Gal3/CD66` + the
inducible PD-L1/IDO1/iNOS/HLA-DR) belong with `PROCESS`. So wiring level 2 does not introduce a parallel
scheme — it extends the pair set (with validation) and adds a within-compartment interpretation layer to
`lineage.py` beneath the broad `compartment` call.
