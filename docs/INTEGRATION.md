# PhenoCycler ↔ Xenium integration

Picks up where `python -m phenocycler.pipeline` stops — at `data/phenotype/broad/`, eight
broad lineages, zero Unassigned — and integrates that with 10x Xenium spatial
transcriptomics processed by [Xenium_Analysis](https://github.com/smith6jt-cop/Xenium_Analysis),
using the curated panel taxonomy from
[XeniumPanelExplorer](https://github.com/smith6jt-cop/XeniumPanelExplorer).

```bash
conda env create -f environment-integration.yml
conda activate phenocycler_integration
pip install -e .

python -m phenocycler.integration.manifest --report      # what pairs with what
python -m phenocycler.integration.pipeline --status      # stage table
python -m phenocycler.integration.pipeline               # run everything, both tissues
python -m phenocycler.integration.pipeline --tissue pancreas    # one tissue
```

---

## Tissue

The cohort is **20 donors × {PhenoCycler, Xenium} × {pancreas, pancreatic lymph node}**.

Tissue is a property of the *dataset*, fixed when it is imported, not a mode of the analysis.
It is derived from the ROI (`panc` → pancreas, `pln_1..3` → pancreatic lymph node), written
into the cell table as a column at creation, and validated: a section whose ROI maps to no
known tissue is refused rather than defaulting to pancreas. Every stage therefore runs on one
tissue or on both, and a bare run covers both:

```bash
python -m phenocycler.integration.pipeline                                # both tissues
python -m phenocycler.integration.pipeline --tissue pancreatic_lymph_node # one tissue
python -m phenocycler.integration.pipeline --roi pln_2                    # one section
```

`--roi` beats `--tissue`; with neither, the scope is every ROI of every tissue in
`[integration] tissues`. Stage runners take that scope and iterate the manifest themselves —
they are *not* called once per ROI, because their output paths are stage-level
(`niche_profiles.csv`, `islets.parquet`) and a per-ROI call would have each ROI overwrite the
last.

The two tissues diverge in exactly two places, both biological. Everything else — export,
import, vocab, registration, niches, the grid, donor composition, pseudo-cell linking — is
identical.

**1. Structures.** Islets exist in pancreas. A lymph node has none, so hormone-seeded DBSCAN
has nothing to seed on and islet matching (S6) has nothing to match. The structure level is
simply not run there; `--status` marks it `n/a`, and a combined run still matches islets in
the pancreas half. B-cell follicles are the obvious substitute landmark and the code would
support them (`call_landmarks` is generic over a seed lineage), but a follicle called from a
5K RNA panel on one side and a 55-plex protein panel on the other is a far weaker
correspondence than an islet, and shipping it would invite structure-level claims the data
cannot support. `tissues.STRUCTURE_SEEDS` is where it would go.

This propagates into QC. The islet-RMSE gate is *undefined* without islets, not missing, so
for those tissues the verdict rests on tissue Dice alone. Counting the absent metric as a
failed gate would cap every lymph-node donor at `WARN` forever and permanently withhold
`crossmodal` from half the cohort — an artefact of the gate design, not a property of the
registration.

**2. Identity panels.** XeniumPanelExplorer ships `pancreas` and `thymus`, no lymph node. But
the cores a lymph node needs are exactly the *shared* ones — T cell, B/plasma, NK/ILC,
myeloid, granulocyte, endothelial, pericyte, neural, epithelial — and the three
pancreas-specific cores (endocrine, exocrine, fibroblast/stellate) do not apply. So pLN
resolves to the shared pool: 8 cores against the pancreas's 11. That is the correct answer,
not a degraded fallback. Scoring a lymph-node cell against `01_pancreas_endocrine` would not
be harmlessly uninformative — the broad call is an argmax over those axes, so an axis that
cannot exist in the tissue is a competitor that can win.

The same logic narrows donor-level concordance: `Endocrine` and `Exocrine` are dropped from
the comparable-lineage set for pLN, because a near-zero fraction on both sides is two
instruments correctly finding nothing, and letting it into the Spearman would inflate the
concordance with a free correct answer.

`phenocycler/integration/tissues.py` is the single place that knows any of this.

### One caveat on the PhenoCycler side

The core pipeline derives `donor_id` as the first digit-run of the qptiff name
(`cells_parquet.py`), so donor 6539's pancreas image and its lymph-node image both land in
`data/cells/donor_id=6539` and cannot be told apart afterwards. Each tissue therefore needs
its own core-pipeline data dir:

```ini
[integration]
pheno_tissue_dirs = pancreas=/path/panc/data,pancreatic_lymph_node=/path/pln/data
```

Left unset, both tissues resolve to the same directory and `export_pheno` exports the first
and **refuses** the rest with an actionable message. That refusal is the point: exporting one
source under two ROIs would write a pancreas section into the lymph-node partition — real
cells, correctly labelled, attributed to the wrong organ, which every downstream comparison
would take at face value. The Xenium side has no such problem; its sections are keyed
`{serial}__{roi}` and one slide serial legitimately carries both.

---

## The two modes, and why the distinction is enforced

```ini
[integration]
mode = sequential      # or same_slide
```

**`sequential`** — serial sections cut from one FFPE block. The two modalities measure
*different cells*. No registration, however good, changes that. Integration therefore
happens at four levels, none of which claims a cell-to-cell correspondence:

| level | pairing unit | needs registration? |
|---|---|---|
| donor | `donor_id` | no |
| niche | kNN composition class | no |
| structure | matched islet | yes |
| pseudo-cell | inferred link | yes, and QC-gated |

**`same_slide`** — one physical section imaged twice. Cells *are* the same cells, so a
genuine paired protein+RNA matrix is recoverable and `sameslide.py` produces it.

`pair_donor` assigns by **mutual-nearest centroid** within a 5 µm cap. `match_by_iou` does
the better polygon-overlap assignment and is callable if you can supply GeoSeries, but is not
wired into `pair_donor`: the missing piece is a loader for PhenoCycler's GeoJSON and Xenium's
zarr boundaries, warped through the transform. This cohort is serial-section, so there is no
same-slide data to validate such a loader against — it is a follow-up, not a gap being
papered over.

Every mode-specific module refuses to run in the wrong mode, as a hard error rather than a
warning. The validity of every downstream claim rests on which situation holds, and a
warning is too easy to scroll past.

---

## Pipeline

```
manifest      S0   donor ↔ Xenium-run pairing; stamps donor_id/roi into the data
export_pheno  S1a  phenotype/broad + cells + restore_gated  →  cells_pheno/
import_xenium S1b  h5ad | SpatialData zarr | raw bundle     →  cells_xen/
vocab         S2   harmonised lineage vocabulary + protein↔gene crosswalk
structures    S3   tissue / islets / ducts / vessels, same algorithm both sides [pancreas]
register      S4   coarse rigid → affine → non-rigid; image + point-set tracks
transform     S5   warp coordinates and masks into the fixed frame
match         S6   islet ↔ islet assignment          [sequential, pancreas]
grid          S7   joint niches (no registration) + registered hex grid
donor         S8   donor-level concordance                     [no registration]
crossmodal    S9   pseudo-cell linking                         [sequential, INFERENCE]
sameslide     S11  cell ↔ cell pairing                         [same_slide]
qc            S10  gates, nulls, report
figures       S10  overlays and concordance plots
```

Everything after S1 speaks one modality-agnostic parquet schema (`contract.py`), which is
what lets registration, matching and QC be written once for both modalities and both modes —
and tested on synthetic frames with no imaging data at all.

---

## What the data actually looks like

### The join key had to be built, not just read

`donor_id` is a clean key on paper: both repos already share `donor_metadata_panc.xlsx`
verbatim, same columns, same 26 nPOD donors. But:

- **No Xenium artifact carries it.** Every Xenium output is keyed by the XETG slide serial
  (`0041323`), which encodes no donor, no region and no disease status.
- **`xenium_paths.csv` — the one file that maps donors to runs — is read by nothing** in the
  Xenium repo. It is hand-maintained and unenforced. It is vendored here at
  `data/integration/xenium_paths.csv` and made load-bearing.
- **The processed samples are not in it.** `0041323` and `0041326` have finished h5ads and
  appear in no manifest, so their donor is recorded nowhere. They enter the manifest as
  `pair_status = donor_unknown` and are excluded from paired analyses until someone fills in
  `data/integration/donor_overrides.csv`. This is visible in `--report`, not silent.
- **A slide serial is not a section.** One serial can carry both a pancreas and a
  pancreatic-lymph-node region (donor 6539 → `0059865` with `__Panc__` *and* `__pLN__`
  bundles), so `data/raw/{serial}.zarr` collides. The section key here is
  `{serial}__{roi}`.

`donor_id` is normalised to a string everywhere — the workbook stores it numerically,
`xenium_paths.csv` as text, PhenoCycler derives it by regex. Left alone, `6539` and `6539.0`
are different keys and the merge yields zero pairs.

### Two eight-class vocabularies that are not the same partition

```
PhenoCycler   Epithelial  Fibroblast  Muscle  Immune              Endocrine  Endothelial  Neural  Neutrophil
Xenium        Exocrine    Stromal ──────────  Immune_T/B/Myeloid  Endocrine  Vascular     Neural  (in Myeloid)
```

- **Immune.** PhenoCycler collapses T/B/myeloid into one class, but `CD3e_pos`, `CD20_pos`
  and `CD163_pos` survive in `restore_gated_redsea/` — so `export_pheno` recovers the split
  as `immune_subclass` with no re-running of anything.
- **Fibroblast vs Muscle.** Xenium cannot split them: `ACTA2` is absent from the 5K
  reference (it arrives only via the custom-100 add-on) and `12_pericyte_smooth_muscle`
  contributes exactly one SMC-hinted gene. The common vocabulary keeps `Stromal`.
- **Neutrophil.** PhenoCycler has an MPO gate; Xenium folds granulocytes into `Myeloid`.
  Reachable on the Xenium side only as a fine label or via the granulocyte identity core.
- **Epithelial.** Xenium deliberately has no epithelial mask — `16_epithelial_general` is
  `embedding_extra`, because KRT/EPCAM are shared by ductal, acinar and endocrine cells.

The mapping is *derived* from the pinned XeniumPanelExplorer submodule
(`external/XeniumPanelExplorer` @ `076d1192`), which ships `panel_roles.csv` and
`identity_core/*.csv` specifically as a machine-readable contract for Python consumers. A
drift guard fails loudly if an upstream lineage stops being covered.

### The `evidence_positive` correction

PhenoCycler's `Epithelial` is a **default sink**: `lineage.py` assigns it to every cell that
fails all the positive gates, flagged `epi_default`. Xenium has no sink — `broad_lowsignal`
cells still land in their best-scoring lineage.

Comparing raw fractions therefore measures a methodological difference, not biology. Every
composition figure reports `frac_all`, `frac_evidence_positive` and `frac_fallback`, and the
headline concordance uses the evidence-positive number on both sides.

### Off-panel markers

46 of 55 PhenoCycler markers have a gene on the deployed 5,081-gene panel. But **INS, GCG,
SST and Vimentin are off-panel** — precisely the markers that define PhenoCycler's Endocrine
and Fibroblast classes. Those comparisons route through Xenium's curated surrogate identity
panels (PDX1, ISL1, NEUROD1, ABCC8 … for endocrine) instead of 1:1 pairs.

That is orthogonal evidence rather than circular confirmation — scientifically a feature —
but it thins the anchor set available for pseudo-cell linking, so `crossmodal.py` refuses to
run below `crossmodal_min_anchors` rather than producing links from three genes.

---

## Registration

Two tracks, both run, best-by-QC wins, choice recorded.

**Track A — point set.** RANSAC over candidate islet correspondences (gated on size and area
ratio before geometry) → ICP over mutual nearest neighbours. Needs no images.

**Track B — image.** Coarse orientation search → opencv ECC (Euclidean then affine) →
SimpleITK B-spline. Uses DAPI↔DAPI where available, otherwise cell-density rasters.

Three things that matter:

- **Rotation and mirroring are not optional.** Sections are mounted by hand at arbitrary
  rotation and are sometimes flipped. The coarse stage searches rotation exhaustively and
  evaluates the mirror hypothesis explicitly; `fit_rigid_umeyama(allow_reflection=True)` can
  represent the flip, which textbook Umeyama cannot.
- **Non-rigid deformation is real but bounded.** The B-spline is coarse (control points
  hundreds of microns apart) and hard-capped at `reg_max_disp_um`. Given enough freedom a
  B-spline will align *any* two images, and a "successful" registration of unrelated
  sections is worse than an obvious failure.
- **Resolution is managed explicitly.** The coarse search runs at ~25 µm/px, ECC and the
  B-spline are capped by pixel count, and QC runs at 10 µm. Because `Transform` holds its
  affine in **microns**, a transform estimated at any resolution drops into the next stage
  unchanged. Registering a 10–20 mm section at the 2 µm working resolution would otherwise
  be intractable.

**Reading morphology images.** The Xenium ingest sets `morphology_focus=False`,
`morphology_mip=False`, `aligned_images=False` ("skip the 35GB TIFFs"), so the zarr does not
contain them. `rasterize.read_xenium_morphology()` reads `morphology_focus/` directly from
the original bundle, at a pyramid level near the target resolution.

### QC gates

| verdict | condition | permits |
|---|---|---|
| PASS | tissue Dice ≥ 0.80 **and** islet RMSE ≤ 150 µm | everything |
| WARN | one gate met | all but pseudo-cell linking |
| FAIL | neither | donor-level and niche only |

In a tissue with no structural unit the islet-RMSE gate does not exist, so the verdict is
PASS on Dice alone and FAIL otherwise — there is no WARN band, because there is no second
gate to half-meet. The dropped gate is stated in the `reasons` column rather than being
silently absent.

Verdicts are per **(donor, ROI)**, not per donor: a donor's pancreas can register cleanly
while its lymph node does not, and `crossmodal` gates on the section it is about to link.

A FAIL donor is **not** dropped from the study — only from the spatial analyses. Losing good
donor-level data because two sections would not align would be the wrong trade.

---

## Outputs

```
data/integration/
├── manifest.csv                     pairing table + pair_status + tissue [tracked]
├── vocab_crosswalk.csv              lineage + marker crosswalk, all tissues [tracked]
├── vocab_crosswalk_<tissue>.csv     per-tissue crosswalk (cores differ)  [tracked]
├── xenium_paths.csv                 vendored donor → run map           [tracked]
├── donor_overrides.csv              donor ids not recorded upstream    [tracked]
├── cells_pheno/donor_id=*/roi=*/    contract cell tables
├── cells_xen/donor_id=*/roi=*/
├── structures/{phenocycler,xenium}/donor_id=*/roi=*/{islets,ducts,vessels}.parquet
├── registration/donor_id=*/roi=*/   transform.json, displacement.npz, qc.json
├── paired/islets.parquet            matched islet pairs      [sequential]
├── paired/cells.parquet             paired protein+RNA cells [same_slide]
├── niches/                          niche_profiles.csv, grid_bins.parquet
├── crossmodal/pseudocell_links.parquet   INFERRED links
├── qc/                              registration_qc.csv, composition.csv,
│                                    disease_trend.csv, integration_report.txt
├── figures/                         registration overlays, concordance, trends
└── export/                          QuPath + Xenium Explorer round-trips
```

The Xenium Explorer export uses the **`cell_id`-matched CSV**, never the zarr path: the zarr
writer is position-indexed and the Xenium pipeline drops ~7% of cells at QC, so positions no
longer align with `cells.parquet`.

---

## Verifying it works

```bash
pytest tests/                                        # no data needed
python -m phenocycler.integration.pipeline --status
python -m phenocycler.integration.pipeline --status --tissue pancreatic_lymph_node
python -m phenocycler.integration.manifest --report
python -m phenocycler.integration.pipeline --donor 6539 --mode sequential
```

The test suite asserts recovery of *known* transforms — a synthetic warp is applied and the
registration must recover it, including the mirrored case and the non-rigid displacement
sign (which fails silently by doubling the misalignment rather than removing it).

**Scientific acceptance**, in order of how much they demand:

1. Donor composition should correlate across modalities (evidence-positive fractions).
2. The ND → AAB → T1D β-loss gradient should appear in **both** modalities independently.
3. On matched islets, protein INS⁺ fraction should track Xenium Beta fraction, and insulitis
   grades should agree.
4. Permuted-pairing and rotation nulls must destroy all of the above.

## Notes

- The core pipeline is untouched: `python -m phenocycler.pipeline` still runs its eight
  stages in the lean `phenocycler_analysis` environment with no new dependencies.
- `Xenium_Analysis/notebooks/06_xenium_phenocycler_integration.ipynb` is superseded. It
  applied a 40×/20× "magnification" scale factor to coordinates that are already in microns,
  synthesised a disc segmentation instead of using the real polygons in the zarr, registered
  with translation only, and was never wired to the real cohort.
- `Xenium_Analysis` hardcodes `ROOT = /blue/maigan/smith6jt/Xenium_Analysis` in 11+ files, so
  nothing here imports it. Its *artifacts* are read; the handful of shared constants (DBSCAN
  parameters, insulitis thresholds) are vendored with provenance comments.
