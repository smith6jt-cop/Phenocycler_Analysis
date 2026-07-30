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

### Section keys on the PhenoCycler side

Both of a donor's images are named for the same donor:

```
6539_Scan1.er.qptiff        pancreas
6539pLN_Scan1.er.qptiff     pancreatic lymph node
```

The core pipeline used to key partitions on `regexp_extract("Image", '[0-9]+')` — the first
digit-run — so both landed in `data/cells/donor_id=6539`. Two sections in one partition is
not a labelling nuisance. Every core stage treats a partition as one image: `redsea.py` reads
`image.iloc[0]` and masks *every* cell in the partition with that one section's GeoJSON,
`restore.py` fits one threshold set per partition, and each section's micron coordinates are
measured from its own slide origin. A shared partition would mask one section with the
other's polygons, normalise two images as one, and interleave two coordinate frames into a
single overlapping point cloud — with nothing raised anywhere.

So the partition key is the **section key**: the donor's digits plus whatever region token
the filename carries, lowercased.

| image | partition | donor | roi |
|---|---|---|---|
| `6539_Scan1.er.qptiff` | `6539` | 6539 | `panc` |
| `6539pLN_Scan1.er.qptiff` | `6539pln` | 6539 | `pln_1` |
| `6539pLN2_Scan1.er.qptiff` | `6539pln2` | 6539 | `pln_2` |

`phenocycler/sections.py` is the single definition, shared by the DuckDB expression that
*builds* the partitions and the Python parser that *takes them apart* — two languages, one
rule, checked against each other in the test suite. An unrecognised token like `6539spleen` is
refused rather than folded into the pancreas, and `redsea.donor_image` asserts its partition
holds exactly one image as a backstop.

A bare `pLN` resolves to `pln_1`. That is an **inference, not a measurement**: `pln_1` is the
region every row of the vendored `xenium_paths.csv` carries (26 of the 26 donors listed there;
three have a second, one a third). But that CSV is the upstream manifest, not the run cohort —
the analysis cohort is 20 donors and **which 20 is not recorded anywhere in this repo**, so the
statistic is drawn from a superset. The mapping is wrong for any donor whose single PhenoCycler
lymph-node section corresponds to that donor's `pln_2` or `pln_3`. It would surface as a
registration that will not align, not as silently wrong output; the fix is a per-donor override
in `sections.REGION_TO_ROI`.

### The post-Xenium re-stain — protein and RNA on the same cells

After the Xenium run, that section is stained for INS and GCG **on the PhenoCycler**. The
output is an ordinary qptiff, so it flows through the whole core pipeline unchanged. What
differs is the slide underneath — it images the *Xenium* section:

```
PhenoCycler section  ──── serial ────▶  Xenium section
 (55-plex protein)                       (5K RNA)
                                             │
                                        same slide
                                             ▼
                                   post-Xenium re-stain
                                    (INS / GCG protein)
```

This is the only place in the cohort where protein and RNA are measured on **one cell**. The
section key therefore carries a second dimension — which physical slide was imaged — so a
re-stain never lands in the PhenoCycler section's partition:

| image | partition | donor | roi | slide |
|---|---|---|---|---|
| `6539_Scan1.er.qptiff` | `6539` | 6539 | `panc` | pheno |
| `6539xen_Scan1.er.qptiff` | `6539xen` | 6539 | `panc` | **xenium** |
| `6539xenpLN_Scan1.er.qptiff` | `6539xenpln` | 6539 | `pln_1` | **xenium** |

Re-stains export to `cells_pheno_xif/` rather than `cells_pheno/` — same `(donor, roi)`, but
a different piece of tissue, and one partition must stay one section.

`postxen` (S1c) pairs them against the transcriptome. It runs in **both** modes: the
correspondence is a property of that section pair, not of the cohort, so gating it on
`[integration] mode` would answer the wrong question. No registration is applied — a residual
translation between two exports of one section is a coordinate-frame difference, estimated as
a median offset and removed.

**What this unlocks.** INS/GCG/SST are off the 5K panel, so Xenium normally calls Beta/Alpha
from surrogate cores (PDX1, ISL1, NEUROD1, ABCC8). With the re-stain those calls are *scored
against hormone protein on the same cell* — `qc/postxen_qc.csv` reports the agreement and
per-subtype recall. The hormone call is then stamped onto the Xenium table, so S6b types both
sides of the serial pair by protein rather than protein-versus-surrogate-RNA.

**One caveat, enforced by scope.** Post-Xenium tissue has been through protease digestion and
decrosslinking, so epitopes are degraded and intensities are *not* comparable to a fresh
PhenoCycler stain. Only identity travels (`CARRY_COLUMNS`), never raw marker values. CD3e is
excluded on the same grounds — a lower-abundance surface marker does not survive that
chemistry, which is what the reported unreliability predicts.

If your re-stains use a filename token other than `xen` / `xpanc` / `xenpln` / `xpln`, `parse`
**refuses** it and names the file; add one line to `sections.REGION_TO_ROI`. That is
deliberately louder than defaulting, which would merge two sections.

This also makes pairing exact rather than assumed. The manifest now knows *which* PhenoCycler
sections exist, so a donor whose lymph node was never scanned is `xenium_only` for that ROI
instead of being reported as a pairing that does not exist.

Two consequences worth knowing:

- The core pipeline's `donor_id=*` partitions are **sections**, so a 20-donor two-tissue
  cohort has ~40 of them. `cfg.discover_sections()` is the work-unit iterator every core
  stage uses; `cfg.discover_donors()` returns the 20 unique donors, which is what the manifest
  and the donor workbook join want. Conflating the two makes `6539pln` a phantom donor.
- `[integration] pheno_tissue_dirs` survives only as an override for a tissue processed
  through a *separate* core-pipeline run. The normal case needs one run and no config.

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

`pair_donor` takes a `method`:

| `method` | criterion |
|---|---|
| `iou` | **polygon overlap** — two segmentations of one cell overlap substantially, two neighbours do not |
| `centroid` | mutual-nearest centroid within a 5 µm cap |
| `auto` (default) | IoU, falling back to centroids when polygons are unavailable |

IoU is preferred wherever it can run. In packed islet tissue an adjacent nucleus is often
nearer than the partner outline's centroid, which is precisely the case centroid distance
gets wrong and overlap does not.

Polygons come from `boundaries.py`, which reads QuPath GeoJSON (`data/redsea_scratch/geojson/`)
and Xenium `cell_boundaries.parquet`, and warps the moving side through the transform
**vertex-wise** so the non-rigid field is applied, not just its affine part.

The two formats disagree about units — the GeoJSON is in full-resolution qptiff **pixels**,
every contract table is in **microns** — and getting that wrong does not error, it puts one
section's polygons a thousand microns from their own cells. Rather than read the qptiff's
resolution tags (which needs the image mounted and returns nothing when the tags are absent,
so a missing scale silently becomes 1.0), the conversion is **fitted from the data**: both
files describe the same cells under the same ids, so a per-axis linear fit of polygon centroid
against `x_um`/`y_um` recovers it exactly. A fitted scale outside 0.05–2.0 µm/px, or a median
residual above 5 µm, means the two files are not the same section and is an error rather than
a silently wrong overlay.

`auto` records which matcher actually ran in `stats['method']`. A silent fallback would make a
polygon-loading failure look like a successful centroid run, and the two have different error
characteristics. `method='iou'` re-raises instead of falling back, so an explicit request
cannot be quietly downgraded.

The same matching is used by `postxen.py` (S1c) for the post-Xenium re-stain, where it matters
most: that really is one section imaged twice, so overlap is the correct criterion by
construction.

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
verbatim, with the same columns and the same nPOD donors (26 appear in the vendored
`xenium_paths.csv`; the run cohort is a 20-donor subset of them). But:

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
