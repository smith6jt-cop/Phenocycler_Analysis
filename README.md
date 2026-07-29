# Phenocycler / CODEX Analysis — REDSEA → RESTORE → Broad Lineage

A notebook- and Python-based pipeline for **Phenocycler (CODEX) spatial-proteomics**
analysis, from raw imaging output through **broad lineage classification**. It is a
faithful, optimized port of the "Senior phenotyping pipeline" from
[Islet-Explorer-Senior](https://github.com/smith6jt-cop/Islet-Explorer-Senior),
restructured into a reusable, testable package with per-donor parallelism and an
optional GPU backend.

The pipeline corrects three dominant artifacts of dense multiplexed tissue imaging —
**lateral marker spillover** (REDSEA), **per-image intensity / autofluorescence drift**
(RESTORE), and **false-endocrine over-calling** (a threshold-relative hormone floor) —
and then assigns every cell to one of **seven** mutually-exclusive lineages with
**zero "Unassigned"** cells.

```
 raw .qptiff  ─┐
               ├─(QuPath segmentation, prerequisite)─┐
 GeoJSON ──────┘                                     │
                                                     ▼
 Cellmeasurements.csv ──▶ [1] cells_parquet (DuckDB) ──▶ data/cells/
                                                     │
 qptiff + GeoJSON ─────────▶ [2] REDSEA spillover ───▶ data/cells_redsea/
                                                     │
                            [3] RESTORE normalize ───▶ data/restore_gated_redsea/
                                (SSC, per-image)      │   + restore_gated_redsea_extra/ (CD99/B3TUBB/MPO)
                                                     │
                            [4] hormone floor (K=5) ─▶ data/restore_gated_redsea/  (false-endocrine fix)
                                                     │
                            [5] broad lineage ────────▶ data/phenotype/broad/   (7 lineages, 0 Unassigned)
                                                     │
                            [6] QuPath export ────────▶ data/phenotype/qupath_class/*.csv  (optional QC)
                            [7] identity figures ─────▶ data/phenotype/celltype_marker_*.png
                                                     │
 Xenium run ──────────────▶ [8] INTEGRATION ─────────▶ data/integration/
 (10x spatial RNA)              (optional, separate env)   paired islets / niches / QC
```

## Scope

**In scope (this repo):** raw QuPath outputs → cells parquet → REDSEA → RESTORE →
broad lineage → optional QuPath re-import for visual QC.

**Also in this repo (optional, separate environment):** PhenoCycler ↔ Xenium integration —
`phenocycler/integration/`, see **[docs/INTEGRATION.md](docs/INTEGRATION.md)**. It starts
from `data/phenotype/broad/` and integrates it with 10x Xenium spatial transcriptomics in two
modes, `sequential` (serial sections) and `same_slide`. The core pipeline is unaffected: it
gains no dependencies and still runs in the lean environment.

**Out of scope / downstream:** image processing before segmentation (illumination
correction, stitching, registration, autofluorescence removal — use
[KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI)), scVI embedding, trajectory /
pseudotime, per-lineage subclustering, and the interactive R Shiny app.

## What REDSEA and RESTORE do

- **REDSEA** (Bai et al., *Front. Immunol.* 2021) removes signal that bleeds across
  shared cell boundaries. `phenocycler/redsea.py` is a scalable pixel-level
  reimplementation (the reference `redseapy` does not scale): per donor it rasterizes
  the QuPath GeoJSON to an int32 mask, sums each qptiff channel per cell and per 1-px
  boundary band with `np.bincount`, builds an 8-connected sparse contact graph, and
  applies `corrected = clip(data − α·(F @ edge))` (subtract-only, α=1, 1-px band by
  default).
- **RESTORE** (Chang et al., *Commun. Biol.* 2020) normalizes per-image intensity using
  mutually-exclusive reference markers. `phenocycler/restore.py` drives the **vendored**
  RESTORE (`external/RESTORE`, git submodule) headlessly, fitting KMeans/GMM/**SSC**
  per image × marker pair (SSC is the default), with a robust cohort-median guard on
  degenerate thresholds. Pan-Cytokeratin is the universal negative reference (except
  `CD3e ← CD163`).
- **Hormone floor** rewrites `{INS,GCG,SST}_pos = (_norm ≥ K)` (K=5) *before* typing:
  RESTORE's per-image hormone threshold lands in the noise for β-loss donors, so this
  threshold-relative floor rejects the noise-tail false positives while keeping real,
  separated β/α/δ cells (`hormone_floor.py`). CD99/B3TUBB/MPO are gated in a separate
  RESTORE pass (`restore --extra`) so the validated 10-marker gates stay byte-identical.
- **Broad lineage** assigns each cell by a strict hierarchy — Endocrine (INS/GCG/SST or
  bright CD99) → Immune (CD3e/CD20/CD163/MPO — neutrophils are immune) → Endothelial →
  Neural (B3TUBB) → structural argmax (Epithelial / Fibroblast / Muscle) → Epithelial
  default — producing **seven** lineages and no "Unassigned" bucket.

## Installation

```bash
git clone --recurse-submodules https://github.com/smith6jt-cop/Phenocycler_Analysis.git
cd Phenocycler_Analysis
bash setup.sh          # inits the RESTORE submodule, creates the conda env, verifies imports
conda activate phenocycler_analysis
pip install spams-bin  # SSC model for RESTORE (if not already pulled by environment.yml)
```

If you cloned without `--recurse-submodules`: `git submodule update --init --recursive`.

**Consuming this as a dependency (e.g. from Islet-Explorer-Senior).** The package is
config-driven, so it can be added as a git submodule and installed *editable* into an
existing analysis environment instead of its own conda env:

```bash
pip install -e .    # into your existing env (e.g. Senior's scvi-env, Python 3.13 / numpy 2.2)
```

`environment.yml` pins a self-contained reference env (Python 3.11); when embedding in
another project, run under that project's interpreter and point `config.ini` (or the
`PHENOCYCLER_*` env vars) at its data.

## Quick start

1. Edit `config.ini` `[paths]` to point at your raw `.qptiff` images, the QuPath
   `Cellmeasurements.csv`, and the donor-metadata workbook (see `DATA_README.md`).
2. Export per-cell GeoJSON from QuPath (`scripts/groovy/export_cells_geojson.groovy`).
3. Run the whole pipeline (idempotent — skips completed steps):

```bash
python -m phenocycler.pipeline            # cells → redsea → restore(+extra) → hormone_floor → lineage → qupath → figures
python -m phenocycler.pipeline --status   # just show what exists
```

Or run steps individually / interactively via the notebooks:

```
notebooks/00_prerequisites_qupath.ipynb   # upstream + GeoJSON export
notebooks/01_build_cells_parquet.ipynb
notebooks/02_redsea_spillover.ipynb
notebooks/03_restore_normalize.ipynb
notebooks/04_broad_lineage.ipynb
notebooks/05_qupath_export.ipynb
notebooks/00_run_full_pipeline.ipynb      # thin orchestrator

notebooks/06_integration_prep.ipynb       # PhenoCycler <-> Xenium (optional)
notebooks/07_sequential_registration.ipynb
notebooks/08_sequential_integration.ipynb
notebooks/09_same_slide_integration.ipynb
```

Every step is also a module CLI, e.g. `python -m phenocycler.redsea --all --jobs 4`.

## Performance & acceleration

Upstream, only the DuckDB build and scVI were parallel; REDSEA, RESTORE and lineage
ran strictly serially. **Every stage is per-donor embarrassingly parallel**, which is
the primary optimization here.

| Step | Engine | Parallelism | GPU |
|------|--------|-------------|-----|
| 1. cells_parquet | DuckDB (out-of-core) | multithreaded (`[compute] duckdb_threads`) | n/a (I/O bound) |
| 2. REDSEA | NumPy + SciPy-sparse | per-donor pool (`--jobs`) **or** SLURM array (one donor/task) | optional CuPy (`--gpu`): `bincount` + sparse `F@edge` |
| 3. RESTORE | vendored RESTORE (SSC) | per-donor apply pool (`--jobs`); subsample cap keeps SSC tractable | — |
| 4. hormone floor | pandas | per-donor pool (`--jobs`) | — |
| 5. broad lineage | vectorized NumPy | per-donor pool (`--jobs`) | — |
| 6. QuPath export | pandas | serial (I/O-light) | — |
| 7. identity figures | pandas + matplotlib | streaming | — |

Set `[compute] n_jobs` in `config.ini` (or `--jobs N`) for the process pools, and
`use_gpu = true` (+ `pip install cupy-cuda12x`) for the REDSEA CuPy path. For
cohort-scale HPC runs use the SLURM chain:

```bash
bash scripts/slurm/run_full_pipeline.sh 40   # 40 SECTIONS: cells → redsea[array 0-39] → restore(+extra) → floor/lineage/qupath/figures
```

The array is sized in **sections**, not donors — the work unit is one image, so a donor with
a pancreas and a lymph node contributes two (`6539` and `6539pln`). Full setup, the
positivity-refinement loop, and the silent-failure checklist:
**[docs/HIPERGATOR.md](docs/HIPERGATOR.md)**.

## Repository structure

```
phenocycler/            installable package
  config.py             config.ini loader + derived paths (no hardcoded paths)
  cells_parquet.py      Step 1 — DuckDB CSV -> partitioned parquet
  redsea.py             Step 2 — pixel REDSEA (+ parallel driver + CuPy backend)
  restore.py            Step 3 — RESTORE driver (headless stubs, robust LUT, parallel apply)
  hormone_floor.py      Step 4 — threshold-relative endocrine floor (false-endocrine fix)
  lineage.py            Step 5 — deterministic 7-class broad-lineage hierarchy
  marker_taxonomy.py    TYPE / PROCESS / EXCLUDED marker split (single source of truth)
  qupath_export.py      Step 6 — per-image UUID-keyed CSVs for QuPath
  figures.py            Step 7 — cell-type × marker identity dotplot + heatmap
  reassess_diag.py      read-only acceptance yardstick (endocrine + guardrail metrics)
  pipeline.py           idempotent orchestrator + status table
  parallel.py           per-donor process pool
  gpu.py                optional CuPy/RAPIDS backend with CPU fallback
phenocycler/integration/  PhenoCycler <-> Xenium integration (optional extra)
  contract.py           modality-agnostic cell-table schema
  manifest.py           S0  donor <-> Xenium-run pairing (+ donor_id stamping)
  export_pheno.py       S1a phenotype/broad + cells + gates -> contract
  import_xenium.py      S1b h5ad | SpatialData zarr | raw bundle -> contract
  vocab.py              S2  harmonised lineages + protein<->gene crosswalk
  structures.py         S3  islets / ducts / vessels (hormone-seeded DBSCAN)
  rasterize.py          S3b common micron grid; qptiff + morphology_focus readers
  register.py           S4  rigid -> affine -> B-spline; image + point-set tracks
  transform.py          S5  micron-space affine + displacement field
  match.py              S6  islet <-> islet Hungarian assignment [sequential]
  grid.py               S7  joint niches + registered hex grid
  donor.py              S8  donor-level concordance (no registration)
  crossmodal.py         S9  pseudo-cell linking [sequential, inference]
  sameslide.py          S11 cell <-> cell pairing [same_slide]
  qc.py / figures.py    S10 gates, nulls, report, overlays
  pipeline.py           idempotent orchestrator (--mode, --status)
external/RESTORE/       vendored RESTORE (git submodule, commit 38df59b)
external/XeniumPanelExplorer/  vendored panel taxonomy (submodule, commit 076d1192)
notebooks/              thin step notebooks + orchestrator
scripts/groovy/         QuPath export/import Groovy scripts
scripts/slurm/          HiPerGator SLURM (per-donor array) scripts
tests/                  pytest: REDSEA math, RESTORE guard, 7-class lineage invariants
config.ini              paths + tunables (+ [integration])
environment-integration.yml  union env for the integration extra
docs/INTEGRATION.md     integration design + data caveats
pyproject.toml          editable-install metadata (pip install -e .)
```

## Testing

```bash
pytest tests/            # REDSEA compensation, contact matrix, rasterize; RESTORE robust
                         # guard; lineage hierarchy (zero Unassigned); config; parallel;
                         # integration: contract, manifest crosswalk, vocabulary totality,
                         # structure calling, registration recovery, matching, mode guards
```

The tests need no imaging data (they use synthetic masks/frames) and run in ~15 s. The
integration tests apply *known* transforms to synthetic sections and assert that
registration recovers them — including the mirrored case and the non-rigid displacement
sign, which fails silently by doubling the misalignment rather than removing it.

## Prerequisites (upstream of this repo)

1. **Image processing** — [KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI)
   (illumination correction → stitching → deconvolution → EDF → registration →
   autofluorescence removal).
2. **Segmentation** — QuPath cell detection, then export the per-cell measurement CSV
   and full-resolution GeoJSON (`scripts/groovy/export_cells_geojson.groovy`).

## Citations

- **REDSEA** — Bai, Y. et al. *Adjacent Cell Marker Lateral Spillover Compensation and
  Reinforcement for Multiplexed Images.* Front. Immunol. 12:652631 (2021).
- **RESTORE** — Chang, Y.H. et al. *RESTORE: Robust intEnSiTy nORmalization mEthod for
  multiplexed imaging.* Commun. Biol. 3:111 (2020).
- **KINTSUGI** — Smith, J.A. et al. *Protocol for processing and analyzing multiplexed
  images...* STAR Protocols 6:103976 (2025).

## Acknowledgments

University of Florida Research Computing (HiPerGator). Pipeline logic ported from
`smith6jt-cop/Islet-Explorer-Senior`; RESTORE vendored from `smith6jt-cop/RESTORE`.
