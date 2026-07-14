# Phenocycler / CODEX Analysis — REDSEA → RESTORE → Broad Lineage

A notebook- and Python-based pipeline for **Phenocycler (CODEX) spatial-proteomics**
analysis, from raw imaging output through **broad lineage classification**. It is a
faithful, optimized port of the "Senior phenotyping pipeline" from
[Islet-Explorer-Senior](https://github.com/smith6jt-cop/Islet-Explorer-Senior),
restructured into a reusable, testable package with per-donor parallelism and an
optional GPU backend.

The pipeline corrects the two dominant artifacts of dense multiplexed tissue imaging —
**lateral marker spillover** (REDSEA) and **per-image intensity / autofluorescence drift**
(RESTORE, over a curated web of directional mutually-exclusive marker pairs) — and then
assigns every cell to one of **5 compartments** by an ordered residual gating tree, with an
explicit **`Other`** terminal for cells that fail every gate (real panel gaps, not force-assigned).

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
                                (SSC, per-image;      │   {m}_pos / _norm / _log2r — ONE pass, all markers
                                 curated pairs)       │
                                                     │
                            [4] broad phenotyping ────▶ data/phenotype/broad/  (5 compartments + Other;
                                (ordered residual)    │                         compartment + cell_type)
                                                     │
                            [5] QuPath export ────────▶ data/phenotype/qupath_class/*.csv  (optional QC)
                            [6] identity figures ─────▶ data/phenotype/celltype_marker_*.png
```

## Scope

**In scope (this repo):** raw QuPath outputs → cells parquet → REDSEA → RESTORE →
broad lineage → optional QuPath re-import for visual QC.

**Out of scope / downstream:** image processing before segmentation (illumination
correction, stitching, registration, autofluorescence removal — use
[KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI)), scVI embedding, trajectory /
pseudotime, islet aggregation, spatial neighborhood analysis, per-lineage
subclustering, and the interactive R Shiny app.

## What REDSEA and RESTORE do

- **REDSEA** (Bai et al., *Front. Immunol.* 2021) removes signal that bleeds across
  shared cell boundaries. `phenocycler/redsea.py` is a scalable pixel-level
  reimplementation (the reference `redseapy` does not scale): per donor it rasterizes
  the QuPath GeoJSON to an int32 mask, sums each qptiff channel per cell and per 1-px
  boundary band with `np.bincount`, builds an 8-connected sparse contact graph, and
  applies `corrected = clip(data − α·(F @ edge))` (subtract-only, α=1, 1-px band by
  default).
- **RESTORE** (Chang et al., *Commun. Biol.* 2020) normalizes per-image intensity using
  mutually-exclusive marker pairs. `phenocycler/restore.py` drives the **vendored**
  RESTORE (`external/RESTORE`, git submodule) headlessly, fitting KMeans/GMM/**SSC**
  per image × pair (SSC is the default), with a robust cohort-median guard on degenerate
  thresholds. The pairs are a **curated directional web** (`MARKER_PAIRS` in `config.py`):
  each `[target, counterpart]` thresholds the target against a biologically mutually-exclusive
  counterpart (e.g. `E_cadherin ← CD31`, `B3TUBB ← CD3e`, intra-islet `INS ↔ GCG`) — not a
  single universal reference. All markers are thresholded in ONE pass.
- **Broad phenotyping** (`lineage.py`) assigns each cell by an **ordered residual gating tree**:
  Epithelial (E-cadherin) → Endothelial (CD31) → Neural (B3TUBB) → Immune (marker union + NK) →
  Mesenchymal (SMA then Vimentin), each gate on the residual of the prior; cells failing every
  gate → **`Other`**. Endocrine (β/α/δ) and Exocrine are the hormone⁺/⁻ sub-branches of
  Epithelial. Output: `compartment` (5 + Other) + finer `cell_type`.

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
python -m phenocycler.pipeline            # cells → redsea → restore → lineage → qupath → figures
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
| 4. broad phenotyping | vectorized NumPy | per-donor pool (`--jobs`) | — |
| 5. QuPath export | pandas | serial (I/O-light) | — |
| 6. identity figures | pandas + matplotlib | streaming | — |

Set `[compute] n_jobs` in `config.ini` (or `--jobs N`) for the process pools, and
`use_gpu = true` (+ `pip install cupy-cuda12x`) for the REDSEA CuPy path. For
cohort-scale HPC runs use the SLURM chain:

```bash
bash scripts/slurm/run_full_pipeline.sh 15   # 15 donors: cells → redsea[array 0-14] → restore → lineage
```

## Repository structure

```
phenocycler/            installable package
  config.py             config.ini loader + derived paths (no hardcoded paths)
  cells_parquet.py      Step 1 — DuckDB CSV -> partitioned parquet
  redsea.py             Step 2 — pixel REDSEA (+ parallel driver + CuPy backend)
  restore.py            Step 3 — RESTORE driver (headless stubs, robust LUT, parallel apply)
  lineage.py            Step 4 — ordered residual gating tree (5 compartments + Other)
  marker_taxonomy.py    TYPE / PROCESS / EXCLUDED marker split (single source of truth)
  qupath_export.py      Step 5 — per-image UUID-keyed CSVs for QuPath
  figures.py            Step 6 — cell-type × marker identity dotplot + heatmap
  pipeline.py           idempotent orchestrator + status table
  parallel.py           per-donor process pool
  gpu.py                optional CuPy/RAPIDS backend with CPU fallback
external/RESTORE/       vendored RESTORE (git submodule, commit 38df59b)
notebooks/              thin step notebooks + orchestrator
scripts/groovy/         QuPath export/import Groovy scripts
scripts/R/              RESTORE efficacy diagnostic (mxnorm) + its installer (R 4.5.1)
scripts/slurm/          HiPerGator SLURM (per-donor array) scripts
tests/                  pytest: REDSEA math, RESTORE guard, ordered-residual lineage invariants
config.ini              paths + tunables
pyproject.toml          editable-install metadata (pip install -e .)
```

## Testing

```bash
pytest tests/            # REDSEA compensation, contact matrix, rasterize; RESTORE robust
                         # guard; lineage hierarchy (zero Unassigned); config; parallel
```

The tests need no imaging data (they use synthetic masks/frames) and run in ~2 s.

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
