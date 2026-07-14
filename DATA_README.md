# Data layout & inputs

All pipeline data lives under `data/` (git-ignored). Paths are configured in
`config.ini` `[paths]`; nothing is hardcoded in the code.

## Inputs you provide

| Input | config.ini key | What it is |
|-------|----------------|------------|
| Raw images | `images_dir` | Directory of multichannel `.qptiff` files (one per donor image). REDSEA reads these one channel at a time. |
| Per-cell measurements | `cells_csv` | The full QuPath detection export (`Cellmeasurements.csv`): one row per cell, columns include `Image`, `Object ID`, `Centroid X/Y µm`, `Parent`, `Cell: Area µm^2`, and `Cell: <marker>: Mean` for every marker. |
| Per-cell boundaries | (derived, in `data/redsea_scratch/geojson/`) | Full-resolution GeoJSON exported from QuPath, **feature `id` == detection UUID**. Produced by `scripts/groovy/export_cells_geojson.groovy`. |
| Donor metadata | `donor_metadata` | Workbook with a donor-id column and a disease-status column (default `donor_id`, `disease.status`; e.g. ND / AAB / T1D). Configurable via `[metadata]`. |

The `object_id` (QuPath `PathObject.getID()`, a UUID) is the join key that threads
through every stage, so REDSEA, RESTORE, lineage and the QuPath re-import all align
exactly with no centroid rounding.

### Producing the inputs (upstream)

1. Process raw images with [KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI)
   (illumination correction → stitching → deconvolution → EDF → registration →
   autofluorescence removal).
2. In QuPath: run cell detection/segmentation, then
   - `Measure → Export measurements` → `Cellmeasurements.csv`, and
   - run `scripts/groovy/export_cells_geojson.groovy` per image (GUI closed) to write
     `data/redsea_scratch/geojson/cells__<image>.geojson`.

The `<image>` tag in the GeoJSON filename is the QuPath image name with
`[^A-Za-z0-9._-]` replaced by `_` (e.g. `6539_Scan1.er.qptiff - resolution #1` →
`cells__6539_Scan1.er.qptiff_-_resolution__1.geojson`). `phenocycler.redsea`
reconstructs exactly this tag from each donor's `image` column, so do not rename.

## Outputs (hive-partitioned by donor)

```
data/
├── cells/donor_id=<id>/*.parquet                 # Step 1: raw MFI + morphology + region
├── cells_redsea/donor_id=<id>/data_0.parquet     # Step 2: REDSEA spillover-corrected means
├── restore_redsea/                               # Step 3: threshs.pkl, positive_fractions.csv, qc/
├── restore_thresholds_redsea.csv                 # Step 3: per-image KMeans/GMM/SSC thresholds (chosen flag)
├── restore_gated_redsea/donor_id=<id>/data_0.parquet         # Step 3: {m}_pos/_norm/_log2r (all markers, one pass)
├── phenotype/
│   ├── broad/donor_id=<id>/data_0.parquet        # Step 4: compartment (5 + Other), cell_type, object_id, donor_id
│   ├── broad_lineage_composition.png             # Step 4: composition + disease-trend figure
│   ├── celltype_marker_dotplot.png / _heatmap.png# Step 6: identity QC figures
│   └── qupath_class/pheno_class_<id>.csv         # Step 5: object_id, compartment, cell_type, image (for QuPath)
├── redsea_reassess/                              # read-only acceptance-yardstick CSVs + figures
└── redsea_scratch/
    ├── geojson/cells__<image>.geojson            # QuPath boundaries (input to REDSEA)
    ├── masks/<id>.tif                            # int32 instance masks (--keep-mask)
    └── intermediates/<id>.npz                    # data/edge/sizes/contact (--save-intermediates)
```

### The broad-lineage gating markers → 7 lineages

Main pass (`restore_gated_redsea`, 10 markers):
`Pan_Cytokeratin, Vimentin, SMA, CD31, INS, GCG, SST, CD3e, CD20, CD163`. Extra pass
(`restore_gated_redsea_extra`, via `restore --extra`): `CD99, B3TUBB, MPO`. Together they
map to the **seven** lineages Epithelial, Fibroblast, Muscle, Endothelial, Endocrine
(INS/GCG/SST or bright CD99), Immune (CD3e/CD20/CD163/**MPO** — neutrophils are immune),
**Neural** (B3TUBB). Before typing, the norm floor rewrites `{INS,GCG,SST}_pos = (_norm ≥ 5)`
and `{CD3e,CD20,CD163}_pos = (_norm ≥ 2)` (MPO gated at 2 in lineage). The full imaging panel
is larger (~59-plex); only these gate the broad lineage.

## Region scheme

`build_cells_parquet` derives `cell_region` from the QuPath `Parent` annotation:
`Islet_N` → `core`, `Islet_N_exp20um` → `peri`, `Annotation (Tissue)` / `Root object`
→ `tissue`, else `other`; and `islet_num` from `Islet_<N>`.

## Storage & reproducibility

- The full cohort (~15–22 donors, ~23M cells) is tens of GB; `data/` is git-ignored.
- Keep raw `.qptiff` and `Cellmeasurements.csv` as the source of truth — every
  `data/` artifact is regenerable from them via `python -m phenocycler.pipeline`.
- For human tissue: ensure IRB approval, de-identification, and secure storage.

## Smoke-testing without real data

`build_cells_parquet(cfg, limit=200_000)` builds a small subset; a single-donor REDSEA
run (`python -m phenocycler.redsea --donor <id>`) is the cheapest end-to-end check. The
unit tests (`pytest tests/`) exercise the REDSEA math, RESTORE guard and lineage
hierarchy with synthetic inputs and need no imaging data.
