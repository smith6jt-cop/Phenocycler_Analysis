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

## Outputs (hive-partitioned by section)

`donor_id=<key>` is a **section** key, not a bare donor id: the donor's digits plus the
region token its qptiff carries, lowercased. A donor with two images gets two partitions.

| image | partition |
|---|---|
| `6539_Scan1.er.qptiff` | `donor_id=6539` (pancreas) |
| `6539pLN_Scan1.er.qptiff` | `donor_id=6539pln` (pancreatic lymph node) |

One partition is one image throughout: one GeoJSON, one RESTORE threshold set, one coordinate
frame. `phenocycler/sections.py` defines the key and `cfg.discover_sections()` iterates them;
`cfg.discover_donors()` returns the unique donors behind them.

```
data/
├── cells/donor_id=<key>/*.parquet                # Step 1: raw MFI + morphology + region
├── cells_redsea/donor_id=<id>/data_0.parquet     # Step 2: REDSEA spillover-corrected means
├── restore_redsea/                               # Step 3: threshs.pkl, positive_fractions.csv, qc/
├── restore_thresholds_redsea.csv                 # Step 3: per-image KMeans/GMM/SSC thresholds (chosen flag)
├── restore_gated_redsea/donor_id=<id>/data_0.parquet         # Step 3: {m}_pos/_norm/_log2r (10 markers)
├── restore_gated_redsea_extra/donor_id=<id>/data_0.parquet   # Step 3 (extra): CD99/B3TUBB/MPO _pos/_norm
├── restore_gated_redsea.pre_hormonefloor/                    # Step 4: un-floored backup, auto-created by hormone_floor on first run (rollback point)
├── phenotype/
│   ├── broad/donor_id=<id>/data_0.parquet        # Step 5: broad_lineage, assign_margin, epi_default, score_* (8)
│   ├── broad_lineage_composition.png             # Step 5: composition + disease-trend figure
│   ├── celltype_marker_dotplot.png / _heatmap.png# Step 7: identity QC figures
│   └── qupath_class/pheno_class_<id>.csv         # Step 6: object_id, broad_lineage, image (for QuPath)
├── redsea_reassess/                              # read-only acceptance-yardstick CSVs + figures
├── integration/                                  # PhenoCycler <-> Xenium (optional; see docs/INTEGRATION.md)
│   ├── manifest.csv                              # donor <-> Xenium-run pairing + pair_status   [tracked]
│   ├── vocab_crosswalk.csv                       # harmonised lineages + protein<->gene pairs   [tracked]
│   ├── xenium_paths.csv                          # vendored donor -> Xenium bundle map          [tracked]
│   ├── donor_overrides.csv                       # donor ids not recorded upstream              [tracked]
│   ├── cells_pheno/donor_id=*/roi=*/data_0.parquet   # contract cell table (protein side)
│   ├── cells_xen/donor_id=*/roi=*/data_0.parquet     # contract cell table (RNA side)
│   ├── structures/{phenocycler,xenium}/donor_id=*/roi=*/{islets,ducts,vessels}.parquet
│   ├── registration/donor_id=*/roi=*/            # transform.json, *_displacement.npz, qc.json
│   ├── paired/islets.parquet                     # matched islet pairs      [sequential]
│   ├── paired/cells.parquet                      # paired protein+RNA cells [same_slide]
│   ├── niches/                                   # niche_profiles.csv, grid_bins.parquet
│   ├── crossmodal/pseudocell_links.parquet       # INFERRED links, not measurements
│   ├── qc/                                       # registration_qc.csv, composition.csv,
│   │                                             # disease_trend.csv, integration_report.txt
│   ├── figures/                                  # registration overlays, concordance, trends
│   └── export/                                   # QuPath + Xenium Explorer round-trips
└── redsea_scratch/
    ├── geojson/cells__<image>.geojson            # QuPath boundaries (input to REDSEA)
    ├── masks/<id>.tif                            # int32 instance masks (--keep-mask)
    └── intermediates/<id>.npz                    # data/edge/sizes/contact (--save-intermediates)
```

### The broad-lineage gating markers → 8 lineages

Main pass (`restore_gated_redsea`, 10 markers):
`Pan_Cytokeratin, Vimentin, SMA, CD31, INS, GCG, SST, CD3e, CD20, CD163`. Extra pass
(`restore_gated_redsea_extra`, via `restore --extra`): `CD99, B3TUBB, MPO`. Together they
map to the **eight** lineages Epithelial, Fibroblast, Muscle, Endothelial, Endocrine
(INS/GCG/SST or bright CD99), Immune, **Neural** (B3TUBB), **Neutrophil** (MPO). Before
typing, the hormone floor rewrites `{INS,GCG,SST}_pos = (_norm ≥ 5)`. The full imaging
panel is larger (~59-plex); only these gate the broad lineage.

## Region scheme

`build_cells_parquet` derives `cell_region` from the QuPath `Parent` annotation:
`Islet_N` → `core`, `Islet_N_exp20um` → `peri`, `Annotation (Tissue)` / `Root object`
→ `tissue`, else `other`; and `islet_num` from `Islet_<N>`.

## Integration inputs (optional)

| Input | config.ini key | What it is |
|-------|----------------|------------|
| Xenium run map | `[integration] xenium_paths_csv` | `donor_id, roi, batch, source_path` — vendored from Xenium_Analysis, where nothing reads it. Made load-bearing here. |
| Donor overrides | `[integration] donor_overrides_csv` | `xenium_sample_key, donor_id, roi[, section_gap_um, block_id]`. For Xenium runs whose donor is recorded nowhere upstream (e.g. `0041323`, `0041326`). |
| Panel taxonomy | `[integration] panel_explorer` | XeniumPanelExplorer submodule — `panel_roles.csv` + `identity_core/*.csv`, shipped upstream as a machine-readable contract for Python consumers. |
| Xenium outputs | resolved from the manifest | `{sample}_phenotyped.h5ad`, or the SpatialData zarr, or the raw `output-XETG...` bundle. Morphology images are read from the **bundle** — the Xenium ingest excludes them from the zarr. |

The cross-modality join key is `donor_id`, normalised to a string on both sides (the donor
workbook stores it numerically, PhenoCycler derives it by regex). The Xenium *section* key is
`{serial}__{roi}`, never the bare serial: one slide serial can carry both a pancreas and a
lymph-node region.

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
