# Data contracts

The workflow treats raw QuPath exports, biological registries, and every
derived dataset as explicit contracts. A directory's existence is never proof
that a stage completed.

For the scientific rationale and stage algorithms, see
[the workflow guide](docs/WORKFLOW.md). For installation and commands, see
[the project README](README.md).

## Source files

The current cohort contract contains exactly one image record per donor.
Different image records may point to the same batched measurement CSV.

| Source | Required contract |
|---|---|
| QuPath measurement CSV | Exact `Image` and `Object ID` values; centroids, parent, cell area/solidity, cell marker means, and nucleus means for registered nuclear/state markers |
| qptiff | Original channel stack used for REDSEA; channel order and names must agree with the declared map |
| cell GeoJSON | One feature per detection, with feature `id` equal to the QuPath UUID and a usable cell geometry |
| nucleus geometry | Either `nucleusGeometry` embedded in the cell GeoJSON or a separate UUID-matched nucleus GeoJSON |
| image metadata | Positive X/Y pixel size, panel ID, channel map, and a non-empty segmentation version; the current REDSEA path requires isotropic X/Y pixels |

The canonical key is the QuPath `Object ID`, represented downstream as
`object_id`. Centroid matching and filename parsing are not substitutes for
this UUID.

The measurement CSV is expected to contain these QuPath fields:

```text
Image
Object ID
Centroid X µm
Centroid Y µm
Parent
Cell: Area µm^2
Cell: Solidity
Cell: <marker>: Mean
Nucleus: <marker>: Mean       # where required by the registry
```

The ingest stage derives `cell_region` and `islet_num` from `Parent`, but these
context fields do not define calibration controls or cell identity. Other
QuPath morphology, classification, and distance columns may remain in the raw
CSV but are intentionally not copied into the canonical cell table. Nucleus
and cell raster areas used by geometry QC come from the contracted GeoJSON.

Use
[`scripts/groovy/export_cells_geojson.groovy`](scripts/groovy/export_cells_geojson.groovy)
for UUID-preserving geometry export. Always pass the intended output directory
as its first script argument (`--args "<output-dir>"`); the exporter refuses
to run without it.

## Cohort specification

Create a starter specification:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
```

Paths are resolved relative to the specification file. A minimal specification
is:

```json
{
  "spec_version": 1,
  "cohort_name": "islet-phenocycler",
  "expected_donors": ["D001"],
  "defaults": {
    "measurement_csv": "raw/CellMeasurements.csv",
    "segmentation_version": "qupath-project-commit-or-model-version"
  },
  "images": [
    {
      "donor_id": "D001",
      "image_id": "exact QuPath Image column value",
      "qptiff": "raw/D001.qptiff",
      "cell_geojson": "raw/D001_cells.geojson",
      "nucleus_geojson": null,
      "pixel_size_um_x": 0.325,
      "pixel_size_um_y": 0.325,
      "panel_id": "panel-v1",
      "channel_map": null
    }
  ]
}
```

`nucleus_geojson: null` declares that each cell feature carries
`nucleusGeometry`; a path declares a separate nucleus FeatureCollection.
`channel_map: null` asks the manifest builder to derive canonical marker names
from the qptiff channel metadata. An explicit map uses:

```json
[
  {
    "marker": "DAPI",
    "channel_name": "DAPI",
    "channel_index": 0
  },
  {
    "marker": "E_cadherin",
    "channel_name": "E-cadherin",
    "channel_index": 1
  }
]
```

The array must enumerate every qptiff channel in contiguous zero-based index
order. `channel_name` must exactly match qptiff metadata; `marker` is the
canonical registry name. Every marker, channel name, and channel index must be
unique within an image. Images with the same `panel_id` must have equivalent
mappings. Different panel IDs may intentionally differ.

Build the immutable cohort manifest:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

The builder validates geometry and panel consistency, captures source
fingerprints, and assigns a cohort `content_id`. The default hash modes are:

| Source | Default | Meaning |
|---|---|---|
| measurement CSV | `full` | SHA-256 over every byte |
| qptiff | `sampled` | deterministic windows plus size and sampling layout; explicitly weaker than a full hash |
| geometry | `full` | exact file fingerprint plus UUID/geometry summary |

Select stronger or weaker modes explicitly when creating the manifest:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json \
  --measurement-hash-mode full \
  --qptiff-hash-mode full
```

Do not hand-edit the generated cohort manifest. Edit the small specification
and recreate it so its identity can be verified.

## Panel availability

Ingest forms the union, not the intersection, of marker columns across declared
measurement exports. A marker that was not acquired for a panel is retained as
null for that donor. `panel_availability.json` records each image, panel, and
observed marker set.

Availability has distinct meanings:

| Condition | Representation | Downstream consequence |
|---|---|---|
| Marker not acquired in donor's panel | null authoritative value; `marker_not_acquired_in_panel` audit reason | per-cell marker state `unavailable` |
| Marker acquired, required declared source present, cell value missing | null for that cell | that cell is `unavailable` for the marker |
| Marker acquired, required declared source column absent | fatal contract error | stage stops |
| Marker measurable, donor-marker calibration invalid | measurement retained | state `indeterminate`, not `unavailable` |
| Marker measurable, within bootstrap decision interval | measurement and evidence retained | state `indeterminate` |

Missing markers never become negatives, never make a cell `Other`, and never
transfer their fixed error allocation to another marker.

## Ingested cell schema

Ingest writes one donor-partitioned Parquet dataset and preserves the union of
panel marker columns. Core columns are:

| Column | Meaning |
|---|---|
| `donor_id` | exact donor ID from the cohort manifest |
| `image` | exact QuPath `Image` value |
| `object_id` | canonical QuPath UUID |
| `X_centroid`, `Y_centroid` | centroid in micrometres |
| `parent_raw` | QuPath parent context retained for audit only |
| `cell_region`, `islet_num` | context derived from `Parent` |
| `cell_area`, `cell_solidity` | structural measurements retained from the CSV |
| `<marker>` | sanitized whole-cell mean |
| `Nucleus__<marker>` | retained nucleus mean where available |

The ingest floor removes only degenerate objects below
`[cells] min_cell_area` (5 µm² by default). Biological and geometry decisions
are not made at ingest. CSV parser rejects are written to
`ingest_rejects.parquet`; any rejected row makes ingest fail.

## Geometry QC schema

Geometry QC is an immutable per-cell decision made before REDSEA. Its principal
fields are:

| Column | Meaning |
|---|---|
| `analysis_eligible` | cell may receive an accepted identity; ineligible cells remain present but type as unavailable |
| `estimation_eligible` | cell may help estimate donor-local control models |
| `spillover_context_eligible` | polygon may provide physical REDSEA adjacency context; false labels are removed before contact construction |
| reason fields | explicit failed rules and audit context |
| raster/cell area ratio | agreement between rasterized polygon and QuPath area |
| nucleus/cell area ratio | identifies cells with no meaningful cytoplasm |
| duplicate metadata | UUID/centroid/area duplicate checks |

Analysis eligibility and estimation eligibility are intentionally different. A
real cell with no measurable cytoplasm, for example, may remain in the analysis
universe but must not estimate a cytoplasmic background model. Missing geometry
or prohibited duplicates is fatal.

## Authoritative expression

[`phenocycler/marker_registry.json`](phenocycler/marker_registry.json) declares
the only accepted source for each marker:

- biological kind: identity marker or process-state marker;
- active, deferred, or standalone use;
- whole-cell or nucleus compartment;
- REDSEA-corrected, REDSEA passthrough, or uncompensated policy;
- incompatible reference, marker family, and broad/subtype use.

The `expression` dataset contains exactly one selected column named
`<marker>` for every registered marker, plus:

```text
donor_id
object_id
image
X_centroid
Y_centroid
cell_region
islet_num
qc_analysis_eligible
qc_estimation_eligible
```

Keeping one authoritative value avoids ambiguous mixtures of raw,
whole-cell, nuclear, and corrected columns downstream.
`audit/expression_availability.parquet` explains every donor-marker selection
or absence.

## Reference-control schema

`reference_controls` stores donor-local masks in a wide donor partition. The
canonical long representation records:

```text
mask_schema_version
method_version
registry_version
registry_fingerprint
source
donor_id
reference
object_id
reference_positive
seed_rank
```

Controls are selected from finite, positive reference-marker values among
estimation-eligible cells with finite coordinates. Ranking uses only the
reference intensity, deterministic spatial round-robin selection, and UUID
tie-breaking. Target intensity, target prevalence, QuPath `Classification`,
phenotype labels, and annotation-derived groups are forbidden inputs.

`audit/reference_control_models.parquet` records candidate counts, selected
counts, status, and provenance.

## Marker-evidence schema

`marker_evidence` is wide: one row per UUID and a repeated set of columns for
each active marker `<m>`.

| Column | Meaning |
|---|---|
| `qc_analysis_eligible` | frozen geometry-QC analysis decision propagated to typing |
| `<m>__corrected_intensity` | authoritative selected measurement passed to calibration |
| `<m>__state` | `positive`, `negative`, `indeterminate`, or `unavailable` |
| `<m>__log2_threshold_ratio` | `(log1p(value) - log1p(threshold)) / log(2)`; zero at the donor-local threshold |
| `<m>__empirical_tail_probability` | finite-sample upper-tail p-value under clean controls |
| `<m>__expression_probability` | `1 - empirical_tail_probability` |
| `<m>__empirical_tail_evidence` | `-log10(empirical_tail_probability)` |
| `<m>__model_valid` | valid model and stable positive/negative cell call |
| `<m>__calibration_status` | donor-marker model disposition |
| `<m>__reference` | incompatible marker used to construct the target-independent control mask |
| `<m>__marker_alpha` | fixed registry-derived marker error allocation |
| `<m>__threshold` | empirical threshold on the original intensity scale |
| `<m>__threshold_ci_low`, `<m>__threshold_ci_high` | 95% bootstrap interval |
| `<m>__n_controls` | clean control count |
| `<m>__control_contamination_fraction` | robust high-outlier fraction removed from controls |

Valid calibration statuses are:

```text
valid
insufficient_controls
insufficient_tail_resolution
excess_control_contamination
invalid_background
```

An invalid donor-marker model does not erase its measurements. Measurable cells
become `indeterminate`, preserving the distinction from panel or cell-level
unavailability. Model details are also written to
`audit/marker_calibration_models.parquet`.

## Assignment schema

`assignments` retains both the broad and subtype decisions. Stable handoff
columns are:

| Column | Meaning |
|---|---|
| `broad_type` | accepted broad class, or explicit ambiguity/unavailability label |
| `specific_type` | accepted parent-constrained subtype, terminal broad type, or `<Parent>-unclassified` |
| `assignment_status` | final decision status |
| `confidence` | leading broad probability; for an accepted subtype, the minimum of accepted broad and subtype probabilities |
| `best_broad_type` | leading broad candidate even when it is not accepted |
| `best_specific_type` | leading candidate at the most specific level evaluated |
| `broad_assignment_status` | broad-level decision status |
| `specific_assignment_status` | specific/final decision status |
| `ranked_broad_probabilities` | ordered JSON broad alternatives |
| `ranked_specific_probabilities` | ordered JSON alternatives at the most specific level evaluated |
| `ranked_probabilities` | compatibility alias for the specific-level ranking |

Detailed `broad_*` and `subtype_*` columns preserve model probability, margin,
stability, reason, authoritative markers, conflicts, and rules provenance.
Statuses include:

```text
anchor
inferred
ambiguous
other
unavailable
```

`Other` means every required broad model was valid and negative. It does not
mean “could not classify.”

## State schema

`cell_states` is keyed by donor and UUID but remains separate from
`assignments`. The current state stage reports Ki67 and PCNA:

| Column | Meaning |
|---|---|
| `qc_analysis_eligible` | frozen geometry-QC decision controlling whether a state may be called |
| `<m>__value` | authoritative state-marker value |
| `<m>__state` | `positive`, `negative`, or `unavailable` |
| `<m>__threshold` | donor-local state threshold |
| `<m>__model_status` | `valid`, `marker_unavailable`, or `insufficient_measurements` |
| `<m>__method` | versioned state method |

These calls describe proliferation state and never alter broad or specific
identity. The threshold is fitted only on estimation-eligible cells; a state is
called only for analysis-eligible cells, and geometry-ineligible cells remain
`unavailable`. State thresholds do not claim the same frequency-neutral
mutually-exclusive-control interpretation as identity-marker calibration.

## QuPath export schema

The export stage writes one immutable CSV per donor/image:

```text
qupath_export/cell_assignments_<donor>.csv
```

Columns are:

```text
object_id
image
broad_type
specific_type
assignment_status
confidence
best_broad_type
best_specific_type
broad_assignment_status
specific_assignment_status
ranked_broad_probabilities
ranked_specific_probabilities
ranked_probabilities
compartment
cell_type
```

`compartment` and `cell_type` are compatibility aliases for `broad_type` and
`specific_type`; they do not replace the uncertainty fields. Ambiguous and
unavailable cells must not be imported as confident named classes. The current
QuPath handoff exports identity; process-state annotations remain in the
separate `cell_states` artifact.

Use
[`scripts/groovy/import_cell_assignments.groovy`](scripts/groovy/import_cell_assignments.groovy).
First require a successful pipeline `status`, which verifies the fingerprinted
export. The script then validates the exact image and every row in the
canonical assignment CSV before changing a QuPath class. Detections outside
the canonical ingest universe, including fragments below the ingest area
floor, are intentionally left unchanged. At the selected broad or specific
level, accepted statuses
(`anchor`, `inferred`, and `other`) receive the corresponding label;
`ambiguous` and `unavailable` receive explicit uncertainty classes rather than
the best candidate.

## Run layout

Given a base `data_dir = data`, the content-addressed layout is:

```text
data/
├── qupath_manifest.json
└── runs/
    ├── LATEST
    └── <run_id>/
        ├── run.json
        ├── manifests/
        │   ├── ingest.json
        │   ├── geometry.json
        │   ├── redsea.json
        │   ├── expression.json
        │   ├── controls.json
        │   ├── calibrate.json
        │   ├── type.json
        │   ├── states.json
        │   └── qupath_export.json                  # when exported
        ├── cells/donor_id=<id>/*.parquet
        ├── ingest_rejects.parquet
        ├── panel_availability.json
        ├── geometry_qc/donor_id=<id>/*.parquet
        ├── redsea/donor_id=<id>/*.parquet
        ├── expression/donor_id=<id>/*.parquet
        ├── reference_controls/donor_id=<id>/*.parquet
        ├── marker_evidence/donor_id=<id>/*.parquet
        ├── assignments/donor_id=<id>/*.parquet
        ├── cell_states/donor_id=<id>/*.parquet
        ├── qupath_export/cell_assignments_<id>.csv # when exported
        ├── audit/
        │   ├── expression_availability.parquet
        │   ├── reference_control_models.parquet
        │   ├── marker_calibration_models.parquet
        │   └── state_models.parquet
        └── scratch/redsea/
```

`run.json` exists only after all analytical stage manifests are current.
`runs/LATEST` contains the completed run ID. The QuPath export manifest is
separate and may be created after the analytical run.

## Manifest contents and immutability

Each stage manifest records:

- stage and method version;
- cohort, registry, rules, and configuration identities;
- relevant source-code hash and repository metadata;
- exact upstream artifact identifiers;
- exact expected donor set;
- output Parquet schema, row counts, file fingerprints, and UUID-universe
  hashes;
- exact fingerprints for stage-owned availability and audit sidecars.

A stage directory is immutable within a run. Non-empty output without a valid
manifest is a partial artifact and fails. Changed inputs, configuration, code,
schemas, donors, output bytes, or sidecar bytes make a manifest stale and also
fail. Correct the source contract and run again; the changed content resolves
to a new `run_id`.

## Storage, privacy, and retention

Multiplexed images and per-cell artifacts can be large. Keep `data/` outside
version control, retain the raw source files and generated cohort manifest
together, and archive the registries/rules used by published analyses. For
human tissue, follow the governing IRB, de-identification, access-control, and
retention requirements.
