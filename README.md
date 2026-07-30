# Phenocycler analysis

This package turns an explicitly declared QuPath cohort into donor-local marker
evidence, hierarchical cell identities, and UUID-keyed tables for review in
QuPath.

The production workflow is:

```text
QuPath cohort manifest
  -> measurement ingest
  -> geometry QC
  -> REDSEA boundary-spillover correction
  -> authoritative per-marker expression
  -> deterministic donor-local reference controls
  -> marker calibration with uncertainty
  -> broad and parent-constrained specific typing
  -> separate process-state annotation
  -> uncertainty-preserving QuPath export
```

The command-line pipeline, its production APIs, and its versioned registries
are authoritative. Notebooks are review-oriented views that may invoke those
same immutable stage APIs; they are not an alternative implementation and must
not introduce thresholds, donor exclusions, or typing rules of their own.

## What the workflow guarantees

- Every source image is declared in an immutable cohort manifest. Donors are
  never inferred from filenames.
- The QuPath `object_id` UUID is the cell key from ingest through export.
- Geometry QC precedes REDSEA and uses geometry only. Analysis eligibility,
  model-estimation eligibility, and physical spillover context are separate
  decisions. Spillover-ineligible polygons are removed before REDSEA contact
  construction; analysis-ineligible cells remain in the UUID universe but
  export as unavailable.
- Each registered marker has one authoritative intensity source and
  compartment. An acquired marker with a missing required source is an error;
  a marker absent from a panel remains explicitly unavailable.
- Reference controls are donor-local, deterministic, spatially balanced, and
  selected without looking at the target marker, phenotype annotations, or the
  frequency of target-positive cells.
- Calibration preserves the corrected intensity, a zero-centered log2
  threshold ratio, an empirical upper-tail probability, the threshold
  uncertainty interval, and a four-state call:
  `positive`, `negative`, `indeterminate`, or `unavailable`.
- Broad types are evaluated simultaneously. There are no ordered residual
  gates, prevalence-derived cutoffs, or forced assignments.
- Specific types are constrained by an accepted broad parent. State markers
  such as Ki67 and PCNA are reported separately and cannot change cell
  identity.
- Every run and stage is content addressed. Partial, stale, or donor-incomplete
  output fails closed instead of being silently reused.

The detailed scientific and operational contract is in
[the workflow guide](docs/WORKFLOW.md). Input, output, and column contracts are
in [the data guide](DATA_README.md).

## Install

Python 3.10 or newer is required. In an activated analysis environment:

```bash
python -m pip install -e .
```

For development:

```bash
python -m pip install -e '.[test]'
pytest -q
```

## 1. Export the QuPath inputs

For every image, retain:

- the QuPath cell-measurement CSV;
- the matching `.qptiff`;
- full-resolution cell and nucleus geometry;
- the exact value in the CSV `Image` column;
- pixel calibration, panel identity, channel mapping, and segmentation version.

Use [`scripts/groovy/export_cells_geojson.groovy`](scripts/groovy/export_cells_geojson.groovy)
to export geometry with the QuPath detection UUID intact. Nucleus geometry may
be embedded as `nucleusGeometry` in the cell GeoJSON or supplied as a separate
UUID-matched GeoJSON. Always pass the intended output directory as the first
script argument (`--args "<output-dir>"`); the exporter refuses to run without
it.

## 2. Create the cohort manifest

Generate a versioned specification template:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
```

Edit the specification, then fingerprint and validate its sources:

```bash
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

Measurement CSVs receive a full SHA-256 fingerprint by default. Large qptiffs
use an explicitly recorded deterministic sampled fingerprint by default; use
`--qptiff-hash-mode full` when full-file hashing is practical.

The specification schema and a complete example are in
[DATA_README.md](DATA_README.md#cohort-specification).

## 3. Configure and run

At minimum, the configuration identifies the data root and cohort manifest:

```ini
[paths]
data_dir = data
qupath_manifest = data/qupath_manifest.json
```

The package defaults supply the marker registry and typing rules. A project can
pin explicit paths to reviewed copies:

```ini
marker_registry = phenocycler/marker_registry.json
typing_rules = phenocycler/typing_rules.json
```

Run and inspect the complete workflow:

```bash
python -m phenocycler.pipeline run --config config.ini --jobs 4
python -m phenocycler.pipeline status --config config.ini
```

Useful bounded runs are:

```bash
# Run the dependency prefix through expression selection.
python -m phenocycler.pipeline run --config config.ini --through expression

# Run named stages only; their prerequisite manifests must already be current.
python -m phenocycler.pipeline run --config config.ini \
  --only controls calibrate type states

# Build all analytical stages without the QuPath handoff.
python -m phenocycler.pipeline run --config config.ini --no-export

# Export a current typing result later.
python -m phenocycler.pipeline export --config config.ini
```

Run `status` immediately before a standalone export and proceed only when all
analytical stages are `CURRENT` (the QuPath line may still be `MISSING`). Run
it again after export and require success before importing the CSVs into
QuPath. This validates the typing artifact and fingerprinted handoff together.

`status` checks the QuPath handoff as well as analytical stages. After an
intentional `--no-export` run it reports that export as missing and returns
nonzero until `export` is run.

The stage order is fixed:

```text
ingest -> geometry -> redsea -> expression -> controls
       -> calibrate -> type -> states
```

`--only` does not infer or rebuild dependencies. This makes every requested
transition explicit and prevents an interactive notebook from quietly
substituting a different upstream result.

## Review notebooks

The notebooks follow the same four conceptual checkpoints:

1. [`01_qupath_ingest_and_geometry_qc.ipynb`](notebooks/01_qupath_ingest_and_geometry_qc.ipynb)
   reviews the cohort contract, ingest, and geometry decisions.
2. [`02_redsea_spillover.ipynb`](notebooks/02_redsea_spillover.ipynb)
   reviews REDSEA alignment, compartments, and authoritative expression.
3. [`03_marker_calibration.ipynb`](notebooks/03_marker_calibration.ipynb)
   reviews reference-control and donor-marker calibration audits.
4. [`04_hierarchical_cell_typing.ipynb`](notebooks/04_hierarchical_cell_typing.ipynb)
   reviews hierarchy, state, and the optional QuPath handoff.

Their optional execution cells call `run_stage` or `export_qupath`, which apply
the same manifest validation and immutability checks as the CLI. Notebook cells
must not write stage files directly.

## Marker evidence, not frequency gates

Calibration is performed independently for each donor and active marker. The
marker registry declares a biologically incompatible counterpart used only to
form a target-independent reference mask. Within that mask, a deterministic
spatial rank selects controls. The target background and tail are then
estimated robustly and a fixed empirical upper-tail decision is bootstrapped.

The key outputs answer different questions:

| Value | Interpretation |
|---|---|
| authoritative intensity | selected corrected signal in the marker's declared compartment |
| `log2_threshold_ratio` | signed distance from the donor-local threshold; zero is the threshold |
| `empirical_tail_probability` | finite-sample upper-tail probability under the donor-local controls |
| threshold interval | bootstrap uncertainty in the decision boundary |
| expression state | stable positive/negative, indeterminate, or unavailable |

This representation addresses donor-to-donor scale and background without
using positive-cell frequency as evidence. A rare population does not weaken
its own cutoff, and a common population does not move its cutoff upward.

For mutually exclusive marker families, the registry preallocates a 1%
family-wise error rate using a fixed Bonferroni allocation. Family membership
is declared before data are inspected; absent panel markers do not donate their
error budget to the markers that remain. This is a 1% family bound, not a 1%
per-marker threshold.

Annotations and existing QuPath classifications are audit context only. They
are forbidden as reference-control inputs and cannot define calibration or
identity.

## Identity and state

Broad identity evaluates evidence for Immune, Endothelial, Epithelial, Neural,
and Mesenchymal classes simultaneously. It retains the best candidate, ranked
class probabilities, assignment confidence, and one of:
`anchor`, `inferred`, `ambiguous`, `other`, or `unavailable`.

`Other` is permitted only when every required broad model is valid and
negative. Missing or invalid evidence cannot be converted into `Other`.
Specific typing is parent constrained, so unresolved children retain a
`<Parent>-unclassified` label rather than crossing broad compartments.

Two panel-specific rules are deliberate:

- SST is subtype-only evidence for Delta after an Epithelial parent is
  accepted. SST never creates broad Epithelial evidence.
- CD11b contributes to Immune evidence and is calibrated against
  E-cadherin, its declared incompatible reference.

Process-state calls are stored in a separate `cell_states` artifact. They do not
rewrite broad or specific identity.

## Runs and failure behavior

The run identifier is derived from the cohort content ID, marker-registry
fingerprint, typing-rules fingerprint, identity-affecting configuration, and
relevant source-code hash. Outputs live below:

```text
data/runs/<run_id>/
```

Each stage writes a manifest containing its exact inputs, configuration, source
hash, expected donors, output schema, row counts, object-ID universe, and
fingerprints for its audit/availability sidecars. A complete analytical run
writes `run.json`; `data/runs/LATEST` records its run ID. QuPath export has its
own manifest because it is a separate handoff.

Existing valid content is reused. The pipeline stops when it finds any of the
following:

- source data, registry, rules, configuration, or code changed;
- a donor or expected partition is missing;
- an output schema or UUID universe changed;
- an audit or availability sidecar changed or disappeared;
- an acquired marker's required source column is absent;
- geometry is missing or object IDs are duplicated;
- a non-empty output exists without a valid manifest.

Do not repair these failures by deleting manifests or mixing files between run
directories. Correct the declared input or configuration and create the new
content-addressed run.

## Outputs

The principal artifacts under a run are:

- `expression/`: one authoritative intensity per registered marker;
- `reference_controls/`: donor/reference membership and deterministic ranks;
- `marker_evidence/`: calibrated per-cell expression evidence and uncertainty;
- `assignments/`: broad and specific identity with ranked alternatives;
- `cell_states/`: process-state annotations kept separate from identity;
- `qupath_export/`: UUID-keyed identity CSVs that retain uncertainty;
- `audit/`: expression availability, control models, calibration models, and
  state models;
- `manifests/`: immutable stage and export contracts.

See [DATA_README.md](DATA_README.md#run-layout) for the full tree and schemas.
Use
[`scripts/groovy/import_cell_assignments.groovy`](scripts/groovy/import_cell_assignments.groovy)
for a validation-first UUID import that keeps broad, specific, and uncertainty
decisions separate. Run `status` first: the importer validates every CSV row
but intentionally leaves QuPath detections outside the canonical ingest
universe (for example, fragments below the structural area floor) unchanged.

## Migration note

The former RESTORE and ordered-residual lineage path is retired. Its thresholds,
positive fractions, normalized columns, and cell labels are not valid inputs
to this workflow and must not be mixed with a content-addressed run.

## Scope

This repository begins with QuPath measurement, geometry, and qptiff outputs.
Image registration, illumination correction, segmentation-model development,
and downstream spatial-neighborhood or trajectory analyses are separate
workflows. The goal here is a reproducible per-cell evidence layer that supports
eventual assignment of every cell while preserving honest ambiguity whenever
the panel or calibration cannot support a specific call.
