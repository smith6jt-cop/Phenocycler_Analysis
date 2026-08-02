# Current QuPath-to-cell-type workflow

This document is the operational and scientific contract for the production
workflow. It explains how a raw QuPath measurement becomes donor-local marker
evidence and then a hierarchical identity without using observed positive-cell
frequency to define the answer.

Start with [the project README](../README.md) for commands and
[the data guide](../DATA_README.md) for schemas.

## Intended outcome

The immediate product is not a single binary gate per marker. It is a
comparable evidence record for every donor, marker, and cell:

```text
authoritative corrected intensity
  + donor-local background-relative scale
  + empirical probability/evidence
  + threshold uncertainty
  + explicit positive/negative/indeterminate/unavailable state
```

Those records support the larger objective of assigning all cells to broad and
specific cell types. The workflow preserves ambiguity when the measurement,
panel, control population, or classifier cannot support an assignment; it does
not manufacture completeness by converting missing evidence to negative or
`Other`.

## Authorities and boundaries

Four versioned inputs always govern a run:

1. The QuPath cohort manifest owns raw-data identity, donor membership, image
   mapping, geometry, panel, channel, pixel calibration, and segmentation
   version.
2. [`marker_registry.json`](../phenocycler/marker_registry.json) owns marker
   meaning, compartment selection, spillover policy, calibration status,
   incompatible references, uses, and fixed testing families.
3. [`typing_rules.json`](../phenocycler/typing_rules.json) owns broad and
   subtype anchors, support markers, exclusions, and parent relationships.
4. The configuration owns fixed geometry, REDSEA, calibration, and compute
   parameters.

Probabilistic production typing adds one inseparable authorization pair: an
immutable threshold-frozen model bundle and its separately signed promotion
manifest. Both paths must be configured together; the manifest must authorize
that exact bundle, donor split, registries, calibration source, validation and
locked-test reports, replay code, and numerical runtime. Leaving both paths
blank is the rules-only default. See
[`PROBABILISTIC_RESCUE.md`](PROBABILISTIC_RESCUE.md) for the fitting and release
contract.

The command-line pipeline and its stage APIs are the production writers.
Notebooks may read and visualize immutable artifacts, help review audits, and
optionally invoke `run_stage` or `export_qupath`. Those calls follow the same
manifest contracts as the CLI. Notebooks must not:

- infer donors or panels from filenames;
- write files directly into `data/runs/<run_id>/`;
- create their own marker pairs, thresholds, or lineage ordering;
- replace unavailable values with zero;
- select controls from phenotype annotations or target-positive frequency;
- overwrite a stage manifest or reuse an output from another run.

## Stage graph

```text
                 ┌──────────── immutable biological registries ────────────┐
                 │                                                         │
QuPath manifest ─┴─> ingest -> geometry -> REDSEA -> expression            │
                                                  │                        │
                                                  ├-> controls -> calibrate├-> type
                                                  │                        │
                                                  └────────────────────────┴-> states
                                                                            │
type ----------------------------------------------------------------> QuPath export
```

The fixed analytical order is:

```text
ingest -> geometry -> redsea -> expression -> controls
       -> calibrate -> type -> states
```

State annotation depends on authoritative expression, not on identity.
QuPath export depends on typing and has a separate manifest.

## Stage 0: lock the QuPath cohort

The small cohort specification is the only place raw source paths are entered.
The manifest builder resolves those paths and captures:

- one donor and exact QuPath `Image` value per image record;
- full measurement-CSV fingerprint;
- qptiff fingerprint and exact marker/channel mapping;
- cell and nucleus GeoJSON fingerprints, feature counts, and UUID universes;
- X/Y pixel calibration;
- panel ID and segmentation version;
- exact expected donor set and cohort content ID.

One image record per donor is the current contract. Measurement CSVs may be
shared across image records when QuPath exported a cohort batch.

The resolver also enforces the central donor exclusions in
[`cohort.py`](../phenocycler/cohort.py). The current exclusions are donor 6457
(failed registration associated with poor DAPI staining) and donor 6579
(maintainer-designated pancreatitis outlier). A manifest containing either
donor fails; exclusions are not silently applied after ingest.

The current REDSEA raster path requires isotropic pixels. A declared X/Y
calibration mismatch, or disagreement between the manifest and qptiff
calibration, fails before correction.

Create it with:

```bash
python -m phenocycler.pipeline manifest template --out cohort-spec.json
python -m phenocycler.pipeline manifest create \
  --spec cohort-spec.json \
  --out data/qupath_manifest.json
```

The manifest performs mapping; the pipeline never extracts a donor ID from an
image or file name. A CSV row whose `Image` value is not declared is fatal.

## Stage 1: ingest the measurement union

Ingest reads every unique measurement CSV named by the cohort manifest using
DuckDB. It:

- unions columns by name across exports;
- retains batch- or panel-specific markers as null where unacquired;
- maps the exact QuPath `Image` value to the manifest donor;
- preserves the UUID as `object_id`;
- retains whole-cell means and required nucleus means;
- retains context and morphology;
- applies only a small structural area floor for degenerate fragments;
- partitions the result by donor.

Panel union matters. Taking the intersection would silently discard useful
markers from donors that acquired them; filling an absent marker with zero
would falsely turn “not measured” into “measured negative.”

The CSV parser stores rejects in `ingest_rejects.parquet`. Any rejected row
stops ingest after the audit is written. `panel_availability.json` records the
observed union, intersection, and per-panel acquisitions.

## Stage 2: geometry-only QC

Geometry QC is deliberately before REDSEA and does not read marker intensity.
For each UUID it compares QuPath measurements with the declared cell/nucleus
polygons and records separate flags:

- `analysis_eligible`: may receive accepted identity and state calls;
- `estimation_eligible`: may help estimate donor-marker background;
- `spillover_context_eligible`: may participate as physical adjacency context.

The checks include minimum and donor-relative maximum area, solidity,
raster-to-QuPath area agreement, nucleus/cell area ratio, geometry presence,
and spatial duplicate detection. Current defaults include a 20 µm² cell area
minimum, four times the donor median maximum, 0.70 solidity minimum, 0.5–2.0
raster-area ratio, and duplicate checks within 1 µm with 15% area tolerance.
The cohort configuration explicitly uses `duplicate_policy =
deterministic_keep` for the known duplicated islet detections: a QC-viable
nucleus-bearing copy is preferred, then the copy retaining the most raster
pixels, then the lexicographically first UUID. Every non-survivor is excluded
from analysis, estimation, and REDSEA physical context. All identity-affecting
values are captured in the run ID.

The three flags must not be collapsed. A real cell with essentially no
cytoplasm may remain analysis eligible for nuclear measurements while being
unfit to estimate a cytoplasmic background. Conversely, a polygon excluded
from cell-level analysis may still be required to describe a real shared
boundary during spillover correction.

A missing nucleus fails closed at the cell level: the object receives
`missing_nucleus_geometry` and is ineligible for analysis and model
estimation. Its cell polygon is also removed from physical spillover context
before REDSEA. Missing canonical cell geometry and prohibited duplicates
remain stage-fatal. Geometry QC does not delete the canonical cell universe.

## Stage 3: REDSEA

REDSEA corrects lateral signal transferred across adjacent cell boundaries in
the qptiff. It consumes the declared image/channel mapping, UUID-matched
geometry, and pixel calibration. The geometry stage is a required upstream
contract. Before boundary and contact construction, canonical polygons whose
`spillover_context_eligible` flag is false are removed from the cell and
nucleus masks. GeoJSON-only fragments outside the canonical ingest universe
may remain as physical boundary context. The current production method emits
compartment-resolved measurements and preserves exact object alignment.

REDSEA rechecks qptiff channel names and contiguous indices against the
manifest and uses the manifest's canonical marker names. A changed or
misdeclared channel map fails rather than shifting marker identity by position.
Before output, an alignment contract also requires complete canonical UUID
coverage and agreement between rerasterized and QuPath means on available
probe channels. At least two of DAPI, INS, Pan_Cytokeratin, CD3e, Vimentin, and
CD31 must be acquired; each available probe requires Pearson correlation at
least 0.80 and a median raster/QuPath ratio between 0.25 and 4.

Nucleus geometry may be embedded as `nucleusGeometry` in each cell feature or
supplied as a separate UUID-matched GeoJSON. Separate nucleus labels are
relabelled into the cell-mask label space and intersected with their own cell
before correction. An embedded feature without a nucleus is retained only for
audit; its polygon is removed before correction and it cannot receive an
accepted cell identity or state call.

The marker registry decides whether an authoritative marker later uses:

- a REDSEA-corrected compartment;
- a REDSEA passthrough for a declared no-spillover marker; or
- an uncompensated raw QuPath compartment.

REDSEA therefore corrects a physical imaging artifact; it does not decide
whether a marker is positive and does not discover cell types. Nuclear and
other no-spillover policies are explicit rather than being silently treated as
whole-cell corrected signal.

Current method and parameter identities, including compartment mode, alpha,
gap bridging, normalization form, and any no-spillover exclusions, are stored
in the stage and run manifests.

## Stage 4: choose one authoritative intensity

The expression stage joins ingest, geometry QC, and REDSEA by exact UUID and
requires the same object universe. For every registered marker it follows the
registry's compartment and spillover policy and emits one column named for the
marker.

This is the boundary that removes competing meanings of “expression.” A
downstream analysis does not choose opportunistically among raw whole-cell,
raw nucleus, corrected cell, or other compartment columns.

For an acquired marker, a missing required source column is a contract failure.
For a marker absent from the declared panel, the selected column is null and
the availability audit records `marker_not_acquired_in_panel`.

The result also carries UUID, donor, image, coordinates, region context, and
the separate analysis/estimation eligibility flags. Existing QuPath
`Classification` or annotation context is not part of authoritative expression
and cannot be used to estimate calibration.

## Stage 5: deterministic donor-local reference controls

An active marker has a predeclared biologically incompatible reference. The
reference is used to identify cells expected to be negative for the target; it
is not used to algebraically normalize the target and it does not set a target
positive fraction.

For each donor and distinct reference marker, control selection:

1. starts with `estimation_eligible` cells;
2. requires finite positive reference intensity and finite X/Y coordinates;
3. divides the donor image extent into an 8 by 8 equal-width grid;
4. ranks cells within each occupied tile by descending reference intensity,
   with a stable UUID-derived tie break;
5. takes the brightest cell from every occupied tile before the
   second-brightest from any tile, continuing in deterministic rounds;
6. stops at 20,000 controls;
7. fails to an all-false mask when fewer than 10,000 are available.

The selection reads no target intensity. It therefore cannot change when the
target-positive population is rare, common, absent, or biologically shifted.
The mask records method, registry fingerprint, source, UUID, membership, and
rank.

The reference mask is built directly from its authoritative intensity, not
from a previously calibrated reference call. Reciprocal declarations such as
E-cadherin/CD31 therefore do not create an estimation cycle.

External control masks must use the canonical schema and exact UUID universe.
They are validated for registry provenance, ranks, counts, and methods.
Columns containing annotation, classification, phenotype, or cell-type labels
are rejected. Annotations remain audit-only context.

## Stage 6: donor-marker calibration

Calibration operates independently for every donor and active marker. It never
pools donors to create a threshold and never estimates a target-positive
frequency.

### Fit

For target marker \(m\), the estimator:

1. intersects the frozen reference-positive mask with estimation eligibility,
   a finite nonnegative target value, and finite coordinates;
2. uses at most 20,000 spatially balanced controls and requires at least
   10,000;
3. transforms controls with `log1p`;
4. estimates the background median and Qn robust scale;
5. removes controls above robust z = 6 as likely target contamination;
6. invalidates the model if contamination exceeds 5% or fewer than 10,000
   clean controls remain;
7. chooses the fixed empirical upper-tail order statistic at the marker's
   preallocated alpha;
8. obtains a deterministic 95% threshold interval from 200 bootstrap
   replicates.

The possible model dispositions are:

```text
valid
insufficient_controls
insufficient_tail_resolution
excess_control_contamination
invalid_background
```

A failed calibration is data, not an invitation to substitute another
threshold. Its measurements remain available but its measurable cell calls are
`indeterminate`.

### Per-cell evidence

The frozen donor-marker model is applied to every measurable cell for audit,
including cells excluded from estimation. The frozen
`qc_analysis_eligible` flag is carried beside that evidence so typing can
enforce the analysis universe. For intensity \(x\) and threshold \(t\):

```text
log2_threshold_ratio = [log1p(x) - log1p(t)] / log(2)
```

This value is zero at the threshold, positive above it, and negative below it.
It is a signed background-relative coordinate, not a raw-intensity fold
change.

The clean control distribution also gives the finite-sample upper-tail
probability:

```text
p = [1 + number of controls with log1p(control) >= log1p(x)]
    / [number of controls + 1]
```

The artifact stores `p`, `1 - p` as expression probability, and `-log10(p)` as
tail evidence. These are ranking/evidence summaries; the categorical decision
also requires threshold-uncertainty stability:

- `positive`: value is above the upper bootstrap bound and `p <= marker_alpha`;
- `negative`: value is at or below the lower bootstrap bound and
  `p > marker_alpha`;
- `indeterminate`: measurable but between the bounds, or the model is invalid;
- `unavailable`: no finite nonnegative measurement for this cell.

This distinction prevents a panel absence from masquerading as biological
negative expression.

### Fixed 1% family-wise error

Each declared broad or subtype marker family has
`familywise_alpha = 0.01`. For a family with \(K\) predeclared members, each
member receives:

```text
marker_alpha = 0.01 / K
```

For a cell under the family null, the union bound makes the probability of at
least one false-positive member call at most 1% when the member p-values are
valid. This Bonferroni bound does not require independence.

This is a 1% family-wise bound, not 1% for every marker. The family membership
is versioned before any donor is processed. A panel-absent member remains in
the denominator; its unused allocation is not redistributed. INS and GCG
participate in both broad epithelial and endocrine-subtype families and use
the stricter of their two fixed allocations.

## Why this avoids frequency, scale, and background failure

| Failure mode | Preventive design |
|---|---|
| Rare target population disappears under a prevalence rule | target prevalence is never an estimator input |
| Common target population raises its own threshold | the control mask is selected without target values |
| Donor intensity scales differ | each donor receives a separate background and threshold |
| Background varies across the image | deterministic controls are spatially balanced |
| A few incompatible-reference cells contain target signal | robust contamination screen and fail-closed limit |
| A mutually exclusive counterpart is absent or inadequate | calibration becomes indeterminate; it is not replaced by a frequency gate |
| Threshold is unstable | bootstrap band produces indeterminate cells near the boundary |
| Marker was not acquired | explicit unavailable state rather than zero/negative |
| Many markers inflate false positives | fixed 1% family-wise allocation |

No GMM, K-means, Otsu, NNMF, positive-fraction target, quantile normalized to
the observed positive population, or ordered residual gate participates in
production marker calibration.

## Stage 7: simultaneous broad and hierarchical specific typing

Typing consumes calibrated marker evidence, not raw intensity columns. The
versioned rules evaluate five broad classes simultaneously:

```text
Immune
Endothelial
Epithelial
Neural
Mesenchymal
```

Each class has curated anchors, support markers, and exclusions. The additive
multinomial model has no intercept, so absence of positive evidence does not
create a class through a donor-specific prevalence prior. Coefficient signs
are constrained by the biological roles.

The production default is `typing_mode=rules_only`: the type stage constructs
a transparent, unfitted classifier from fixed registry weights and never fits
the donor being classified. This lane accepts authoritative anchors, resolved
`Other`, and explicit uncertainty states, but it cannot issue a non-anchor
`inferred` call because it has no promoted donor-bootstrap stability model.

The optional `typing_mode=probabilistic_bundle` lane is enabled only when
`[typing] model_bundle` and `promotion_manifest` are both configured. Its
bundle was fitted outside production on development/calibration donors, owns
frozen probability/margin/stability gates, and supplies whole-development-
donor bootstrap stability. The public production API rejects an unsigned or
unpromoted bundle. Candidate evaluation uses the visibly distinct
`typing_mode=probabilistic_candidate` lane and cannot authorize production.
The ordinary type stage never fits or retunes a model.

Rules-only softmax values are explicitly uncalibrated scores. Promoted-bundle
values are donor/class-balanced, held-out-temperature-scaled probabilities,
not expected population prevalence. Empirical marker tail probabilities and
p-values remain audit quantities and are never substituted for model support.
Every output row stamps its score semantics and model provenance.

### Broad decision

For each cell:

- geometry-QC-ineligible cells retain their ranked audit alternatives but
  receive `Unavailable` broad/specific labels, `unavailable` statuses, and no
  accepted confidence;
- one non-conflicted authoritative anchor accepts its broad class as
  `anchor`;
- conflicting anchors or an anchor/exclusion conflict become `ambiguous`;
- `Other` is emitted only when every required broad marker model is valid and
  negative;
- missing or invalid required models without positive evidence produce
  `unavailable`;
- a non-anchor inference is available only through a promoted bundle and
  requires valid positive support, a unique winner, and all frozen
  probability, margin, and donor-bootstrap stability gates;
- otherwise the result remains `ambiguous`.

The type stage does not invent stability from the cells it is classifying.
Without the exact verified bundle/promotion pair, non-anchor cells remain
`ambiguous` or `unavailable`; with it, `inferred` is emitted only after every
frozen gate passes. A reviewed label never overwrites a production cell.

The output retains the best broad candidate and ranked probabilities even when
acceptance fails. This makes ambiguity useful for review and future rules
rather than a terminal information loss.

### Parent-constrained subtype decision

Subtype rules run only within an accepted `anchor` or `inferred` broad parent.
A subtype cannot pull a cell across broad compartments. If a parent has subtype
rules but none can be accepted, the cell becomes
`<Parent>-unclassified`; a parent with no subtype model retains its broad name.

The current overlay declares:

| Broad parent | Specific candidates |
|---|---|
| Immune | `T_cell`, `B_lineage`, `Macrophage`, `Myeloid_APC`, `Neutrophil` |
| Epithelial | `Beta`, `Alpha`, `Delta` |
| Endothelial | `Lymphatic_endothelial`, `CD34_endothelial` |
| Mesenchymal | `SMA_positive_mesenchymal` |
| Neural | terminal broad label; no specific model |

Two rules are locked in validation:

- SST is active only for the Delta subtype after Epithelial is accepted. It is
  calibrated against Pan_Cytokeratin and cannot create broad Epithelial
  evidence.
- CD11b is active Immune support/subtype evidence and is calibrated against
  E-cadherin.

The full output preserves both broad and subtype status, reason, confidence,
margin, stability, per-gate pass/failure fields, anchor/conflict evidence,
alternatives, score semantics, typing mode, model provenance, registry
fingerprint, and typing-rules fingerprint. The run and type-stage manifests
separately bind the verified promotion artifact.

## Stage 8: keep process state orthogonal

Identity answers what a cell is; state answers what it is doing. Combining the
two creates contradictory labels and lets transient biology change lineage.

The current state stage reports Ki67 and PCNA proliferation separately. Each
available marker receives a donor-local `log10-mean-plus-k-SD` state threshold
with current defaults `k = 3` and at least 200 positive finite measurements.
Only estimation-eligible cells fit the threshold, and only analysis-eligible
cells receive positive/negative state calls; other cells are unavailable. This
method is explicitly not presented as the mutually-exclusive,
frequency-neutral identity calibration. It never changes `broad_type` or
`specific_type`.

State output lives in `cell_states`; model summaries live in
`audit/state_models.parquet`.

## Stage 9: uncertainty-preserving QuPath handoff

The export stage joins assignments to image context by UUID and writes:

```text
qupath_export/cell_assignments_<donor>.csv
```

The handoff contains broad and specific identity, separate broad/specific
statuses, confidence, best broad/specific candidates, and separate ranked
alternatives for both levels. Compatibility aliases mirror the authoritative
broad/specific fields.

An importer or reviewer must inspect `assignment_status`; it must not turn
`Ambiguous` or `Unavailable` into a confident named class. State is not folded
into this identity export and remains a separate artifact. Missing image
context, duplicate UUIDs, or an altered export file is fatal.

Use
[`import_cell_assignments.groovy`](../scripts/groovy/import_cell_assignments.groovy)
for the QuPath import. First require a successful pipeline `status` so the
export file fingerprint and typing stage are current. The script validates the
exact image, UUID uniqueness, and every CSV assignment in a first pass, then
imports either the broad or specific decision. QuPath detections outside the
canonical ingest universe are intentionally left unchanged. Ambiguous and
unavailable results receive explicit uncertainty classes rather than the best
candidate.

## Content-addressed execution

Before writing output, the pipeline resolves a run ID from:

- cohort content ID;
- marker-registry fingerprint;
- typing-rules fingerprint;
- identity-affecting configuration;
- relevant source-code byte hash.

The ID is the first 20 hexadecimal characters of the resulting SHA-256. All
artifacts are written below:

```text
data/runs/<run_id>/
```

Every stage manifest records its exact upstream manifest IDs, method,
configuration, expected donors, relevant source hash, output schema, file
fingerprints, row counts, UUID universe, and fingerprints for its
availability/model-audit sidecars. Stage files are immutable.

If all eight analytical manifests are current, the pipeline writes `run.json`
and updates `data/runs/LATEST`. QuPath export has its own manifest and does not
change the analytical run identity.

This design gives three useful outcomes:

- an identical invocation validates and reports `[current]`;
- a legitimate scientific or source change resolves to another run directory;
- a partial or silently altered directory raises instead of being accepted.

Routine CLI status uses fast validation: unchanged path, size, and modification
time reuse the captured digest, while changed metadata triggers a rehash.
Archival audits can call manifest/cohort validation with `mode="content"` (or
`validation_mode="content"`) to recompute the recorded full or explicitly
sampled digest.

## Running, resuming, and checking

Run the complete workflow:

```bash
python -m phenocycler.pipeline run --config config.ini --pipelined
```

Run a dependency prefix:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --pipelined --through calibrate
```

Run exact named stages whose dependencies already have current manifests:

```bash
python -m phenocycler.pipeline run --config config.ini \
  --only type states
```

For the narrow case where `duplicate_policy = fatal` stopped geometry at the
first donor containing duplicated islet detections, reuse the verified ingest
and completed duplicate-free geometry prefix in the corrected
`deterministic_keep` run:

```bash
python -m phenocycler.pipeline resume \
  --config config.ini \
  --from-run <failed_run_id> \
  --from-donor <first_unfinished_donor>
```

This is not an unchecked partial-output override. It requires identical ingest
code, cohort input, ingest filter, donor set, and ingest content; requires the
only geometry configuration change to be `fatal -> deterministic_keep`;
verifies each prefix donor's exact row/UUID universe, nucleus-exclusion
semantics, and absence of duplicate groups; and records both reuse operations
in fingerprinted receipts. Ingest is hard-linked, not copied. Prefix geometry
values are preserved while its two all-null duplicate fields are re-encoded
with stable Parquet types. The default is to continue through the donor task
queue, so reusable-prefix REDSEA work can overlap geometry beginning at the
unfinished donor. Use `--through geometry` for a bounded resume.

For a stopped run that predates donor receipts but has a complete, current
geometry manifest, use the separate recovery contract:

```bash
python -m phenocycler.pipeline recover \
  --config config.ini \
  --from-run <stopped_run_id>
```

The source process must be fully stopped first. Recovery validates the source
ingest and geometry manifests, producing code, complete scientific
configuration, donor set, row counts, schemas, and UUID universes. It then
hard-links the complete geometry dataset and publishes one current geometry
receipt per donor. An unmanifested REDSEA donor is accepted only when its
single final `data_0.parquet` opens cleanly, passes a full-column scan, has the
exact current schema and donor identity, and matches that donor's ingest
row/UUID universe. Source scratch and incomplete donor directories are
ignored. The accepted REDSEA donors receive receipts and can enter expression
as soon as the queue starts; unfinished donors return to REDSEA.

Add `--prepare-only` to perform and record the handoff without launching work.
Recovery writes `geometry_recovery.json` and, when applicable,
`redsea_recovery.json`; both are fingerprinted by the donor receipts.

### Donor receipts and task resources

Ingest remains a cohort barrier because the shared measurement exports are
scanned together. Every later task is donor-local:

```text
geometry[d] -> redsea[d] -> expression[d] -> controls[d]
                                  │              └─ calibrate[d] -> type[d]
                                  └─ states[d]
```

A downstream task becomes ready only after each required upstream donor
receipt validates. A receipt records its donor, method, scientific
configuration, relevant code hash, upstream receipt or manifest identifiers,
output schema, row count, UUID-universe hash, file fingerprints, and any
donor-local audit sidecar. Writes use a temporary file and atomic rename; the
receipt is published last. Cohort stage manifests remain the archival
completion contract. They are sealed from the complete receipt set and
fingerprint every receipt. During an active queue, the calibration review
notebook may instead read a completed donor only after validating that donor's
receipt tree; it uses donor-local audit sidecars until the sealed cohort audit
is available. This read path does not authorize notebook-side stage execution
or make a partial cohort complete.

The queue has three dedicated process pools and enforces both a total task
limit and group limits:

```ini
[compute]
n_jobs = 4
geometry_workers = 2
redsea_workers = 1
downstream_workers = 4
use_gpu = false
gpu_device = 0
```

At most `n_jobs` tasks run at once. The group settings are ceilings within that
total, so this example permits two memory-heavy geometry tasks, one REDSEA
task, and any remaining slot for downstream work. When `use_gpu = true`, the
REDSEA pool is forced to one worker for the configured device even if
`redsea_workers` is larger. Separate pools keep the GPU context in the REDSEA
worker and keep geometry memory pressure bounded.

`status` reports `IN PROGRESS n/N donor receipts` until a cohort stage is
sealed. Re-running the same pipelined command validates and skips current
receipts. If all receipts exist but execution stopped before the aggregate
manifest was written, the next invocation seals that manifest without
recomputing donors.

Check every stage and the QuPath export:

```bash
python -m phenocycler.pipeline status --config config.ini
```

`status` returns nonzero if anything is missing or stale. Build the handoff
separately with:

```bash
python -m phenocycler.pipeline export --config config.ini
```

Run `status` immediately before a standalone export and proceed only if every
analytical stage is `CURRENT`; the QuPath line may be `MISSING`. Run it again
after export and require success before QuPath import. This validates the
typing stage and fingerprinted export as one handoff review.

This means an intentional `--no-export` analytical run is complete at
`run.json`, but `status` still reports the optional handoff missing and returns
nonzero until it is exported.

`--only` is intentionally literal; it does not guess or rerun prerequisites.
A missing required manifest stops the selected stage.

## Failure and recovery contract

| Failure | Meaning | Correct response |
|---|---|---|
| source fingerprint changed | raw export no longer matches the cohort manifest | verify the intended source and recreate the cohort manifest |
| cohort donor set differs | incomplete, extra, excluded, or duplicated donor record | correct the cohort specification |
| unmapped `Image` value | CSV contains an image not owned by the manifest | correct the exact `image_id` or source export |
| CSV parser rejects | one or more source rows could not be parsed | inspect `ingest_rejects.parquet` and fix/re-export the source |
| missing/duplicate geometry UUID | canonical geometry universe is invalid | re-export or correct segmentation/geometry |
| anisotropic or mismatched pixel calibration | current REDSEA raster units are not valid | correct the manifest/source calibration; do not rescale silently |
| acquired source column missing | panel says marker exists but its required compartment does not | correct channel mapping, registry, or upstream export |
| qptiff channel absent from marker registry | a channel has no declared biological/measurement policy | review and add a versioned registry entry or correct the channel map |
| REDSEA/QuPath alignment failure | UUID, channel, scale, or raster units do not agree | correct the manifest/export; do not accept the corrected values |
| insufficient reference controls | donor/reference cannot support the fixed estimator | retain indeterminate results; review panel/reference design |
| insufficient tail resolution | control count cannot resolve the fixed marker alpha | retain indeterminate results; do not loosen alpha post hoc |
| excess contamination | incompatible-reference controls contain too much target signal | audit biology/pairing and retain indeterminate results |
| marker panel absent | marker was not measured | retain unavailable state |
| required typing evidence invalid | broad absence cannot be established | retain unavailable/ambiguous; do not call `Other` |
| cohort output exists without manifest | interrupted or foreign stage-mode partial artifact | inspect it; never bless it by fabricating a manifest |
| donor output exists without donor receipt | interruption occurred after the atomic data write but before provenance publication | quarantine that exact donor partition and its donor audit, then rerun the pipelined command |
| fatal duplicate policy stopped at the first duplicate-bearing donor | ingest and a duplicate-free geometry prefix may be reusable | use the verified `pipeline resume` contract with the source run and first unfinished donor |
| stopped pre-receipt run has a complete geometry manifest and some complete REDSEA donors | the old stage-barrier process cannot attach to the donor queue | stop every old worker, then use the verified `pipeline recover` contract |
| manifest or output stale | bytes, schema, UUID universe, config, or code differ | correct the cause and create the newly addressed run |
| audit/availability sidecar stale | a stage-owned diagnostic changed or disappeared | treat the stage as stale; do not regenerate the sidecar independently |

Never repair a contract failure by copying a donor partition from another run,
editing a generated manifest, reallocating alpha to available markers, or
adding a notebook-only cutoff. Cross-run reuse is permitted only through the
narrow verified resume command above, which creates a newly addressed run and
fingerprinted provenance receipts.

For a verified transient interruption that left output before its manifest,
quarantine the exact unmanifested stage directory and its stage-owned
sidecars, then rerun the same stage. Do not remove current upstream manifests
or broadly clear the run root. Under pipelined execution, a current donor
receipt is resumable and must not be removed; only an exact orphan partition
without a receipt is quarantined. The aggregate sidecar ownership is:

| Stage | Sidecars |
|---|---|
| ingest | `ingest_rejects.parquet`, `panel_availability.json` |
| expression | `audit/expression_availability.parquet` |
| controls | `audit/reference_control_models.parquet` |
| calibrate | `audit/marker_calibration_models.parquet` |
| states | `audit/state_models.parquet` |

## Registry and rule change control

A proposed marker or typing change should be made in the owning versioned file,
not scattered across notebooks or stage code.

For a marker-registry change, review:

- biological kind and broad/subtype/state use;
- authoritative compartment and spillover policy;
- active versus deferred calibration;
- incompatible reference and its feasibility in each panel;
- fixed family membership and its effect on marker alpha;
- panel-absent behavior;
- SST subtype-only and CD11b/E-cadherin validation invariants.

For a typing-rule change, review:

- parent class;
- anchors versus support evidence;
- exclusions and likely co-expression conflicts;
- required marker availability;
- whether the rule is broad or parent-constrained subtype;
- validation evidence and rules-version change.

Both file fingerprints enter the run ID, so a reviewed change cannot silently
reuse the old assignments.

## Review checklist

Before using a run for biological inference:

- `status` reports every required stage current.
- `run.json` and all eight analytical manifests exist.
- expected donor counts and UUID universes agree across stages.
- `panel_availability.json` matches the intended acquisitions.
- expression-availability audit has no unexplained acquired-source failures.
- geometry exclusions and estimation-only differences have been reviewed.
- reference-control counts and spatial coverage are plausible for every donor.
- calibration-status and contamination summaries are reviewed per donor and
  marker.
- indeterminate and unavailable fractions are reported, not discarded.
- broad anchor conflicts and ranked alternatives are reviewed.
- `Other` cells meet the all-required-valid-and-negative contract.
- subtype labels occur only inside accepted parents.
- state results have not been merged into identity.
- the QuPath review uses UUID and preserves assignment status.

## Migration

The former RESTORE and ordered-residual lineage outputs are retired. They are
not upstream inputs, compatibility fallbacks, or valid threshold sources for
this workflow. Rebuild from the declared QuPath cohort rather than mixing old
normalized or gated files into a content-addressed run.
