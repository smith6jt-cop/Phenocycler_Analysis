# Error-aware refinement workflow for multiplex proteomic images

## Outcome

This workflow turns raw multiplex images, masks, and per-cell measurements into
cell identities through a sequence of **frozen, plotted decisions**. It is
designed for whole-slide data with millions of cells, for which manual review
is valuable but cannot be the computational engine.

The central operating rule is:

> Optimize a stage only against evidence that can measure that stage's error,
> using development donors; preserve uncertainty when the evidence cannot
> distinguish the alternatives; and test the frozen workflow once on donors
> that were never used for tuning.

`Ambiguous` is therefore not a temporary synonym for “pick the highest score.”
It is an abstention. A cell receives a named final class only after the source
of uncertainty is corrected and the **entire frozen pipeline is rerun**. A
manual review label remains detached reference data and never overwrites a
production assignment.

This design does not require QuPath. The same contracts can begin with an
OME-TIFF/qptiff/OME-Zarr image, a labeled mask image or polygons, and a table
containing stable `image_id` and `object_id` values. Napari, Mantis, or
cytomapper can provide image-linked review.

## Start with the existing repositories

The most useful pieces in `smith6jt-cop` should be reused in this order:

1. [Phenocycler_Analysis](https://github.com/smith6jt-cop/Phenocycler_Analysis)
   supplies immutable manifests, stable object identity, donor-local controls,
   REDSEA integration, calibrated marker states, hierarchical typing, and the
   staged review notebooks. It remains the provenance backbone.
2. [KINTSUGI](https://github.com/smith6jt-cop/KINTSUGI) supplies the strongest
   local examples of segmentation overlays, pre/post-REDSEA biaxial plots,
   single-positive preservation, morphology QC, cluster heatmaps, UMAPs, and
   spatial recoloring. Its diagnostic Otsu/quantile gates must not replace the
   registered donor-control calibration.
3. [CyLinter](https://github.com/smith6jt-cop/cylinter) supplies the cached,
   restartable human-in-the-loop pattern for finding image, segmentation, and
   intensity artifacts.
4. [REDSEA](https://github.com/smith6jt-cop/REDSEA) supplies boundary-aware
   lateral spillover correction. Its acceptance test must measure both the
   removal of implausible neighbor signal and retention of real dim positives.
5. [ark-analysis](https://github.com/smith6jt-cop/ark-analysis) supplies Mesmer
   and Pixie notebooks, pixel/cell cluster heatmaps, overlays, and interactive
   cluster merging. Use these as discovery and error-detection comparators.
6. [scimap](https://github.com/smith6jt-cop/scimap) supplies scalable AnnData
   handling, visual marker gating, rule-table phenotyping, clustering, and
   spatial analysis. It is a good transparent comparator to the local rules.
7. [InstanSeg](https://github.com/smith6jt-cop/instanseg) and
   [MCMICRO](https://github.com/smith6jt-cop/mcmicro) are candidate segmentation
   and image-processing implementations, not sources of phenotype truth.
8. [UTAG](https://github.com/smith6jt-cop/utag),
   [BANKSY](https://github.com/smith6jt-cop/Banksy_py),
   [mistyR](https://github.com/smith6jt-cop/mistyR), and
   [Squidpy](https://github.com/smith6jt-cop/squidpy) can test spatial
   plausibility after a spatial-free identity call. They must not erase a real
   isolated rare cell merely because its neighbors differ.

Do not revive the retired RESTORE workflow or reuse desired abundance,
cross-donor uniformity, disease status, an old cluster, or an old cell label as
an implicit gate. RESTORE is useful only for its mutually exclusive-control
diagnostic ideas.

## Evidence base

The design combines the modular/provenance principles of
[MCMICRO](https://doi.org/10.1038/s41592-021-01308-y) and
[MITI](https://doi.org/10.1038/s41592-022-01415-4); the image and cell QC model
of [CyLinter](https://doi.org/10.1038/s41592-024-02328-0); the whole-cell
segmentation baseline in [Mesmer](https://doi.org/10.1038/s41587-021-01094-0);
and the spillover preservation tests in
[REDSEA](https://doi.org/10.3389/fimmu.2021.652631).

Cell-typing candidates occupy different roles rather than competing as a
single undifferentiated leaderboard:

- [scimap](https://doi.org/10.21105/joss.06604) provides transparent visual
  calibration and hierarchical rules at very large scale.
- [GammaGateR](https://doi.org/10.1093/bioinformatics/btae356) provides a
  reproducible soft-gating comparator when its distributional assumptions fit.
- [Nimbus](https://doi.org/10.1038/s41592-025-02826-9) is an image-aware
  marker-positivity candidate; it does not define the final ontology.
- [Astir](https://doi.org/10.1016/j.cels.2021.08.012) is a scalable,
  prior-marker probabilistic baseline with an unknown class.
- [CELESTA](https://doi.org/10.1038/s41592-022-01498-z) adds spatial context
  without requiring a labeled reference. It is a comparator, not its own
  validation, because the same neighborhood cannot both make and prove a call.
- [STELLAR](https://doi.org/10.1038/s41592-022-01651-8) transfers labels over a
  spatial graph and can discover novel states, but requires representative
  donor-disjoint reference labels.
- [CellSighter](https://doi.org/10.1038/s41467-023-40066-7) uses image crops and
  masks and is useful when enough local image-level labels exist.
- [STARLING](https://doi.org/10.1038/s41467-024-55214-w) is valuable for
  segmentation-aware discovery and mixed-profile flags; it does not itself
  name canonical populations.
- [Pixie](https://doi.org/10.1038/s41467-023-40068-5) uses pixel composition to
  reduce failures of integrated-mean expression and is a strong orthogonal
  discovery comparator.

Published head-to-head results are study-specific. The
[PHLEX/TYPEx comparison](https://doi.org/10.1038/s41467-024-48870-5) reinforces
the need to optimize and validate on the local panel and specimens rather than
selecting a method from its headline F1 score.

For pancreas, the cell hierarchy and review policy should also reflect the
[Damond et al. type 1 diabetes IMC study](https://doi.org/10.1016/j.cmet.2018.11.014):
hormone loss can precede beta-cell loss. Low INS alone must not automatically
make a cell “not beta,” while apparent hormone co-positivity in a dense islet
must be checked for merged masks and lateral spillover before it is described
as polyhormonal biology.

## Architecture and information flow

```text
immutable source images + masks + panel metadata
  -> development refinement namespace (never data/runs/<run_id>)
  -> donor-disjoint schema-v2 split manifest with preregistered review and release governance
  -> stage-specific candidate comparisons and authorized labels on development donors
  -> preliminary fit + content-validated calibration assignments
  -> prediction-blinded calibration probability review + detached challenge audit
  -> reviewed, calibration-only threshold provenance + freeze without refitting
  -> standalone frozen validation, then explicitly opened locked-donor evaluation
  -> structured broad/subtype reports + policy-derived promotion manifest
  -> new immutable production run
```

The refinement namespace should contain only derived diagnostics, proposed
models, and detached annotations. It must never modify a source stage or its
manifest. The existing `gold_standard_v1/v2` bundles remain untouched audit
artifacts. New development labels and each review split require separately
authorized, content-addressed provenance before annotations may fit, tune, or
release a candidate.

Plotting and exploratory comparison helpers live in top-level `refinement.py`,
deliberately outside the production `phenocycler/` package. The release path is
stricter: `candidate_evaluation` creates complete, non-production assignment
artifacts; `threshold_selection` derives the calibration operating grid from
their exact reviewed rows; and `evaluation_reporting` derives validation and
locked-test metrics from exact assignment/sample/ledger artifacts. The
fail-closed promotion mechanism is described in
[PROBABILISTIC_RESCUE.md](PROBABILISTIC_RESCUE.md). Creating and promoting a
real learned artifact remains a later, explicit decision that correctly creates
a new content-addressed run. Production requires the final threshold-provenanced
bundle and its separate `TypingModelPromotionManifest` together. The `type`
stage also requires the new run's sealed cohort `calibrate` output to have
exactly the same semantic content hash as the calibration artifact recorded
during fitting; the donor scheduler waits at that cohort barrier.

### Implemented probabilistic release sequence

The preliminary model and the final frozen bundle are different immutable
artifacts. First fit once with explicit seed gates and the non-promotable flag:

```bash
python -m phenocycler.pipeline probabilistic fit \
  --evidence data/runs/<source-run>/marker_evidence \
  --labels development-labels.parquet \
  --splits donor-splits.json \
  --source-run-manifest data/runs/<source-run>/run.json \
  --label-provenance development-label-provenance.json \
  --evidence-provenance development-evidence-provenance.json \
  --bundle-version pancreas-preliminary-v1 \
  --unfrozen-candidate \
  --probability-threshold 0.70 \
  --margin-threshold 0.10 \
  --stability-threshold 0.80 \
  --out model-bundles/pancreas-preliminary-v1.json
```

The three CLI values above set the inferred probability, margin, and stability
seed gates; the preliminary bundle also records the fixed
`anchor_probability=0.95` and `negative_probability=0.20` defaults. All five
preliminary thresholds are bound by schema-v3 threshold provenance. The fit
entry point creates only an unfrozen, non-promotable candidate; it cannot attach
threshold provenance or create a final bundle. Generate a
complete, non-production calibration assignment artifact next:

```bash
python -m phenocycler.pipeline probabilistic evaluate \
  --model-bundle model-bundles/pancreas-preliminary-v1.json \
  --evidence data/runs/<source-run>/marker_evidence \
  --source-run-manifest data/runs/<source-run>/run.json \
  --output-root refinement/calibration-assignments-v1 \
  --manifest refinement/calibration-assignments-v1.json \
  --split calibration
```

After prediction-blinded review, call
`build_threshold_selection_provenance()` with that manifest, the current
marker/rule registries and source-run manifest, and target-keyed samples,
ledgers, and `ReviewLabelProvenance` objects for broad plus every fitted
subtype. Pass assignment-run-bound private sampling audits as
`calibration_sample_frames`, the exact projections shown to reviewers as
`calibration_reviewer_frames`, and the committed review key from its protected
file; neither frame can substitute for the other. Its schema-v3 artifact
records all five preliminary gates and derives and embeds the complete
operating grid and selected metrics from those exact rows; callers cannot
provide a result table or assert a selected metric. Each grid minimum must be
at least as strict as its corresponding inferred seed gate. Freeze it onto the
same scientific fit:

```bash
python -m phenocycler.pipeline probabilistic freeze \
  --model-bundle model-bundles/pancreas-preliminary-v1.json \
  --threshold-provenance model-bundles/pancreas-threshold-selection-v3.json \
  --bundle-version pancreas-typing-v1 \
  --out model-bundles/pancreas-typing-v1.json
```

This changes the full bundle content ID while preserving the
threshold-independent `model_fit_content_id`; coefficients, temperatures, and
bootstrap members are not refit. Use `probabilistic evaluate` again with the
frozen bundle and `--split validation`, then build a
`SplitEvaluationReport` with `build_split_evaluation_report()`. Repeat for
`--split locked_test --allow-locked-test` only after validation is frozen, and
pass `allow_locked_test=True` to the report builder. Finally,
`TypingModelPromotionManifest.create()` accepts a prespecified
`PromotionGatePolicy`, both structured reports, exact calibration source-replay
inputs, current registries, and protected review/release keys. It replays the
frozen calibration decision and both held-out sources; approval is derived from
confidence-bound total/rescue and per-class FP/FN gates rather than supplied by
the caller.

The locked-test option is an explicit governance acknowledgement, not a
cryptographic or stateful one-use lock. The study protocol and access controls
must enforce the single opening, archive its immutable artifacts, and prohibit
candidate regeneration after the result is inspected. See
[PROBABILISTIC_RESCUE.md](PROBABILISTIC_RESCUE.md) for the exact Python builder
calls and complete lineage checks.

Every notebook begins by validating the source content IDs and split
fingerprint. Every decision writes one compact row to a decision ledger:

| Field | Meaning |
|---|---|
| `checkpoint` | image QC, segmentation, spillover, calibration, broad type, subtype, or release |
| `candidate_id` | content fingerprint of parameters/model/code |
| `split` | development fold, validation, or locked test |
| `metric` | prespecified quantity; never desired abundance |
| `estimate`, `ci_low`, `ci_high` | donor/ROI-clustered estimate and interval |
| `gate`, `direction` | frozen acceptance criterion |
| `decision` | accept, reject, inspect, or unable-to-decide |
| `rationale` | short root-cause statement |
| `review_bundle_id` | fingerprint of detached labels, when used |

## Efficient notebook suite

Notebooks 01–06 are operator interfaces, not alternate algorithm
implementations. Notebooks 01–04 retain their existing guarded package-stage
calls, notebook 05 is read-only, and notebook 06 defaults to an inert,
output-free template. Notebook 06 may write only new immutable evaluation
reports and a promotion manifest when its explicit release flags are enabled.

| Notebook | Purpose | Existing base |
|---|---|---|
| Proposed `00_contract_and_donor_splits` | freeze IDs, panel completeness, development/calibration/validation/locked-test donors, metrics, and error costs | manifests and artifact hashes; no checked-in notebook yet |
| [01_qupath_ingest_and_geometry_qc](../notebooks/01_qupath_ingest_and_geometry_qc.ipynb) | review raw ingest, geometry availability, segmentation exclusions, and traced ROIs | immutable ingest/geometry artifacts + KINTSUGI/CyLinter patterns |
| [02_redsea_spillover](../notebooks/02_redsea_spillover.ipynb) | compare compartment spillover correction while protecting true positives | REDSEA/KINTSUGI plots |
| [03_marker_calibration](../notebooks/03_marker_calibration.ipynb) | review positive/negative/indeterminate marker evidence and operating-point uncertainty | donor-local calibration + Nimbus/GammaGateR comparators |
| [04_hierarchical_cell_typing](../notebooks/04_hierarchical_cell_typing.ipynb) | inspect missing populations, simultaneous broad typing, parent-constrained subtypes, and QuPath handoff | scimap/Astir/Pixie/STARLING comparisons remain development-only |
| [05_error_aware_refinement](../notebooks/05_error_aware_refinement.ipynb) | build broad-type probability/challenge samples and plot the frozen broad validation point; it does not evaluate subtype or all-cell specific targets | detached `refinement.py` helpers |
| [06_probabilistic_validation_and_promotion](../notebooks/06_probabilistic_validation_and_promotion.ipynb) | freeze and gate the full-target validation report, explicitly open the locked test once, replay the exact calibration and held-out sources, and create the signed candidate-bound promotion manifest | manifests and decision ledger |

For a 31-million-cell cohort, notebooks must not load whole donor tables unless
the plot truly needs cell-level data. Use:

- projected Parquet/DuckDB aggregation for counts and quantiles;
- a deterministic, donor/ROI-balanced sketch for scatterplots and UMAP;
- chunked full-cohort prediction only after a candidate is frozen;
- on-demand pyramid tiles for image galleries;
- content-addressed cached summaries, so changing a plot style does not rerun a
  scientific stage;
- donor or ROI bootstrap units, never millions of cells treated as independent
  replicates.

Notebook output should be stripped before version control. The durable output
is the compact summary/decision table, not a multi-megabyte embedded figure.
Notebook 05 is intentionally broad-only. Create a separate blinded sample,
adjudicated ledger, provenance artifact, plots, and structured report for each
fitted subtype parent; a broad ledger cannot stand in for subtype evidence. It
uses the library HMAC sampler, exports only the reviewer projection, passes the
private audit frame and source-ingest review-context content ID to provenance
validation, and never opens locked-test labels. Promotion also requires a
separate all-cell `level="specific", parent=None` review/report for validation
and locked test, so notebook 05 cannot satisfy promotion by itself.
Notebook 06 is the output-free operator template for that complete release
sequence. Its action flags default to false; validation must pass the exact
preregistered policy before the locked-test flag can be opened, and its
challenge-sample counts remain diagnostic rather than headline metrics. When a
challenge diagnostic is supplied, notebook 06 requires its exact validation
candidate/assignment lineage and an explicit disposition; a catastrophic
finding stops that candidate and starts a new development cycle.

## Plotted decisions by stage

### 0. Contract, panel, and donor splits

**Plots**

- sample-by-marker availability/status heatmap;
- donor × disease/batch/slide mosaic to expose confounding;
- cell yield and exclusion flow by donor;
- train/calibration/locked-test coverage of tissues, batches, and rare classes.

**Decision**

Split by donor before any tuning. Keep all images/ROIs from a donor in one
split. If batch and biology are fully confounded, record that limitation rather
than “correcting” it away. Freeze the intended metrics, review sample sizes,
and per-class FP/FN costs now. The schema-v2 `DonorSplitManifest` must record a
positive `review_n_per_stratum` separately for calibration, validation, and
locked test, and the SHA-256 commitment of a secret key containing at least 16
bytes. Keep that key in a protected file outside the repository; notebooks
receive only a key-file path or environment-variable placeholder and must never
print or embed the secret. The review sizes and commitment participate in the
split content ID and therefore cannot change after predictions are inspected.
Its `release_governance` object preregisters
`threshold_selection_policy_content_id` and `threshold_grid_content_id` before
calibration labels are opened. It also records
`promotion_policy_content_id`, the content ID of the release policy frozen
before evaluation, and `signing_key_sha256`, the commitment to a separate
protected release-signing key. The latter two may be blank for non-production
development, but both must be nonblank for production promotion. Keep the
review and signing keys external and distinct; only commitments belong in the
split, and no notebook, manifest, log, or reviewer export may contain either
secret.

### 1. Image, registration, and segmentation QC

**Plots**

- multichannel raw thumbnails and negative-control channels;
- saturation, background, signal-to-background, focus, tissue-loss, and
  registration-displacement distributions by slide/cycle;
- candidate segmentation overlays on the same stratified ROIs;
- object precision/recall, Panoptic Quality, split/merge counts, nuclei per
  cell, tissue coverage, area, solidity, eccentricity, and cell density;
- galleries sampled from each exclusion reason, not only the worst-looking
  objects.

**False positives** are debris, fragments, duplicate masks, and one cell split
into several objects. **False negatives** are missed cells, over-merged cells,
weak nuclei, and cell types systematically lost by morphology filters.

**Decision**

Choose the segmentation/QC candidate on manually traced, donor-held-out ROIs
that span density and morphology. Require both object-level precision and
recall; a high average IoU alone does not reveal systematic missed rare cells.
Store exclusion flags and reasons instead of deleting rows.

The current source run already makes this stage a priority: 4,737,702 objects
(`15.06%`) have missing nuclei, and total geometry analysis eligibility is
`79.96%`. These exclusions need an image-level, risk-stratified review before
they can be assumed to be harmless debris.

### 2. Quantification, channel spillover, and lateral spillover

Treat detector/isotope crosstalk and lateral neighbor contamination as
different mechanisms. If single-stain controls exist, estimate a channel
compensation matrix separately. REDSEA then addresses boundary-localized
neighbor signal.

**Plots**

- emitter-versus-receiver slopes in single-stain or mutually exclusive controls;
- target-cell border/interior ratio and receiver signal versus source-neighbor
  signal/contact length;
- pre/post biaxial plots for mutually exclusive marker pairs;
- pre/post raw-image and mask overlays for isolated, boundary, and dense cells;
- authentic single-positive retention, implausible double-positive reduction,
  double-negative inflation, clipping fraction, and correction magnitude;
- candidate-versus-baseline downstream marker recall on the labeled audit set.

**False positives** arise when neighbor signal or channel crosstalk produces a
marker in the wrong object. **False negatives** arise when correction subtracts
real dim signal or when a merged mask is treated as a spillover problem.

**Decision**

Select a Pareto point: neighbor-associated signal and biologically impossible
co-positivity must fall, while signal in isolated/known-positive cells and
held-out marker recall are retained. REDSEA's published result is a warning
against optimizing only double-positive reduction: subtraction without the
corresponding reinforcement can inflate double negatives.

Persist alignment correlations, raster/intensity ratios, clamp fractions,
isolated-cell counts, and preservation metrics as audit artifacts. At present
several of these are printed but not retained by the local pipeline.

### 3. Normalization and marker operating points

Compare raw/log or asinh values, donor-control scaling, and—only when justified
by shared controls—more aggressive batch correction. The
[mIF normalization benchmark](https://doi.org/10.1093/bioinformatics/btab877)
shows why the locally best Pareto point matters: a method that removes batch can
also remove biology when the design is confounded.

**Plots**

- control and target densities by donor with thresholds and bootstrap intervals;
- reference-control coefficient of variation and contamination;
- positive-versus-negative control effect size;
- marker PR curves as the primary operating plot for rare positives, plus ROC;
- reliability/Brier plot when a true probability model is claimed;
- threshold-sensitivity curves and crops immediately above/below each candidate;
- donor-marker status heatmap and indeterminate fraction;
- PCA/UMAP colored separately by batch and biology, plus batch-classifier AUC.

**False positives** are controlled by true negative controls, exclusion markers,
and image-localization review. **False negatives** require independently labeled
positive examples; negative controls alone cannot estimate sensitivity.

**Decision**

Each marker returns `positive`, `negative`, `indeterminate`, or `unavailable`.
Choose a class-specific operating point on development donors using a
prespecified FP/FN cost, then lock it. Never select a threshold because it
produces the expected number of positive cells. A marker whose controls or
separation fail remains unavailable.

Nimbus can be benchmarked as an image-aware marker caller; GammaGateR/scimap as
table-based comparators. Their outputs do not become truth. Every marker and
platform is validated locally.

### 4. Discovery and hierarchical typing

Run discovery on a deterministic donor/ROI-balanced subset. Compare seeds,
resolutions, and resamples. Discovery asks “what population might be missing?”;
it does not automatically name cells.

**Plots**

- cluster-marker heatmap/dot plot and spatial recoloring;
- donor composition and cluster stability (Jaccard/ARI) across seeds/resamples;
- STARLING segmentation-error or mixed-profile probability;
- broad then parent-constrained subtype confusion matrices;
- per-class precision, recall, specificity, FPR, FNR, F1, and AUPRC with
  donor-bootstrap intervals;
- top-two margin, calibration, and coverage-versus-error curves;
- method-disagreement UpSet/table for local rules, scimap/Astir, and an
  orthogonal Pixie/STARLING discovery result;
- image galleries stratified by true positive, false positive, false negative,
  ambiguous, unavailable, rare class, donor, and segmentation risk.

**Decision**

Use simultaneous broad classes first. Evaluate subtype only inside an accepted
parent. Require compatible positive anchors and exclusion/conflict checks.
Choose a candidate by nested donor-blocked validation, emphasizing worst-donor
behavior and rare-class PR rather than micro-accuracy. Spatial methods remain
an audit unless they are themselves evaluated on donor-held-out labels; if a
spatial method helps make a call, spatial coherence is no longer independent
validation for that call.

For the implemented probabilistic rescue lane, fitting is additionally
fail-closed:

- training evidence must resolve through the exact source `RunManifest` and
  its sealed `calibrate` `StageManifest`; a free-standing directory fingerprint
  is not sufficient;
- temperature scaling uses calibration donors and fails for uniform,
  boundary-optimal, or locally nonidentifiable solutions;
- exactly the requested number of whole-development-donor bootstrap draws are
  recorded. Class-incomplete and nonconverged draws are marked invalid, kept in
  the stability denominator, and are never retried away; at least 80% must be
  valid. Production rows stamp both
  `typing_broad_bootstrap_replicates` and
  `typing_broad_valid_bootstrap_replicates`, with corresponding parent-specific
  subtype totals and valid counts; and
- a one-class subtype parent remains explicitly structural/rules-only because
  named-versus-alternative probability cannot be estimated from one class.

Calibrate probability, top-two margin, and stability jointly through
`build_threshold_selection_provenance()`. The builder reads one current,
content-validated calibration assignment artifact and verifies it against exact
assignment-run-bound sampling-audit frames, reviewer frames, and adjudicated
ledgers for broad plus every fitted subtype. The HMAC sampler takes no
caller-authored salt or sample size; both derive from the schema-v2 split and
committed secret. The threshold builder derives the operating grid itself; it
does not accept a
caller-authored grid, selected row, or metric values. Hold authoritative
anchors and `Other` fixed, rethreshold only the eligible rescue lane, and
preserve the valid-positive-support and unique-winner gates at every point.
Plot total named-class coverage separately from rescue-lane coverage and
selective error; otherwise, anchor calls can hide a rescue branch that makes
zero calls. The declared total gates (`maximum_selective_error`,
`minimum_coverage`) and rescue-specific gates
(`maximum_rescue_selective_error`, `minimum_rescue_coverage`) must hold for
every modeled target at the same threshold tuple. The selector fails if no
common point meets all constraints; aggregate FP/FN cost breaks ties only among
feasible tuples. Schema-v3 provenance records all five preliminary assignment
thresholds. Its probability, margin, and stability grid minima cannot be more
permissive than the corresponding preliminary seed gates, because subtype
scores only exist after a broad parent passes those upstream gates.

### 5. Dual FP/FN review

One review sample cannot efficiently estimate population error and find rare
failure modes. Build two immutable bundles:

1. **Probability sample** — random or explicitly probability-weighted within
   donor/ROI and predicted class. It supports unbiased, inverse-probability
   weighted estimates and confidence intervals.
2. **Challenge sample** — enriched for likely errors: low margin, method
   disagreement, exclusion-positive calls, impossible co-expression, large
   spillover correction, high segmentation-error probability, dense/border
   cells, rare candidates, invalid/unavailable markers, and batch outliers. It
   discovers root causes but must not estimate cohort prevalence.

Review both predicted positives **and** predicted negative/Other/unassigned
cells. A sensitive union of independent marker, morphology, pixel-model, and
discovery candidates supplies the false-negative challenge pool. Use two
blinded reviewers and adjudication. Use `build_review_sampling_frame()` for the
probability sample. Its HMAC-SHA256 selection and opaque stratum codes bind the
split, path-independent candidate assignment-run ID, and target; exact
candidate manifest/output IDs remain separately bound by review provenance.
There is no caller salt or caller `n`. Keep its private sampling-audit frame
separate from both model outputs and the reviewer view. Its exact fields are:

```text
donor_id, object_id, sampling_stratum, sampling_method, sampling_salt,
stratum_population, stratum_sampled, sampling_weight, selection_probability
```

Do not give that audit frame to reviewers. Create their export only through
`blinded_review_sample_from_dataset(sampling_audit,
source_run_manifest=source_run_manifest, splits=split_manifest)`. The API
requires the source run to match the split and its complete donor universe,
resolves exactly one referenced `ingest` stage, and content-validates that
stage and snapshot before and after its donor-local read. It accepts neither a
free-standing snapshot nor an arbitrary location table. Its view contains only
`donor_id`, `object_id`, the secret-keyed opaque `sampling_stratum`, and image
identifiers/coordinates/region fields. It contains no weights, populations,
method, salt, selection probability, prediction, status, reason, candidate,
model score/probability, margin, stability, confidence, or gate flags. The
private audit contract verifies unique keys, complete stratum counts,
`sampling_weight = stratum_population / stratum_sampled`, and
`selection_probability = 1 / sampling_weight`.

The separately frozen adjudication ledger contains at least:

```text
donor_id, object_id, reference_label, reviewer_1, reviewer_2,
adjudicated, root_cause, evidence_sufficient
```

Schema-v5 `ReviewLabelProvenance` fingerprints that complete ledger, the
private sampling-audit frame, and the exact reviewer frame that was shown;
requires the same donor/object-key coverage in all three; and binds them to the
donor split, split purpose, current assignment run ID, target level/parent,
exact typing-model bundle, candidate manifest/output, and
`review_context_artifact_content_id`. For image review, pass
`ingest_manifest.output.content_sha256`, matching the source-run-bound snapshot
used to build the reviewer view. Calibration reviews bind the preliminary full
bundle ID; validation and locked-test reviews bind the frozen full bundle ID.
Audit weights and predictions are joined back only after this validation and
complete adjudication. A ledger from another run, model, candidate output,
review-context artifact, target, audit, reviewer frame, or split cannot
validate, and locked-test labels require an explicit release-time open.

For validation and locked-test reports, build three kinds of probability
target: broad; one conditional subtype target for every fitted probabilistic
parent; and a mandatory end-to-end `level="specific", parent=None` target.
Build the specific audit from the complete assignment frame, stratifying on
`specific_assignment_status` and `specific_type`; do not restrict it to cells
predicted into a broad parent. Its classes are the frozen final ontology
leaves. This all-cell population detects broad misrouting and false negatives
from structural/rules-only or fallback specific assignments that a conditional
subtype sample cannot see. It complements rather than replaces broad and
per-parent subtype review.

Root-cause codes should distinguish `image_artifact`, `registration`,
`segmentation_split`, `segmentation_merge`, `cell_missed`, `spillover`,
`marker_threshold`, `marker_unavailable`, `ontology_missing`,
`broad_model`, `subtype_model`, and `insufficient_evidence`. Mitigation occurs
at the earliest responsible stage, followed by a complete rerun; it is not a
per-cell override.

### 6. Locked validation and release

Freeze code, registry, image masks, segmenter, correction parameters, marker
operating points, ontology, model, and abstention thresholds before opening the
locked donor labels.

**Release plots**

- weighted confusion matrix;
- precision/recall/FPR/FNR/NPV by class with two-stage survey-bootstrap
  simultaneous intervals;
- coverage and ambiguity by class and donor;
- selective error versus coverage;
- pipeline sensitivity matrix across reasonable upstream alternatives;
- drift versus development distributions;
- abundance plots only after accuracy is frozen, as biological results rather
  than calibration targets.

Release requires all prespecified class gates and no unresolved catastrophic
root cause. If the test fails, report the failure; do not retune on the locked
donors. A new development/validation cycle needs a new future test set.

Approval is represented only by a separate immutable
`TypingModelPromotionManifest`; the strict model-bundle schema has no embedded
promotion-report field. Its schema-v3 validation and locked-test reports each
cover broad, the end-to-end all-cell specific target, and every probabilistic
conditional subtype model. Every target is bound to its own complete
`ReviewLabelProvenance`, sampling audit, reviewer frame, source-run review
context, candidate-assignment manifest/output, exact donor split, and
donor-plus-stratified-cell survey-bootstrap intervals. The fixed policy derives
pass/fail from confidence-bound total coverage, rescue coverage, selective
error, rescue selective error, and per-class FP/FN rates. Caller-supplied pass
flags and hash-only metric attestations are rejected.

`build_split_evaluation_report()` accepts the frozen candidate manifest and
one `EvaluationTargetReview(sample_frame, reviewer_frame, review_labels,
provenance)` per required target, not precomputed results. It
validates the complete source/assignment object universe before and after the
read and computes sampling-weighted metrics and simultaneous two-stage
survey-bootstrap intervals. `TypingModelPromotionManifest.create()` then embeds
both structured reports and the policy, replays the exact threshold-selection
and held-out sources, and signs the complete package/runtime-bound receipt. It
refuses missing targets, wrong-split donors, insufficient reviewed support, too
few bootstrap replicates, or any failed confidence-bound gate.

Production configuration must name both artifacts:

```ini
[typing]
model_bundle = model_bundles/pancreas-typing-v1.json
promotion_manifest = model_bundles/pancreas-typing-v1-promotion.json
```

Configuring only one fails. Before any donor is typed, the cohort calibration
manifest must be sealed and its output semantic content hash must exactly match
the source calibration content recorded in the model's evidence provenance.
The pipelined scheduler waits for that sealed manifest; preprocessing drift is
an error, not a reason to reuse the old model.

## How `Ambiguous` cells receive a final classification

```text
Ambiguous production cell
  -> assign uncertainty/root-cause flags
  -> enter probability sample and/or challenge sample
  -> blinded image + marker evidence review
  -> correct the responsible upstream rule/model/artifact for a class of cells
  -> refit only on development donors
  -> validate on donor-disjoint data
  -> one-time locked test and exact promotion manifest
  -> configure bundle + promotion and rerun every cell
  -> named class only if its acceptance rules pass
```

The outcomes are deliberately asymmetric:

- valid broad evidence plus valid subtype evidence → named subtype;
- valid broad evidence but unresolved subtype → `<Parent>-unclassified`;
- conflicting valid evidence or inadequate margin/stability → `Ambiguous`;
- required measurements/models missing or invalid → `Unavailable`;
- all required broad models valid and negative → `Other`;
- reviewer believes a class is likely but evidence remains inadequate → the
  review label is retained as reference, while production remains ambiguous.

This gives the user step-by-step control through review and versioned promotion
without turning classification into an unauditable sequence of manual edits.

## Historical preflight for the latest completed source run

Do not optimize its `confidence` or margin thresholds. The findings below are a
read-only record of completed pre-refinement run `c8412757d12de36ca056` and the
typing implementation that produced it. They explain why the new safeguards
were required; they do **not** describe the current v2 typing code:

- 31,454,041 cells;
- broad status: 12.60% anchor, 50.20% Other, 9.14% Ambiguous, 28.06%
  Unavailable, and **0 inferred calls**;
- that zero counted only rows with `broad_assignment_status == "inferred"`.
  The QuPath export still populated broad and specific display labels from
  authoritative anchors, resolved `Other`, explicit uncertainty states,
  terminal broad parents, and `<Parent>-unclassified` fallbacks. Visible
  QuPath classifications therefore do not imply that probabilistic inference
  occurred;
- the then-current `type_calibrated_cells()` supplied neither a fitted seed
  model nor grouped donor stability, so its documented inference branch was
  disconnected;
- `expression_probability` is `1 -` an empirical background-tail p-value, not
  a calibrated expression posterior. E-cadherin-negative cells have a median
  value of 0.954;
- the then-current logits consumed every finite score without masking categorical
  negative or invalid/indeterminate evidence. In a minimal fixture, an
  indeterminate marker with invalid evidence can become an inferred call if a
  stability value is supplied;
- the resulting broad “confidence” is internally inconsistent: median scores
  are 0.978 for `Other`, 0.000026 for Neural anchors, and 0.0083 for
  Endothelial anchors.

Those conditions are now enforced by the implemented
[`state-safe-threshold-support-v1` and probabilistic bundle contract](PROBABILISTIC_RESCUE.md):
background-tail scores are audit-only, invalid/indeterminate evidence is masked,
inference requires valid positive support, donor-bootstrap stability is supplied
only by an immutable opt-in bundle, and score semantics and gate failures are
stamped on every assignment. The historical run remains immutable and must not
be reinterpreted under the new semantics.

There is also a governance boundary to resolve explicitly: current repository
policy treats existing classifications and gold-standard reviews as audit-only
and forbids their use for tuning or release gates. Development and production
APIs now enforce that boundary, but they do not authorize labels themselves.
Keep both `model_bundle` and `promotion_manifest` blank until a deliberate
policy change creates a new, separately fingerprinted development, validation,
and locked-test label namespace and the candidate passes release.

## Minimum promotion criteria

A candidate can be proposed for production only when all are true:

1. object IDs and the full source universe are preserved;
2. every exclusion has an explicit reason and an audited sample;
3. segmentation precision **and** recall pass on held-out ROIs;
4. spillover correction reduces neighbor-associated false signal without
   unacceptable loss of isolated/known positives;
5. marker FP and FN operating points are validated separately;
6. invalid and indeterminate evidence cannot create an inferred identity;
7. broad, conditional subtype, and end-to-end all-cell specific weighted
   coverage/error and per-class FP/FN bounds pass; rare-class precision-recall
   curves are inspected diagnostically (or gated only if separately
   preregistered in the study SOP);
8. total named-class coverage and rescue-lane coverage are reported separately
   beside conditional accuracy;
9. the probability-sample review is complete and its software-enforced
   provenance binds the exact assignment-run-bound sampling audit, reviewer
   frame, model, candidate manifest/output, and source-run stage output used as
   the review-context artifact. If the SOP uses an enriched challenge review,
   retain its detached audit ledger with exact candidate/assignment lineage,
   explicit evidence sufficiency, root cause, catastrophic flag, and
   disposition; it cannot supply headline
   metrics or software promotion pass/fail, nor can it estimate population
   error; a catastrophic finding requires a new development cycle;
10. the locked test is run only after freezing, and its result is not used for
    another tuning pass;
11. a separate promotion manifest binds the exact candidate/split to structured
    validation and locked-test reports for broad, all-cell specific, and every
    fitted subtype target, every confidence-bound policy gate passes, and its
    policy/signing key match the split's nonblank preregistered release
    commitments; and
12. the production run's sealed calibration content exactly matches the
    evidence artifact used to train the candidate.

If no candidate meets these criteria, retaining more `Ambiguous` cells is the
correct result.
