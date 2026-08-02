# Probabilistic rescue for ambiguous cell identities

## Status and scope

The repository now contains a complete, fail-closed path for promoting a
donor-trained probabilistic model into the existing hierarchical `type` stage.
It is a rescue lane for cells that do not receive an authoritative rule anchor;
it does not replace marker calibration, the biological hierarchy, or the
`Ambiguous`/`Unavailable` abstention states.

The production default remains **rules only**. The checked-in `config.ini` has
empty `[typing] model_bundle` and `promotion_manifest` values, and this
repository does not currently ship a biologically fitted model bundle. The APIs
have synthetic contract tests, but no synthetic fixture is a production model.

A real bundle is intentionally not produced from the historical QuPath labels
or `gold_standard_v1/v2` reviews. Repository policy declares those artifacts
audit-only. Model development requires a new, explicitly authorized,
content-fingerprinted label set collected under a frozen donor split. Until
that governance step occurs, leave the optional configuration blank.

The broader image, segmentation, marker, and FP/FN review workflow is described
in [REFINEMENT_WORKFLOW.md](REFINEMENT_WORKFLOW.md). This document covers the
implemented probabilistic fitting and production boundary.

## Implemented components

| Component | Responsibility |
|---|---|
| `phenocycler.hierarchical_typing` | versioned state-safe features, constrained classifiers, hierarchical decisions, and per-gate diagnostics |
| `phenocycler.refinement_contracts` | immutable donor splits and purpose-restricted label provenance |
| `phenocycler.probabilistic_training` | donor-balanced fitting, held-out temperature scaling, whole-donor bootstrap fitting, and production-aware grid primitives |
| `phenocycler.candidate_evaluation` | content-validated calibration, validation, and explicitly opened locked-test assignments without a production promotion |
| `phenocycler.threshold_selection` | calibration-only reviewed operating-grid derivation and immutable threshold provenance |
| `phenocycler.evaluation_reporting` | weighted broad/subtype/all-cell-specific FP/FN reports with two-stage donor + stratified-cell survey bootstrap inference |
| `phenocycler.typing_model_bundle` | strict JSON schema, simultaneous release bounds, signed replay receipts, stable fit/full-bundle identities, no-refit threshold freezing, compatibility checks, model materialization, and chunked ensemble stability |
| `phenocycler.typing_stage` | optional bundle inference, abstention logic, geometry exclusion, and row-level provenance |
| `phenocycler.pipeline` / `donor_pipeline` | run identity, promotion and calibration-content gates, immutable type-stage inputs, and donor receipts |

Fitting is a development API, not a pipeline stage. This separation prevents an
ordinary production run from fitting on the cells it is about to classify.

## State-safe feature contract

The feature contract is
`state-safe-threshold-support-v1` (`FEATURE_SCHEMA_VERSION`). The same
`evidence_feature_matrix()` transformer is used for model fitting, probability
prediction, bootstrap stability, and production assignment. It accepts either
long cell-marker evidence or the wide production marker-evidence table while
preserving marker order and donor/object keys.

For each marker, evidence is resolved in this order:

1. Prefer the donor-calibrated signed `log2_threshold_ratio`.
2. Convert that coordinate to bounded positive support:

   ```text
   support = max(0, 2 * sigmoid(log(2) * log2_threshold_ratio) - 1)
   ```

   Support is zero at and below the calibrated threshold, rises monotonically
   above it, and approaches one. It is deliberately not called a posterior.
3. If the calibrated coordinate is absent, an explicit categorical `positive`
   or `negative` state supplies support one or zero.
4. `evidence_probability` is accepted only as an explicit numeric pilot input
   when neither calibrated coordinate nor categorical decision is present.
5. `empirical_tail_probability`, default `tail_probability`, and
   `expression_probability` are audit quantities. A background-tail p-value is
   never converted to `1 - p` and fed to the classifier.

The categorical and model-validity gates then take precedence over every
numeric value:

- a valid categorical negative has zero support;
- indeterminate or explicitly invalid evidence has zero model input and cannot
  support an inferred class;
- unavailable or absent evidence remains unavailable in the assignment
  contract and becomes zero only in the finite matrix used for multiplication;
- nonfinite values, malformed booleans, unknown states, duplicate cell-marker
  rows, and probability-like values outside `[0, 1]` are rejected;
- inference requires valid positive anchor/support evidence for the winning
  class in addition to probability, margin, and stability gates.

This directly closes the historical failure mode in which a high `1 - p` value
or invalid/indeterminate measurement could influence a logit. Regression tests
cover long/wide parity, categorical-negative precedence, invalid high-valued
evidence, tail-probability audit-only behavior, and the valid-support inference
gate.

## Donor split and label governance

### Freeze the split first

`DonorSplitManifest` defines four mutually exclusive sets and gives the entire
split a canonical content ID:

| Split | Permitted use |
|---|---|
| `development` | fit base coefficients and whole-donor bootstrap members |
| `calibration` | fit temperature scaling and choose frozen acceptance thresholds |
| `validation` | evaluate the frozen candidate; never read by the fitting API |
| `locked_test` | one-time release evaluation; never read by the fitting API |

Development and calibration must be nonempty. Every class at each fitted level
requires development labels from at least two donors, and both development and
calibration labels must cover that level's complete class vocabulary. This
prevents a retained bootstrap ensemble from being determined by one rare-class
donor. All ROIs and sections from one donor remain in one split.

Schema-v2 split manifests also preregister one positive
`review_n_per_stratum` for each of `calibration`, `validation`, and
`locked_test`, plus the lowercase SHA-256 commitment of a secret review key
containing at least 16 bytes. Store that key in a protected file outside the
repository and expose only its path (for example through
`PHENOCYCLER_REVIEW_BLINDING_KEY_FILE`). Never put the key itself in a
notebook, manifest, command line, log, or reviewer export. The commitment and
sample sizes are part of the split content ID, so neither can be chosen after
seeing candidate assignments.

The same schema has a `release_governance` object. Its
`threshold_selection_policy_content_id` and `threshold_grid_content_id`
preregister calibration choices before labels are opened. Its
`promotion_policy_content_id` preregisters the content-addressed
`PromotionGatePolicy`, and `signing_key_sha256` commits only the SHA-256 digest
of the separate release signing key. The promotion-policy and signing-key
fields may remain blank for a non-production development split, but **both
must be nonblank for production promotion**. Keep the review key and release
signing key as separate protected files outside the repository; only their
paths enter a process, and neither secret is printed, embedded in a notebook,
or serialized into a split, model, report, or promotion artifact.

```python
import hashlib
import os
from pathlib import Path

from phenocycler.refinement_contracts import DonorSplitManifest

review_key_path = Path(os.environ["PHENOCYCLER_REVIEW_BLINDING_KEY_FILE"])
release_signing_key_path = Path(
    os.environ["PHENOCYCLER_RELEASE_SIGNING_KEY_FILE"]
)
review_key = review_key_path.read_bytes()
release_signing_key = release_signing_key_path.read_bytes()
if len(review_key) < 16:
    raise ValueError("review blinding key must contain at least 16 bytes")
if len(release_signing_key) < 16:
    raise ValueError("release signing key must contain at least 16 bytes")

split_manifest = DonorSplitManifest(
    split_version="pancreas-donor-split-v2",
    source_run_id=SOURCE_RUN_ID,
    cohort_content_id=COHORT_CONTENT_ID,
    development=DEVELOPMENT_DONORS,
    calibration=CALIBRATION_DONORS,
    validation=VALIDATION_DONORS,
    locked_test=LOCKED_TEST_DONORS,
    review_n_per_stratum={
        "calibration": 40,
        "validation": 40,
        "locked_test": 40,
    },
    review_blinding_key_sha256=hashlib.sha256(review_key).hexdigest(),
    threshold_selection_policy_content_id=(
        PRECOMMITTED_THRESHOLD_POLICY.content_id
    ),
    threshold_grid_content_id=PRECOMMITTED_THRESHOLD_GRID_CONTENT_ID,
    promotion_policy_content_id=PRECOMMITTED_RELEASE_POLICY.content_id,
    release_signing_key_sha256=hashlib.sha256(release_signing_key).hexdigest(),
)
del review_key, release_signing_key
```

The Python constructor calls the second digest
`release_signing_key_sha256`; serialized schema-v2 JSON exposes it as
`release_governance.signing_key_sha256`. The precommitted policy must be the
same policy later supplied to promotion; changing its gates requires a new
split rather than a post-result edit.

Before production use, `resolve_run_context()` additionally requires:

- the split's `cohort_content_id` to equal the current QuPath cohort content ID;
- the four split sets together to be an exact partition of the current donor
  set—no missing or foreign donor.

The split's historical `source_run_id` remains provenance. It is not required
to equal the new production run ID because adding a promoted bundle correctly
creates a new content-addressed run.

### Authorize new labels explicitly

Development labels use stable `donor_id`/`object_id` keys and declare `level`,
`parent`, `label`, explicit `evidence_sufficient=true`, `audit_only=false`,
`label_purpose=probabilistic_typing_development`, and
`label_source=independent_blinded_review` on every row. They also require two
distinct nonblank reviewers, `adjudicated=true`, and a nonblank root-cause
code. All governance fields are included in the fingerprint; omitting them
fails.
`DevelopmentLabelProvenance.from_labels()` fingerprints
their canonical content and accepts only the explicit purpose
`probabilistic_typing_development`. Its content identity also binds the labels
to the exact donor-split content ID. Moving a donor between development and
calibration invalidates that authorization even when the label rows are
unchanged. The fitting CLI reads this pre-existing immutable provenance JSON;
it cannot mint authorization for the file it is about to fit.

`validate_development_labels()` rejects:

- any change from the recorded label fingerprint;
- labels from a different source run;
- donors outside the split;
- all validation or locked-test labels;
- duplicate or conflicting cell/level keys;
- unresolved, insufficient-evidence, or audit-only training rows.

For blinded image review, `build_review_sampling_frame()` is the only supported
probability sampler. Its HMAC-SHA256 draw is bound to the split content ID,
split name, path-independent candidate assignment-run ID, and target
`level`/`parent`. Exact candidate manifest/output IDs are review-lineage fields
but do not change the selected cells. The caller supplies neither a salt nor
`n`; the function uses the preregistered split value. Even the stratum code is
secret-keyed, so an ordinal code cannot reveal a predicted class or
disposition.

Keep the returned sampling-audit frame private. It contains only keys and
opaque strata plus `sampling_method`, `sampling_salt`, `stratum_population`,
`stratum_sampled`, `sampling_weight`, and `selection_probability`.
`review_sample_frame_fingerprint()` validates and fingerprints this audit
frame. `blinded_review_sample_from_dataset()` creates the distinct
reviewer-facing view from `source_run_manifest` and the exact donor split: only
donor/object keys, opaque stratum codes, and requested image-location columns.
The API requires the source run to match the split, requires its donor universe
to equal the split, resolves exactly one referenced `ingest` stage, and
content-validates that stage and snapshot before and after the donor-local
read. It does not accept a free-standing snapshot or caller-authored location
table and rejects prediction and sampling-design leaks. Do not hand weights,
population counts, method, or salt to reviewers.

Schema-v5 `ReviewLabelProvenance` fingerprints both the private audit frame and
the exact reviewer frame that was shown, and binds the ledger to the donor
split, assignment run, model bundle, candidate manifest/output, target, and
`review_context_artifact_content_id`. For image review, that context ID is the
content ID of the same source-run `ingest` output resolved by
`blinded_review_sample_from_dataset()`. Labels, reviewer frame, and audit must
cover exactly the same donor/object keys. An evidence-sufficient adjudicated
row requires two distinct reviewers, a resolved reference label, and a
root-cause code. Only after the ledger validates and adjudication is complete
may the audit weights and candidate predictions be rejoined for metrics.
Locked-test labels remain closed unless the release-time caller explicitly
opens them.

```python
from phenocycler.refinement_contracts import (
    ReviewLabelProvenance,
    blinded_review_sample_from_dataset,
    build_review_sampling_frame,
    validate_review_sampling_against_assignments,
    validate_review_labels,
)

review_key = review_key_path.read_bytes()
validation_sampling_audit = build_review_sampling_frame(
    validation_assignments,
    splits=split_manifest,
    split_name="validation",
    candidate_assignment_run_id=assignment_run_id,
    level="broad",
    parent=None,
    blinding_key=review_key,
)
validate_review_sampling_against_assignments(
    validation_sampling_audit,
    validation_assignments,
    splits=split_manifest,
    split_name="validation",
    candidate_assignment_run_id=assignment_run_id,
    level="broad",
    parent=None,
    blinding_key=review_key,
)
validation_reviewer_view = blinded_review_sample_from_dataset(
    validation_sampling_audit,
    source_run_manifest=source_run_manifest,
    splits=split_manifest,
)
assert validation_reviewer_view.attrs[
    "review_context_artifact_content_id"
] == ingest_manifest.output.content_sha256
del review_key

validation_review_provenance = ReviewLabelProvenance.from_labels(
    validation_labels,
    bundle_id="pancreas-validation-review-v1",
    splits=split_manifest,
    split_name="validation",
    sample_frame=validation_sampling_audit,
    reviewer_frame=validation_reviewer_view,
    assignment_run_id=assignment_run_id,
    typing_model_bundle_content_id=frozen_bundle.content_id,
    candidate_evaluation_manifest_content_id=validation_candidate.content_id,
    candidate_assignment_output_content_id=(
        validation_candidate.output_content_sha256
    ),
    level="broad",
    parent=None,
)
if not validation_review_provenance.complete_evidence_sufficient:
    raise ValueError("every sampled cell must be adjudicated with sufficient evidence")
validated_labels = validate_review_labels(
    validation_labels,
    provenance=validation_review_provenance,
    splits=split_manifest,
    sample_frame=validation_sampling_audit,
    reviewer_frame=validation_reviewer_view,
    assignment_run_id=assignment_run_id,
    typing_model_bundle_content_id=frozen_bundle.content_id,
    candidate_evaluation_manifest_content_id=validation_candidate.content_id,
    candidate_assignment_output_content_id=(
        validation_candidate.output_content_sha256
    ),
)
# Only now, after complete adjudication, join audit weights and predictions.
reviewed_validation = (
    validation_sampling_audit.merge(
        validated_labels, on=["donor_id", "object_id"], validate="1:1"
    ).merge(
        validation_assignments, on=["donor_id", "object_id"], validate="1:1"
    )
)
```

Historical labels do not acquire development authorization merely by being
copied or renamed. A deliberate review/governance decision must create the new
label namespace and provenance artifact.

`DevelopmentEvidenceProvenance` separately binds fitting to the exact marker-
evidence source artifact, source run, donor split, feature schema, row count,
and content fingerprint of the labeled projection. The source is not an
arbitrary path hash: it must be the output recorded by the source
`RunManifest`'s single `calibrate` `StageManifest`. The run name must equal the
split's `source_run_id`; run/stage identities, donor sets, output path, and the
current dataset content are verified. The recorded
`source_run_manifest_content_id`, `source_stage_manifest_content_id`, and
`source_artifact_content_id` identify that exact chain; the last is the
calibration stage output's semantic content hash. Modified, substituted, or
different-run evidence therefore cannot retain an old bundle's provenance.

After label governance approval, create these two immutable manifests in a
separate preparation step (never inside the fitting command):

```python
from phenocycler.artifacts import RunManifest
from phenocycler.hierarchical_typing import FEATURE_SCHEMA_VERSION
from phenocycler.probabilistic_training import load_labeled_evidence
from phenocycler.refinement_contracts import (
    DEVELOPMENT_LABEL_PURPOSE,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
)

label_provenance = DevelopmentLabelProvenance.from_labels(
    labels,
    bundle_id="approved-pancreas-development-v1",
    splits=split_manifest,
    purpose=DEVELOPMENT_LABEL_PURPOSE,
)
label_provenance.write_json("development-label-provenance.json")

source_run_manifest = RunManifest.read_json(SOURCE_RUN_MANIFEST_PATH)
labeled_evidence = load_labeled_evidence(EVIDENCE_PATH, labels)
evidence_provenance = DevelopmentEvidenceProvenance.from_evidence(
    labeled_evidence,
    splits=split_manifest,
    evidence_path=EVIDENCE_PATH,
    source_run_manifest=source_run_manifest,
    feature_schema_version=FEATURE_SCHEMA_VERSION,
)
evidence_provenance.write_json("development-evidence-provenance.json")
```

## Preliminary fitting, calibration selection, and freezing

The implemented release path deliberately separates the scientific fit from
the acceptance decision:

```text
preliminary/unfrozen fit
  -> complete calibration candidate assignments
  -> exact prediction-blinded reviews for broad + every fitted subtype
  -> replayable threshold-selection provenance (schema v3)
  -> freeze those thresholds onto the same fit, without refitting
  -> complete validation reviews/report for broad + all-cell specific + fitted subtypes
  -> explicitly opened locked-test reviews/report for the same required targets
  -> policy-derived promotion manifest
```

### 1. Fit one preliminary bundle

The preliminary fit needs explicit *seed* gates because it must execute the
same hierarchical assignment contract when it creates calibration candidates.
Those values are not a production selection and the resulting bundle cannot be
promoted. The stable `model_fit_content_id` excludes these gates, bundle
version, creation time, and later threshold provenance.

```python
from phenocycler.hierarchical_typing import AssignmentThresholds
from phenocycler.probabilistic_training import fit_probabilistic_typing_bundle

preliminary = fit_probabilistic_typing_bundle(
    evidence,
    development_and_calibration_labels,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    splits=split_manifest,
    label_provenance=label_provenance,
    evidence_provenance=evidence_provenance,
    bundle_version="pancreas-preliminary-v1",
    thresholds=AssignmentThresholds(
        inferred_probability=0.70,  # seed gate, not the release threshold
        inferred_margin=0.10,
        inferred_stability=0.80,
        anchor_probability=0.95,
        negative_probability=0.20,
    ),
    n_bootstrap=100,
    random_state=20260801,
)
assert preliminary.threshold_selection_provenance is None
preliminary.write_json("model-bundles/pancreas-preliminary-v1.json")
```

The equivalent CLI makes the non-promotable state explicit:

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

The CLI flags set the three inferred seed gates; it also records the
`AssignmentThresholds` defaults `anchor_probability=0.95` and
`negative_probability=0.20`. Schema-v3 threshold provenance binds all five
preliminary values. The fitting entry point creates only this unfrozen
candidate; it cannot attach threshold-selection provenance or create a final
bundle. The separate freeze operation does that without refitting.

This call performs the following for the broad model and each subtype parent
with authorized labels:

1. Build the class/marker layout from the versioned typing registry.
2. Require every class in at least two development donors, then fit a
   no-intercept additive multinomial model on development donors.
   Class mass is balanced first across classes and then across donors, so a
   common class or heavily reviewed donor cannot dominate.
3. Constrain anchor/support coefficients to be nonnegative, exclusion
   coefficients to be nonpositive, and every unregistered coefficient to zero.
   Fits that do not converge fail.
4. Fit one positive softmax temperature on held-out calibration-donor labels.
   Because calibration also gives equal total weight to every class and then
   to every donor within a class, the resulting values estimate the declared
   donor/class-balanced review distribution. They are not cohort-prevalence
   posteriors and must not be interpreted as expected class abundance. The
   fitted temperature is restricted to `[0.05, 20]`; uniform probabilities, a
   boundary optimum, or a locally flat/nonidentifiable objective fail fitting.
5. Make exactly `n_bootstrap` whole-development-donor draws with replacement;
   there are no retries or discarded draws. A class-incomplete or nonconverged
   draw is retained with its failure reason and `bootstrap_valid=false`.
   Stability counts agreement only from valid, unique-winner members but still
   divides by all requested draws, so failures lower stability. At least 80%
   of the requested draws must be valid or fitting aborts.
6. Record every donor draw, coefficient matrix, label count, fit parameter,
   validity/failure flag, temperature, and split/label/evidence fingerprint in
   the preliminary bundle.

Set `require_all_subtype_parents=True` when release policy requires a fitted
rescue model for every learnable configured subtype parent. A parent with only
one configured subtype has no estimable named-versus-alternative softmax risk;
it is recorded in `structural_rules_only_parents` and deliberately receives no
probabilistic subtype model. Its registered structural rules and authoritative
anchors remain available. For a multi-class parent with no fitted subtype
ensemble, production still permits authoritative anchors, but probabilistic
subtype inference has no stability and therefore abstains (or fitting fails
when `require_all_subtype_parents=True`).

### 2. Generate a content-validated calibration candidate

Candidate evaluation is isolated from production promotion and preserves the
complete donor/object universe of the source calibration artifact. The
manifest binds the preliminary bundle, exact source run/stage/artifact,
registry and rules fingerprints, split donors, assignment run, output schema,
row count, object IDs, and output content. Every scientific consumer invokes
`validate_current()` before using it; a dropped, added, duplicated,
substituted, or modified row fails.

```bash
python -m phenocycler.pipeline probabilistic evaluate \
  --model-bundle model-bundles/pancreas-preliminary-v1.json \
  --evidence data/runs/<source-run>/marker_evidence \
  --source-run-manifest data/runs/<source-run>/run.json \
  --output-root refinement/calibration-assignments-v1 \
  --manifest refinement/calibration-assignments-v1.json \
  --split calibration
```

Candidate inference is an explicit nonproduction API and stamps
`typing_mode=probabilistic_candidate`. The public production typing entry point
will not accept a bundle unless it is paired with an in-memory promotion
manifest whose external-key signature has already been verified; only that
path stamps `typing_mode=probabilistic_bundle`.

Build a private `build_review_sampling_frame()` audit, its reviewer-only
`blinded_review_sample_from_dataset()` projection from the exact source-run
`ingest` output, and an adjudicated ledger for the broad target and separately
for every probabilistic subtype target. Record that ingest output's content ID
as each schema-v5 provenance object's `review_context_artifact_content_id`, and
pass the exact reviewer projection to `ReviewLabelProvenance.from_labels()` so
its fingerprint is bound alongside the private audit. Pass the private audits
as `calibration_sample_frames` and their corresponding exact projections as
`calibration_reviewer_frames`; the two contracts are not interchangeable. The
mapping key is `("broad", "")` for broad and `("subtype", parent)` for a
subtype. Each modeled class needs reviewed positive and negative support.

### 3. Derive replayable threshold provenance from exact reviews

`build_threshold_selection_provenance()` accepts no caller-authored operating
table, selected row, or result metrics. It validates the complete candidate
artifact before and after reading; verifies its assignments against the exact
assignment-run-bound sampling audits, reviewer frames, and adjudicated ledgers;
derives the full probability × margin × stability grid; and deterministically
selects a point. The declared total and
rescue coverage/error gates must pass independently for broad and every fitted
subtype. Aggregate FP/FN cost is only the tie-break among threshold tuples that
are feasible for every target.

```python
from phenocycler.candidate_evaluation import CandidateEvaluationManifest
from phenocycler.threshold_selection import (
    ThresholdSelectionPolicy,
    build_threshold_selection_provenance,
)

calibration_candidate = CandidateEvaluationManifest.read_json(
    "refinement/calibration-assignments-v1.json"
)
threshold_policy = ThresholdSelectionPolicy(
    policy_version="pancreas-calibration-policy-v1",
    false_positive_cost=2.0,
    false_negative_cost=1.0,
    maximum_selective_error=0.10,
    minimum_coverage=0.50,
    maximum_rescue_selective_error=0.15,
    minimum_rescue_coverage=0.20,
)
threshold_provenance = build_threshold_selection_provenance(
    selection_version="pancreas-threshold-selection-v3",
    model_bundle=preliminary,
    calibration_candidate_manifest=calibration_candidate,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    source_run_manifest=source_run_manifest,
    calibration_sample_frames=calibration_sample_frames,
    calibration_reviewer_frames=calibration_reviewer_frames,
    calibration_review_ledgers=calibration_review_ledgers,
    calibration_review_provenance=calibration_review_provenances,
    review_blinding_key=review_key_path.read_bytes(),
    policy=threshold_policy,
    probability_thresholds=(0.70, 0.75, 0.80, 0.85, 0.90),
    margin_thresholds=(0.10, 0.20, 0.25, 0.30),
    stability_thresholds=(0.80, 0.85, 0.90, 0.95),
)
threshold_provenance.write_json(
    "model-bundles/pancreas-threshold-selection-v3.json"
)
```

Schema v3 embeds all five preliminary assignment thresholds, the normalized
reviewed scientific rows, and the complete
derived grid. Its strict loader rejects duplicate/nonfinite JSON and rederives
the grid, selected thresholds, and selected metrics, so changing those values
and recomputing only the outer content hash is insufficient. It binds both the
preliminary bundle `content_id` and the threshold-independent
`model_fit_content_id`. The minimum probability, margin, and stability values
in the calibration grid must each be at least as strict as its corresponding
preliminary seed gate. A more-permissive grid is rejected because subtype
scores only exist downstream of broad parents accepted by those seed gates.

### 4. Freeze the selected thresholds without refitting

```bash
python -m phenocycler.pipeline probabilistic freeze \
  --model-bundle model-bundles/pancreas-preliminary-v1.json \
  --threshold-provenance model-bundles/pancreas-threshold-selection-v3.json \
  --bundle-version pancreas-typing-v1 \
  --out model-bundles/pancreas-typing-v1.json
```

`TypingModelBundle.freeze_threshold_selection()` replaces only the three
selected inferred gates; the schema-v3 record proves the complete five-gate
preliminary state, while anchor and negative gates remain unchanged. Freeze
embeds the exact provenance and creates a new full bundle `content_id`. It
asserts that the `model_fit_content_id` is unchanged. There is no coefficient,
temperature, or bootstrap refit in this step.

## Data-driven FP/FN and plotted decisions

Use an unbiased probability sample to estimate performance and a separate
challenge sample to find rare failures. The challenge sample should enrich
low-margin cells, model disagreements, exclusion conflicts, large spillover
corrections, segmentation risks, rare candidates, `Other`, `Ambiguous`, and
`Unavailable` cells. Only the probability sample and its validated sampling
weights may determine threshold or held-out report metrics. The challenge
sample remains a detached root-cause audit and must not estimate prevalence or
headline error rates. If the study supplies one at release, retain its exact
validation candidate/assignment lineage plus explicit evidence-sufficiency,
root-cause, catastrophic, and disposition fields. A catastrophic finding is a
categorical stop requiring a new development cycle; it is not converted into a
challenge-sample error rate or a cell-level override.

At minimum, inspect these plotted decisions by class and donor:

| Plot | Decision it supports | FP/FN protection |
|---|---|---|
| reliability curve plus Brier/log-loss | whether temperature-scaled values behave as probabilities | detects overconfident false calls and underconfident true calls |
| precision-recall curve | probability operating region, especially for rare classes | exposes FP burden hidden by ROC curves |
| weighted confusion matrix | class-specific error direction | separates wrong accepted calls from abstentions |
| coverage versus selective error | acceptable abstention/error tradeoff | prevents accuracy gains obtained only by rejecting nearly everything |
| probability × margin × stability heatmap | frozen joint gate | shows whether one permissive gate defeats the others |
| per-donor metric/interval plot | worst-donor robustness | prevents a large donor from masking donor-specific FN or FP failures |
| bootstrap-stability distribution | ensemble agreement cutoff | identifies calls dependent on a small donor subset |
| stratified image galleries | root cause | distinguishes model errors from segmentation, spillover, or unavailable evidence |

The threshold-provenance builder derives the declared grid directly from the
content-validated assignments and reviewed reference rows. Fixed authoritative
anchors and `Other` dispositions are retained; only the structurally eligible
rescue lane is rethresholded. `Other` is a resolved reference negative and a
non-call, not a modeled named-class success. Valid positive support and a
strict unique winner (`margin > 1e-12`) remain invariant gates at every grid
point and cannot be weakened by choosing zero thresholds.

Optional sample weights support inverse-probability estimates from the
probability sample. The table reports total named-class coverage/selective
error separately from rescue-eligible weight, rescue-call weight, rescue
coverage, and rescue selective error. This prevents fixed anchor calls from
concealing a rescue lane that accepts no cells. Wrong accepted named calls
count as both a false positive for the called class and a false negative for a
named reference class; abstentions count as false negatives but not false
positives. `false_positive_cost` and `false_negative_cost` make the intended
asymmetry explicit. No desired class prevalence enters the score.

The prespecified `ThresholdSelectionPolicy` sets total and rescue-specific
coverage/error gates plus FP/FN costs. The selector requires one common
threshold tuple to satisfy every modeled target and raises rather than
weakening the criteria. Threshold choice belongs to calibration donors. After
freezing, use validation donors to estimate generalization and open locked-test
labels only for the release decision; never tune after either evaluation.

Bootstrap confidence intervals and comparisons should use donors or ROIs as
the independent unit, not millions of cells as independent replicates. Report
coverage beside conditional accuracy and inspect both predicted positives and
predicted-negative/Other cells so false negatives remain observable.
One-vs-rest reliability retains resolved `Other` references as negatives for
every modeled class and verifies that the complete softmax row sums to one; it
must not discard the very cells that expose forced-class overconfidence.

Notebooks are the plotted decision interface, not an alternative model
implementation. They should call these package APIs, read only their declared
split, cache content-addressed summaries, and write the selected operating row
and rationale to the decision ledger. Development notebooks must not open the
locked-test ledger, and changing plot styling must not refit a candidate.
Notebook 05 reads a frozen validation `CandidateEvaluationManifest` and is
intentionally **broad-only**. Each probabilistic subtype needs its own blinded
sample, reviewer-frame-bound adjudicated ledger, provenance, plots, and report
outside notebook 05. Validation and locked-test release reports also require a
separate all-cell `level="specific", parent=None` target. Therefore notebook
05's broad results cannot satisfy promotion on their own. Threshold sweeps
belong to calibration review before the final bundle is frozen.

## Standalone held-out reports and promotion

Create a new complete assignment artifact for the frozen bundle; do not reuse
the preliminary calibration output:

```bash
python -m phenocycler.pipeline probabilistic evaluate \
  --model-bundle model-bundles/pancreas-typing-v1.json \
  --evidence data/runs/<source-run>/marker_evidence \
  --source-run-manifest data/runs/<source-run>/run.json \
  --output-root refinement/validation-assignments-v1 \
  --manifest refinement/validation-assignments-v1.json \
  --split validation
```

After freezing one complete review provenance per target, derive metrics from
the assignment/sample/ledger artifacts rather than supplying metric values:

```python
from phenocycler.evaluation_reporting import (
    EvaluationTargetReview,
    build_split_evaluation_report,
)
from phenocycler.typing_model_bundle import (
    PromotionGatePolicy,
    load_typing_model_bundle,
)

frozen_bundle = load_typing_model_bundle(
    "model-bundles/pancreas-typing-v1.json",
    marker_registry=marker_registry,
    typing_registry=typing_registry,
)
release_policy = PromotionGatePolicy(
    minimum_donors=5,
    minimum_rows=200,
    minimum_coverage_lower=0.50,
    maximum_selective_error_upper=0.10,
    minimum_rescue_coverage_lower=0.20,
    maximum_rescue_selective_error_upper=0.15,
    maximum_class_false_positive_rate_upper=0.10,
    maximum_class_false_negative_rate_upper=0.20,
    minimum_class_positive_rows=20,
    minimum_class_negative_rows=20,
    minimum_class_positive_donors=3,
    minimum_class_negative_donors=3,
    confidence_level=0.95,
    # For eight required targets at 95% study-wise confidence, the finite-tail
    # floor is 12,800; preregister a modest buffer above that boundary.
    bootstrap_replicates=16000,
)
if release_policy.content_id != frozen_bundle.split_manifest.promotion_policy_content_id:
    raise ValueError("release policy differs from the donor-split precommitment")

def evaluation_reviews(
    sampling_audits, reviewer_frames, review_ledgers, provenances
):
    def key(provenance):
        parent = provenance.parent if provenance.level == "subtype" else None
        return provenance.level, parent

    return tuple(
        EvaluationTargetReview(
            sample_frame=sampling_audits[key(provenance)],
            reviewer_frame=reviewer_frames[key(provenance)],
            review_labels=review_ledgers[key(provenance)],
            provenance=provenance,
        )
        for provenance in provenances
    )

validation_report = build_split_evaluation_report(
    candidate_evaluation="refinement/validation-assignments-v1.json",
    model_bundle=frozen_bundle,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    source_run_manifest=source_run_manifest,
    target_reviews=evaluation_reviews(
        validation_sampling_audits,
        validation_reviewer_frames,
        validation_review_ledgers,
        validation_review_provenances,
    ),
    promotion_policy=release_policy,
    review_blinding_key=review_key_path.read_bytes(),
)
```

`build_split_evaluation_report()` rejects preliminary bundles and calibration
assignments. It requires exactly these distinct, complete targets: broad;
`level="specific", parent=None`; and one conditional subtype target for every
fitted probabilistic subtype parent. Construct the specific sampling audit from
the complete split assignment frame with `build_review_sampling_frame(...,
level="specific", parent=None)`, then create its source-run-bound reviewer
frame and schema-v5 provenance exactly as above. Its strata use
`specific_assignment_status` and `specific_type`, but its sampling population
is **all cells**, not only cells predicted into a particular broad parent.

The specific target reports end-to-end FP/FN rates over the frozen final
ontology leaves. It catches a true leaf sent down the wrong broad branch and
false negatives from structural/rules-only or fallback specific assignments—
errors that a conditional per-parent subtype review cannot observe. Broad and
conditional subtype targets remain mandatory because they localize different
failure mechanisms. The report computes sampling-weighted total/rescue
coverage and selective error plus per-class FP and FN rates. Intervals use donor
resampling plus [Rao-Wu-style](https://doi.org/10.1080/01621459.1988.10478591)
rescaling inside non-census review strata, finite-population correction,
target-wise simultaneous max-deviation bands, and a conservative clustered
zero-event guard, informed by conservative complex-survey proportion guidance
from [Korn and Graubard](https://publications.gc.ca/collections/collection_2016/statcan/12-001/CS12-001-24-2-eng.pdf)
and [Dean and Pagano](https://doi.org/10.1093/jssam/smv024). It is an explicit
cluster-Kish/Clopper-Pearson approximation, not a claim of an exact finite-sample
survey interval. Zero-denominator bootstrap draws remain conservative rather
than being discarded. The policy confidence `C` is allocated across the
validation and locked-test reports first; each report then allocates error
between its simultaneous band and boundary guard and across its `T` required
targets. The preregistered replicate count must therefore provide at least 20
draws in the smallest simultaneous tail: use at least `80*T/(1-C)` replicates
(12,800 for `T=8`, `C=0.95`). Preregistering 16,000 in the example leaves a
finite-tail buffer; a smaller count can deliberately yield `[0, 1]` fail-closed
intervals rather than an apparently precise release claim.

Generate locked-test assignments only after validation is frozen:

```bash
python -m phenocycler.pipeline probabilistic evaluate \
  --model-bundle model-bundles/pancreas-typing-v1.json \
  --evidence data/runs/<source-run>/marker_evidence \
  --source-run-manifest data/runs/<source-run>/run.json \
  --output-root refinement/locked-test-assignments-v1 \
  --manifest refinement/locked-test-assignments-v1.json \
  --split locked_test \
  --allow-locked-test
```

```python
locked_test_report = build_split_evaluation_report(
    candidate_evaluation="refinement/locked-test-assignments-v1.json",
    model_bundle=frozen_bundle,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    source_run_manifest=source_run_manifest,
    target_reviews=evaluation_reviews(
        locked_test_sampling_audits,
        locked_test_reviewer_frames,
        locked_test_review_ledgers,
        locked_test_review_provenances,
    ),
    promotion_policy=release_policy,
    review_blinding_key=review_key_path.read_bytes(),
    allow_locked_test=True,
)
```

The `allow_locked_test` flag is an explicit governance acknowledgement,
**not** a cryptographic or stateful one-use lock. Access control and the study
protocol must predeclare the single opening; archive that immutable
assignment/review result and do not regenerate candidates after inspecting it.

## Immutable bundle and promotion checks

`TypingModelBundle` is strict, content-addressed JSON—not a pickle. Its full
content ID covers schema/method/feature versions, score semantics,
registry/rules fingerprints, split, label and training-evidence provenance,
thresholds and their selection provenance, base and bootstrap coefficients,
donor draws and validity, temperature, and exact fitting parameters. The
separate stable `model_fit_content_id` covers the scientific fit but excludes
release gates and metadata. `write_json()` refuses to overwrite an existing
bundle.

The bundle contains no embedded promotion-report field; its strict loader
rejects such unknown metadata. Production authorization exists only in a
separate `TypingModelPromotionManifest`, created after frozen validation and
locked-test evaluation. Production promotion first requires the bundle's
schema-v2 split to contain nonblank
`release_governance.promotion_policy_content_id` and
`release_governance.signing_key_sha256`. The supplied policy must have the
precommitted content ID, and the protected external signing key must match the
precommitted digest. A blank field or either mismatch fails closed. The secret
key remains outside every JSON artifact and notebook.

The signed replay receipt hashes the complete installed `phenocycler` package,
project/environment specifications, Python executable, platform ABI, and the
actual installed NumPy, pandas, SciPy, and PyArrow files—not package version
strings alone. A deployment with a solved environment or immutable container
should also set `PHENOCYCLER_RUNTIME_IMAGE_DIGEST` to its lowercase SHA-256
digest before both promotion and production loading. Runtime drift then
invalidates the receipt instead of silently rerunning threshold-edge decisions
under different numerical kernels.

The promotion manifest embeds and binds:

- the exact candidate bundle content ID and donor-split content ID;
- structured validation and locked-test `SplitEvaluationReport` objects;
- each report's exact candidate manifest/output, assignment run, donors,
  sampling-audit and reviewer-frame fingerprints, target-level
  `ReviewLabelProvenance`, weighted metrics, source-run-bound review-context
  artifact, and two-stage survey-bootstrap confidence intervals;
- the prespecified confidence-bound `PromotionGatePolicy`; and
- the fixed decision `approved_for_production`.

Every target review must identify this candidate model, the correct split,
source run, assignment run, exact candidate manifest/output, private
sampling-audit frame, exact reviewer frame, and review-context artifact. Each
report must contain broad, the all-cell end-to-end specific target, and every
fitted conditional subtype target. The policy derives the decision; callers
cannot provide pass flags. A promotion manifest cannot be reused for another
bundle or split. Example study-specific gates (which must be prespecified, not
copied blindly) are:

```python
from datetime import datetime, timezone

from phenocycler.threshold_selection import ThresholdSelectionReplayInputs
from phenocycler.typing_model_bundle import TypingModelPromotionManifest

threshold_source_replay_inputs = ThresholdSelectionReplayInputs(
    preliminary_bundle=preliminary,
    calibration_candidate_manifest=calibration_candidate,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    source_run_manifest=source_run_manifest,
    calibration_sample_frames=calibration_sample_frames,
    calibration_reviewer_frames=calibration_reviewer_frames,
    calibration_review_ledgers=calibration_review_ledgers,
    review_blinding_key=review_key_path.read_bytes(),
)

review_key = review_key_path.read_bytes()
release_key = release_signing_key_path.read_bytes()
if len(review_key) < 16 or len(release_key) < 16:
    raise ValueError("protected review and release keys must contain at least 16 bytes")

promotion = TypingModelPromotionManifest.create(
    promotion_version="pancreas-typing-v1-release",
    model_bundle=frozen_bundle,
    policy=release_policy,
    validation_report=validation_report,
    locked_test_report=locked_test_report,
    threshold_source_replay_inputs=threshold_source_replay_inputs,
    marker_registry=marker_registry,
    typing_registry=typing_registry,
    review_blinding_key=review_key,
    release_signing_key=release_key,
    issuer="pancreas-typing-release-operator",
    issued_at_utc=datetime.now(timezone.utc).isoformat(),
)
promotion.write_json("model-bundles/pancreas-typing-v1-promotion.json")
del review_key, release_key
```

`issued_at_utc` is canonical, signed UTC issuance metadata; it is not an
expiration time or a substitute for an external trusted timestamp. Retain the
immutable artifact in the study audit store if independent proof of issuance
time is required.

The policy confidence is study-wise across both reports; each serialized report
uses the deterministic adjusted level `1 - (1 - C) / 2` and exactly the
preregistered bootstrap count. The manifest checks every confidence-bound
total/rescue and per-class FP/FN gate for broad, the mandatory all-cell
specific target, and every fitted subtype target. It also replays the exact
calibration candidate, source evidence, sampling audits, reviewer contexts,
and ledgers that produced the frozen thresholds, without refitting. Hash-only
report summaries and caller-declared approvals are not accepted.

The bundle/promotion loaders and production context resolution reject:

- an unknown schema, method, feature contract, score semantics, or bundle/model
  field;
- a changed or internally inconsistent content ID;
- marker-registry or typing-rules fingerprint drift;
- class/marker order or coefficient-shape differences;
- NaN/Inf values, sign violations, nonzero unregistered coefficients, or
  coefficients above the production cap;
- an invalid temperature, empty or wrong-sized bootstrap ensemble, a donor-
  draw/member/validity count mismatch, less than 80% valid draws, too few
  development donors, overlapping training/calibration donors, or a
  nonconverged base fit;
- unknown subtype parents;
- a changed bundle/promotion file, a promotion for another model or split, a
  different cohort content ID, or a nonexact donor partition.

The bundle and promotion file fingerprints are captured when the run context is
resolved, recorded as inputs to the cohort `type` manifest and every donor
`type` receipt, and checked again by content immediately before inference.
Their content IDs also enter the semantic run identity. Changing either
artifact creates a new run; changing one in place makes the old input stale.

Before probabilistic typing, the current run's cohort `calibrate` manifest must
be sealed and content-valid. Its marker-evidence semantic content hash must
exactly equal the calibration-output content hash recorded in the candidate's
training-evidence provenance. The donor-pipelined scheduler holds every `type`
task at this cohort barrier; it does not infer from donor-local calibration
receipts alone. Any calibration/preprocessing drift aborts typing and requires
a newly trained and evaluated candidate.

There is no fallback from a configured bundle to rules-only typing. A missing,
changed, incompatible, unpromoted, or computationally invalid bundle aborts the
type stage before a valid donor receipt can be written.

## Opt-in production configuration

Leave the default blank for rules-only typing:

```ini
[typing]
model_bundle =
promotion_manifest =
```

After a candidate passes the declared validation and promotion process, point
to both immutable JSON files. They must be configured together; setting only
one is an error. Relative paths are resolved from the configuration file:

```ini
[typing]
model_bundle = model_bundles/pancreas-typing-v1.json
promotion_manifest = model_bundles/pancreas-typing-v1-promotion.json
```

Then run the normal pipeline. There is no special inference command:

```bash
python -m phenocycler.pipeline run --config config.ini --pipelined
python -m phenocycler.pipeline status --config config.ini
```

Rules-only output is stamped with `typing_mode=rules_only` and the explicit
score semantics `uncalibrated_rule_weight_softmax_score`. The stable provenance
schema uses empty version/content values and zero total/valid bootstrap counts
when no bundle is configured. Bundle output is stamped on every row with:

- `typing_mode=probabilistic_bundle`;
- `typing_model_bundle_content_id` and `typing_model_bundle_version`;
- `typing_model_method_version` and `typing_feature_schema_version`;
- `typing_score_semantics=donor_class_balanced_temperature_scaled_softmax_probability`;
- `typing_broad_score_semantics` plus parent-specific
  `typing_subtype_score_semantics` and `typing_subtype_model_source` (partial
  subtype bundles explicitly stamp an uncalibrated rules fallback);
- `typing_broad_bootstrap_replicates` and
  `typing_broad_valid_bootstrap_replicates`, plus the applicable parent-
  specific total and valid subtype bootstrap counts.

The ordinary registry/rules fingerprints remain present in both modes.

## Production inference decision order

Final assignment is deterministic and uncertainty preserving. Geometry-QC
ineligibility has final precedence and produces `Unavailable`. For each
otherwise eligible broad or applicable subtype decision, the order is:

1. No usable identity evidence: `Unavailable`.
2. Exactly one valid authoritative anchor with no exclusion conflict: accept
   the anchor. Anchor calls do not require a fitted probabilistic rescue model.
3. Multiple anchors, or an anchor/exclusion conflict: `Ambiguous`.
4. All required models valid and negative: `Other` at that level.
5. Otherwise, accept the model's best class as `inferred` only when all five
   gates pass:
   - valid positive anchor/support evidence for that winning class;
   - a unique top logit (ties abstain even when numeric thresholds are zero);
   - frozen temperature-scaled probability threshold;
   - frozen top-two margin threshold;
   - frozen whole-donor ensemble-stability threshold.
6. If required evidence is missing and there is no valid positive evidence:
   `Unavailable`; otherwise: `Ambiguous`.

Subtype evaluation occurs only after a broad parent is accepted by anchor or
inference. A terminal parent remains its broad label. A parent with configured
subtypes but no accepted subtype becomes `<Parent>-unclassified`; it is never
reassigned across broad compartments.

The output retains explicit `*_probability_pass`, `*_margin_pass`,
`*_stability_pass`, `*_valid_evidence_pass`, `*_unique_winner_pass`, and
`*_inference_failure_reasons` fields. These fields determine why an ambiguous
cell abstained; they are not permission to edit that cell manually. Correct the
responsible rule, evidence stage, label set, or model for a class of cells,
freeze a new artifact, and rerun the full immutable type stage.

## Why no real bundle is generated yet

The software can fit, serialize, validate, and consume a model, but it cannot
manufacture biological ground truth. The historical assignments and review
bundles were created for audit under a different policy and cannot estimate an
independent false-negative rate, authorize coefficient fitting, or serve as a
locked release test.

Producing a real bundle therefore requires, in order:

1. approve and fingerprint a donor split—including review, calibration-policy/
   grid, promotion-policy, and external-key commitments—before inspecting
   held-out labels;
2. collect new blinded probability samples with exact reviewer frames and
   source-run review-context lineage, plus any SOP-prespecified detached
   challenge samples used for root-cause discovery;
3. adjudicate image-linked labels with evidence-sufficiency and root-cause
   fields;
4. authorize a separate development-label artifact with the required purpose;
5. fit a preliminary bundle using only development/calibration labels;
6. create and content-validate complete calibration assignments;
7. derive schema-v3 threshold provenance from the exact broad and subtype
   calibration samples, ledgers, and candidate output;
8. freeze those thresholds onto the same `model_fit_content_id` without
   refitting;
9. create standalone frozen validation assignments and the structured report
   containing broad, all-cell specific, and every fitted subtype target;
10. explicitly open the locked test, create the same complete target report,
    and do not tune or regenerate a candidate after inspecting it;
11. create the separate promotion manifest whose policy derives approval from
    both exact reports; and
12. configure both accepted immutable artifacts and create a new production
    run whose calibration content exactly matches training.

Until those steps are complete, the scientifically correct behavior is the
current rules-only mode with explicit abstention—not a placeholder bundle and
not a forced final label for every ambiguous cell.
