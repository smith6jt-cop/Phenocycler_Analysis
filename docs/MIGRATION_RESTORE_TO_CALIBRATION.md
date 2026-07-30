# What the RESTORE → donor-calibration migration retired, and what is unresolved

This branch replaces the RESTORE + ordered-residual lineage path with donor-local marker
calibration and hierarchical typing. Merging `main` into it deletes several modules that
carried deliberate scientific decisions, and those decisions do not all have an equivalent
on this side yet.

This file exists so that none of them is lost by omission. It records what was removed and
what the replacement does instead — it does **not** claim the questions are settled.

## Removed, with a direct replacement

| removed | replaced by |
|---|---|
| `phenocycler/restore.py` | `marker_calibration.py` + `calibration_stage.py` (donor-local, control-anchored) |
| `phenocycler/lineage.py` | `hierarchical_typing.py` + `typing_stage.py`, driven by `typing_rules.json` |
| `phenocycler/figures.py` | no replacement in-tree (see "Unresolved" below) |
| `phenocycler/reassess_diag.py` | no replacement in-tree (see "Unresolved" below) |
| `scripts/slurm/0{1..4}_*.sh` | one manifest-driven `scripts/slurm/run_full_pipeline.sh` |
| `redsea.donor_image` parquet guard | `artifacts._validate_cohort_images` (one image per donor, enforced at manifest construction — strictly stronger) |

## Removed with NO equivalent — open questions

These came from the `immune-floor-7class` line of work (PR #4) and its tests
(`tests/test_lineage.py`, `tests/test_hormone_floor.py`, both deleted here).

**1. The `_norm` separation floors.** `hormone_floor` rewrote `{INS,GCG,SST}_pos := _norm >= 5`
and `{CD3e,CD20,CD163}_pos := _norm >= 2`, with MPO gated at the immune K inside `lineage.py`
rather than in the stage. `_norm` was threshold-relative, so K was a *signal-separation*
margin in multiples of the per-image RESTORE threshold.

This branch's `stable_positive` is a *sampling-uncertainty* margin on where the threshold
sits — a different quantity, and a much narrower one. Whether it removes the need for a
separation floor is plausible but unverified in-tree, and `reassess_diag.py`, the only
instrument that could adjudicate it, is deleted.

The concrete evidence the floors were tuned against: donor 6476's CD3e-positive fraction
(40.3% → 2.65% under K=2, judged an over-call) versus donor 6579's (12.2% → 10.6%, judged a
real pancreatitis T-cell mode and deliberately preserved).

To settle it, either produce per-donor CD3e+/INS+ fractions from the new calibration showing
the 6476-class over-call is absent without a floor, or add a `min_log2_threshold_ratio`
predicate to the `stable_positive` conjunction — `marker_calibration.py` already computes
`log2_threshold_ratio` per cell, so the hook exists.

**2. Donor 6579 is a direct contradiction, not a merge conflict.** `phenocycler/cohort.py`
hard-excludes it as a pancreatitis outlier and `ensure_eligible_donors` raises if it appears.
The K=2 immune floor was validated *on* 6579 — it is the counterexample that proved the floor
did not simply delete immune signal. Excluding the donor and removing the floor at the same
time removes both the guard and its evidence base.

**3. The taxonomy was re-cut, implicitly.** `typing_rules.json` dissolves `Endocrine` (INS and
GCG become Epithelial anchors; Beta/Alpha/Delta become Epithelial subtypes) and collapses
`Fibroblast` + `Muscle` into `Mesenchymal`, dropping the argmax-over-three structural
tiebreak. Seven broad classes become five. This is defensible as a simplification, but it
arrives as a JSON file rather than as a stated decision, and it changes the meaning of every
downstream composition statistic.

**4. The immune floor cannot be reinstated on its original marker set.** CD20 and CD163 are
`calibration_status: "deferred"` in `marker_registry.json`, with B lineage anchoring on CD79a
instead — so two of the three markers the K=2 floor was validated against are not calibrated
here.

**5. The CD99-bright endocrine gate is retired.** `cd99_bright = 3.0` (validated cohort-wide
at 36% islet-core, 11× enrichment) has no home: CD99 is `kind: "state"`, `deferred`. The
config field survives only so `integration/export_pheno.py` imports.

## Unresolved: no diagnostic output

This branch ships no plotting or acceptance-check module for the core pipeline. The
diagnostic gap is the more serious half of that: without an equivalent of `reassess_diag`'s
objective (endocrine keratin/Vimentin load, impossible co-expression) and its guardrails
(acinar keratin retention, vessel Vimentin, ND→T1D β-loss and T-gain), question 1 above
cannot be answered empirically.
