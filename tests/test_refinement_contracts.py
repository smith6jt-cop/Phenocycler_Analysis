from __future__ import annotations

import hashlib
from dataclasses import replace

import pandas as pd
import pytest

from phenocycler.artifacts import (
    CodeMetadata,
    InputArtifact,
    RunManifest,
    StageManifest,
)
from phenocycler.refinement_contracts import (
    DEVELOPMENT_EVIDENCE_PURPOSE,
    DEVELOPMENT_LABEL_PURPOSE,
    DEVELOPMENT_LABEL_SOURCE,
    REVIEW_SAMPLING_METHOD,
    DevelopmentEvidenceProvenance,
    DevelopmentLabelProvenance,
    DonorSplitManifest,
    ReviewLabelProvenance,
    blinded_review_sample_frame,
    blinded_review_sample_from_dataset,
    build_review_sampling_frame,
    evidence_frame_fingerprint,
    label_frame_fingerprint,
    review_sample_frame_fingerprint,
    validate_development_labels,
    validate_evidence_source,
    validate_review_labels,
    validate_review_sampling_against_assignments,
    validate_training_evidence,
)

REVIEW_BLINDING_KEY = "refinement-contract-tests-key-v1"


def _splits() -> DonorSplitManifest:
    return DonorSplitManifest(
        split_version="donor-split-v1",
        source_run_id="source-run",
        cohort_content_id="cohort-id",
        development=("D2", "D1"),
        calibration=("D3",),
        validation=("D4",),
        locked_test=("D5",),
        review_n_per_stratum={
            "calibration": 1,
            "validation": 1,
            "locked_test": 1,
        },
        review_blinding_key_sha256=hashlib.sha256(
            REVIEW_BLINDING_KEY.encode()
        ).hexdigest(),
    )


def _labels() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "donor_id": ["D1", "D2", "D3", "D3"],
            "object_id": ["a", "b", "c", "c"],
            "level": ["broad", "broad", "broad", "subtype"],
            "parent": ["", "", "", "Immune"],
            "label": ["Immune", "Epithelial", "Immune", "T_cell"],
            "evidence_sufficient": [True, True, True, True],
            "audit_only": [False, False, False, False],
            "label_purpose": [DEVELOPMENT_LABEL_PURPOSE] * 4,
            "label_source": [DEVELOPMENT_LABEL_SOURCE] * 4,
            "reviewer_1": ["reviewer-A"] * 4,
            "reviewer_2": ["reviewer-B"] * 4,
            "adjudicated": [True] * 4,
            "root_cause": ["independent_marker_review"] * 4,
        }
    )


def test_split_identity_and_label_fingerprint_are_order_invariant(tmp_path):
    first = _splits()
    second = DonorSplitManifest.from_dict(first.to_dict())
    assert first.content_id == second.content_id
    first.validate_cohort(["D5", "D4", "D3", "D2", "D1"])

    labels = _labels()
    assert label_frame_fingerprint(labels) == label_frame_fingerprint(
        labels.sample(frac=1, random_state=7)
    )

    path = tmp_path / "splits.json"
    first.write_json(path)
    assert DonorSplitManifest.read_json(path) == first
    with pytest.raises(FileExistsError, match="immutable"):
        first.write_json(path)


def test_split_rejects_overlap_and_nonexact_cohort():
    with pytest.raises(ValueError, match="overlap"):
        DonorSplitManifest(
            split_version="v",
            source_run_id="run",
            cohort_content_id="cohort",
            development=("D1", "D2"),
            calibration=("D2",),
        )
    with pytest.raises(ValueError, match="exact cohort"):
        _splits().validate_cohort(["D1", "D2", "D3", "D4", "D6"])


def test_only_separate_development_labels_can_enter_fitting():
    labels = _labels()
    provenance = DevelopmentLabelProvenance.from_labels(
        labels,
        bundle_id="new-development-review-v1",
        splits=_splits(),
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    observed = validate_development_labels(
        labels,
        provenance=provenance,
        splits=_splits(),
    )
    assert set(observed["donor_id"]) == {"D1", "D2", "D3"}

    with pytest.raises(ValueError, match="explicitly authorized"):
        DevelopmentLabelProvenance.from_labels(
            labels,
            bundle_id="gold_standard_v2",
            splits=_splits(),
            purpose="audit_only",
        )
    with pytest.raises(KeyError, match="label_purpose"):
        label_frame_fingerprint(labels.drop(columns="label_purpose"))

    leaked = pd.concat(
        [
            labels,
            pd.DataFrame(
                {
                    "donor_id": ["D5"],
                    "object_id": ["locked"],
                    "level": ["broad"],
                    "parent": [""],
                    "label": ["Immune"],
                    "evidence_sufficient": [True],
                    "audit_only": [False],
                    "label_purpose": [DEVELOPMENT_LABEL_PURPOSE],
                    "label_source": [DEVELOPMENT_LABEL_SOURCE],
                    "reviewer_1": ["reviewer-A"],
                    "reviewer_2": ["reviewer-B"],
                    "adjudicated": [True],
                    "root_cause": ["independent_marker_review"],
                }
            ),
        ],
        ignore_index=True,
    )
    leaked_provenance = DevelopmentLabelProvenance.from_labels(
        leaked,
        bundle_id="new-development-review-v2",
        splits=_splits(),
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    with pytest.raises(ValueError, match="locked-test"):
        validate_development_labels(
            leaked,
            provenance=leaked_provenance,
            splits=_splits(),
        )


def test_label_fingerprint_detects_changes_and_audit_only_rows_fail():
    labels = _labels()
    provenance = DevelopmentLabelProvenance.from_labels(
        labels,
        bundle_id="dev",
        splits=_splits(),
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    changed = labels.copy()
    changed.loc[0, "label"] = "Epithelial"
    with pytest.raises(ValueError, match="fingerprint"):
        validate_development_labels(
            changed,
            provenance=provenance,
            splits=_splits(),
        )

    audit = labels.assign(audit_only=True)
    with pytest.raises(ValueError, match="audit-only"):
        label_frame_fingerprint(audit)

    inconsistent = labels.copy()
    inconsistent.loc[
        inconsistent["level"].eq("broad") & inconsistent["object_id"].eq("c"), "label"
    ] = "Epithelial"
    with pytest.raises(ValueError, match="matching broad parent"):
        label_frame_fingerprint(inconsistent)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("label_source", "model_prediction", "label_source"),
        ("reviewer_1", "", "non-empty reviewers"),
        ("reviewer_2", "reviewer-A", "distinct reviewers"),
        ("adjudicated", False, "must be adjudicated"),
        ("root_cause", "", "root_cause"),
    ],
)
def test_development_labels_require_independent_blinded_adjudication(column, value, message):
    labels = _labels()
    labels.loc[0, column] = value
    with pytest.raises(ValueError, match=message):
        label_frame_fingerprint(labels)

    with pytest.raises(KeyError, match="reviewer_1"):
        label_frame_fingerprint(_labels().drop(columns="reviewer_1"))


def test_label_and_evidence_provenance_are_bound_to_the_exact_split(tmp_path):
    splits = _splits()
    labels = _labels()
    label_provenance = DevelopmentLabelProvenance.from_labels(
        labels,
        bundle_id="dev-v1",
        splits=splits,
        purpose=DEVELOPMENT_LABEL_PURPOSE,
    )
    label_path = tmp_path / "labels.json"
    label_provenance.write_json(label_path)
    assert DevelopmentLabelProvenance.read_json(label_path) == label_provenance

    rearranged = DonorSplitManifest(
        split_version="donor-split-v2",
        source_run_id="source-run",
        cohort_content_id="cohort-id",
        development=("D1", "D3"),
        calibration=("D2",),
        validation=("D4",),
        locked_test=("D5",),
    )
    with pytest.raises(ValueError, match="different donor split"):
        validate_development_labels(
            labels,
            provenance=label_provenance,
            splits=rearranged,
        )


def test_evidence_provenance_resolves_the_exact_calibration_stage(tmp_path):
    splits = _splits()
    output = tmp_path / "marker_evidence"
    for donor in splits.donors:
        partition = output / f"donor_id={donor}"
        partition.mkdir(parents=True)
        pd.DataFrame({"object_id": [f"{donor}-cell"], "A__state": ["positive"]}).to_parquet(
            partition / "part.parquet", index=False
        )
    source = tmp_path / "producer.py"
    source.write_text("METHOD = 'calibrate-v1'\n")
    raw = tmp_path / "raw.txt"
    raw.write_text("immutable input\n")
    code = CodeMetadata.capture(tmp_path, source_paths=(source,))
    stage = StageManifest.create(
        stage="calibrate",
        method_version="calibrate-v1",
        expected_donors=splits.donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=output,
        config={"contract": "test"},
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "calibrate.json")
    run = RunManifest.create(
        run_name=splits.source_run_id,
        expected_donors=splits.donors,
        stage_manifest_paths=(stage_path,),
        config={"contract": "test"},
        code=code,
    )
    evidence = pd.DataFrame(
        {
            "donor_id": ["D1", "D2", "D3"],
            "object_id": ["D1-cell", "D2-cell", "D3-cell"],
            "A__state": ["positive", "positive", "positive"],
        }
    )
    provenance = DevelopmentEvidenceProvenance.from_evidence(
        evidence,
        splits=splits,
        evidence_path=output,
        source_run_manifest=run,
        feature_schema_version="features-v1",
    )
    assert provenance.source_run_manifest_content_id == run.content_id
    assert provenance.source_stage_manifest_content_id == stage.content_id
    assert provenance.source_artifact_content_id == stage.output.content_sha256
    validate_evidence_source(
        output,
        provenance=provenance,
        splits=splits,
        source_run_manifest=run,
    )

    copied = tmp_path / "copied_evidence"
    copied.mkdir()
    with pytest.raises(ValueError, match="not the calibration output"):
        validate_evidence_source(
            copied,
            provenance=provenance,
            splits=splits,
            source_run_manifest=run,
        )

    evidence = pd.DataFrame(
        {
            "donor_id": ["D1", "D2", "D3"],
            "object_id": ["a", "b", "c"],
            "A__state": ["positive", "negative", "positive"],
        }
    )
    evidence_provenance = DevelopmentEvidenceProvenance(
        source_run_id=splits.source_run_id,
        split_manifest_content_id=splits.content_id,
        source_run_manifest_content_id="r" * 64,
        source_stage_manifest_content_id="s" * 64,
        source_artifact_content_id="e" * 64,
        labeled_frame_fingerprint=evidence_frame_fingerprint(evidence),
        feature_schema_version="features-v1",
        row_count=len(evidence),
        purpose=DEVELOPMENT_EVIDENCE_PURPOSE,
    )
    evidence_path = tmp_path / "evidence.json"
    evidence_provenance.write_json(evidence_path)
    assert DevelopmentEvidenceProvenance.read_json(evidence_path) == evidence_provenance
    changed = evidence.copy()
    changed.loc[0, "A__state"] = "negative"
    with pytest.raises(ValueError, match="content fingerprint"):
        validate_training_evidence(
            changed,
            provenance=evidence_provenance,
            splits=splits,
            feature_schema_version="features-v1",
        )


def _review_labels(donor="D4"):
    return pd.DataFrame(
        {
            "donor_id": [donor],
            "object_id": ["reviewed"],
            "reference_label": ["Immune"],
            "reviewer_1": ["reviewer-A"],
            "reviewer_2": ["reviewer-B"],
            "adjudicated": [True],
            "root_cause": ["broad_model"],
            "evidence_sufficient": [True],
        }
    )


def _review_sample(donor="D4"):
    return pd.DataFrame(
        {
            "donor_id": [donor],
            "object_id": ["reviewed"],
            "image": ["image-1"],
            "X_centroid": [10.0],
            "Y_centroid": [20.0],
            "cell_region": ["tissue"],
            "sampling_stratum": ["S0001"],
            "sampling_method": [REVIEW_SAMPLING_METHOD],
            "sampling_salt": ["review-sample-v1"],
            "stratum_population": [10],
            "stratum_sampled": [1],
            "sampling_weight": [10.0],
            "selection_probability": [0.1],
            "review_bundle": ["probability_sample"],
        }
    )


def _reviewer_frame(sample, *, context_content_id="e" * 64):
    frame = blinded_review_sample_frame(sample)
    frame.attrs["review_context_artifact_content_id"] = context_content_id
    return frame


def _review_manifest(
    labels,
    sample,
    *,
    splits,
    split_name,
    bundle_id,
    level="broad",
    parent=None,
):
    return ReviewLabelProvenance.from_labels(
        labels,
        bundle_id=bundle_id,
        splits=splits,
        split_name=split_name,
        sample_frame=sample,
        reviewer_frame=_reviewer_frame(sample),
        assignment_run_id="assignment-run",
        typing_model_bundle_content_id="b" * 64,
        candidate_evaluation_manifest_content_id="c" * 64,
        candidate_assignment_output_content_id="d" * 64,
        level=level,
        parent=parent,
    )


def test_review_manifest_is_split_restricted_and_locked_test_is_sealed(tmp_path):
    splits = _splits()
    labels = _review_labels()
    sample = _review_sample()
    provenance = _review_manifest(
        labels,
        sample,
        splits=splits,
        split_name="validation",
        bundle_id="validation-review-v1",
    )
    observed = validate_review_labels(
        labels,
        provenance=provenance,
        splits=splits,
        sample_frame=sample,
        reviewer_frame=_reviewer_frame(sample),
        assignment_run_id="assignment-run",
        typing_model_bundle_content_id="b" * 64,
        candidate_evaluation_manifest_content_id="c" * 64,
        candidate_assignment_output_content_id="d" * 64,
    )
    assert observed.loc[0, "reference_label"] == "Immune"
    assert provenance.complete_evidence_sufficient is True
    assert provenance.level == "broad"
    assert provenance.parent == ""
    path = tmp_path / "review-manifest.json"
    provenance.write_json(path)
    assert ReviewLabelProvenance.read_json(path) == provenance

    escaped = _review_labels("D3")
    escaped_sample = _review_sample("D3")
    with pytest.raises(ValueError, match="escape"):
        _review_manifest(
            escaped,
            escaped_sample,
            splits=splits,
            split_name="validation",
            bundle_id="bad-validation",
        )

    locked = _review_labels("D5")
    locked_sample = _review_sample("D5")
    locked_provenance = _review_manifest(
        locked,
        locked_sample,
        splits=splits,
        split_name="locked_test",
        bundle_id="locked-v1",
    )
    with pytest.raises(ValueError, match="cannot be opened"):
        validate_review_labels(
            locked,
            provenance=locked_provenance,
            splits=splits,
            sample_frame=locked_sample,
            reviewer_frame=_reviewer_frame(locked_sample),
            assignment_run_id="assignment-run",
            typing_model_bundle_content_id="b" * 64,
            candidate_evaluation_manifest_content_id="c" * 64,
            candidate_assignment_output_content_id="d" * 64,
        )


def test_review_sample_is_blinded_and_bound_to_the_manifest():
    splits = _splits()
    labels = _review_labels()
    sample = _review_sample()
    provenance = _review_manifest(
        labels, sample, splits=splits, split_name="validation", bundle_id="review"
    )
    assert provenance.sample_frame_fingerprint == review_sample_frame_fingerprint(sample)
    assert provenance.review_context_artifact_content_id == "e" * 64

    with pytest.raises(ValueError, match="must carry.*review_context"):
        ReviewLabelProvenance.from_labels(
            labels,
            bundle_id="unstamped-review",
            splits=splits,
            split_name="validation",
            sample_frame=sample,
            reviewer_frame=blinded_review_sample_frame(sample),
            assignment_run_id="assignment-run",
            typing_model_bundle_content_id="b" * 64,
            candidate_evaluation_manifest_content_id="c" * 64,
            candidate_assignment_output_content_id="d" * 64,
        )

    exposed = sample.assign(broad_best_probability=0.9)
    with pytest.raises(ValueError, match="not prediction-blinded"):
        review_sample_frame_fingerprint(exposed)

    changed = sample.assign(X_centroid=11.0)
    with pytest.raises(ValueError, match="blinded content fingerprint"):
        validate_review_labels(
            labels,
            provenance=provenance,
            splits=splits,
            sample_frame=changed,
            reviewer_frame=_reviewer_frame(sample),
            assignment_run_id="assignment-run",
            typing_model_bundle_content_id="b" * 64,
            candidate_evaluation_manifest_content_id="c" * 64,
            candidate_assignment_output_content_id="d" * 64,
        )


def test_review_sampling_is_replayed_against_the_complete_candidate_population():
    assignments = pd.DataFrame(
        {
            "donor_id": ["D4", "D4", "D4"],
            "object_id": ["a", "b", "c"],
            "broad_assignment_status": ["inferred", "inferred", "other"],
            "broad_label": ["Immune", "Immune", "Other"],
        }
    )
    sample = build_review_sampling_frame(
        assignments,
        splits=_splits(),
        split_name="validation",
        candidate_assignment_run_id="c" * 64,
        level="broad",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )
    validate_review_sampling_against_assignments(
        sample,
        assignments,
        splits=_splits(),
        split_name="validation",
        candidate_assignment_run_id="c" * 64,
        level="broad",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )

    reviewer_frame = blinded_review_sample_frame(sample)
    assert list(reviewer_frame) == ["donor_id", "object_id", "sampling_stratum"]
    audit_only = set(sample) - {"donor_id", "object_id", "sampling_stratum"}
    assert audit_only.isdisjoint(reviewer_frame.columns)

    reweighted = sample.copy()
    stratum = reweighted.loc[0, "sampling_stratum"]
    selected = reweighted["sampling_stratum"].eq(stratum)
    reweighted.loc[selected, "stratum_population"] += 100
    sampled = reweighted.loc[selected, "stratum_sampled"].astype(float)
    population = reweighted.loc[selected, "stratum_population"].astype(float)
    reweighted.loc[selected, "sampling_weight"] = population / sampled
    reweighted.loc[selected, "selection_probability"] = sampled / population
    # The sample-only arithmetic is internally consistent, but it is false
    # relative to the full candidate assignment population.
    review_sample_frame_fingerprint(reweighted)
    with pytest.raises(ValueError, match="populations differ"):
        validate_review_sampling_against_assignments(
            reweighted,
            assignments,
            splits=_splits(),
            split_name="validation",
            candidate_assignment_run_id="c" * 64,
            level="broad",
            parent=None,
            blinding_key=REVIEW_BLINDING_KEY,
        )


def test_reviewer_locations_are_read_from_a_content_validated_dataset(tmp_path):
    assignments = pd.DataFrame(
        {
            "donor_id": ["D4", "D4"],
            "object_id": ["a", "b"],
            "broad_assignment_status": ["inferred", "other"],
            "broad_label": ["Immune", "Other"],
        }
    )
    sample = build_review_sampling_frame(
        assignments,
        splits=_splits(),
        split_name="validation",
        candidate_assignment_run_id="c" * 64,
        level="broad",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )
    root = tmp_path / "cells"
    for donor_id in _splits().donors:
        partition = root / f"donor_id={donor_id}"
        partition.mkdir(parents=True)
        object_ids = ["a", "b"] if donor_id == "D4" else [f"{donor_id}-cell"]
        pd.DataFrame(
            {
                "donor_id": [donor_id] * len(object_ids),
                "object_id": object_ids,
                "image": [f"image-{donor_id}"] * len(object_ids),
                "X_centroid": list(range(10, 10 + len(object_ids))),
                "Y_centroid": list(range(30, 30 + len(object_ids))),
                "cell_region": ["tissue"] * len(object_ids),
            }
        ).to_parquet(partition / "cells.parquet", index=False)
    path = root / "donor_id=D4" / "cells.parquet"
    source = tmp_path / "ingest_producer.py"
    source.write_text("METHOD = 'ingest-review-context-v1'\n")
    raw = tmp_path / "raw-review-context.txt"
    raw.write_text("immutable input\n")
    code = CodeMetadata.capture(tmp_path, source_paths=(source,))
    stage = StageManifest.create(
        stage="ingest",
        method_version="ingest-review-context-v1",
        expected_donors=_splits().donors,
        inputs=(InputArtifact.from_path("raw", raw),),
        output_root=root,
        config={"contract": "review-context-test"},
        code=code,
    )
    stage_path = stage.write_json(tmp_path / "ingest-review-context.json")
    run = RunManifest.create(
        run_name=_splits().source_run_id,
        expected_donors=_splits().donors,
        stage_manifest_paths=(stage_path,),
        config={"contract": "review-context-test"},
        code=code,
    )
    reviewer = blinded_review_sample_from_dataset(
        sample, source_run_manifest=run, splits=_splits()
    )
    assert set(reviewer) == {
        "donor_id",
        "object_id",
        "sampling_stratum",
        "image",
        "X_centroid",
        "Y_centroid",
        "cell_region",
    }
    assert reviewer.attrs["review_context_artifact_content_id"] == stage.output.content_sha256

    changed = pd.read_parquet(path)
    changed.loc[0, "X_centroid"] = 999.0
    changed.to_parquet(path, index=False)
    with pytest.raises(Exception, match="changed|content"):
        blinded_review_sample_from_dataset(
            sample, source_run_manifest=run, splits=_splits()
        )


def test_specific_review_sampling_covers_all_cells_not_only_predicted_parent():
    assignments = pd.DataFrame(
        {
            "donor_id": ["D4", "D4", "D4"],
            "object_id": ["correct-beta", "misrouted-beta", "abstained-beta"],
            "specific_assignment_status": ["anchor", "anchor", "ambiguous"],
            "specific_type": ["Beta", "Immune-T", "Epithelial-unclassified"],
        }
    )
    sample = build_review_sampling_frame(
        assignments,
        splits=_splits(),
        split_name="validation",
        candidate_assignment_run_id="c" * 64,
        level="specific",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )
    assert set(sample["object_id"]) == set(assignments["object_id"])
    validate_review_sampling_against_assignments(
        sample,
        assignments,
        splits=_splits(),
        split_name="validation",
        candidate_assignment_run_id="c" * 64,
        level="specific",
        parent=None,
        blinding_key=REVIEW_BLINDING_KEY,
    )


def test_review_manifest_completeness_is_derived_and_revalidated():
    splits = _splits()
    labels = _review_labels()
    labels.loc[0, "reference_label"] = "Ambiguous"
    labels.loc[0, "evidence_sufficient"] = False
    sample = _review_sample()
    provenance = _review_manifest(
        labels,
        sample,
        splits=splits,
        split_name="validation",
        bundle_id="incomplete-review",
    )
    assert provenance.complete_evidence_sufficient is False
    validate_review_labels(
        labels,
        provenance=provenance,
        splits=splits,
        sample_frame=sample,
        reviewer_frame=_reviewer_frame(sample),
        assignment_run_id="assignment-run",
        typing_model_bundle_content_id="b" * 64,
        candidate_evaluation_manifest_content_id="c" * 64,
        candidate_assignment_output_content_id="d" * 64,
    )

    forged = replace(
        provenance,
        complete_evidence_sufficient=True,
        content_id="",
    )
    with pytest.raises(ValueError, match="completeness/evidence-sufficiency"):
        validate_review_labels(
            labels,
            provenance=forged,
            splits=splits,
            sample_frame=sample,
            reviewer_frame=_reviewer_frame(sample),
            assignment_run_id="assignment-run",
            typing_model_bundle_content_id="b" * 64,
            candidate_evaluation_manifest_content_id="c" * 64,
            candidate_assignment_output_content_id="d" * 64,
        )

    missing_field = provenance.to_dict()
    del missing_field["complete_evidence_sufficient"]
    with pytest.raises(ValueError, match="missing required fields"):
        ReviewLabelProvenance.from_dict(missing_field)

    wrong_type = provenance.to_dict()
    wrong_type["complete_evidence_sufficient"] = 1
    with pytest.raises(TypeError, match="must be a boolean"):
        ReviewLabelProvenance.from_dict(wrong_type)


def test_review_manifest_enforces_subtype_parent_contract():
    splits = _splits()
    labels = _review_labels()
    sample = _review_sample()
    with pytest.raises(ValueError, match="subtype reviews require parent"):
        _review_manifest(
            labels,
            sample,
            splits=splits,
            split_name="validation",
            bundle_id="missing-parent",
            level="subtype",
        )
    with pytest.raises(ValueError, match="broad/specific reviews require blank parent"):
        _review_manifest(
            labels,
            sample,
            splits=splits,
            split_name="validation",
            bundle_id="unexpected-parent",
            parent="Immune",
        )

    subtype = _review_manifest(
        labels,
        sample,
        splits=splits,
        split_name="validation",
        bundle_id="subtype-review",
        level="subtype",
        parent="Immune",
    )
    assert subtype.level == "subtype"
    assert subtype.parent == "Immune"
