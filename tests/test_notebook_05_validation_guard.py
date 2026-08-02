import json
from pathlib import Path

NOTEBOOK = Path(__file__).resolve().parents[1] / "notebooks" / "05_error_aware_refinement.ipynb"


def _code_cells() -> list[str]:
    payload = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return [
        "".join(cell.get("source", ()))
        for cell in payload["cells"]
        if cell.get("cell_type") == "code"
    ]


def _all_cell_source() -> str:
    payload = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    return "\n".join("".join(cell.get("source", ())) for cell in payload["cells"])


def test_notebook_05_never_scans_assignments_outside_validation_donors():
    cells = _code_cells()
    source = "\n".join(cells)
    path_guard = (
        "if any(path is None for path in (CONFIG_PATH, MODEL_BUNDLE_PATH, "
        "SOURCE_RUN_MANIFEST_PATH, CANDIDATE_EVALUATION_MANIFEST_PATH)):"
    )
    assert path_guard in source
    assert "candidate_evaluation = load_candidate_evaluation_manifest(" in source
    assert "allow_locked_test=False" in source
    assert "if candidate_evaluation.split_name != 'validation':" in source
    assert "assignments_dir = Path(candidate_evaluation.output.root)" in source
    assert "review_donors = tuple(candidate_evaluation.donors)" in source
    assert "review_donors = tuple(model_bundle.split_manifest.validation)" not in source
    assert "review_donors = tuple(context.donors)" not in source
    assert "RunContext" not in source
    assert "resolve_run_context" not in source
    assert "context." not in source
    assert "PROMOTION_MANIFEST_PATH" not in source
    assert "model_bundle.threshold_selection_provenance is not None" in source
    assert source.index(path_guard) < source.index(
        "candidate_evaluation = load_candidate_evaluation_manifest("
    )
    assert source.index("review_donors = tuple(candidate_evaluation.donors)") < source.index(
        "assignment_glob ="
    )

    assignment_query_cells = [
        cell for cell in cells if "read_parquet(?) a" in cell
    ]
    assert assignment_query_cells
    for cell in assignment_query_cells:
        assert cell.index("connection.register('allowed_review_donors'") < cell.index(
            "read_parquet(?) a"
        )
        assert cell.count("read_parquet(?) a") == cell.count(
            "JOIN allowed_review_donors d"
        )


def test_notebook_05_uses_the_exact_source_ingest_snapshot_for_image_locations():
    source = "\n".join(_code_cells())
    assert "source_run_manifest = RunManifest.read_json(SOURCE_RUN_MANIFEST_PATH)" in source
    assert "reference.stage == 'ingest'" in source
    assert "ingest_manifest = ingest_references[0].validate_current(" in source
    assert "validate_stage=True, validation_mode='content'" in source
    assert "ingest_manifest.expected_donors != source_run_manifest.expected_donors" in source
    assert "blinded_review_sample_from_dataset(" in source
    assert "source_run_manifest=source_run_manifest" in source
    assert "splits=model_bundle.split_manifest" in source
    assert "location_snapshot=" not in source
    assert "locations=" not in source
    assert (
        "blinded_probability_review.attrs.get('review_context_artifact_content_id')"
        in source
    )
    assert "!= ingest_manifest.output.content_sha256" in source
    assert "cells_glob =" not in source
    assert "validate_dataset_snapshot_current(ingest_manifest.output, mode='content')" in source
    assert (
        "(Path(ingest_manifest.output.root) / 'donor_id=*' / "
        "'*.parquet').as_posix()" in source
    )
    assert "[cells_glob]" not in source
    assert "sample_locations" not in source
    assert source.index("ingest_references[0].validate_current(") < source.index(
        "blinded_review_sample_from_dataset("
    )


def test_notebook_05_exports_only_the_blinded_probability_review_view():
    source = "\n".join(_code_cells())
    assert "PHENOCYCLER_REVIEW_BLINDING_KEY_FILE" in source
    assert "review_blinding_key_path.read_bytes()" in source
    assert "build_review_sampling_frame(" in source
    assert "validate_review_sampling_against_assignments(" in source
    assert "split_name='validation'" in source
    assert "level='broad', parent=None" in source
    assert "blinding_key=review_blinding_key" in source
    assert "candidate_assignment_run_id=assignment_run_id" in source
    assert "BASELINE_PER_STRATUM" not in source
    assert "SAMPLE_SALT" not in source
    assert "blinded_stratum_codes" not in source
    assert "REVIEW_SAMPLING_METHOD" not in source
    assert "display(probability_sampling_audit" not in source
    assert "display(blinded_probability_review.head" in source
    assert "blinded_review_sample_from_dataset(" in source
    assert "probability_sampling_audit, source_run_manifest=source_run_manifest" in source
    assert "blinded_review_sample_frame(" not in source
    assert "review_sample_frame_fingerprint(\n    probability_sampling_audit" in source
    assert "sample_frame=probability_sampling_audit" in source
    assert "sample_frame=blinded_probability_review" not in source
    assert "reviewer_frame=blinded_probability_review" in source
    assert "assignment_run_id = candidate_evaluation.assignment_run_id" in source
    assert "assignment_run_id=assignment_run_id" in source
    assert "typing_model_bundle_content_id=model_bundle.content_id" in source
    assert (
        "candidate_evaluation_manifest_content_id=candidate_evaluation.content_id"
        in source
    )
    assert (
        "candidate_assignment_output_content_id="
        "candidate_evaluation.output_content_sha256" in source
    )
    assert "review_context_artifact_content_id=" not in source
    validation_index = source.index("reference_labels = validate_review_labels(")
    audit_join_index = source.index("probability_sampling_audit.merge(")
    assignment_join_index = source.index("sample_assignments, on=['donor_id', 'object_id']")
    assert validation_index < audit_join_index < assignment_join_index


def test_notebook_05_requires_schema_v3_frozen_thresholds_and_stays_broad_only():
    source = "\n".join(_code_cells())
    all_source = _all_cell_source()
    assert "model_bundle.threshold_selection_provenance.schema_version == 3" in source
    assert "candidate_evaluation.split_name != 'validation'" in source
    assert "review_provenance.split_name != 'validation'" in source
    assert source.count("allow_locked_test=False") >= 2
    assert "level='subtype'" not in source
    assert "level='specific'" not in source
    assert "schema-v5 `ReviewLabelProvenance`" in all_source
    assert "mandatory all-cell `level='specific', parent=None`" in all_source
    assert "cannot satisfy promotion" in all_source


def test_notebook_05_documents_preregistered_release_governance():
    source = _all_cell_source()
    assert "promotion_policy_content_id" in source
    assert "signing_key_sha256" in source
    assert "both must be nonblank for production promotion" in source
    assert "external signing key" in source


def test_notebook_05_reports_zero_rescue_calls_separately_from_total_coverage():
    source = "\n".join(_code_cells())
    assert "'scope': 'total named calls'" in source
    assert "'scope': 'probabilistic rescue lane'" in source
    assert "['Ambiguous', 'Unavailable', 'Other']" in source
    assert "structural_rescue_lane = reviewed['broad_threshold_eligible'].astype(bool)" in source
    assert "ZERO PROBABILISTIC-RESCUE CALLS" in source
    assert "Total versus probabilistic-rescue coverage" in source
    assert (
        "challenge_margin_threshold = float(model_bundle.thresholds.inferred_margin)"
        in source
    )
    assert (
        "challenge_probability_threshold = "
        "float(model_bundle.thresholds.inferred_probability)" in source
    )
    assert "CHALLENGE_RANK_CONTEXT = f'{assignment_run_id}|" in source
    assert "CHALLENGE_RANK_CONTEXT = f'{candidate_evaluation.content_id}|" not in source
    assert "a.broad_margin < 0.25" not in source
    assert "a.broad_best_probability >= 0.80" not in source


def test_notebook_05_code_cells_compile():
    for index, source in enumerate(_code_cells()):
        compile(source, f"{NOTEBOOK}:code-cell-{index}", "exec")
