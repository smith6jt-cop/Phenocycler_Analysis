import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "06_probabilistic_validation_and_promotion.ipynb"
WORKFLOW = ROOT / "docs" / "REFINEMENT_WORKFLOW.md"


def _payload() -> dict:
    return json.loads(NOTEBOOK.read_text(encoding="utf-8"))


def _code_cells() -> list[str]:
    return [
        "".join(cell.get("source", ()))
        for cell in _payload()["cells"]
        if cell.get("cell_type") == "code"
    ]


def _all_source() -> str:
    return "\n".join(
        "".join(cell.get("source", ())) for cell in _payload()["cells"]
    )


def test_notebook_06_is_valid_output_free_and_code_cells_compile():
    payload = _payload()
    assert payload["nbformat"] == 4
    assert payload["cells"]
    for index, cell in enumerate(payload["cells"]):
        if cell["cell_type"] != "code":
            continue
        assert cell.get("execution_count") is None
        assert cell.get("outputs") == []
        compile("".join(cell.get("source", ())), f"{NOTEBOOK}:cell-{index}", "exec")


def test_validation_is_frozen_and_passes_before_explicit_locked_open():
    source = "\n".join(_code_cells())
    validation_build = source.index("validation_report = build_split_evaluation_report(")
    validation_gate = source.index("validation_report.validate_for_promotion(")
    validation_write = source.index("validation_report.write_json(VALIDATION_REPORT_PATH)")
    locked_guard = source.index("OPEN_LOCKED_TEST = False")
    locked_build = source.index("locked_test_report = build_split_evaluation_report(")
    assert validation_build < validation_write < validation_gate < locked_guard < locked_build
    assert source.count("OPEN_LOCKED_TEST = False") == 1
    assert "OPEN_LOCKED_TEST = True" not in source
    assert "allow_locked_test=False" in source[validation_build:locked_guard]
    assert "validation_report.split_name != 'validation'" in source
    assert "if not VALIDATION_GATE_PASSED or validation_report is None:" in source
    archive_reload = source.index(
        "archived_validation_report = SplitEvaluationReport.read_json(", locked_guard
    )
    boundary_replay = source.index(
        "archived_validation_report.validate_candidate_source_replay(", locked_guard
    )
    boundary_gate = source.index(
        "archived_validation_report.validate_for_promotion(", locked_guard
    )
    locked_candidate_read = source.index(
        "locked_test_candidate = load_candidate_evaluation_manifest(", locked_guard
    )
    assert locked_guard < archive_reload < boundary_replay < boundary_gate < locked_candidate_read
    assert "archived_validation_report.content_id != validation_report.content_id" in source
    assert "if (BUILD_LOCKED_TEST_REPORT or LOAD_EXISTING_LOCKED_TEST_REPORT) and not OPEN_LOCKED_TEST:" in source
    assert "allow_locked_test=True" in source[locked_guard:locked_build + 2000]
    assert "locked_test_report.split_name != 'locked_test'" in source
    assert source.count("validate_candidate_source_replay(") >= 4
    assert "locked_test_report.validate_for_promotion(" in source[locked_build:]
    assert "locked_test_report.write_json(LOCKED_TEST_REPORT_PATH)" in source[locked_build:]


def test_promotion_uses_exact_calibration_source_replay_and_current_create_api():
    cells = _code_cells()
    promotion_cell = next(
        cell for cell in cells if "TypingModelPromotionManifest.create(" in cell
    )
    required_replay_fields = (
        "preliminary_bundle=preliminary_bundle",
        "calibration_candidate_manifest=calibration_candidate",
        "marker_registry=marker_registry",
        "typing_registry=typing_registry",
        "source_run_manifest=source_run_manifest",
        "calibration_sample_frames=",
        "calibration_reviewer_frames=",
        "calibration_review_ledgers=",
        "review_blinding_key=review_blinding_key",
    )
    assert "threshold_source_replay_inputs = ThresholdSelectionReplayInputs(" in promotion_cell
    for field in required_replay_fields:
        assert field in promotion_cell
    required_create_arguments = (
        "promotion_version=PROMOTION_VERSION",
        "model_bundle=frozen_bundle",
        "policy=promotion_policy",
        "validation_report=validation_report",
        "locked_test_report=locked_test_report",
        "threshold_source_replay_inputs=threshold_source_replay_inputs",
        "marker_registry=marker_registry",
        "typing_registry=typing_registry",
        "review_blinding_key=review_blinding_key",
        "release_signing_key=release_signing_key",
        "issuer=PROMOTION_ISSUER",
        "issued_at_utc=datetime.now(timezone.utc).isoformat()",
    )
    for argument in required_create_arguments:
        assert argument in promotion_cell
    assert promotion_cell.index("if not OPEN_LOCKED_TEST:") < promotion_cell.index(
        "TypingModelPromotionManifest.create("
    )
    assert "promotion_manifest.write_json(PROMOTION_MANIFEST_PATH)" in promotion_cell


def test_keys_are_loaded_from_protected_files_and_never_rendered():
    source = "\n".join(_code_cells())
    assert "PHENOCYCLER_REVIEW_BLINDING_KEY_FILE" in source
    assert "PHENOCYCLER_RELEASE_SIGNING_KEY_FILE" in source
    assert "key_path.read_bytes()" in source
    assert "release_signing_key = _read_protected_key(RELEASE_KEY_ENV)" in source
    assert "del release_signing_key" in source
    assert "del review_blinding_key" in source
    forbidden = (
        "print(release_signing_key",
        "print(review_blinding_key",
        "display(release_signing_key",
        "display(review_blinding_key",
        "repr(release_signing_key",
        "repr(review_blinding_key",
    )
    assert not any(item in source for item in forbidden)


def test_plots_use_replayed_intervals_and_explicit_fp_fn_coverage_gates():
    source = "\n".join(_code_cells())
    assert "for metric in ('coverage', 'rescue_coverage', 'selective_error', 'rescue_selective_error')" in source
    assert "result.false_positive_rate" in source
    assert "result.false_negative_rate" in source
    assert "interval.estimate" in source
    assert "interval.lower" in source
    assert "interval.upper" in source
    assert "policy.minimum_coverage_lower" in source
    assert "policy.minimum_rescue_coverage_lower" in source
    assert "policy.maximum_selective_error_upper" in source
    assert "policy.maximum_rescue_selective_error_upper" in source
    assert "policy.maximum_class_false_positive_rate_upper" in source
    assert "policy.maximum_class_false_negative_rate_upper" in source
    assert "yerr=np.vstack" in source


def test_challenge_sample_is_diagnostic_only_and_absent_from_release_calls():
    source = _all_source()
    cells = _code_cells()
    report_cells = [cell for cell in cells if "build_split_evaluation_report(" in cell]
    promotion_cell = next(
        cell for cell in cells if "TypingModelPromotionManifest.create(" in cell
    )
    assert "CHALLENGE_DIAGNOSTIC_ONLY = True" in source
    assert "diagnostic only" in source.lower()
    assert "not a software-enforced headline metric gate" in source
    assert "expected_lineage" in source
    assert "candidate_evaluation_manifest_content_id" in source
    assert "candidate_assignment_output_content_id" in source
    assert "CHALLENGE_DISPOSITION_READY" in source
    assert "catastrophic challenge findings require a new development cycle" in source
    assert all("challenge_diagnostic" not in cell for cell in report_cells)
    assert "challenge_diagnostic" not in promotion_cell


def test_workflow_links_the_exact_notebook_06_name():
    workflow = WORKFLOW.read_text(encoding="utf-8")
    expected = (
        "[06_probabilistic_validation_and_promotion]"
        "(../notebooks/06_probabilistic_validation_and_promotion.ipynb)"
    )
    assert expected in workflow
    assert "`06_locked_validation_and_promotion`" not in workflow
    assert "Notebooks 01–06 are operator interfaces" in workflow
    assert "Notebooks 01–04 retain their existing guarded package-stage" in workflow
    assert "notebook 05 is read-only" in workflow
    assert "probability-sample review is complete and its software-enforced" in workflow
    assert "it cannot supply headline\n   metrics or software promotion pass/fail" in workflow
