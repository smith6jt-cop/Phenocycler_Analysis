"""Unit tests for the Gate-4 RESTORE production apply (restore_apply.py).

Two tiers:
  * Tier 1 — the pure per-cell tri-state mapper ``cell_calls_from_evaluation`` (no I/O).
  * Tier 2 — the ``run_apply`` / ``_apply_one_donor`` orchestration. The FROZEN math
    (``evaluate_locked_pair`` / ``load_validation_sample``) is already covered by
    test_restore_validation.py's 64 tests, so here it is stubbed to isolate the adapter plumbing:
    the accepted-pair loop, hard-fail halting, the tri-state join, atomic writes, and the manifest.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import phenocycler.restore_apply as ra
from phenocycler import lineage, load_config
from phenocycler.restore_validation import ACCEPTED_PAIRS, PairState


# --------------------------------------------------------------------------- #
# Tier 1 — pure mapper
# --------------------------------------------------------------------------- #

def test_valid_calls_by_divisor():
    """VALID: positive where value>divisor, negative where measured-and-not, unavailable off valid_idx."""
    ev = {
        "state": "VALID",
        "n_input": 4,
        "valid_idx": np.array([0, 1, 2]),                     # cell 3 lacks a required measurement
        "full_lineage_positive": np.array([True, False, False, False]),
        "full_normalized": np.array([2.0, 0.5, 0.5, np.nan]),
    }
    pos, norm, cst = ra.cell_calls_from_evaluation(ev)
    assert list(cst) == ["positive", "negative", "negative", "unavailable"]
    assert pos[0] and not pos[1:].any()
    assert norm[0] == 2.0 and np.isnan(norm[3])


def test_confirmed_absence_is_negative_not_unavailable():
    """NO_TARGET_EXPRESSION_CONFIRMED -> every cell a confident NEGATIVE (never 'unavailable')."""
    ev = {
        "state": "NO_TARGET_EXPRESSION_CONFIRMED",
        "n_input": 3,
        "valid_idx": np.array([0, 1, 2]),
        "full_lineage_positive": np.zeros(3, bool),
        "full_normalized": np.full(3, np.nan),
    }
    pos, norm, cst = ra.cell_calls_from_evaluation(ev)
    assert list(cst) == ["negative", "negative", "negative"]
    assert not pos.any()


@pytest.mark.parametrize("state", sorted(s.value for s in ra.HARD_FAIL_STATES))
def test_technical_states_map_to_unavailable(state):
    """Every technical (non-callable) state maps to all-unavailable (the driver halts before writing these)."""
    ev = {
        "state": state,
        "n_input": 3,
        "valid_idx": np.array([0, 1, 2]),
        "full_lineage_positive": np.zeros(3, bool),
        "full_normalized": np.full(3, np.nan),
    }
    pos, _, cst = ra.cell_calls_from_evaluation(ev)
    assert list(cst) == ["unavailable", "unavailable", "unavailable"]
    assert not pos.any()


def test_state_sets_partition_pairstate():
    """Callable vs hard-fail sets partition PairState exactly (guards a newly added state)."""
    assert ra.CALLABLE_STATES | ra.HARD_FAIL_STATES == frozenset(PairState)
    assert ra.CALLABLE_STATES.isdisjoint(ra.HARD_FAIL_STATES)
    assert ra.CALLABLE_STATES == {PairState.VALID, PairState.NO_TARGET_EXPRESSION_CONFIRMED}


def test_required_targets_match_accepted_pairs():
    assert ra.REQUIRED_TARGETS == tuple(t for t, _ in ACCEPTED_PAIRS)


# --------------------------------------------------------------------------- #
# Tier 2 — driver orchestration (frozen math stubbed)
# --------------------------------------------------------------------------- #

def _write_manifest(cfg, drop=()):
    rows = [
        {"target": t, "reference": r, "review": "ACCEPTED", "evidence": "gate-2 sign-off (test)"}
        for t, r in ACCEPTED_PAIRS
        if (t, r) not in set(drop)
    ]
    cfg.restore_pair_reviews_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(cfg.restore_pair_reviews_csv, index=False)


def _stub_frozen(monkeypatch, *, sample_n=6, canonical_n=None, state_for=None):
    """Stub load_validation_sample / evaluate_locked_pair / locked_pair_metrics for the adapter tests."""
    state_for = state_for or {}

    def fake_load(sample_root, raw_cells, redsea_cells, donor):
        df = pd.DataFrame({"object_id": [f"{donor}-{i}" for i in range(sample_n)]})
        audit = {"n": sample_n, "canonical_n": sample_n if canonical_n is None else canonical_n,
                 "image": f"img-{donor}"}
        return df, audit

    def fake_eval(donor_df, donor, target, reference, *, config=None,
                  expression_review=None, pair_review=None):
        n = len(donor_df)
        state = state_for.get(target, PairState.VALID)
        valid = np.arange(n)
        pos = np.zeros(n, bool)
        norm = np.full(n, np.nan)
        if state is PairState.VALID:
            pos[::2] = True                                   # even cells positive for every accepted target
            norm[valid] = np.where(pos[valid], 2.0, 0.5)
        return {"state": state.value, "n_input": n, "valid_idx": valid,
                "full_lineage_positive": pos, "full_normalized": norm,
                "canonical_divisor": 1.0 if state is PairState.VALID else None,
                "donor": str(donor), "target": target, "reference": reference}

    def fake_metrics(ev):
        return {"donor": ev["donor"], "target": ev["target"], "reference": ev["reference"],
                "state": ev["state"], "canonical_divisor": ev["canonical_divisor"],
                "n_input": ev["n_input"], "reference_to_target_ratio": 1.5,
                "divisor_reproducibility": 1.0, "target_low_sbr": False}

    monkeypatch.setattr(ra, "load_validation_sample", fake_load)
    monkeypatch.setattr(ra, "evaluate_locked_pair", fake_eval)
    monkeypatch.setattr(ra, "locked_pair_metrics", fake_metrics)


def _setup(tmp_path, monkeypatch, *, make_sample=True, make_manifest=True, drop=(), **stub):
    cfg = load_config(data_dir=tmp_path)
    # The policy keys are UNSET by default and config_from_pipeline refuses to guess -- see
    # test_unset_policy_refuses_rather_than_defaulting. Every test that actually runs an apply must
    # therefore state the policy, exactly as a real run does via the parent repo's config.ini.
    cfg.restore_control_definition = "reference_only"
    cfg.restore_divisor_statistic = "p99"
    if make_sample:
        cfg.restore_allcells_sample_dir.mkdir(parents=True, exist_ok=True)
    if make_manifest:
        _write_manifest(cfg, drop=drop)
    _stub_frozen(monkeypatch, **stub)
    return cfg


def test_run_apply_writes_gated_and_manifest(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, sample_n=6)
    manifest = ra.run_apply(cfg, donors=["TEST"])

    gated = pd.read_parquet(cfg.restore_gated_dir / "donor_id=TEST" / "data_0.parquet")
    assert list(gated["donor_id"].unique()) == ["TEST"]
    for t, _ in ACCEPTED_PAIRS:
        for suffix in ("_pos", "_norm", "_state"):
            assert f"{t}{suffix}" in gated.columns
        assert set(gated[f"{t}_state"].unique()).issubset({"positive", "negative", "unavailable"})
    # even cells positive, odd cells negative for every target
    assert (gated["CD3e_state"] == "positive").sum() == 3

    # cohort manifest: one row per accepted pair for this donor
    assert len(manifest) == len(ACCEPTED_PAIRS)
    assert set(manifest["state"]) == {"VALID"}
    assert cfg.restore_apply_manifest_csv.exists()
    # per-donor json manifest
    assert (cfg.restore_gated_dir / "donor_id=TEST" / "apply_manifest.json").exists()


def test_apply_output_atomic_and_overwrites(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, sample_n=4)
    ra.run_apply(cfg, donors=["TEST"])
    part = cfg.restore_gated_dir / "donor_id=TEST"
    assert not list(part.glob("*.tmp*"))                      # atomic write left no temp file
    # a second run replaces the partition cleanly (no partial/leftover files)
    ra.run_apply(cfg, donors=["TEST"])
    assert not list(part.glob("*.tmp*"))
    assert (part / "data_0.parquet").exists()


def test_end_to_end_apply_then_lineage(tmp_path, monkeypatch):
    """restore_apply output feeds lineage.assign_donor: even cells (all-positive) -> Immune; odd -> Other."""
    cfg = _setup(tmp_path, monkeypatch, sample_n=6)
    ra.run_apply(cfg, donors=["TEST"])
    gated = str(cfg.restore_gated_dir / "donor_id=TEST" / "data_0.parquet")
    out = lineage.assign_donor("TEST", gated, cfg.cells_dir).set_index("object_id")
    assert (out.loc[["TEST-0", "TEST-2", "TEST-4"], "compartment"] == "Immune").all()
    assert (out.loc[["TEST-1", "TEST-3", "TEST-5"], "compartment"] == "Other").all()


def test_missing_sample_dir_fails_loud(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, make_sample=False)
    with pytest.raises(SystemExit, match="sample"):
        ra.run_apply(cfg, donors=["TEST"])


def test_missing_pair_manifest_fails_loud(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, make_manifest=False)
    with pytest.raises(SystemExit, match="manifest"):
        ra.run_apply(cfg, donors=["TEST"])


def test_pair_not_accepted_fails_loud(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, drop=[("CD11b", "EpCAM")])
    with pytest.raises(SystemExit, match="ACCEPTED"):
        ra.run_apply(cfg, donors=["TEST"])


def test_full_coverage_guard(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, sample_n=4, canonical_n=10)  # sample covers 4/10 canonical cells
    with pytest.raises(ValueError, match="canonical"):
        ra.run_apply(cfg, donors=["TEST"])


def test_hard_fail_state_halts(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, state_for={"CD3e": PairState.MODEL_UNSTABLE})
    with pytest.raises(ValueError, match="MODEL_UNSTABLE"):
        ra.run_apply(cfg, donors=["TEST"])


def test_confirmed_absence_state_is_written_not_halted(tmp_path, monkeypatch):
    """A maintainer-confirmed absence is a CALLABLE state: it writes (all-negative), never halts."""
    cfg = _setup(tmp_path, monkeypatch, state_for={"CD11b": PairState.NO_TARGET_EXPRESSION_CONFIRMED})
    ra.run_apply(cfg, donors=["TEST"])
    gated = pd.read_parquet(cfg.restore_gated_dir / "donor_id=TEST" / "data_0.parquet")
    assert (gated["CD11b_state"] == "negative").all()
    assert (gated["CD3e_state"] == "positive").sum() == 3     # other pairs unaffected


# --------------------------------------------------------------------------- #
# Tier 3 — the v10 pair-method policy reaches the frozen math
#
# Before this, `control_definition` / `divisor_statistic` were reachable ONLY from the
# `restore_validation evaluate-locked` CLI, so a production apply always ran the dataclass defaults
# no matter which sweep arm had been chosen. These tests pin the config -> evaluate_locked_pair path.
# --------------------------------------------------------------------------- #

def test_config_from_pipeline_reads_ini_policy(tmp_path):
    cfg = load_config(data_dir=tmp_path)
    cfg.restore_control_definition = "reference_and_target"
    cfg.restore_divisor_statistic = "p99"
    pvc = ra.config_from_pipeline(cfg)
    assert pvc.control_definition == "reference_and_target"
    assert pvc.divisor_statistic == "p99"


def test_unset_policy_refuses_rather_than_defaulting(tmp_path):
    """An UNSET policy must fail loudly and name the file to fix.

    This replaces `test_config_defaults_match_the_frozen_dataclass`, which asserted that PipelineConfig
    mirrored PairValidationConfig's defaults so the indirection would be "behaviour-neutral". That
    reasoning is what created the trap: PairValidationConfig defaults to `divisor_statistic="max"`,
    which restore_apply CANNOT run (MODEL_UNSTABLE on 115 Vimentin<-E_cadherin, 6436 B3TUBB<-EpCAM,
    6591 CD68<-E_cadherin), while the frozen production `p99` lives only in the PARENT repo's
    config.ini. A run launched without `--config config.ini` therefore silently selected the unrunnable
    policy instead of the frozen one. Behaviour-neutral defaults are not neutral when the default is
    wrong.
    """
    cfg = load_config(data_dir=tmp_path)          # submodule config.ini declares neither key
    assert cfg.restore_control_definition == ""
    assert cfg.restore_divisor_statistic == ""
    with pytest.raises(SystemExit, match="is not set"):
        ra.config_from_pipeline(cfg)


def test_each_policy_key_is_named_individually_when_missing(tmp_path):
    """Half-configured is as dangerous as unconfigured, so each key is checked on its own."""
    cfg = load_config(data_dir=tmp_path)
    cfg.restore_control_definition = "reference_only"
    with pytest.raises(SystemExit, match="divisor_statistic is not set"):
        ra.config_from_pipeline(cfg)

    cfg = load_config(data_dir=tmp_path)
    cfg.restore_divisor_statistic = "p99"
    with pytest.raises(SystemExit, match="control_definition is not set"):
        ra.config_from_pipeline(cfg)


def test_the_repo_root_config_ini_carries_the_frozen_policy(tmp_path):
    """The guard is only safe if the parent config.ini really does set both keys."""
    root_ini = Path(__file__).resolve().parents[2] / "config.ini"
    if not root_ini.exists():                      # submodule checked out standalone
        pytest.skip("parent repository config.ini is not present")
    pvc = ra.config_from_pipeline(load_config(root_ini))
    assert pvc.control_definition == "reference_only"
    assert pvc.divisor_statistic == "p99"


@pytest.mark.parametrize(
    "attr,value,match",
    [
        ("restore_control_definition", "reference_and_targe", "control_definition"),
        ("restore_divisor_statistic", "p90", "divisor_statistic"),
    ],
)
def test_invalid_policy_fails_before_any_donor(tmp_path, attr, value, match):
    """A typo in config.ini must fail with the allowed values named, not mid-run inside the frozen math."""
    cfg = load_config(data_dir=tmp_path)
    setattr(cfg, attr, value)
    with pytest.raises(SystemExit, match=match):
        ra.config_from_pipeline(cfg)


def test_run_apply_passes_ini_policy_to_evaluate_locked_pair(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, sample_n=4)
    cfg.restore_control_definition = "reference_and_target"
    cfg.restore_divisor_statistic = "p95"

    seen = []
    inner = ra.evaluate_locked_pair

    def capture(donor_df, donor, target, reference, *, config=None, **kw):
        seen.append(config)
        return inner(donor_df, donor, target, reference, config=config, **kw)

    monkeypatch.setattr(ra, "evaluate_locked_pair", capture)
    ra.run_apply(cfg, donors=["TEST"])

    assert len(seen) == len(ACCEPTED_PAIRS)
    assert {c.control_definition for c in seen} == {"reference_and_target"}
    assert {c.divisor_statistic for c in seen} == {"p95"}


def test_explicit_config_overrides_the_ini(tmp_path, monkeypatch):
    """An explicitly passed PairValidationConfig wins over config.ini (pilots and tests)."""
    from phenocycler.restore_validation import PairValidationConfig

    cfg = _setup(tmp_path, monkeypatch, sample_n=4)
    cfg.restore_divisor_statistic = "p95"

    seen = []
    inner = ra.evaluate_locked_pair

    def capture(donor_df, donor, target, reference, *, config=None, **kw):
        seen.append(config)
        return inner(donor_df, donor, target, reference, config=config, **kw)

    monkeypatch.setattr(ra, "evaluate_locked_pair", capture)
    ra.run_apply(cfg, donors=["TEST"], config=PairValidationConfig(divisor_statistic="p999"))
    assert {c.divisor_statistic for c in seen} == {"p999"}


# --------------------------------------------------------------------------- #
# Technical failure policy (2026-07-27)
#
# At 7 universally-expressed pairs, halting the cohort on one donor-pair was right. At 24 pairs x 20
# donors -- sparse immune markers, hormones that legitimately fail on individual donors -- it blocks
# 19 good donors on one bad donor-marker. The tri-state design already has the honest answer.
# What NEITHER policy may do is absorb the failure as a biological negative.
# --------------------------------------------------------------------------- #

def test_halt_is_the_default(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, state_for={"CD3e": PairState.MODEL_UNSTABLE})
    with pytest.raises(ValueError, match="MODEL_UNSTABLE"):
        ra.run_apply(cfg, donors=["TEST"])


def test_unavailable_policy_blocks_rather_than_negating(tmp_path, monkeypatch):
    cfg = _setup(tmp_path, monkeypatch, sample_n=6, state_for={"CD3e": PairState.MODEL_UNSTABLE})
    ra.run_apply(cfg, donors=["TEST"], on_technical_failure="unavailable")
    gated = pd.read_parquet(cfg.restore_gated_dir / "donor_id=TEST" / "data_0.parquet")
    # the failed marker is UNAVAILABLE for every cell -- never 'negative', which would let a lower
    # gate silently claim the cell as though CD3e had been measured and found absent
    assert (gated["CD3e_state"] == "unavailable").all()
    assert not (gated["CD3e_state"] == "negative").any()
    assert (gated["CD68_state"] == "positive").sum() == 3          # other pairs unaffected
    # and it is recorded, not swallowed
    failures = pd.read_csv(cfg.restore_gated_dir / "technical_failures.csv")
    assert failures.target.tolist() == ["CD3e"]
    assert failures.state.tolist() == ["MODEL_UNSTABLE"]


def test_unavailable_marker_makes_lineage_unresolved_not_other(tmp_path, monkeypatch):
    """The whole point: a blocked cell must be Unresolved, so it is visibly un-typed rather than
    quietly counted as a negative in the residual."""
    cfg = _setup(tmp_path, monkeypatch, sample_n=6, state_for={"CD3e": PairState.MODEL_UNSTABLE})
    ra.run_apply(cfg, donors=["TEST"], on_technical_failure="unavailable")
    gated = str(cfg.restore_gated_dir / "donor_id=TEST" / "data_0.parquet")
    out = lineage.assign_donor("TEST", gated, cfg.cells_dir).set_index("object_id")
    # odd cells are negative for every other gate; CD3e is unavailable, so Immune is indeterminate
    assert (out.loc[["TEST-1", "TEST-3", "TEST-5"], "compartment"] == "Unresolved").all()
    assert (out.loc[["TEST-1", "TEST-3", "TEST-5"], "assignment_reason"] == "blocked_unavailable_gate").all()
