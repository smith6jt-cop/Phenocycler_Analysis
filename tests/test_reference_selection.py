import pandas as pd
import pytest

from phenocycler import reference_selection as rs


def _screen_row(target, reference, ratio, divisor, *, undersized, donor, fold=5.0,
                called=0.1, n_assigned=10_000, target_n=2000):
    return {
        "donor": donor,
        "target": target,
        "reference": reference,
        "reference_n": int(round(ratio * target_n)),
        "target_n": target_n,
        "reference_to_target_ratio": ratio,
        "reference_undersized": bool(undersized),
        "candidate_full_maximum": divisor,
        "reference_fold": fold,
        "threshold_positive_n": int(called * n_assigned),
        "n_assigned": n_assigned,
    }


def _synthetic_screen():
    """20 donors; mirrors the real CD11b / CD68 / E_cadherin frequency structure."""
    rows = []
    for d in range(20):
        donor = f"d{d:02d}"
        # CD11b: baseline E_cadherin passes (marginal -- undersized on 6/20 donors),
        # EpCAM a stronger comparator, CD3e fails hard (undersized on every donor + low divisor).
        rows.append(_screen_row("CD11b", "E_cadherin", 1.5, 13000, undersized=(d < 6), donor=donor))
        rows.append(_screen_row("CD11b", "EpCAM", 3.6, 13000, undersized=(d < 3), donor=donor))
        rows.append(_screen_row("CD11b", "CD3e", 0.04, 9900, undersized=True, donor=donor, called=0.24))
        # CD68: baseline E_cadherin passes; CD3e fails.
        rows.append(_screen_row("CD68", "E_cadherin", 7.4, 300, undersized=False, donor=donor))
        rows.append(_screen_row("CD68", "CD3e", 0.2, 290, undersized=True, donor=donor))
        # E_cadherin: only CD31, undersized on every donor -> open.
        rows.append(_screen_row("E_cadherin", "CD31", 0.19, 1000, undersized=True, donor=donor))
    return pd.DataFrame(rows)


def test_reference_matrix_ratio_is_pure_proportion():
    pe = _synthetic_screen()
    m = rs.build_reference_matrix(pe)
    row = m[(m.target == "CD11b") & (m.reference == "EpCAM")].iloc[0]
    assert row["ratio_median"] == pytest.approx(3.6)
    assert bool(row["passes_frequency_rule"]) is True


def test_undersized_references_are_rejected_and_are_the_divisor_outlier():
    pe = _synthetic_screen()
    m = rs.build_reference_matrix(pe).set_index(["target", "reference"])
    assert m.loc[("CD11b", "CD3e"), "recommendation"] == "rejected"
    assert m.loc[("CD68", "CD3e"), "recommendation"] == "rejected"
    assert m.loc[("E_cadherin", "CD31"), "recommendation"] == "rejected"
    # divisor-consistency cross-check: the undersized CD3e reference is also the low-divisor outlier
    assert m.loc[("CD11b", "CD3e"), "divisor_vs_consensus_median"] < 0.95


def test_recommendation_keeps_passing_baseline_and_opens_when_none_pass():
    pe = _synthetic_screen()
    m = rs.build_reference_matrix(pe)
    decisions = {d.target: d for d in rs.recommend_references(m)}
    # CD11b baseline E_cadherin passes -> kept, with a marginal note naming the stronger EpCAM
    assert decisions["CD11b"].recommended_reference == "E_cadherin"
    assert decisions["CD11b"].status == "accepted"
    assert "EpCAM" in decisions["CD11b"].rationale
    # CD68 baseline passes -> kept
    assert decisions["CD68"].recommended_reference == "E_cadherin"
    # E_cadherin has no passing candidate -> open; no reference is invented
    assert decisions["E_cadherin"].recommended_reference is None
    assert decisions["E_cadherin"].status == "open"


def test_maintainer_overrides_win_including_a_failing_reference():
    pe = _synthetic_screen()
    overrides = {"CD11b": "EpCAM", "E_cadherin": "CD31"}
    m = rs.build_reference_matrix(pe, overrides=overrides).set_index(["target", "reference"])
    # CD11b switched to the stronger EpCAM; the marginal baseline drops to comparator
    assert m.loc[("CD11b", "EpCAM"), "recommendation"] == "accepted"
    assert m.loc[("CD11b", "E_cadherin"), "recommendation"] == "comparator"
    assert m.loc[("CD11b", "CD3e"), "recommendation"] == "rejected"
    # the override accepts CD31 for E_cadherin even though it FAILS the frequency rule
    assert m.loc[("E_cadherin", "CD31"), "recommendation"] == "accepted"
    decisions = {
        d.target: d
        for d in rs.recommend_references(
            rs.build_reference_matrix(pe, overrides=overrides), overrides=overrides
        )
    }
    assert decisions["CD11b"].recommended_reference == "EpCAM"
    assert decisions["E_cadherin"].recommended_reference == "CD31"
    assert decisions["E_cadherin"].status == "accepted"
    assert "advisory" in decisions["E_cadherin"].rationale


def test_frozen_overrides_constant_matches_maintainer_decisions():
    assert rs.ACCEPTED_REFERENCE_OVERRIDES == {"CD11b": "EpCAM", "E_cadherin": "CD31"}


def test_build_reference_matrix_requires_diagnostic_columns():
    pe = _synthetic_screen().drop(columns=["reference_to_target_ratio"])
    with pytest.raises(ValueError, match="reference_to_target_ratio"):
        rs.build_reference_matrix(pe)


# --------------------------------------------------------------------------- #
# Promotion candidates — Pan-CK / Ker8_18 / EpCAM as TARGETS (2026-07-27)
#
# These have no entry in BASELINE_PAIRS: nothing was predeclared for them, so the dossier must not
# report "the baseline fails", which would read as though a reference had been tried and rejected.
# --------------------------------------------------------------------------- #

def _promotion_screen():
    rows = []
    for d in range(20):
        donor = f"d{d:02d}"
        # Pan_Cytokeratin: both candidate references pass; Vimentin is the stronger one.
        rows.append(_screen_row("Pan_Cytokeratin", "CD31", 1.4, 2000, undersized=False, donor=donor))
        rows.append(_screen_row("Pan_Cytokeratin", "Vimentin", 3.1, 2050, undersized=False, donor=donor))
        # EpCAM: neither passes -> open, and still no baseline to blame.
        rows.append(_screen_row("EpCAM", "CD31", 0.3, 800, undersized=True, donor=donor))
        rows.append(_screen_row("EpCAM", "Vimentin", 0.6, 810, undersized=True, donor=donor))
    return pd.DataFrame(rows)


def test_promotion_candidates_are_in_the_target_order():
    """Otherwise build_bundle silently drops them: it iterates _TARGET_ORDER, not the matrix."""
    for target in ("Pan_Cytokeratin", "Ker8_18", "EpCAM"):
        assert target in rs._TARGET_ORDER


def test_promotion_candidate_recommendation_does_not_blame_a_baseline():
    from phenocycler.restore_validation import BASELINE_PAIRS
    assert "Pan_Cytokeratin" not in {t for t, _ in BASELINE_PAIRS}      # premise of this test
    m = rs.build_reference_matrix(_promotion_screen())
    d = {x.target: x for x in rs.recommend_references(m)}["Pan_Cytokeratin"]
    assert d.status == "accepted"
    assert d.recommended_reference == "Vimentin"                        # strongest passing
    assert "no predeclared baseline" in d.rationale
    assert "fails the frequency rule" not in d.rationale
    assert "Gate-1->3" in d.rationale                                   # promotion is not automatic


def test_promotion_candidate_with_no_passing_reference_is_open():
    m = rs.build_reference_matrix(_promotion_screen())
    d = {x.target: x for x in rs.recommend_references(m)}["EpCAM"]
    assert d.status == "open"
    assert d.recommended_reference is None
