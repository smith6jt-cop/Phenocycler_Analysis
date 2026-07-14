"""Unit tests for the extracted diagnostic helpers (no imaging data — synthetic arrays)."""
from __future__ import annotations

import numpy as np
import pytest

from phenocycler import restore
from phenocycler.diagnostics import read_donor, r_otsu  # noqa: F401  (import-smoke of the diagnostics module)


# --------------------------------------------------------------------------- #
# idx_select — behavior-preservation vs the pre-refactor inline patch math
# --------------------------------------------------------------------------- #
def _old_patch_idx_select(T, R, neg_q, ratio_x=0.75, ratio_y=0.5):
    """Verbatim pre-refactor math from patch_restore_idx_floor (restore.py, old lines 177-180)."""
    c0 = np.quantile(T, neg_q); c1 = np.quantile(R, neg_q)
    return (((R > np.quantile(R, ratio_y)) & (T <= c0)) |
            ((T > np.quantile(T, ratio_x)) & (R <= c1)))


def test_idx_select_matches_old_patch_math():
    rng = np.random.default_rng(0)
    T = rng.gamma(2.0, 50.0, 5000); R = rng.gamma(2.0, 50.0, 5000)
    for neg_q in (0.3, 0.5, 0.7):
        mask, c0, c1 = restore.idx_select(T, R, neg_q)
        np.testing.assert_array_equal(mask, _old_patch_idx_select(T, R, neg_q))
        assert c0 == float(np.quantile(T, neg_q))
        assert c1 == float(np.quantile(R, neg_q))
        assert mask.dtype == bool and mask.shape == T.shape


def test_idx_select_escape_hatch_is_copositive_floor_50():
    rng = np.random.default_rng(3)
    T = rng.gamma(2.0, 50.0, 2000); R = rng.gamma(2.0, 50.0, 2000)
    mask, f0, f1 = restore.idx_select(T, R, neg_q=0.0)
    assert f0 == 50.0 and f1 == 50.0
    expected = ((T > 50.0) & (R > np.quantile(R, 0.5))) | ((R > 50.0) & (T > np.quantile(T, 0.75)))
    np.testing.assert_array_equal(mask, expected)


def test_idx_select_nonzero_q_guard_lifts_collapsed_ceiling():
    # Item 3: REDSEA zero-clip. >50% of the target clamped to 0 => quantile(T, 0.5) == 0 collapses the
    # negative arm to exact-zero cells only; the nonzero_q guard floors the ceiling at a nonzero quantile.
    rng = np.random.default_rng(5)
    T = np.clip(rng.normal(30, 40, 5000), 0, None); R = rng.gamma(2.0, 50.0, 5000)
    T[:2700] = 0.0
    m0, c0, _ = restore.idx_select(T, R, 0.5, 0.75, 0.5, 0.0)      # default: no guard
    mg, cg0, _ = restore.idx_select(T, R, 0.5, 0.75, 0.5, 0.25)    # guard on
    assert c0 == 0.0                     # ceiling collapsed under heavy zero-clip
    assert cg0 > 0.0                     # guard floors it at the nonzero-cell quantile
    assert mg.sum() >= m0.sum()          # negative arm no longer restricted to exact-zero cells


def test_idx_select_nonzero_q_is_noop_without_zeros():
    # No zero-clip + small nonzero_q => quantile(T, neg_q) already dominates => identical selection.
    rng = np.random.default_rng(6)
    T = rng.gamma(2.0, 50.0, 3000) + 1.0; R = rng.gamma(2.0, 50.0, 3000) + 1.0
    m0, _, _ = restore.idx_select(T, R, 0.5, 0.75, 0.5, 0.0)
    mg, _, _ = restore.idx_select(T, R, 0.5, 0.75, 0.5, 0.10)
    np.testing.assert_array_equal(m0, mg)


# --------------------------------------------------------------------------- #
# neg_stat — the four threshold statistics
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("stat", ["mean3sd", "mad", "pctile", "logmad"])
def test_neg_stat_variants_finite_positive(stat):
    rng = np.random.default_rng(1)
    v = rng.gamma(2.0, 30.0, 2000)
    t = restore.neg_stat(v, stat)
    assert np.isfinite(t) and t > 0


def test_neg_stat_mean3sd_is_mean_plus_k_std_ddof0():
    v = np.array([0.0, 0.0, 0.0, 10.0])
    assert restore.neg_stat(v, "mean3sd", k=3.0) == float(v.mean() + 3.0 * v.std())


def test_neg_stat_degenerate_is_nan():
    assert np.isnan(restore.neg_stat(np.array([5.0]), "mean3sd"))
    assert np.isnan(restore.neg_stat(np.array([np.nan, np.nan]), "mad"))


def test_neg_stat_unknown_raises():
    with pytest.raises(ValueError):
        restore.neg_stat(np.arange(10.0), "nope")


# --------------------------------------------------------------------------- #
# r_otsu — pure Otsu valley
# --------------------------------------------------------------------------- #
def test_r_otsu_valley_between_two_modes():
    pytest.importorskip("skimage")
    rng = np.random.default_rng(2)
    x = np.concatenate([rng.normal(10.0, 1.0, 1000), rng.normal(100.0, 5.0, 1000)])
    thr = r_otsu(x)
    assert 10.0 < thr < 100.0


def test_r_otsu_degenerate_is_nan():
    pytest.importorskip("skimage")
    assert np.isnan(r_otsu(np.array([3.0])))
    assert np.isnan(r_otsu(np.array([3.0, 3.0, 3.0])))


# --------------------------------------------------------------------------- #
# fit_clusters — 2-cluster negative-population selection (KMeans/GMM; SSC needs spams)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("model", ["KMeans", "GMM"])
def test_fit_clusters_selects_reference_broad_negative(model):
    rng = np.random.default_rng(4)
    # reference-positive negative pop (broad on the reference axis) vs a tight target-positive arm
    neg = np.column_stack([rng.normal(5, 2, 400), rng.normal(200, 60, 400)])   # T low, R high+broad
    sig = np.column_stack([rng.normal(200, 20, 400), rng.normal(5, 2, 400)])    # T high, R low
    cloud = np.vstack([neg, sig])
    out_cloud, labels, negi = restore.fit_clusters(cloud, model=model, subsample=10_000, seed=0)
    assert set(np.unique(labels)) <= {0, 1}
    # the chosen negative cluster should be the reference-broad one (higher mean reference)
    assert out_cloud[labels == negi][:, 1].mean() > out_cloud[labels != negi][:, 1].mean()


# --------------------------------------------------------------------------- #
# partner_low_thresh — threshold the target among reference-LOW cells (+ r_otsu relocation)
# --------------------------------------------------------------------------- #
def test_partner_low_thresh_valley_on_reference_low():
    pytest.importorskip("skimage")
    rng = np.random.default_rng(0)
    n = 4000
    # reference: half low (~5), half high (~300); target among reference-LOW is bimodal (background + β)
    R = np.concatenate([rng.normal(5, 2, n), rng.normal(300, 40, n)])
    T = np.concatenate([rng.normal(5, 2, n // 2), rng.normal(400, 60, n // 2),   # reference-low: bg + signal
                        rng.normal(5, 2, n)])                                      # reference-high: target low
    thr = restore.partner_low_thresh(T, R, 0.5)
    assert np.isfinite(thr) and 5 < thr < 400   # valley lands between the background and signal modes


def test_partner_low_thresh_absent_when_unimodal():
    pytest.importorskip("skimage")
    rng = np.random.default_rng(1)
    n = 4000
    R = np.concatenate([rng.normal(5, 2, n), rng.normal(300, 40, n)])
    T = rng.normal(5, 2, 2 * n)                  # target is background everywhere (β-loss analogue)
    assert not np.isfinite(restore.partner_low_thresh(T, R, 0.5))   # no valley -> +inf (absent)


def test_partner_low_thresh_absent_when_too_few():
    rng = np.random.default_rng(2)
    R = rng.normal(300, 10, 500); T = rng.normal(400, 10, 500)
    # partner_low_q=0.1 -> ~50 reference-low cells, below min_n
    assert not np.isfinite(restore.partner_low_thresh(T, R, 0.1, min_n=200))


def test_r_otsu_and_partner_low_reexported_from_diagnostics():
    # relocated into restore.py; diagnostics must re-export the SAME objects (notebook import unchanged)
    from phenocycler import diagnostics
    assert diagnostics.r_otsu is restore.r_otsu
    assert diagnostics.partner_low_thresh is restore.partner_low_thresh
    assert diagnostics.bimodal_thresh is restore.bimodal_thresh


# --------------------------------------------------------------------------- #
# bimodal_thresh — standalone binarization (proliferation Ki67/PCNA)
# --------------------------------------------------------------------------- #
def test_bimodal_thresh_finds_crossover():
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(0)
    v = np.concatenate([rng.normal(5, 2, 3600), rng.normal(300, 40, 400)])   # negative bulk + 10% positive
    thr = restore.bimodal_thresh(v)
    assert np.isfinite(thr) and 5 < thr < 300


def test_bimodal_thresh_absent_when_unimodal():
    pytest.importorskip("sklearn")
    rng = np.random.default_rng(1)
    v = rng.normal(5, 2, 4000)                                               # single mode -> no proliferation
    assert not np.isfinite(restore.bimodal_thresh(v))


def test_logmean_ksigma_upper_tail():
    # graded marker (proliferation): robust upper-tail cutoff isolates a small positive fraction
    rng = np.random.default_rng(0)
    v = np.clip(rng.lognormal(3.0, 0.6, 20000), 0, None)                     # broad graded, right tail
    t3 = restore.logmean_ksigma(v, 3.0); t2 = restore.logmean_ksigma(v, 2.0)
    assert np.isfinite(t3) and t3 > t2                                       # higher k => higher cutoff
    assert 0 < (v >= t3).mean() < 0.05                                       # small proliferating fraction
    assert not np.isfinite(restore.logmean_ksigma(v[:10]))                   # too few cells -> inf
