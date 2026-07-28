"""Islet matching, niches, and the donor-level comparison — including the epi_default fix."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.contract import CellTable, coerce_cell_table
from phenocycler.integration.donor import COMPARABLE, compare_donor, composition
from phenocycler.integration.grid import (NICHE_LINEAGES, call_niches, composition_vectors,
                                          compare_grids, grid_composition, hex_bin,
                                          majority_smooth, niche_agreement, niche_profiles)
from phenocycler.integration.match import (concordance, match_structures, permuted_null)


def _structures(n, seed, jitter=0.0, drop=0.0, transform=None):
    r = np.random.default_rng(seed)
    xy = r.uniform(300, 3700, (n, 2))
    df = pd.DataFrame({
        "structure_id": [f"s{i}" for i in range(n)],
        "x_um": xy[:, 0], "y_um": xy[:, 1],
        "n_cells": r.integers(30, 180, n).astype(float),
        "area_um2": r.uniform(4e3, 1.1e4, n),
        "eccentricity": r.uniform(0.1, 0.6, n),
        "frac_beta": r.uniform(0.2, 0.8, n),
        "immune_per_100_endo": r.uniform(0, 20, n),
        "insulitis_grade": r.choice(["no_insulitis", "peri_insulitis", "insulitis"], n),
    })
    if transform is not None:
        moved = transform.apply(df.loc[:, ["x_um", "y_um"]].to_numpy())
        df["x_um"], df["y_um"] = moved[:, 0], moved[:, 1]
    if jitter:
        df["x_um"] += r.normal(0, jitter, len(df))
        df["y_um"] += r.normal(0, jitter, len(df))
    if drop:
        df = df[r.random(len(df)) >= drop].reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
# Matching
# --------------------------------------------------------------------------- #

def test_identity_matching_recovers_every_pair():
    fixed = _structures(25, 1)
    pairs = match_structures(fixed.copy(), fixed, max_dist_um=200, area_ratio=2.5)
    assert len(pairs) == 25
    assert (pairs["distance_um"] < 1e-9).all()
    assert pairs["mutual_nn"].all()


def test_matching_under_jitter_and_dropout():
    """Islets appear and vanish between sections; a match fraction of 1.0 would be suspicious."""
    fixed = _structures(30, 2)
    moving = fixed.copy()
    r = np.random.default_rng(3)
    moving["x_um"] += r.normal(0, 20, len(moving))
    moving["y_um"] += r.normal(0, 20, len(moving))
    moving = moving[r.random(len(moving)) >= 0.2].reset_index(drop=True)

    pairs = match_structures(moving, fixed, max_dist_um=200, area_ratio=2.5)
    assert 0.6 * len(moving) <= len(pairs) <= len(moving)
    assert pairs["distance_um"].median() < 60


def test_assignment_is_one_to_one():
    """Hungarian, not greedy: no fixed islet may be claimed twice."""
    fixed = _structures(20, 4)
    moving = pd.concat([fixed, fixed], ignore_index=True)     # every islet duplicated
    moving["structure_id"] = [f"m{i}" for i in range(len(moving))]
    pairs = match_structures(moving, fixed, max_dist_um=200, area_ratio=2.5)
    assert pairs["fixed_id"].is_unique
    assert pairs["moving_id"].is_unique


def test_distance_gate_is_enforced():
    """No returned pair may exceed the cap, and true partners beyond it must be rejected.

    Asserting zero matches would be wrong: after a 500 um shift some islets land near a
    *different* islet by chance, and those spurious pairs are exactly what the permuted null
    exists to expose. The invariant that must hold is the gate itself.
    """
    fixed = _structures(15, 5)
    moving = fixed.copy()
    moving["x_um"] += 500.0                                    # far beyond the cap

    pairs = match_structures(moving, fixed, max_dist_um=200)
    assert (pairs["distance_um"] <= 200).all()
    # None of the true correspondences (same structure_id) survive the shift.
    assert not (pairs["moving_id"] == pairs["fixed_id"]).any()


def test_area_gate_is_enforced():
    fixed = _structures(15, 6)
    moving = fixed.copy()
    moving["area_um2"] *= 10.0                                 # 10x area mismatch
    moving["n_cells"] *= 10.0
    assert len(match_structures(moving, fixed, max_dist_um=200, area_ratio=2.5)) == 0


def test_permuted_null_separates_signal_from_chance():
    """A dense islet field always produces *some* matches within 200 um by chance."""
    fixed = _structures(30, 7)
    real = _structures(30, 7, jitter=15.0)
    stats = permuted_null(real, fixed, max_dist_um=200, area_ratio=2.5, n_perm=60)
    assert stats["n_real"] > stats["null_p95"]
    assert stats["z"] > 3.0

    scrambled = _structures(30, 999)
    null_stats = permuted_null(scrambled, fixed, max_dist_um=200, area_ratio=2.5, n_perm=60)
    assert null_stats["z"] < stats["z"]


def test_concordance_reports_no_correlation_when_there_is_none():
    """The metric must not manufacture agreement from independent draws."""
    a = _structures(40, 11)
    b = _structures(40, 12)
    pairs = pd.DataFrame({
        "xenium_frac_beta": a["frac_beta"].to_numpy(),
        "phenocycler_frac_beta": b["frac_beta"].to_numpy(),
        "xenium_insulitis_grade": a["insulitis_grade"].to_numpy(),
        "phenocycler_insulitis_grade": b["insulitis_grade"].to_numpy(),
    })
    out = concordance(pairs, "xenium", "phenocycler")
    assert abs(out["spearman_frac_beta"]) < 0.5


def test_concordance_detects_real_agreement():
    a = _structures(40, 13)
    pairs = pd.DataFrame({
        "xenium_frac_beta": a["frac_beta"].to_numpy(),
        "phenocycler_frac_beta": a["frac_beta"].to_numpy() + 0.02,
        "xenium_insulitis_grade": a["insulitis_grade"].to_numpy(),
        "phenocycler_insulitis_grade": a["insulitis_grade"].to_numpy(),
    })
    out = concordance(pairs, "xenium", "phenocycler")
    assert out["spearman_frac_beta"] > 0.95
    assert out["insulitis_agreement"] == pytest.approx(1.0)


def test_empty_inputs_return_empty():
    assert len(match_structures(pd.DataFrame(), _structures(5, 1))) == 0
    assert len(match_structures(_structures(5, 1), pd.DataFrame())) == 0


# --------------------------------------------------------------------------- #
# Niches / grid
# --------------------------------------------------------------------------- #

def _table(modality, seed=1, n=1500, extent=3000.0):
    r = np.random.default_rng(seed)
    lin = r.choice(["Endocrine", "Exocrine", "Vascular", "T_NK"], n, p=[.25, .55, .12, .08])
    df = pd.DataFrame({
        "cell_id": [f"{modality[0]}{i}" for i in range(n)],
        "x_um": r.uniform(0, extent, n), "y_um": r.uniform(0, extent, n),
        "lineage_native": lin, "lineage_common": lin, "celltype_fine": lin,
    })
    df = coerce_cell_table(df, modality=modality, donor_id="T", roi="panc")
    return CellTable(df=df, modality=modality, donor_id="T", roi="panc")


def test_composition_vectors_sum_to_one():
    t = _table("phenocycler")
    comp = composition_vectors(t, k=25)
    assert comp.shape == (len(t.df), len(NICHE_LINEAGES))
    assert np.allclose(comp.sum(axis=1), 1.0, atol=1e-5)


def test_composition_is_invariant_to_row_order():
    t = _table("phenocycler", seed=9)
    a = composition_vectors(t, k=20)
    shuffled = t.df.sample(frac=1.0, random_state=0).reset_index(drop=True)
    b = composition_vectors(CellTable(df=shuffled, modality="phenocycler",
                                      donor_id="T", roi="panc"), k=20)
    assert np.allclose(np.sort(a, axis=0), np.sort(b, axis=0), atol=1e-6)


def test_niches_are_fitted_jointly_so_ids_mean_the_same_thing():
    """Clustering each section separately would give two unrelated labellings."""
    p, x = _table("phenocycler", 3), _table("xenium", 4)
    niches = call_niches([p, x], k=25, n_niches=6, smooth_k=10)
    assert set(niches) == {("phenocycler", "T", "panc"), ("xenium", "T", "panc")}
    prof_p = niche_profiles(p, niches[("phenocycler", "T", "panc")])
    prof_x = niche_profiles(x, niches[("xenium", "T", "panc")])
    agree = niche_agreement(prof_p, prof_x)
    # Same generative distribution both sides, so matching niche ids must be similar.
    assert agree["cosine"].median() > 0.8


def test_majority_smoothing_reduces_speckle():
    r = np.random.default_rng(2)
    xy = r.uniform(0, 1000, (600, 2))
    labels = (xy[:, 0] > 500).astype(int)
    noisy = labels.copy()
    flip = r.random(len(labels)) < 0.2
    noisy[flip] = 1 - noisy[flip]
    smoothed = majority_smooth(noisy, xy, k=15, passes=2)
    assert (smoothed == labels).mean() > (noisy == labels).mean()


def test_hex_binning_is_deterministic_and_covers_all_points():
    r = np.random.default_rng(8)
    xy = r.uniform(0, 2000, (400, 2))
    a, b = hex_bin(xy, 100.0), hex_bin(xy, 100.0)
    assert np.array_equal(a, b)
    assert len(np.unique(a, axis=0)) > 1


def test_grid_composition_and_comparison():
    p, x = _table("phenocycler", 5), _table("xenium", 6)
    pb, xb = grid_composition(p, 200.0), grid_composition(x, 200.0)
    merged, stats = compare_grids(pb, xb, min_cells=5)
    assert stats["n_bins"] > 0
    assert 0.0 <= stats["median_cosine"] <= 1.0


# --------------------------------------------------------------------------- #
# Donor level — the epi_default asymmetry
# --------------------------------------------------------------------------- #

def _table_with_sink(modality, sink_frac, seed=1, n=2000):
    """A section where `sink_frac` of the Exocrine calls came from a fallback rule."""
    r = np.random.default_rng(seed)
    lin = r.choice(["Endocrine", "Exocrine", "Vascular", "T_NK"], n, p=[.25, .55, .12, .08])
    ev = np.ones(n, dtype=bool)
    exo = np.flatnonzero(lin == "Exocrine")
    ev[r.choice(exo, int(len(exo) * sink_frac), replace=False)] = False
    df = pd.DataFrame({
        "cell_id": [f"{modality[0]}{i}" for i in range(n)],
        "x_um": r.uniform(0, 3000, n), "y_um": r.uniform(0, 3000, n),
        "lineage_native": lin, "lineage_common": lin, "celltype_fine": lin,
        "evidence_positive": ev,
    })
    df = coerce_cell_table(df, modality=modality, donor_id="T", roi="panc")
    return CellTable(df=df, modality=modality, donor_id="T", roi="panc")


def test_composition_reports_all_three_fractions():
    t = _table_with_sink("phenocycler", 0.4, seed=21)
    comp = composition(t)
    for col in ("frac_all", "frac_evidence_positive", "frac_fallback"):
        assert col in comp.columns
    exo = comp[comp["lineage_common"] == "Exocrine"].iloc[0]
    assert exo["frac_fallback"] == pytest.approx(0.4, abs=0.05)
    # The corrected fraction differs from the raw one — that is the whole point.
    assert exo["frac_evidence_positive"] != pytest.approx(exo["frac_all"], abs=1e-3)


def test_evidence_positive_correction_changes_the_comparison():
    """PhenoCycler's Epithelial is a default sink; Xenium's Exocrine is a positive call.

    Comparing raw fractions measures that methodological difference, not biology.
    """
    p = _table_with_sink("phenocycler", 0.40, seed=31)
    x = _table_with_sink("xenium", 0.05, seed=31)
    merged, stats = compare_donor(p, x)

    assert stats["frac_fallback_pheno"] > stats["frac_fallback_xen"]
    exo = merged[merged["lineage_common"] == "Exocrine"].iloc[0]
    assert abs(exo["delta_frac_ev"]) != pytest.approx(abs(exo["delta_frac_all"]), abs=1e-4)
    assert "spearman_composition_ev" in stats


def test_comparable_lineages_exclude_the_unresolvable_ones():
    assert "Unassigned" not in COMPARABLE
    # Granulocyte has no Xenium broad class, so a difference there is not informative.
    assert "Granulocyte" not in COMPARABLE
    assert {"Endocrine", "Vascular", "T_NK"} <= set(COMPARABLE)
