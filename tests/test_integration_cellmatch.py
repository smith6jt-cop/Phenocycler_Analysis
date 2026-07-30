"""S6b — cell-to-cell matching across serial sections, and the guards that make it honest.

The premise this module rests on was measured, not assumed: with 10-15 um cells and 5 um
sections, ~70% of cells straddle the cut but only ~28% are detectable in both, the shared
population is heavily biased toward large cells (0% of 8 um lymphocytes), and unconstrained
matching tops out at 0.72 precision because cell spacing is close to the matching cap.
Type-constrained matching under a <=2 um local residual reaches 0.95.

Every test here pins one of the guards that keeps that true, because the failure mode is not
a crash — it is a confidently wrong pairing that reads exactly like a right one.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.cellmatch import (
    MATCHABLE_TYPES, fine_type, match_islet_cells, match_typed, refine_local,
)
from phenocycler.integration.contract import CellTable, coerce_cell_table
from phenocycler.integration.transform import Transform, from_rigid


def _islet(n=90, centre=(500.0, 500.0), radius=60.0, seed=0):
    """A synthetic islet with a realistic Beta/Alpha/Delta mix and known cell identities."""
    rng = np.random.default_rng(seed)
    ang = rng.uniform(0, 2 * np.pi, n)
    rad = radius * np.sqrt(rng.uniform(0, 1, n))
    types = rng.choice(["Beta", "Alpha", "Delta"], n, p=[0.6, 0.3, 0.1])
    return pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(n)],
        "x_um": centre[0] + rad * np.cos(ang),
        "y_um": centre[1] + rad * np.sin(ang),
        "lineage_native": ["Endocrine"] * n,
        "lineage_common": ["Endocrine"] * n,
        "endocrine_subtype": types,
        "structure_id": ["islet_1"] * n,
    })


def _table(df, modality, tissue="pancreas"):
    out = coerce_cell_table(df.copy(), modality=modality, donor_id="T1", roi="panc")
    out["structure_id"] = df["structure_id"].to_numpy()
    if "endocrine_subtype" in df.columns:
        out["endocrine_subtype"] = df["endocrine_subtype"].to_numpy()
    return CellTable(df=out, modality=modality, donor_id="T1", roi="panc", tissue=tissue)


def _serial_pair(shared_frac=0.30, jitter_um=1.0, seed=1, n=90):
    """Two sections of one islet: `shared_frac` of cells appear in both, the rest are unique.

    Mirrors the geometry — a minority of cells are detectable in both sections, and the two
    fragments of a shared cell have slightly different centroids.
    """
    rng = np.random.default_rng(seed)
    base = _islet(n=n, seed=seed)
    n_shared = int(n * shared_frac)
    shared = base.iloc[:n_shared]

    a = base.copy()
    b = shared.copy()
    b[["x_um", "y_um"]] += rng.normal(0, jitter_um, (len(b), 2))
    extra = _islet(n=n - n_shared, seed=seed + 100)          # cells unique to section B
    extra["cell_id"] = [f"b{i}" for i in range(len(extra))]
    b = pd.concat([b, extra], ignore_index=True)
    return a, b, set(shared["cell_id"])


# --------------------------------------------------------------------------- #
# Scope
# --------------------------------------------------------------------------- #

def test_only_endocrine_types_are_matchable():
    """Immune cells are geometrically unrecoverable — an 8 um cell cannot put enough nucleus
    on both sides of a 5 um cut. The module must refuse them rather than report a low rate."""
    assert MATCHABLE_TYPES == ("Beta", "Alpha", "Delta")
    for immune in ("T", "B", "Myeloid", "Macrophage", "Neutrophil"):
        assert immune not in MATCHABLE_TYPES


def test_fine_type_prefers_the_protein_derived_call():
    """`endocrine_subtype` is the hormone argmax — protein on both sides once post-Xenium IF
    exists. Preferring it means the two sides are typed by the same kind of evidence."""
    df = pd.DataFrame({"endocrine_subtype": ["Beta", "", "Alpha"],
                       "celltype_fine": ["Alpha", "Delta", "Beta"]})
    assert fine_type(df).tolist() == ["Beta", "Delta", "Alpha"]

    # Anything outside the matchable set is blanked, not passed through.
    df = pd.DataFrame({"celltype_fine": ["Beta", "Macrophage", "Acinar"]})
    assert fine_type(df).tolist() == ["Beta", "", ""]


# --------------------------------------------------------------------------- #
# Local refinement
# --------------------------------------------------------------------------- #

def test_refinement_recovers_a_known_local_offset():
    a, b, _shared = _serial_pair(shared_frac=1.0, jitter_um=0.2)
    offset = np.array([3.0, -2.0])
    moving = a.loc[:, ["x_um", "y_um"]].to_numpy() + offset
    fixed = b.loc[:, ["x_um", "y_um"]].to_numpy()

    _tf, residual, n_in = refine_local(moving, fixed)
    assert n_in > 10
    assert residual < 1.0, f"residual {residual} um — refinement did not recover the offset"


@pytest.mark.parametrize("offset_um,expect", [
    (5.0, "ok"),          # a real local correction — converges tightly
    (30.0, "residual"),   # converges, but to a 4-5 um residual: the gate must catch it
    (120.0, "shift"),     # a fit this large is matching a different islet
    (400.0, "diverge"),   # beyond ICP's capture range entirely
])
def test_every_misalignment_scale_fails_safe(offset_um, expect):
    """Three independent guards, and every failure mode lands in one of them.

    This is the map of which guard fires when — worth pinning, because the dangerous outcome
    is not an error but a fit that looks plausible and is wrong. Note the 30 um case: ICP
    happily converges, and only the residual gate stands between that and a claimed
    measurement at ~0.40 precision.
    """
    a, b, _ = _serial_pair(shared_frac=1.0)
    moving = a.loc[:, ["x_um", "y_um"]].to_numpy() + np.array([offset_um, 0.0])
    fixed = b.loc[:, ["x_um", "y_um"]].to_numpy()

    _tf, residual, n_in = refine_local(moving, fixed, max_shift_um=25.0, max_dist_um=500.0)
    if expect == "ok":
        assert n_in > 10 and residual < 2.0
    elif expect == "residual":
        assert n_in > 0 and residual > 2.0, "must not slip under the residual gate"
    elif expect == "shift":
        assert n_in == -1, "an implausible local shift must be rejected"
    else:
        assert not np.isfinite(residual), "a non-convergent fit must not report a residual"


# --------------------------------------------------------------------------- #
# Type constraint — the thing that makes this a measurement
# --------------------------------------------------------------------------- #

def test_type_constraint_suppresses_the_chance_match_rate():
    """Restricting to one type cuts candidate density, so the null collapses while true
    matches are untouched. That is the entire difference between 0.72 and 0.95 precision."""
    a, b, shared = _serial_pair(shared_frac=0.3, jitter_um=1.0)
    pairs, stats = match_typed(a, b, max_dist_um=3.0, n_perm=25)

    assert stats["n_matched"] > 0
    assert stats["signal_over_null"] > 3.0, (
        f"signal/null {stats['signal_over_null']:.1f} — matching is not above chance")
    # Pairs must be within-type by construction.
    a_type = dict(zip(a["cell_id"], a["endocrine_subtype"]))
    b_type = dict(zip(b["cell_id"], b["endocrine_subtype"]))
    for _, r in pairs.iterrows():
        assert a_type[r["cell_id_moving"]] == b_type[r["cell_id_fixed"]] == r["celltype"]


def test_matched_pairs_are_mostly_the_genuinely_shared_cells():
    """The precision claim, checked against ground truth rather than asserted."""
    a, b, shared = _serial_pair(shared_frac=0.3, jitter_um=1.0)
    pairs, _stats = match_typed(a, b, max_dist_um=3.0, n_perm=10)
    assert len(pairs) >= 10

    correct = sum(1 for _, r in pairs.iterrows()
                  if r["cell_id_moving"] == r["cell_id_fixed"] and r["cell_id_moving"] in shared)
    precision = correct / len(pairs)
    assert precision > 0.7, f"precision {precision:.2f} — below the simulated ceiling"


# --------------------------------------------------------------------------- #
# The islet-level guards
# --------------------------------------------------------------------------- #

@pytest.fixture(scope="module")
def cfg():
    return load_config()


def test_islet_is_skipped_when_the_local_residual_is_too_large():
    """At 5 um residual precision falls to 0.40 — worse than a coin flip on a claimed
    measurement. The gate must refuse rather than emit pairs with a caveat."""
    a, b, _ = _serial_pair(shared_frac=0.3)
    moving = _table(a, "xenium")
    fixed = _table(b, "phenocycler")
    # 30 um out: ICP converges, but to a ~4.5 um residual — the regime where precision is
    # 0.40. Nothing errors; only the gate stands between this and a claimed measurement.
    tf = from_rigid(0.0, (30.0, 0.0))

    pairs, stats = match_islet_cells(moving, fixed, "islet_1", "islet_1", tf,
                                     max_residual_um=2.0)
    assert pairs.empty
    assert stats["skipped"] is True
    assert "residual" in stats["reason"]


def test_islet_is_skipped_when_matches_are_not_above_chance():
    a, b, _ = _serial_pair(shared_frac=0.3)
    moving, fixed = _table(a, "xenium"), _table(b, "phenocycler")
    pairs, stats = match_islet_cells(moving, fixed, "islet_1", "islet_1", Transform(),
                                     min_signal_over_null=1e6)
    assert pairs.empty
    assert stats["skipped"] is True
    assert "chance" in stats["reason"]


def test_islet_is_skipped_when_there_are_too_few_cells():
    small = _islet(n=3)
    t = _table(small, "xenium")
    pairs, stats = match_islet_cells(t, t, "islet_1", "islet_1", Transform())
    assert pairs.empty and stats["skipped"] is True
    assert "too few cells" in stats["reason"]


def test_end_to_end_islet_match_reports_residual_and_null():
    """The happy path: a well-aligned islet yields pairs, a small residual, and a null the
    match count clearly beats."""
    a, b, shared = _serial_pair(shared_frac=0.35, jitter_um=0.8)
    moving, fixed = _table(a, "xenium"), _table(b, "phenocycler")
    tf = from_rigid(0.0, (2.0, -1.5))

    pairs, stats = match_islet_cells(moving, fixed, "islet_1", "islet_1", tf)
    assert stats["skipped"] is False, stats.get("reason")
    assert stats["residual_um"] < 2.0
    assert stats["signal_over_null"] > 3.0
    assert set(pairs.columns) >= {"cell_id_moving", "cell_id_fixed", "celltype", "distance_um"}

    correct = sum(1 for _, r in pairs.iterrows() if r["cell_id_moving"] == r["cell_id_fixed"])
    assert correct / len(pairs) > 0.7


def test_refinement_is_type_blind():
    """Refining on the typed cells it then matches would manufacture the agreement it claims.

    Pinned by construction: `refine_local` takes plain coordinate arrays and has no access to
    a type column at all, so the coupling cannot be reintroduced without changing its
    signature.
    """
    import inspect

    params = set(inspect.signature(refine_local).parameters)
    assert not ({"types", "celltype", "type_col"} & params), (
        "refine_local must not see cell types — it runs before, and independent of, matching")
