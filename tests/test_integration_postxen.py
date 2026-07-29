"""S1c — the post-Xenium re-stain: the one place protein and RNA are on the same cell.

The re-stain is an ordinary PhenoCycler qptiff, so it flows through the whole core pipeline
unchanged. What makes it different is the slide underneath: it images the *Xenium* section.
That is a same-slide relationship inside a serial-section cohort, and everything here pins
one of the consequences.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler import sections as S
from phenocycler.integration.contract import coerce_cell_table, partition_path, write_cell_table
from phenocycler.integration.postxen import (
    CARRY_COLUMNS, pair_restain, stamp_xenium, surrogate_agreement,
)


# --------------------------------------------------------------------------- #
# The slide dimension
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("key,roi,slide", [
    ("6539", "panc", S.SLIDE_PHENO),
    ("6539pln", "pln_1", S.SLIDE_PHENO),
    ("6539xen", "panc", S.SLIDE_XENIUM),
    ("6539xenpln", "pln_1", S.SLIDE_XENIUM),
])
def test_section_key_carries_which_slide_was_imaged(key, roi, slide):
    """A re-stain shares a (donor, roi) with the PhenoCycler section but is different tissue.

    Without the slide in the key, both would land in one partition and two sections' cells
    would be averaged together with nothing to indicate it had happened.
    """
    sec = S.parse(key)
    assert (sec.donor_id, sec.roi, sec.slide) == ("6539", roi, slide)
    assert sec.on_xenium_slide == (slide == S.SLIDE_XENIUM)


def test_pheno_and_restain_sections_never_share_a_partition():
    panc = S.parse("6539")
    restain = S.parse("6539xen")
    assert panc.donor_id == restain.donor_id and panc.roi == restain.roi
    assert panc.key != restain.key, "same donor and roi, different tissue — keys must differ"


def test_an_unlisted_restain_token_is_refused_not_guessed():
    """If the cohort's re-stains use a token that is not registered, the error names it and
    one line fixes it. Defaulting would silently merge two sections."""
    with pytest.raises(S.SectionError, match="unknown region token"):
        S.parse("6539postxen")


def test_key_round_trips_through_the_slide():
    assert S.key_for("6539", "panc", S.SLIDE_XENIUM) == "6539xen"
    assert S.key_for("6539", "pln_1", S.SLIDE_XENIUM) == "6539xenpln"
    assert S.key_for("6539", "panc", S.SLIDE_PHENO) == "6539"
    for roi in ("panc", "pln_1"):
        for slide in (S.SLIDE_PHENO, S.SLIDE_XENIUM):
            sec = S.parse(S.key_for("6539", roi, slide))
            assert (sec.roi, sec.slide) == (roi, slide)


# --------------------------------------------------------------------------- #
# Same-slide pairing
# --------------------------------------------------------------------------- #

def _cells(n, seed, jitter=0.0, offset=(0.0, 0.0), subtypes=None, extent=400.0):
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, extent, (n, 2))
    if jitter:
        xy = xy + rng.normal(0, jitter, xy.shape)
    return pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(n)],
        "x_um": xy[:, 0] + offset[0],
        "y_um": xy[:, 1] + offset[1],
        "lineage_native": ["Endocrine"] * n,
        "lineage_common": ["Endocrine"] * n,
        "endocrine_subtype": subtypes if subtypes is not None else [""] * n,
    })


def _write(cfg, base, df, modality, donor="6539", roi="panc"):
    out = coerce_cell_table(df.copy(), modality=modality, donor_id=donor, roi=roi)
    for col in ("endocrine_subtype", "celltype_fine"):
        if col in df.columns:
            out[col] = df[col].to_numpy()
    write_cell_table(out, partition_path(base, donor, roi))


@pytest.fixture
def cfg(tmp_path):
    c = load_config()
    c.data_dir = tmp_path
    return c


def test_restain_pairs_to_xenium_at_high_rate(cfg):
    """Same section, so most cells should pair — unlike the ~28% ceiling across a 5 um gap."""
    truth = _cells(300, seed=0)
    _write(cfg, cfg.cells_pheno_xif_dir, truth, "phenocycler")
    # Xenium segments the same cells independently: sub-micron centroid disagreement.
    _write(cfg, cfg.cells_xen_dir, _cells(300, seed=0, jitter=0.6), "xenium")

    pairs, stats = pair_restain(cfg, "6539", "panc")
    assert stats["pair_rate_xenium"] > 0.7, stats
    assert stats["signal_over_null"] > 3.0
    correct = (pairs["restain_cell_id"] == pairs["xen_cell_id"]).mean()
    assert correct > 0.9, f"same-slide matching should be near-exact, got {correct:.2f}"


def test_a_frame_offset_is_removed_rather_than_registered(cfg):
    """Both exports are the same physical section, so a constant offset is a coordinate-frame
    difference, not tissue movement. It is estimated and removed here rather than handed to
    the serial-section registration machinery, which solves a problem this does not have."""
    truth = _cells(300, seed=1)
    _write(cfg, cfg.cells_pheno_xif_dir, truth, "phenocycler")
    _write(cfg, cfg.cells_xen_dir, _cells(300, seed=1, jitter=0.6, offset=(120.0, -80.0)),
           "xenium")

    pairs, stats = pair_restain(cfg, "6539", "panc")
    assert stats["frame_offset_um"] > 100, "the offset should be detected and reported"
    assert stats["pair_rate_xenium"] > 0.7, "and then removed, so pairing still works"
    assert (pairs["restain_cell_id"] == pairs["xen_cell_id"]).mean() > 0.9


def test_unrelated_sections_do_not_pair(cfg):
    """The guard that says 'these are not the same section'."""
    _write(cfg, cfg.cells_pheno_xif_dir, _cells(300, seed=2), "phenocycler")
    _write(cfg, cfg.cells_xen_dir, _cells(300, seed=99), "xenium")

    _pairs, stats = pair_restain(cfg, "6539", "panc")
    assert stats["signal_over_null"] < 2.0, (
        f"unrelated cell fields must not look like a match: {stats['signal_over_null']:.1f}x null")


# --------------------------------------------------------------------------- #
# What the re-stain is actually for
# --------------------------------------------------------------------------- #

def test_hormone_subtype_is_carried_onto_the_xenium_table(cfg):
    """This is the point: Xenium endocrine identity stops being a surrogate inference.

    Once stamped, S6b types both sides of the serial pair by hormone protein rather than
    protein on one side and surrogate RNA on the other.
    """
    rng = np.random.default_rng(3)
    subtypes = rng.choice(["Beta", "Alpha", "Delta"], 200, p=[0.6, 0.3, 0.1])
    _write(cfg, cfg.cells_pheno_xif_dir, _cells(200, seed=3, subtypes=list(subtypes)),
           "phenocycler")
    _write(cfg, cfg.cells_xen_dir, _cells(200, seed=3, jitter=0.5), "xenium")

    pairs, _stats = pair_restain(cfg, "6539", "panc")
    assert "endocrine_subtype" in pairs.columns
    assert set(CARRY_COLUMNS) <= set(pairs.columns)

    n = stamp_xenium(cfg, "6539", "panc", pairs)
    assert n > 150

    from phenocycler.integration.contract import read_cell_table
    xen = read_cell_table(cfg.cells_xen_dir, "6539", "panc", modality="xenium")
    stamped = xen.df["endocrine_subtype"].astype(str)
    assert (stamped != "").sum() == n
    assert set(stamped[stamped != ""]) <= {"Beta", "Alpha", "Delta"}


def test_surrogate_agreement_scores_the_rna_inference(cfg):
    """The number the re-stain exists to produce: how good IS the surrogate endocrine call?

    Built so a deliberate disagreement shows up rather than being smoothed away — the RNA
    side here calls every cell Beta, so agreement must land near the true Beta fraction.
    """
    rng = np.random.default_rng(4)
    subtypes = list(rng.choice(["Beta", "Alpha"], 200, p=[0.7, 0.3]))
    _write(cfg, cfg.cells_pheno_xif_dir, _cells(200, seed=4, subtypes=subtypes), "phenocycler")

    xen_df = _cells(200, seed=4, jitter=0.5)
    xen_df["celltype_fine"] = ["Beta"] * 200
    _write(cfg, cfg.cells_xen_dir, xen_df, "xenium")

    pairs, _ = pair_restain(cfg, "6539", "panc")
    out = surrogate_agreement(pairs, cfg, "6539", "panc")
    assert out["n_scored"] > 100
    assert 0.6 < out["surrogate_agreement"] < 0.8, out
    assert out["recall_beta"] == pytest.approx(1.0)     # RNA says Beta for everything
    assert out["recall_alpha"] == pytest.approx(0.0)    # so every true Alpha is missed


def test_stamping_fills_but_never_blanks(cfg):
    """A Xenium cell with no re-stain partner keeps whatever call it already had."""
    _write(cfg, cfg.cells_pheno_xif_dir,
           _cells(50, seed=5, subtypes=["Beta"] * 50), "phenocycler")
    xen_df = _cells(50, seed=5, jitter=0.5)
    xen_df["endocrine_subtype"] = ["Alpha"] * 50
    # Move half the Xenium cells far away so they cannot pair.
    xen_df.loc[25:, "x_um"] += 5000.0
    _write(cfg, cfg.cells_xen_dir, xen_df, "xenium")

    pairs, _ = pair_restain(cfg, "6539", "panc")
    stamp_xenium(cfg, "6539", "panc", pairs)

    from phenocycler.integration.contract import read_cell_table
    xen = read_cell_table(cfg.cells_xen_dir, "6539", "panc", modality="xenium")
    assert (xen.df["endocrine_subtype"].astype(str) == "").sum() == 0, "nothing may be blanked"


# --------------------------------------------------------------------------- #
# Mode
# --------------------------------------------------------------------------- #

def test_postxen_runs_in_both_modes():
    """The correspondence is a property of the section pair, not the cohort. Gating it on
    `[integration] mode` would answer the wrong question."""
    from phenocycler.integration.pipeline import ORDER_SAME_SLIDE, ORDER_SEQUENTIAL

    assert "postxen" in ORDER_SEQUENTIAL
    assert "postxen" in ORDER_SAME_SLIDE
    # And it must not import the cohort-mode guard.
    import inspect

    from phenocycler.integration import postxen

    assert "require_same_slide" not in inspect.getsource(postxen)


def test_postxen_precedes_the_stages_that_consume_its_output():
    """Stamping must happen before S6b types cells, or S6b falls back to surrogate RNA."""
    from phenocycler.integration.pipeline import ORDER_SEQUENTIAL

    order = ORDER_SEQUENTIAL
    assert order.index("postxen") < order.index("cellmatch")
    assert order.index("postxen") < order.index("structures")
