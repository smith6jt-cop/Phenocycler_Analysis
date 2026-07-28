"""Structure extraction: hormone-seeded islet calling, symmetric across modalities."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.contract import CellTable, coerce_cell_table
from phenocycler.integration.rasterize import (GridSpec, grid_for, rasterize_points,
                                               rasterize_structures, tissue_mask)
from phenocycler.integration.structures import (CT_INSULITIS_THRESH, CT_PERI_THRESH,
                                                StructureParams, call_islets, call_landmarks,
                                                hormone_seed_mask, insulitis_grade,
                                                size_class)


def _field(modality, n_islets=12, seed=7, islet_sizes=None, insulitic=(),
           n_exocrine=3000, extent=4000.0):
    """A synthetic section with a known number of islets of known size."""
    r = np.random.default_rng(seed)
    centres = r.uniform(400, extent - 400, (n_islets, 2))
    sizes = islet_sizes if islet_sizes is not None else r.integers(25, 160, n_islets)
    rows = []
    for k, (c, n) in enumerate(zip(centres, sizes)):
        sub = r.choice(["Beta", "Alpha", "Delta"], n, p=[0.6, 0.3, 0.1])
        for (x, y), s in zip(c + r.normal(0, 22, (n, 2)), sub):
            rows.append((x, y, "Endocrine", s))
        if k in insulitic:
            for x, y in c + r.normal(0, 25, (40, 2)):
                rows.append((x, y, "T_NK", "T"))
    for _ in range(n_exocrine):
        rows.append((r.uniform(0, extent), r.uniform(0, extent), "Exocrine", "Acinar"))
    for _ in range(300):
        rows.append((r.uniform(0, extent), r.uniform(0, extent), "Vascular", "Endothelial"))

    a = np.array(rows, dtype=object)
    df = pd.DataFrame({
        "cell_id": [f"c{i}" for i in range(len(a))],
        "x_um": a[:, 0].astype(float), "y_um": a[:, 1].astype(float),
        "lineage_native": a[:, 2], "lineage_common": a[:, 2], "celltype_fine": a[:, 3],
    })
    if modality == "phenocycler":
        df["endocrine_subtype"] = np.where(df["lineage_common"] == "Endocrine",
                                           df["celltype_fine"], "")
        for m, s in (("INS", "Beta"), ("GCG", "Alpha"), ("SST", "Delta")):
            df[f"feat__{m}"] = np.where(df["celltype_fine"] == s, 12.0, 0.2)
    df = coerce_cell_table(df, modality=modality, donor_id="T1", roi="panc")
    return CellTable(df=df, modality=modality, donor_id="T1", roi="panc"), centres, sizes


@pytest.fixture(scope="module")
def cfg():
    return load_config()


def test_planted_islets_are_recovered(cfg):
    table, _c, sizes = _field("phenocycler", n_islets=12)
    islets, ids = call_islets(table, StructureParams.from_config(cfg),
                              hormone_min_norm=cfg.hormone_min_norm)
    assert len(islets) == 12
    n_endo = int((table.df["lineage_common"] == "Endocrine").sum())
    assert int((ids != "").sum()) == n_endo
    assert sorted(islets["n_cells"]) == sorted(int(s) for s in sizes)


def test_both_modalities_use_the_same_algorithm(cfg):
    """Symmetry is the point: a difference between modalities must not be a method artifact."""
    p, _, _ = _field("phenocycler", n_islets=10, seed=3)
    x, _, _ = _field("xenium", n_islets=10, seed=3)
    params = StructureParams.from_config(cfg)
    pi, _ = call_islets(p, params, hormone_min_norm=cfg.hormone_min_norm)
    xi, _ = call_islets(x, params, hormone_min_norm=cfg.hormone_min_norm)
    assert len(pi) == len(xi) == 10
    assert sorted(pi["n_cells"]) == sorted(xi["n_cells"])


def test_hormone_seeding_beats_label_seeding(cfg):
    """The empirical finding this design is built on.

    Seeding on the broad endocrine LABEL merges islets through scattered label-positive
    cells; seeding on hormone-positive cells does not. Upstream this was 2-9 islets vs 112+.
    """
    r = np.random.default_rng(5)
    table, _c, _s = _field("phenocycler", n_islets=10, seed=5)
    # Scatter generic endocrine-labelled but hormone-negative cells across the section, as
    # CD99-bright and generic endocrine calls really do. Density is chosen so a disc of
    # radius eps=50 um holds more than min_samples=10 of them, which is the regime real
    # tissue is in — below that threshold DBSCAN cannot percolate and the effect does not
    # appear at all.
    n_scatter = 25_000
    extra = pd.DataFrame({
        "cell_id": [f"g{i}" for i in range(n_scatter)],
        "x_um": r.uniform(0, 4000, n_scatter), "y_um": r.uniform(0, 4000, n_scatter),
        "lineage_native": "Endocrine", "lineage_common": "Endocrine",
        "celltype_fine": "Endocrine", "endocrine_subtype": "",
        "feat__INS": 0.3, "feat__GCG": 0.3, "feat__SST": 0.3,
    })
    df = pd.concat([table.df, coerce_cell_table(extra, modality="phenocycler",
                                                donor_id="T1", roi="panc")],
                   ignore_index=True)
    mixed = CellTable(df=df, modality="phenocycler", donor_id="T1", roi="panc")

    params = StructureParams.from_config(cfg)
    hormone_seeded, _ = call_islets(mixed, params, hormone_min_norm=cfg.hormone_min_norm)

    from phenocycler.integration.structures import dbscan_labels
    label_mask = (mixed.df["lineage_common"] == "Endocrine").to_numpy()
    lab = dbscan_labels(mixed.xy()[label_mask], params.eps_um, params.min_samples)
    label_seeded = int(lab.max()) + 1 if (lab >= 0).any() else 0

    assert len(hormone_seeded) == 10
    assert label_seeded < len(hormone_seeded), (
        f"label seeding found {label_seeded} islets, hormone seeding {len(hormone_seeded)} — "
        "the scattered generic-endocrine cells should have merged the label-seeded clusters")


def test_hormone_seed_mask_per_modality(cfg):
    p, _, _ = _field("phenocycler", n_islets=5, seed=11)
    x, _, _ = _field("xenium", n_islets=5, seed=11)
    assert hormone_seed_mask(p, hormone_min_norm=5.0).sum() > 0
    assert hormone_seed_mask(x).sum() > 0
    # Xenium seeds on the specific fine subtypes, not the generic Endocrine label.
    assert hormone_seed_mask(x).sum() == int(
        x.df["celltype_fine"].isin(["Beta", "Alpha", "Delta"]).sum())


def test_insulitis_grade_and_size_class():
    assert insulitis_grade(0.0) == "no_insulitis"
    assert insulitis_grade(CT_PERI_THRESH) == "peri_insulitis"
    assert insulitis_grade(CT_INSULITIS_THRESH) == "insulitis"
    assert size_class(10) == "tiny"
    assert size_class(75) == "med"
    assert size_class(5000) == "xl"


def test_infiltrated_islet_is_graded_insulitic(cfg):
    table, _c, _s = _field("phenocycler", n_islets=6, seed=13, insulitic=(0,))
    islets, _ = call_islets(table, StructureParams.from_config(cfg),
                            hormone_min_norm=cfg.hormone_min_norm)
    assert (islets["insulitis_grade"] == "insulitis").sum() >= 1
    assert islets["immune_per_100_endo"].max() > CT_INSULITIS_THRESH


def test_per_islet_composition_matches_the_planted_ratio(cfg):
    table, _c, _s = _field("phenocycler", n_islets=8, seed=17)
    islets, _ = call_islets(table, StructureParams.from_config(cfg),
                            hormone_min_norm=cfg.hormone_min_norm)
    # Planted at 60/30/10; sampling noise per islet, so check the cohort mean.
    assert islets["frac_beta"].mean() == pytest.approx(0.6, abs=0.08)
    assert islets["frac_alpha"].mean() == pytest.approx(0.3, abs=0.08)


def test_empty_and_degenerate_inputs(cfg):
    df = coerce_cell_table(
        pd.DataFrame({"cell_id": ["a"], "x_um": [0.0], "y_um": [0.0],
                      "lineage_native": ["Exocrine"], "lineage_common": ["Exocrine"]}),
        modality="phenocycler", donor_id="T", roi="panc")
    t = CellTable(df=df, modality="phenocycler", donor_id="T", roi="panc")
    islets, ids = call_islets(t, StructureParams.from_config(cfg))
    assert len(islets) == 0
    assert (ids == "").all()
    assert len(call_landmarks(t, StructureParams.from_config(cfg), "duct")) == 0


# --------------------------------------------------------------------------- #
# Rasterisation
# --------------------------------------------------------------------------- #

def test_grid_roundtrip_is_exact():
    rng = np.random.default_rng(2)
    xy = rng.uniform(0, 4000, (200, 2))
    grid = grid_for(xy, 2.5)
    assert np.abs(grid.to_microns(grid.to_pixels(xy)) - xy).max() < 1e-9


def test_grid_covers_both_point_sets():
    a = np.array([[0.0, 0.0], [100.0, 100.0]])
    b = np.array([[-500.0, 900.0]])
    grid = grid_for(a, 10.0, other=b, margin_um=50)
    px = grid.to_pixels(np.vstack([a, b]))
    assert (px >= 0).all() and (px[:, 0] <= grid.width).all() and (px[:, 1] <= grid.height).all()


def test_rasterize_conserves_point_count():
    rng = np.random.default_rng(4)
    xy = rng.uniform(100, 900, (500, 2))
    grid = grid_for(xy, 5.0, margin_um=50)
    assert rasterize_points(xy, grid).sum() == pytest.approx(500.0)


def test_tissue_mask_is_contiguous_and_bounded():
    rng = np.random.default_rng(6)
    xy = rng.uniform(1000, 3000, (4000, 2))
    grid = grid_for(xy, 10.0, margin_um=300)
    mask = tissue_mask(xy, grid, min_area_um2=1000)
    assert mask.any() and not mask.all()
    # Everything inside the point cloud's bounding box should be covered after closing.
    px = grid.to_pixels(np.array([[1500.0, 1500.0], [2500.0, 2500.0]]))
    for c, r in px:
        assert mask[int(r), int(c)]


def test_rasterize_structures_weights_by_size():
    islets = pd.DataFrame({"x_um": [100.0, 500.0], "y_um": [100.0, 500.0],
                           "n_cells": [10.0, 90.0]})
    grid = GridSpec(0.0, 0.0, 10.0, 80, 80)
    img = rasterize_structures(islets, grid, weight_by="n_cells")
    assert img.sum() == pytest.approx(100.0)
    assert img.max() == pytest.approx(90.0)
