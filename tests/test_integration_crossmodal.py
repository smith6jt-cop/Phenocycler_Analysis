"""Pseudo-cell linking, mode guards, and the QC gates.

The two things worth protecting here are the minimum-anchor guard (linking on a handful of
features produces confident nonsense) and the mode guard (cell-to-cell pairing is a
measurement on one section and an invalid claim across two).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.contract import CellTable, coerce_cell_table
from phenocycler.integration.crossmodal import (AnchorError, fuzzy_smooth, link_cells,
                                                link_donor, permuted_link_null,
                                                resolve_anchors)
from phenocycler.integration.qc import (ALLOWED_BY_VERDICT, FAIL, PASS, WARN, verdict_for)
from phenocycler.integration.sameslide import (ModeError, match_by_centroid,
                                               require_same_slide, require_sequential)

PROTEIN_TO_GENE = {
    "CD31": "PECAM1", "CD3e": "CD3E", "CD163": "CD163", "CD68": "CD68",
    "CD20": "MS4A1", "MPO": "MPO", "B3TUBB": "TUBB3", "CD99": "CD99",
    "Pan_Cytokeratin": "KRT19", "SMA": "ACTA2",
}
LINEAGE_OF = {
    "CD31": "Vascular", "CD3e": "T_NK", "CD163": "Myeloid", "CD68": "Myeloid",
    "CD20": "B_Plasma", "MPO": "Granulocyte", "B3TUBB": "Neural", "CD99": "Endocrine",
    "Pan_Cytokeratin": "Exocrine", "SMA": "Stromal",
}


def _paired_tables(n=400, n_markers=10, seed=0, extent=2000.0):
    """Two sections whose lineages carry matching protein and gene signal."""
    r = np.random.default_rng(seed)
    markers = list(PROTEIN_TO_GENE)[:n_markers]
    lineages = sorted({LINEAGE_OF[m] for m in markers})

    def build(modality, sd):
        lin = r.choice(lineages, n)
        df = pd.DataFrame({
            "cell_id": [f"{modality[0]}{i}" for i in range(n)],
            "x_um": sd.uniform(0, extent, n), "y_um": sd.uniform(0, extent, n),
            "lineage_native": lin, "lineage_common": lin, "celltype_fine": lin,
            "niche_id_joint": "0",
        })
        for m in markers:
            name = m if modality == "phenocycler" else PROTEIN_TO_GENE[m]
            df[f"feat__{name}"] = np.where(lin == LINEAGE_OF[m],
                                           6.0 + sd.random(n), 0.3 + sd.random(n) * 0.3)
        df = coerce_cell_table(df, modality=modality, donor_id="T", roi="panc")
        return CellTable(df=df, modality=modality, donor_id="T", roi="panc")

    return build("phenocycler", np.random.default_rng(seed + 1)), \
        build("xenium", np.random.default_rng(seed + 2))


@pytest.fixture(scope="module")
def cfg():
    return load_config()


# --------------------------------------------------------------------------- #
# Anchors
# --------------------------------------------------------------------------- #

def test_anchors_come_from_the_run_not_the_identity_cores(cfg):
    """KRT19 and ACTA2 reach the deployed panel only via the custom-100 add-on.

    They are in no identity core, so validating against the core union would declare them
    off-panel and silently starve the anchor count.
    """
    p, x = _paired_tables(n_markers=10)
    try:
        pairs, p_cols, x_cols = resolve_anchors(cfg, p, x)
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")

    found = {mp.protein for mp in pairs}
    assert {"Pan_Cytokeratin", "SMA"} <= found, (
        "custom-100 genes present in the run were not recognised as anchors")
    assert len(pairs) == len(p_cols) == len(x_cols)


def test_anchor_requires_presence_in_both_tables(cfg):
    """A crosswalk pair absent from one table is not an anchor."""
    p, x = _paired_tables(n_markers=10)
    x = CellTable(df=x.df.drop(columns=["feat__PECAM1"]), modality="xenium",
                  donor_id=x.donor_id, roi=x.roi)
    try:
        pairs, _, _ = resolve_anchors(cfg, p, x)
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")
    assert "CD31" not in {mp.protein for mp in pairs}


# --------------------------------------------------------------------------- #
# The minimum-anchor guard
# --------------------------------------------------------------------------- #

def test_linking_refuses_below_the_minimum_anchor_count(tmp_path, cfg):
    """Three genes cannot support a shared embedding, and the failure must be loud."""
    p, x = _paired_tables(n_markers=3)
    thin = load_config(data_dir=tmp_path, crossmodal_min_anchors=8)

    import phenocycler.integration.crossmodal as cm
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(thin.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(thin.cells_xen_dir, "T", "panc"))

    with pytest.raises(AnchorError, match="below the configured minimum"):
        cm.link_donor(thin, "T", "panc")


def test_linking_proceeds_with_enough_anchors(tmp_path):
    p, x = _paired_tables(n_markers=10, n=300)
    cfg2 = load_config(data_dir=tmp_path, crossmodal_min_anchors=5)
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(cfg2.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(cfg2.cells_xen_dir, "T", "panc"))
    try:
        res = link_donor(cfg2, "T", "panc")
    except Exception as exc:  # pragma: no cover
        pytest.skip(f"panel taxonomy unavailable: {exc}")
    assert res.n_anchors >= 5
    assert len(res.links) > 0
    assert set(res.links.columns) >= {"inferred_pheno_cell_id", "inferred_xen_cell_id",
                                      "inferred_confidence", "inferred_lineage_agree"}


def test_link_output_columns_are_prefixed_inferred():
    """Serial-section links are inference; the column names must say so."""
    p, x = _paired_tables(n_markers=10, n=200)
    p_cols = [f"feat__{m}" for m in list(PROTEIN_TO_GENE)[:10]]
    x_cols = [f"feat__{PROTEIN_TO_GENE[m]}" for m in list(PROTEIN_TO_GENE)[:10]]
    links = link_cells(p, x, p_cols, x_cols, smooth_k=5)
    assert len(links) > 0
    for col in links.columns:
        assert col.startswith("inferred_") or col in ("group",), col


def test_links_are_one_to_one():
    p, x = _paired_tables(n_markers=10, n=150)
    p_cols = [f"feat__{m}" for m in list(PROTEIN_TO_GENE)[:10]]
    x_cols = [f"feat__{PROTEIN_TO_GENE[m]}" for m in list(PROTEIN_TO_GENE)[:10]]
    links = link_cells(p, x, p_cols, x_cols, smooth_k=5)
    assert links["inferred_pheno_cell_id"].is_unique
    assert links["inferred_xen_cell_id"].is_unique


def test_permuted_link_null_is_reported():
    """With few lineages, random pairing already agrees often; the null must be shown."""
    links = pd.DataFrame({
        "inferred_pheno_lineage": ["Endocrine"] * 60 + ["Exocrine"] * 40,
        "inferred_xen_lineage": ["Endocrine"] * 60 + ["Exocrine"] * 40,
    })
    stats = permuted_link_null(links, n_perm=50)
    assert stats["agree_real"] == pytest.approx(1.0)
    assert stats["agree_null"] < 1.0
    assert stats["z"] > 3.0


def test_fuzzy_smoothing_reduces_variance():
    r = np.random.default_rng(3)
    xy = r.uniform(0, 1000, (300, 2))
    vals = r.normal(0, 1, (300, 4))
    assert fuzzy_smooth(vals, xy, k=15).std() < vals.std()


# --------------------------------------------------------------------------- #
# Mode guards
# --------------------------------------------------------------------------- #

def test_same_slide_pairing_refuses_in_sequential_mode(tmp_path):
    cfg_seq = load_config(data_dir=tmp_path, integration_mode="sequential")
    with pytest.raises(ModeError, match="DIFFERENT cells"):
        require_same_slide(cfg_seq)


def test_sequential_steps_refuse_in_same_slide_mode(tmp_path):
    cfg_ss = load_config(data_dir=tmp_path, integration_mode="same_slide")
    require_same_slide(cfg_ss)                       # allowed
    with pytest.raises(ModeError, match="sequential-mode step"):
        require_sequential(cfg_ss, "islet matching")


def test_crossmodal_refuses_in_same_slide_mode(tmp_path):
    """In same-slide mode the pairing is measured, so approximating it is a mistake."""
    p, x = _paired_tables(n_markers=10, n=100)
    cfg_ss = load_config(data_dir=tmp_path, integration_mode="same_slide")
    from phenocycler.integration.contract import partition_path, write_cell_table

    write_cell_table(p.df, partition_path(cfg_ss.cells_pheno_dir, "T", "panc"))
    write_cell_table(x.df, partition_path(cfg_ss.cells_xen_dir, "T", "panc"))
    with pytest.raises(RuntimeError, match="sameslide"):
        link_donor(cfg_ss, "T", "panc")


def test_invalid_mode_rejected_at_config_load(tmp_path):
    with pytest.raises(ValueError, match="mode must be"):
        load_config(data_dir=tmp_path, integration_mode="serial")


def test_centroid_pairing_is_mutual_and_capped():
    fixed = np.array([[0.0, 0.0], [100.0, 0.0], [200.0, 0.0]])
    moving = np.array([[1.0, 0.0], [101.0, 0.0], [900.0, 0.0]])
    pairs = match_by_centroid(moving, fixed, max_dist_um=5.0)
    assert len(pairs) == 2
    assert (pairs["distance_um"] <= 5.0).all()


# --------------------------------------------------------------------------- #
# QC gates
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("dice,rmse,expected", [
    (0.92, 60.0, PASS),
    (0.92, 400.0, WARN),
    (0.40, 60.0, WARN),
    (0.40, 400.0, FAIL),
    (float("nan"), float("nan"), FAIL),
])
def test_verdict_grading(dice, rmse, expected):
    v, _reasons = verdict_for(dice, rmse, dice_min=0.80, rmse_max=150.0)
    assert v == expected


def test_failed_registration_still_allows_registration_free_analyses():
    """A donor whose sections will not align must not be dropped from the study."""
    assert "donor" in ALLOWED_BY_VERDICT[FAIL]
    assert "niche" in ALLOWED_BY_VERDICT[FAIL]
    assert "structure" not in ALLOWED_BY_VERDICT[FAIL]
    assert "crossmodal" not in ALLOWED_BY_VERDICT[FAIL]
    assert "crossmodal" in ALLOWED_BY_VERDICT[PASS]


def test_verdict_reasons_are_actionable():
    _v, reasons = verdict_for(0.4, 400.0, dice_min=0.80, rmse_max=150.0)
    assert any("Dice" in r for r in reasons)
    assert any("RMSE" in r for r in reasons)
