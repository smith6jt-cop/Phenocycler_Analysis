"""Unit tests for the deterministic 8-class broad-lineage hierarchy (zero Unassigned)."""

from __future__ import annotations

import pandas as pd

from phenocycler import lineage
from phenocycler.config import LINEAGES, EXTRA_MARKERS

# Markers gated in the MAIN restore_gated_redsea pass (the validated 10) vs the EXTRA pass
# (CD99/B3TUBB/MPO, merged by object_id at assignment time).
MAIN_MARKERS = [m for ms in LINEAGES.values() for m in ms if m not in EXTRA_MARKERS]


def _write(tmp_path, rows):
    """Write main-gated + extra-gated + cells partitions for donor TEST from per-cell dicts.

    Each row: {"id": str, <marker>: (pos: bool, norm: float), ...}. A marker in EXTRA_MARKERS is
    routed to the extra-gated frame; all others to the main-gated frame. Unset -> (False, 0.0).
    Returns (main_gated_file, cells_dir, extra_gated_dir).
    """
    ids = [r["id"] for r in rows]

    def _frame(markers):
        d = {"object_id": ids}
        for m in markers:
            d[f"{m}_pos"] = [bool(r.get(m, (False, 0.0))[0]) for r in rows]
            d[f"{m}_norm"] = [float(r.get(m, (False, 0.0))[1]) for r in rows]
        return pd.DataFrame(d)

    gdir = tmp_path / "restore_gated_redsea" / "donor_id=TEST"; gdir.mkdir(parents=True)
    gf = gdir / "data_0.parquet"; _frame(MAIN_MARKERS).to_parquet(gf, index=False)

    xdir = tmp_path / "restore_gated_redsea_extra" / "donor_id=TEST"; xdir.mkdir(parents=True)
    _frame(EXTRA_MARKERS).to_parquet(xdir / "data_0.parquet", index=False)

    cdir = tmp_path / "cells" / "donor_id=TEST"; cdir.mkdir(parents=True)
    pd.DataFrame({"object_id": [str(i) for i in ids], "cell_region": "tissue",
                  "islet_num": ""}).to_parquet(cdir / "data_0.parquet", index=False)
    return str(gf), tmp_path / "cells", tmp_path / "restore_gated_redsea_extra"


def _assign(tmp_path, rows, cd99_bright=3.0):
    gf, cells_dir, extra_dir = _write(tmp_path, rows)
    out, _ = lineage.assign_donor("TEST", gf, cells_dir, extra_dir, cd99_bright=cd99_bright)
    return out


def test_hierarchy_and_zero_unassigned(tmp_path):
    rows = [
        {"id": "e",   "INS": (True, 5.0)},                                  # Endocrine
        {"id": "i",   "CD3e": (True, 4.0)},                                 # Immune
        {"id": "v",   "CD31": (True, 3.0)},                                 # Endothelial
        {"id": "neu", "MPO": (True, 4.0)},                                  # Neutrophil
        {"id": "nrl", "B3TUBB": (True, 4.0)},                               # Neural
        {"id": "cd99", "CD99": (False, 5.0)},                              # Endocrine (bright CD99)
        {"id": "epi", "Pan_Cytokeratin": (False, 0.9),
                      "Vimentin": (False, 0.1), "SMA": (False, 0.1)},       # Epithelial by default
        {"id": "fib", "Vimentin": (True, 2.0),
                      "Pan_Cytokeratin": (False, 0.5), "SMA": (False, 0.3)},# Fibroblast (struct argmax)
        {"id": "mus", "SMA": (True, 3.0),
                      "Pan_Cytokeratin": (False, 0.5), "Vimentin": (False, 0.3)},  # Muscle
        {"id": "e_over_i", "INS": (True, 5.0), "CD3e": (True, 9.0)},        # Endocrine > Immune
        {"id": "i_over_v", "CD3e": (True, 4.0), "CD31": (True, 9.0)},       # Immune > Endothelial
        {"id": "neut_over_nrl", "MPO": (True, 3.0), "B3TUBB": (True, 9.0)}, # Neutrophil > Neural
        {"id": "endo_over_neut", "INS": (True, 5.0), "MPO": (True, 9.0)},   # Endocrine > Neutrophil
        {"id": "imm_over_nrl", "CD3e": (True, 3.0), "B3TUBB": (True, 9.0)}, # Immune > Neural
    ]
    out = _assign(tmp_path, rows)
    call = dict(zip(out["object_id"], out["broad_lineage"]))
    assert call["e"] == "Endocrine"
    assert call["i"] == "Immune"
    assert call["v"] == "Endothelial"
    assert call["neu"] == "Neutrophil"
    assert call["nrl"] == "Neural"
    assert call["cd99"] == "Endocrine"
    assert call["epi"] == "Epithelial"
    assert call["fib"] == "Fibroblast"
    assert call["mus"] == "Muscle"
    assert call["e_over_i"] == "Endocrine"
    assert call["i_over_v"] == "Immune"
    assert call["neut_over_nrl"] == "Neutrophil"
    assert call["endo_over_neut"] == "Endocrine"
    assert call["imm_over_nrl"] == "Immune"

    # core invariant: every cell typed into one of the eight lineages, zero Unassigned
    assert set(out["broad_lineage"]).issubset(set(LINEAGES))
    assert (out["broad_lineage"] == "Unassigned").sum() == 0
    assert len(out) == len(rows)


def test_cd99_bright_gate(tmp_path):
    """Only bright CD99 (_norm >= cd99_bright) is Endocrine; dim CD99 falls to the structural call."""
    rows = [
        {"id": "bright", "CD99": (False, 5.0)},                             # Endocrine
        {"id": "dim", "CD99": (False, 2.0), "Vimentin": (True, 2.0)},       # Fibroblast (CD99 too dim)
    ]
    out = _assign(tmp_path, rows, cd99_bright=3.0).set_index("object_id")
    assert out.loc["bright", "broad_lineage"] == "Endocrine"
    assert out.loc["dim", "broad_lineage"] == "Fibroblast"


def test_epi_default_flag(tmp_path):
    """Cells positive for no structural marker are Epithelial-by-default and flagged."""
    rows = [
        {"id": "weak", "Pan_Cytokeratin": (False, 0.2),
                       "Vimentin": (False, 0.1), "SMA": (False, 0.05)},
        {"id": "strongvim", "Vimentin": (True, 2.0)},
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert bool(out.loc["weak", "epi_default"]) is True
    assert out.loc["weak", "broad_lineage"] == "Epithelial"
    assert bool(out.loc["strongvim", "epi_default"]) is False


def test_long_labels_not_truncated(tmp_path):
    """Guard the <U16 dtype bug: 'Endothelial'/'Neutrophil' must survive intact."""
    rows = [{"id": "v", "CD31": (True, 3.0)}, {"id": "n", "MPO": (True, 3.0)}]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert out.loc["v", "broad_lineage"] == "Endothelial"    # not "Endothelia"
    assert out.loc["n", "broad_lineage"] == "Neutrophil"


def test_score_columns_present(tmp_path):
    """All 8 score_ columns exist (the pipeline uses score_Neural/Neutrophil to detect stale 6-class)."""
    out = _assign(tmp_path, [{"id": "e", "INS": (True, 5.0)}])
    for L in LINEAGES:
        assert f"score_{L}" in out.columns
    assert "score_Neural" in out.columns and "score_Neutrophil" in out.columns
