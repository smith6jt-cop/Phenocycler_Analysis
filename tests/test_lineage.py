"""Unit tests for the deterministic broad-lineage hierarchy (zero Unassigned)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from phenocycler import lineage
from phenocycler.config import LINEAGES

MARKERS = [m for ms in LINEAGES.values() for m in ms]


def _gated_frame(rows):
    """Build a gated parquet-style frame from a list of per-cell dicts.

    Each row dict maps a subset of markers to (pos: bool, norm: float); unset
    markers default to (False, 0.0).
    """
    data = {"object_id": [r["id"] for r in rows]}
    for m in MARKERS:
        data[f"{m}_pos"] = [bool(r.get(m, (False, 0.0))[0]) for r in rows]
        data[f"{m}_norm"] = [float(r.get(m, (False, 0.0))[1]) for r in rows]
    return pd.DataFrame(data)


def _write(tmp_path, gated_df):
    gdir = tmp_path / "restore_gated_redsea" / "donor_id=TEST"
    gdir.mkdir(parents=True)
    gf = gdir / "data_0.parquet"
    gated_df.to_parquet(gf, index=False)
    cdir = tmp_path / "cells" / "donor_id=TEST"
    cdir.mkdir(parents=True)
    pd.DataFrame({"object_id": gated_df["object_id"].astype(str),
                  "cell_region": "tissue", "islet_num": ""}).to_parquet(
        cdir / "data_0.parquet", index=False)
    return gf, tmp_path / "cells"


def test_hierarchy_and_zero_unassigned(tmp_path):
    rows = [
        {"id": "e",   "INS": (True, 5.0)},                          # Endocrine
        {"id": "i",   "CD3e": (True, 4.0)},                         # Immune
        {"id": "v",   "CD31": (True, 3.0)},                         # Endothelial
        {"id": "epi", "Pan_Cytokeratin": (False, 0.9),             # Epithelial by default
                      "Vimentin": (False, 0.1), "SMA": (False, 0.1)},
        {"id": "fib", "Vimentin": (True, 2.0),                      # Fibroblast (struct argmax)
                      "Pan_Cytokeratin": (False, 0.5), "SMA": (False, 0.3)},
        {"id": "mus", "SMA": (True, 3.0),                           # Muscle (struct argmax)
                      "Pan_Cytokeratin": (False, 0.5), "Vimentin": (False, 0.3)},
        {"id": "e_over_i", "INS": (True, 5.0), "CD3e": (True, 9.0)},   # Endocrine wins
        {"id": "i_over_v", "CD3e": (True, 4.0), "CD31": (True, 9.0)},  # Immune wins
    ]
    gf, cells_dir = _write(tmp_path, _gated_frame(rows))
    out = lineage.assign_donor("TEST", str(gf), cells_dir)

    call = dict(zip(out["object_id"], out["broad_lineage"]))
    assert call["e"] == "Endocrine"
    assert call["i"] == "Immune"
    assert call["v"] == "Endothelial"
    assert call["epi"] == "Epithelial"
    assert call["fib"] == "Fibroblast"
    assert call["mus"] == "Muscle"
    assert call["e_over_i"] == "Endocrine"   # endocrine precedence over immune
    assert call["i_over_v"] == "Immune"      # immune precedence over endothelial

    # the core invariant: every cell is typed into one of the six lineages
    assert set(out["broad_lineage"]).issubset(set(LINEAGES))
    assert (out["broad_lineage"] == "Unassigned").sum() == 0
    assert len(out) == len(rows)


def test_epi_default_flag(tmp_path):
    """Cells positive for no structural marker are Epithelial-by-default and flagged."""
    rows = [
        {"id": "weak", "Pan_Cytokeratin": (False, 0.2),
                       "Vimentin": (False, 0.1), "SMA": (False, 0.05)},
        {"id": "strongvim", "Vimentin": (True, 2.0)},
    ]
    gf, cells_dir = _write(tmp_path, _gated_frame(rows))
    out = lineage.assign_donor("TEST", str(gf), cells_dir).set_index("object_id")
    assert bool(out.loc["weak", "epi_default"]) is True
    assert out.loc["weak", "broad_lineage"] == "Epithelial"
    assert bool(out.loc["strongvim", "epi_default"]) is False


def test_endothelial_label_not_truncated(tmp_path):
    """Guard against the <U16 dtype bug: 'Endothelial' (11 chars) must survive."""
    rows = [{"id": "v", "CD31": (True, 3.0)}]
    gf, cells_dir = _write(tmp_path, _gated_frame(rows))
    out = lineage.assign_donor("TEST", str(gf), cells_dir)
    assert out["broad_lineage"].iloc[0] == "Endothelial"   # not "Endothelia"
