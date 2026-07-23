"""Unit tests for the 5-compartment ordered-residual phenotyping (+ explicit Other)."""

from __future__ import annotations

import pandas as pd

from phenocycler import lineage
from phenocycler.config import COMPARTMENT_ORDER, OTHER_LABEL

CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL]


def _write(tmp_path, rows):
    """Write ONE restore_gated parquet (all markers) + a cells partition for donor TEST.
    row = {"id": str, <marker>: (pos: bool, norm: float), ...}; an unset marker -> (False, 0.0)."""
    ids = [r["id"] for r in rows]
    markers = sorted({m for r in rows for m in r if m != "id"})
    d = {"object_id": [str(i) for i in ids]}
    for m in markers:
        d[f"{m}_pos"] = [bool(r.get(m, (False, 0.0))[0]) for r in rows]
        d[f"{m}_norm"] = [float(r.get(m, (False, 0.0))[1]) for r in rows]
    gdir = tmp_path / "restore_gated_redsea" / "donor_id=TEST"; gdir.mkdir(parents=True)
    gf = gdir / "data_0.parquet"; pd.DataFrame(d).to_parquet(gf, index=False)
    cdir = tmp_path / "cells" / "donor_id=TEST"; cdir.mkdir(parents=True)
    pd.DataFrame({"object_id": [str(i) for i in ids], "cell_region": "tissue",
                  "islet_num": ""}).to_parquet(cdir / "data_0.parquet", index=False)
    return str(gf), tmp_path / "cells"


def _assign(tmp_path, rows):
    gf, cells_dir = _write(tmp_path, rows)
    out, _ = lineage.assign_donor("TEST", gf, cells_dir)
    return out


def test_priority_order_and_other(tmp_path):
    rows = [
        {"id": "epi", "E_cadherin": (True, 3)},
        {"id": "eth", "CD31": (True, 3)},
        {"id": "nrl", "B3TUBB": (True, 3)},
        {"id": "imm", "CD3e": (True, 3)},
        {"id": "mus", "SMA": (True, 3)},
        {"id": "fib", "Vimentin": (True, 3)},
        {"id": "other"},                                                 # nothing positive -> Other
        {"id": "epi>eth", "E_cadherin": (True, 3), "CD31": (True, 9)},    # Epithelial wins
        {"id": "eth>mes", "CD31": (True, 3), "SMA": (True, 9)},           # Endothelial > Mesenchymal
        {"id": "nrl>mes", "B3TUBB": (True, 3), "Vimentin": (True, 9)},    # Neural > Mesenchymal
        {"id": "imm>mes", "CD3e": (True, 3), "Vimentin": (True, 9)},      # Immune > Mesenchymal
    ]
    out = _assign(tmp_path, rows)
    c = dict(zip(out["object_id"], out["compartment"]))
    assert c["epi"] == "Epithelial"
    assert c["eth"] == "Endothelial"
    assert c["nrl"] == "Neural"
    assert c["imm"] == "Immune"
    assert c["mus"] == "Mesenchymal" and c["fib"] == "Mesenchymal"
    assert c["other"] == "Other"
    assert c["epi>eth"] == "Epithelial"
    assert c["eth>mes"] == "Endothelial"
    assert c["nrl>mes"] == "Neural"
    assert c["imm>mes"] == "Immune"
    assert set(c.values()).issubset(set(CLASSES))


def test_other_terminal_not_force_assigned(tmp_path):
    out = _assign(tmp_path, [{"id": "x", "CD44": (True, 5)},              # state marker, no gate
                             {"id": "y"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Other" and out.loc["x", "cell_type"] == "Other"
    assert out.loc["y", "compartment"] == "Other"


def test_endocrine_vs_exocrine_subsplit(tmp_path):
    rows = [
        {"id": "beta", "E_cadherin": (True, 3), "INS": (True, 8)},
        {"id": "alpha", "E_cadherin": (True, 3), "GCG": (True, 8)},
        {"id": "delta", "E_cadherin": (True, 3), "SST": (True, 8)},
        {"id": "iapp_only", "E_cadherin": (True, 3), "IAPP": (True, 8)},   # IAPP dropped 2026-07-10 -> hormone-neg
        {"id": "exo", "E_cadherin": (True, 3)},
        {"id": "ductal", "E_cadherin": (True, 3), "Keratin_5": (True, 4)},
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert (out["compartment"] == "Epithelial").all()
    assert out.loc["beta", "cell_type"] == "Endocrine-beta"
    assert out.loc["alpha", "cell_type"] == "Endocrine-alpha"
    assert out.loc["delta", "cell_type"] == "Endocrine-delta"
    assert out.loc["iapp_only", "cell_type"] == "Exocrine"   # IAPP no longer an endocrine marker -> hormone-neg Epithelial
    assert out.loc["exo", "cell_type"] == "Exocrine"
    assert out.loc["ductal", "cell_type"] == "Exocrine-ductal"


def test_immune_subtypes_incl_nk(tmp_path):
    rows = [
        {"id": "t", "CD3e": (True, 3)},
        {"id": "b", "CD20": (True, 3)},
        {"id": "mac", "CD68": (True, 3)},
        {"id": "neut", "MPO": (True, 3)},
        {"id": "dc", "CD11c": (True, 3)},
        {"id": "nk", "CD56": (True, 3), "CD57": (True, 3)},              # CD3e- -> NK, must ENTER Immune
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert (out["compartment"] == "Immune").all()
    assert out.loc["t", "cell_type"] == "T"
    assert out.loc["b", "cell_type"] == "B"
    assert out.loc["mac", "cell_type"] == "Macrophage"
    assert out.loc["neut", "cell_type"] == "Neutrophil"
    assert out.loc["dc", "cell_type"] == "DC"
    assert out.loc["nk", "cell_type"] == "NK"


def test_mesenchymal_muscle_vs_fibroblast(tmp_path):
    rows = [
        {"id": "mus", "SMA": (True, 3)},
        {"id": "fib", "Vimentin": (True, 3)},
        {"id": "both", "SMA": (True, 3), "Vimentin": (True, 3)},         # SMA claimed first -> Muscle
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert (out["compartment"] == "Mesenchymal").all()
    assert out.loc["mus", "cell_type"] == "Muscle"
    assert out.loc["fib", "cell_type"] == "Fibroblast"
    assert out.loc["both", "cell_type"] == "Muscle"


def test_long_labels_not_truncated(tmp_path):
    """Guard the fixed-width dtype: 'Endothelial'(11) compartment + 'Endothelial-lymphatic'(21) cell_type."""
    out = _assign(tmp_path, [{"id": "l", "CD31": (True, 3), "Podoplanin": (True, 3)}]).set_index("object_id")
    assert out.loc["l", "compartment"] == "Endothelial"
    assert out.loc["l", "cell_type"] == "Endothelial-lymphatic"


def test_output_schema(tmp_path):
    """compartment + cell_type present; the old broad_lineage / score_ / epi_default are gone."""
    out = _assign(tmp_path, [{"id": "e", "E_cadherin": (True, 3), "INS": (True, 5)}])
    assert "compartment" in out.columns and "cell_type" in out.columns
    assert "broad_lineage" not in out.columns
    assert not any(c.startswith("score_") for c in out.columns)
    assert "epi_default" not in out.columns


def test_frequent_immune_types_outrank_sparse_b_cells(tmp_path):
    """CD20 is the least reliable immune marker and B cells are sparse in pancreas, so a more frequent
    co-positive identity must win. Previously B was written last and overwrote macrophages."""
    rows = [
        {"id": "b_only", "CD20": (True, 3)},
        {"id": "b_cd79a", "CD79a": (True, 3)},
        {"id": "b_vs_mac", "CD20": (True, 9), "CD68": (True, 3)},     # -> Macrophage, not B
        {"id": "b_vs_dc", "CD20": (True, 9), "CD11c": (True, 3)},     # -> DC, not B
        {"id": "b_vs_neut", "CD20": (True, 9), "MPO": (True, 3)},     # -> Neutrophil, not B
        {"id": "b_vs_t", "CD20": (True, 9), "CD3e": (True, 3)},       # -> T, not B
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    assert (out["compartment"] == "Immune").all()
    # A cell with no competing immune marker is still B.
    assert out.loc["b_only", "cell_type"] == "B"
    assert out.loc["b_cd79a", "cell_type"] == "B"
    # Every more frequent immune identity now outranks it, regardless of relative intensity.
    assert out.loc["b_vs_mac", "cell_type"] == "Macrophage"
    assert out.loc["b_vs_dc", "cell_type"] == "DC"
    assert out.loc["b_vs_neut", "cell_type"] == "Neutrophil"
    assert out.loc["b_vs_t", "cell_type"] == "T"
