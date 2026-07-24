"""Unit tests for the tri-state ordered-residual phenotyping (5 compartments + Other + Unresolved)."""

from __future__ import annotations

import pandas as pd
import pytest

from phenocycler import lineage
from phenocycler.config import COMPARTMENT_ORDER, COMPARTMENT_GATES, OTHER_LABEL, UNRESOLVED_LABEL

CLASSES = COMPARTMENT_ORDER + [OTHER_LABEL, UNRESOLVED_LABEL]
GATE_MARKERS = sorted({m for gate in COMPARTMENT_GATES.values() for m in gate})


def _write(tmp_path, rows):
    """Write ONE tri-state restore_gated parquet + a cells partition for donor TEST.
    row = {"id": str, <gate_marker>: state, ...} with state in {"positive","negative","unavailable"};
    an unset gate marker defaults to "negative" (measured, available). Every gate marker gets a
    ``{m}_state`` column (restore_apply always emits all 7)."""
    ids = [r["id"] for r in rows]
    d = {"object_id": [str(i) for i in ids]}
    for m in GATE_MARKERS:
        d[f"{m}_state"] = [str(r.get(m, "negative")) for r in rows]
    gdir = tmp_path / "restore_gated_redsea" / "donor_id=TEST"; gdir.mkdir(parents=True)
    gf = gdir / "data_0.parquet"; pd.DataFrame(d).to_parquet(gf, index=False)
    cdir = tmp_path / "cells" / "donor_id=TEST"; cdir.mkdir(parents=True)
    pd.DataFrame({"object_id": [str(i) for i in ids], "cell_region": "tissue",
                  "islet_num": ""}).to_parquet(cdir / "data_0.parquet", index=False)
    return str(gf), tmp_path / "cells"


def _assign(tmp_path, rows):
    gf, cells_dir = _write(tmp_path, rows)
    return lineage.assign_donor("TEST", gf, cells_dir)


def test_priority_order_and_other(tmp_path):
    """Immune is FIRST now; the first POSITIVE gate claims; all-negative -> Other."""
    rows = [
        {"id": "imm", "CD3e": "positive"},
        {"id": "epi", "E_cadherin": "positive"},
        {"id": "eth", "CD31": "positive"},
        {"id": "nrl", "B3TUBB": "positive"},
        {"id": "mes", "Vimentin": "positive"},
        {"id": "other"},                                                    # all gates negative -> Other
        {"id": "imm>epi", "CD3e": "positive", "E_cadherin": "positive"},    # Immune wins (first in order)
        {"id": "epi>eth", "E_cadherin": "positive", "CD31": "positive"},    # Epithelial > Endothelial
        {"id": "eth>mes", "CD31": "positive", "Vimentin": "positive"},      # Endothelial > Mesenchymal
        {"id": "nrl>mes", "B3TUBB": "positive", "Vimentin": "positive"},    # Neural > Mesenchymal
    ]
    out = _assign(tmp_path, rows).set_index("object_id")
    c = out["compartment"].to_dict()
    assert c["imm"] == "Immune"
    assert c["epi"] == "Epithelial"
    assert c["eth"] == "Endothelial"
    assert c["nrl"] == "Neural"
    assert c["mes"] == "Mesenchymal"
    assert c["other"] == "Other"
    assert c["imm>epi"] == "Immune"
    assert c["epi>eth"] == "Epithelial"
    assert c["eth>mes"] == "Endothelial"
    assert c["nrl>mes"] == "Neural"
    assert set(c.values()).issubset(set(CLASSES))
    assert out.loc["imm", "assignment_reason"] == "positive"
    assert out.loc["other", "assignment_reason"] == "all_negative"
    assert out.loc["other", "blocked_by"] == ""


def test_immune_outranks_epithelial(tmp_path):
    """Multi-positive across compartments: the highest-priority POSITIVE gate wins (Immune > Epithelial)."""
    out = _assign(tmp_path, [{"id": "x", "CD3e": "positive", "E_cadherin": "positive"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Immune"
    assert out.loc["x", "assignment_reason"] == "positive"


def test_positive_with_unavailable_sibling_claims(tmp_path):
    """Within a gate, a POSITIVE marker dominates an UNAVAILABLE sibling (positive > unavailable)."""
    out = _assign(tmp_path, [{"id": "x", "CD3e": "positive", "CD68": "unavailable"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Immune"
    assert out.loc["x", "assignment_reason"] == "positive"


def test_unavailable_higher_priority_blocks_positive_lower(tmp_path):
    """An UNAVAILABLE higher-priority gate blocks a lower positive claim -> Unresolved, never a wrong class."""
    out = _assign(tmp_path, [{"id": "x", "CD3e": "unavailable", "E_cadherin": "positive"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Unresolved"
    assert out.loc["x", "assignment_reason"] == "blocked_unavailable_gate"
    assert out.loc["x", "blocked_by"] == "Immune"


def test_all_unavailable_is_unresolved(tmp_path):
    """Every gate indeterminate -> Unresolved, blocked by the highest-priority (first) gate."""
    row = {"id": "x"}
    for m in GATE_MARKERS:
        row[m] = "unavailable"
    out = _assign(tmp_path, [row]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Unresolved"
    assert out.loc["x", "blocked_by"] == "Immune"


def test_all_negative_is_other(tmp_path):
    """Valid-and-negative for every gate -> Other (a real panel gap, NOT Unresolved / force-assigned)."""
    out = _assign(tmp_path, [{"id": "x"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Other"
    assert out.loc["x", "assignment_reason"] == "all_negative"
    assert out.loc["x", "cell_type"] == "Other"


def test_confirmed_absence_negative_does_not_block(tmp_path):
    """restore_apply maps a confirmed absence to NEGATIVE (not unavailable); a negative gate never blocks,
    so the cell falls through to a lower gate rather than becoming Unresolved."""
    out = _assign(tmp_path, [{"id": "x", "CD3e": "negative", "CD68": "negative", "CD11b": "negative",
                              "E_cadherin": "positive"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Epithelial"


def test_unavailable_earlier_gate_blocks_even_when_later_gate_would_be_other(tmp_path):
    """Other requires EVERY gate valid-and-negative: an unavailable earlier gate -> Unresolved, not Other."""
    # Immune all-negative(available); Epithelial unavailable; the rest negative -> blocked at Epithelial.
    out = _assign(tmp_path, [{"id": "x", "E_cadherin": "unavailable"}]).set_index("object_id")
    assert out.loc["x", "compartment"] == "Unresolved"
    assert out.loc["x", "blocked_by"] == "Epithelial"


def test_missing_state_column_fails_loud(tmp_path):
    """A required gate-marker {m}_state absent from the gated parquet is a loud error, never an all-negative
    silent default (Gate-4 removal of the old absent->negative fallback)."""
    gf, cells_dir = _write(tmp_path, [{"id": "x", "CD3e": "positive"}])
    pd.read_parquet(gf).drop(columns=["CD11b_state"]).to_parquet(gf, index=False)
    with pytest.raises(KeyError, match="CD11b_state"):
        lineage.assign_donor("TEST", gf, cells_dir)


def test_provenance_columns_present(tmp_path):
    """The broad parquet carries the assignment provenance + per-marker tri-state passthrough."""
    out = _assign(tmp_path, [{"id": "x", "CD3e": "positive"}])
    for col in ("object_id", "donor_id", "compartment", "cell_type", "assignment_reason", "blocked_by"):
        assert col in out.columns
    for m in GATE_MARKERS:
        assert f"{m}_state" in out.columns                     # tri-state passthrough (provenance)
    assert (out["cell_type"] == out["compartment"]).all()      # subtypes deferred -> cell_type == compartment


def test_output_schema(tmp_path):
    """compartment + cell_type present; the old broad_lineage / score_ / epi_default columns are gone."""
    out = _assign(tmp_path, [{"id": "x", "E_cadherin": "positive"}])
    assert "compartment" in out.columns and "cell_type" in out.columns
    assert "broad_lineage" not in out.columns
    assert not any(c.startswith("score_") for c in out.columns)
    assert "epi_default" not in out.columns


def test_assign_and_write_atomic(tmp_path):
    """_assign_and_write writes the partition atomically (no .tmp leftover) and returns the summary."""
    _write(tmp_path, [{"id": "a", "CD3e": "positive"}, {"id": "b"}, {"id": "u", "CD3e": "unavailable"}])
    out_dir = tmp_path / "broad"
    summary = lineage._assign_and_write("TEST", str(tmp_path / "restore_gated_redsea"),
                                        str(tmp_path / "cells"), str(out_dir))
    assert summary["n"] == 3
    assert summary["other"] == 1 and summary["unresolved"] == 1
    part = out_dir / "donor_id=TEST"
    assert (part / "data_0.parquet").exists()
    assert not list(part.glob("*.tmp*"))                       # atomic write left no temp file
