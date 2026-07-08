"""Unit tests for the positivity _norm floor (hormones @ K5 + immune @ K2)."""

from __future__ import annotations

import pandas as pd

from phenocycler import hormone_floor, load_config
from phenocycler.config import HORMONE_MARKERS, IMMUNE_FLOOR_MARKERS


def test_floor_donor_file_rewrites_pos_by_norm(tmp_path):
    """floor_donor_file rewrites {m}_pos = (_norm >= K) per marker; _norm + other columns pass through."""
    df = pd.DataFrame({
        "object_id": ["a", "b", "c", "d"],
        "INS_pos":  [True, True, True, True], "INS_norm":  [1.1, 6.0, 0.5, 20.0],  # hormone K=5
        "CD3e_pos": [True, True, True, True], "CD3e_norm": [1.2, 3.0, 1.9, 8.0],   # immune  K=2
        "other_col": [1, 2, 3, 4],
    })
    src = tmp_path / "in.parquet"; df.to_parquet(src, index=False)
    out = tmp_path / "out.parquet"
    n, report = hormone_floor.floor_donor_file(str(src), str(out), {"INS": 5.0, "CD3e": 2.0})
    r = pd.read_parquet(out)

    assert n == 4
    assert list(r["INS_pos"]) == [False, True, False, True]      # only _norm >= 5 survive
    assert list(r["CD3e_pos"]) == [False, True, False, True]     # only _norm >= 2 survive
    assert list(r["INS_norm"]) == [1.1, 6.0, 0.5, 20.0]          # _norm untouched
    assert list(r["other_col"]) == [1, 2, 3, 4]                  # non-floored columns pass through


def test_default_floors_covers_hormones_and_immune():
    cfg = load_config()
    floors = hormone_floor.default_floors(cfg)
    assert all(floors[m] == cfg.hormone_min_norm for m in HORMONE_MARKERS)
    assert all(floors[m] == cfg.immune_min_norm for m in IMMUNE_FLOOR_MARKERS)
    assert "MPO" not in floors     # MPO (extra dir) is gated in lineage.py, not the stage floor
