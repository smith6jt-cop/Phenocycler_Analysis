"""Unit tests for single-cell QC (no image data required)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from phenocycler import cell_qc


PX2 = 0.5077663810243286 ** 2


def _cfg(tmp_path):
    return SimpleNamespace(
        data_dir=tmp_path,
        cells_dir=tmp_path / "cells",
        cells_redsea_dir=tmp_path / "cells_redsea",
    )


def _write(tmp_path, donor, cells, redsea=None):
    d = tmp_path / "cells" / f"donor_id={donor}"
    d.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cells).to_parquet(d / "data_0.parquet", index=False)
    if redsea is not None:
        r = tmp_path / "cells_redsea" / f"donor_id={donor}"
        r.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(redsea).to_parquet(r / "data_0.parquet", index=False)


def _cells(**over):
    """Six clean 60 um^2 cells, 10 um apart; override any column to inject a defect."""
    n = 6
    base = {
        "object_id": [f"c{i}" for i in range(n)],
        "X_centroid": [10.0 * i for i in range(n)],
        "Y_centroid": [0.0] * n,
        "cell_area": [60.0] * n,
        "nucleus_area": [30.0] * n,
        "nc_area_ratio": [0.5] * n,
        "cell_solidity": [0.92] * n,
        "cell_circularity": [0.52] * n,
        "cell_maxdiam": [10.0] * n,
        "cell_mindiam": [7.0] * n,
        "DAPI": [500.0] * n,
    }
    base.update(over)
    return base


def _redsea(areas):
    return {"object_id": [f"c{i}" for i in range(len(areas))],
            "cell_area_px": [a / PX2 for a in areas]}


def test_clean_cells_all_pass(tmp_path):
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_pass.all()
    assert (qc.qc_exclude_reason == "").all()
    assert qc.qc_estimation_eligible.all()


def test_no_nucleus_is_excluded(tmp_path):
    """The DNA gate: an object without a nucleus is not a cell, whatever its area."""
    cfg = _cfg(tmp_path)
    nuc = [30.0, np.nan, 0.0, 30.0, 30.0, 30.0]
    _write(tmp_path, "1", _cells(nucleus_area=nuc), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == ["", "no_nucleus", "no_nucleus", "", "", ""]
    # Area alone would never have caught these — they are full-sized objects.
    assert (qc.loc[~qc.qc_pass, "cell_area"] == 60.0).all()


def test_area_floor_and_donor_relative_ceiling(tmp_path):
    cfg = _cfg(tmp_path)
    areas = [10.0, 60.0, 60.0, 60.0, 60.0, 400.0]     # one tiny, one 6.7x the median
    _write(tmp_path, "1", _cells(cell_area=areas), _redsea(areas))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == [
        "area_below_min", "", "", "", "", "area_above_max"]
    # The ceiling is 4x the donor's OWN median of the size-independent survivors, not an absolute.
    assert qc.attrs["donor_median_area"] == 60.0
    assert qc.attrs["donor_max_area"] == 240.0


def test_area_ceiling_is_donor_local_not_absolute(tmp_path):
    """A donor with uniformly larger cells gets a proportionally larger ceiling."""
    cfg = _cfg(tmp_path)
    for donor, scale in (("small", 1.0), ("large", 2.0)):
        areas = [60.0 * scale] * 5 + [200.0 * scale]
        _write(tmp_path, donor, _cells(cell_area=areas), _redsea(areas))
    cuts = {
        d: cell_qc.compute_donor_qc(_cfg(tmp_path), d, cell_qc.CellQCThresholds()).attrs[
            "donor_max_area"]
        for d in ("small", "large")
    }
    assert cuts["large"] == pytest.approx(2 * cuts["small"])


def test_low_solidity_is_excluded_but_circularity_is_not(tmp_path):
    """Solidity gates; circularity is recorded only — it is compartment-biased on real data."""
    cfg = _cfg(tmp_path)
    sol = [0.92, 0.50, 0.92, 0.92, 0.92, 0.92]
    circ = [0.52, 0.52, 0.05, 0.52, 0.52, 0.52]      # c2 is extremely non-circular
    _write(tmp_path, "1", _cells(cell_solidity=sol, cell_circularity=circ), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == ["", "low_solidity", "", "", "", ""]
    assert qc.qc_pass.iloc[2]          # the near-zero-circularity cell survives, by design


def test_degenerate_rasterization_is_excluded(tmp_path):
    """The 6450/6523 defect: the polygon reports 60 um^2 but rasterizes to ~2 px."""
    cfg = _cfg(tmp_path)
    px = [60.0 / PX2, 2.0, 60.0 / PX2, 0.0, 60.0 / PX2, 60.0 / PX2]
    _write(tmp_path, "1", _cells(),
           {"object_id": [f"c{i}" for i in range(6)], "cell_area_px": px})
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == [
        "", "raster_degenerate", "", "raster_empty", "", ""]


def test_missing_redsea_skips_raster_predicates_without_failing(tmp_path):
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells())          # no REDSEA output at all
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_pass.all()
    assert qc.raster_area_ratio.isna().all()


def test_no_cytoplasm_is_tier2_called_but_not_estimation_eligible(tmp_path):
    """These are real cells — they keep their call and are only barred from setting thresholds."""
    cfg = _cfg(tmp_path)
    nc = [0.5, 1.0, 0.5, 0.5, 0.5, 0.5]
    _write(tmp_path, "1", _cells(nc_area_ratio=nc), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_pass.all()                          # retained
    assert qc.qc_no_cytoplasm.tolist() == [False, True, False, False, False, False]
    assert qc.qc_estimation_eligible.tolist() == [True, False, True, True, True, True]


def test_duplicate_detections_keep_exactly_one_copy(tmp_path):
    """The islet double-segmentation: one nucleus, two UUIDs ~1 px apart. Drop one, keep one."""
    cfg = _cfg(tmp_path)
    x = [0.0, 0.5, 10.0, 20.0, 30.0, 40.0]      # c0/c1 are the same cell (0.5 um apart)
    area = [60.0, 59.7, 60.0, 60.0, 60.0, 60.0]
    px = [60.0, 12.0, 60.0, 60.0, 60.0, 60.0]   # c1 lost the raster fight -> c0 is the survivor
    _write(tmp_path, "1", _cells(X_centroid=x, cell_area=area),
           {"object_id": [f"c{i}" for i in range(6)], "cell_area_px": [a / PX2 for a in px]})
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == ["", "duplicate_detection", "", "", "", ""]
    assert qc.qc_pass.sum() == 5


def test_duplicate_survivor_is_the_copy_that_won_the_raster(tmp_path):
    """Keep the copy whose polygon won the mask — its REDSEA values are the intact ones."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(X_centroid=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0]),
           {"object_id": [f"c{i}" for i in range(6)],
            "cell_area_px": [v / PX2 for v in [12.0, 60.0, 60.0, 60.0, 60.0, 60.0]]})
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    # c0 is the raster loser this time, so c1 survives instead — the rule follows the evidence.
    assert qc.qc_exclude_reason.tolist() == ["duplicate_detection", "", "", "", "", ""]


def test_duplicate_removal_is_deterministic_without_raster_data(tmp_path):
    """With no REDSEA output the tie-break is the lexicographically smallest object_id."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(X_centroid=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0]))
    a = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    b = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert a.qc_exclude_reason.tolist() == b.qc_exclude_reason.tolist()
    assert a.qc_exclude_reason.tolist() == ["", "duplicate_detection", "", "", "", ""]


def test_a_triple_overlap_keeps_one_not_two(tmp_path):
    """Connected components, so a 3-way cluster leaves a single survivor."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(X_centroid=[0.0, 0.5, 0.9, 20.0, 30.0, 40.0]),
           _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert (qc.qc_exclude_reason == "duplicate_detection").sum() == 2
    assert qc.qc_pass.sum() == 4


def test_different_areas_are_neighbours_not_duplicates(tmp_path):
    """Adjacent cells of clearly different size are two real cells — never merged."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1",
           _cells(X_centroid=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
                  cell_area=[60.0, 25.0, 60.0, 60.0, 60.0, 60.0]),
           _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert (qc.qc_exclude_reason == "duplicate_detection").sum() == 0


def test_duplicate_whose_partner_already_failed_qc_needs_no_choice(tmp_path):
    """Dedup runs on the Tier-1 survivors, so a pair with one bad member resolves for free."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1",
           _cells(X_centroid=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0],
                  nucleus_area=[30.0, np.nan, 30.0, 30.0, 30.0, 30.0]),
           _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    assert qc.qc_exclude_reason.tolist() == ["", "no_nucleus", "", "", "", ""]
    assert (qc.qc_exclude_reason == "duplicate_detection").sum() == 0   # c0 kept, no arbitrary drop


def test_deduplicate_can_be_disabled(tmp_path):
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(X_centroid=[0.0, 0.5, 10.0, 20.0, 30.0, 40.0]),
           _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds(deduplicate=False))
    assert (qc.qc_exclude_reason == "duplicate_detection").sum() == 0


def test_exclusion_reasons_partition_the_excluded_set(tmp_path):
    """First reason wins, so per-reason counts sum exactly to the number excluded."""
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1",
           _cells(cell_area=[10.0, 60.0, 60.0, 60.0, 60.0, 60.0],
                  nucleus_area=[np.nan, 30.0, 30.0, 30.0, 30.0, 30.0],   # c0 fails BOTH
                  cell_solidity=[0.92, 0.50, 0.92, 0.92, 0.92, 0.92]),
           _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    row = cell_qc.donor_summary("1", qc)
    per_reason = sum(row[f"n_{r}"] for r in cell_qc.EXCLUSION_ORDER)
    assert per_reason == row["excluded"] == 2
    assert qc.qc_exclude_reason.iloc[0] == "no_nucleus"   # earlier in EXCLUSION_ORDER wins


def test_apply_filters_canonical_and_keeps_a_reversible_backup(tmp_path):
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(cell_area=[10.0] + [60.0] * 5), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(cfg, "1", cell_qc.CellQCThresholds())
    n = cell_qc.apply_qc_to_canonical(cfg, "1", qc)

    path = tmp_path / "cells" / "donor_id=1" / "data_0.parquet"
    backup = tmp_path / "cells_preqc" / "donor_id=1" / "data_0.parquet"
    assert n == 5
    assert len(pd.read_parquet(path)) == 5
    assert "c0" not in set(pd.read_parquet(path).object_id)
    assert len(pd.read_parquet(backup)) == 6          # pre-QC state preserved

    # The backup must live OUTSIDE cells/. Several stages glob cells/donor_id=*/*.parquet and
    # require exactly one match — redsea.reconcile_existing_donor raises "expected one raw and one
    # REDSEA parquet" the moment a second file appears there.
    assert len(list((tmp_path / "cells" / "donor_id=1").glob("*.parquet"))) == 1


def test_require_nucleus_can_be_disabled(tmp_path):
    cfg = _cfg(tmp_path)
    _write(tmp_path, "1", _cells(nucleus_area=[np.nan] * 6), _redsea([60.0] * 6))
    qc = cell_qc.compute_donor_qc(
        cfg, "1", cell_qc.CellQCThresholds(require_nucleus=False))
    assert qc.qc_pass.all()


@pytest.mark.parametrize(
    "kwargs,match",
    [
        ({"min_cell_area": -1.0}, "min_cell_area"),
        ({"max_area_median_multiple": 1.0}, "max_area_median_multiple"),
        ({"min_cell_solidity": 1.5}, "min_cell_solidity"),
        ({"min_raster_area_ratio": -0.1}, "min_raster_area_ratio"),
        ({"no_cytoplasm_nc_ratio": 0.0}, "no_cytoplasm_nc_ratio"),
    ],
)
def test_invalid_thresholds_are_rejected(kwargs, match):
    with pytest.raises(ValueError, match=match):
        cell_qc.CellQCThresholds(**kwargs)


def test_thresholds_from_config_match_the_shipped_defaults():
    from phenocycler import load_config

    t = cell_qc.CellQCThresholds.from_config(load_config())
    assert t.min_cell_area == 20.0
    assert t.require_nucleus is True
    assert t.max_area_median_multiple == 4.0
    assert t.min_cell_solidity == 0.70
    assert t.min_raster_area_ratio == 0.5
    assert t.no_cytoplasm_nc_ratio == 0.999
