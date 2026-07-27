"""Unit tests for the canonical cells-parquet build.

`build_sql` had no test at all before 2026-07-26, despite owning the pipeline's only cell filter and
defining the canonical cell universe that REDSEA, the RESTORE compartment extraction and
`load_validation_sample` all derive from by join.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
import pytest

from phenocycler import cells_parquet as cp


HEADER = [
    "Image", "Object ID", "Centroid X µm", "Centroid Y µm", "Parent", "Classification",
    "Cell: Area µm^2", "Nucleus: Area µm^2", "Cell: Length µm", "Cell: Circularity",
    "Cell: Solidity", "Cell: Max diameter µm", "Cell: Min diameter µm",
    "Nucleus: Length µm", "Nucleus: Circularity", "Nucleus: Solidity",
    "Nucleus: Max diameter µm", "Nucleus: Min diameter µm", "Nucleus/Cell area ratio",
    "Signed distance to annotation Islet µm", "Signed distance to annotation Tissue µm",
    "Cell: DAPI: Mean", "Cell: E-cadherin: Mean", "Nucleus: DAPI: Mean",
]


def _row(oid, area, parent="Annotation (Tissue)", image="6374_Scan1.er.qptiff", nuc=30.0):
    try:
        ratio = nuc / float(area)
    except (TypeError, ValueError, ZeroDivisionError):
        ratio = ""            # unparseable area -> unparseable ratio, as QuPath would emit
    return {
        **{c: "" for c in HEADER},
        "Image": image, "Object ID": oid, "Centroid X µm": "1.0", "Centroid Y µm": "2.0",
        "Parent": parent, "Classification": "",
        "Cell: Area µm^2": str(area), "Nucleus: Area µm^2": str(nuc),
        "Cell: Length µm": "10", "Cell: Circularity": "0.5", "Cell: Solidity": "0.9",
        "Cell: Max diameter µm": "10", "Cell: Min diameter µm": "7",
        "Nucleus: Length µm": "6", "Nucleus: Circularity": "0.8", "Nucleus: Solidity": "0.95",
        "Nucleus: Max diameter µm": "6", "Nucleus: Min diameter µm": "5",
        "Nucleus/Cell area ratio": str(ratio),
        "Signed distance to annotation Islet µm": "-5",
        "Signed distance to annotation Tissue µm": "3",
        "Cell: DAPI: Mean": "500", "Cell: E-cadherin: Mean": "700", "Nucleus: DAPI: Mean": "600",
    }


def _csv(tmp_path, rows) -> Path:
    p = tmp_path / "Cellmeasurements.csv"
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    return p


def _run(tmp_path, rows, min_cell_area):
    duckdb = pytest.importorskip("duckdb")
    csv_path = _csv(tmp_path, rows)
    out = tmp_path / "out"
    (out / "cells").mkdir(parents=True, exist_ok=True)
    sql, mm = cp.build_sql(csv_path, out, None, min_cell_area)
    duckdb.connect().execute(sql)
    # donor_id is a hive PARTITION key (a `donor_id=<id>` directory), not a stored column, so it has
    # to come back off the path — exactly as `config.discover_donors` reads it in production.
    files = sorted((out / "cells").rglob("*.parquet"))
    parts = []
    for f in files:
        part = pd.read_parquet(f)
        part["donor_id"] = f.parent.name.split("=", 1)[1]
        parts.append(part)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    return df, mm


def test_marker_means_maps_only_whole_cell_means():
    mm = cp.marker_means(HEADER)
    assert mm == {"Cell: DAPI: Mean": "DAPI", "Cell: E-cadherin: Mean": "E-cadherin"}
    assert "Nucleus: DAPI: Mean" not in mm      # nuclear compartment is not carried forward


def test_sanitize_matches_the_redsea_column_convention():
    assert cp.sanitize("E-cadherin") == "E_cadherin"
    assert cp.sanitize("Pan-Cytokeratin") == "Pan_Cytokeratin"
    assert cp.sanitize("Ker8/18") == "Ker8_18"


def test_min_cell_area_predicate_is_inclusive_at_the_threshold(tmp_path):
    """The filter is `>=`, so a cell exactly at the floor is retained."""
    rows = [_row("a", 19.999), _row("b", 20.0), _row("c", 20.001), _row("d", 100.0)]
    df, _ = _run(tmp_path, rows, 20.0)
    assert sorted(df.object_id) == ["b", "c", "d"]


def test_min_cell_area_zero_disables_the_filter(tmp_path):
    rows = [_row("a", 1.0), _row("b", 100.0)]
    df, _ = _run(tmp_path, rows, 0.0)
    assert sorted(df.object_id) == ["a", "b"]


def test_non_finite_area_is_always_dropped(tmp_path):
    """`isfinite` rides inside the filter, so unparseable or infinite areas never survive."""
    rows = [_row("ok", 60.0), _row("blank", ""), _row("junk", "not-a-number"), _row("inf", "inf")]
    df, _ = _run(tmp_path, rows, 20.0)
    assert sorted(df.object_id) == ["ok"]


def test_donor_id_is_the_first_digit_run_of_the_image_name(tmp_path):
    rows = [
        _row("a", 60.0, image="6374_Scan1.er.qptiff"),
        _row("b", 60.0, image="HDL115pancLN_Scan1.er.qptiff"),
        _row("c", 60.0, image="6450panc_Scan1.er.qptiff"),
    ]
    df, _ = _run(tmp_path, rows, 20.0)
    assert dict(zip(df.object_id, df.donor_id)) == {"a": "6374", "b": "115", "c": "6450"}


def test_cell_region_scheme(tmp_path):
    rows = [
        _row("core", 60.0, parent="Islet_7"),
        _row("peri", 60.0, parent="Islet_7_exp20um"),
        _row("tissue", 60.0, parent="Annotation (Tissue)"),
        _row("root", 60.0, parent="Root object (Image)"),
        _row("other", 60.0, parent="Something Else"),
    ]
    df, _ = _run(tmp_path, rows, 20.0)
    got = dict(zip(df.object_id, df.cell_region))
    assert got == {"core": "core", "peri": "peri", "tissue": "tissue",
                   "root": "tissue", "other": "other"}
    assert dict(zip(df.object_id, df.islet_num))["core"] == "7"


def test_morphology_columns_needed_by_cell_qc_are_all_present(tmp_path):
    """cell_qc reads these directly off the canonical parquet — none may go missing."""
    df, _ = _run(tmp_path, [_row("a", 60.0)], 20.0)
    from phenocycler.cell_qc import QC_SOURCE_COLUMNS

    missing = set(QC_SOURCE_COLUMNS) - set(df.columns)
    assert not missing, f"cells parquet is missing QC inputs: {sorted(missing)}"


def test_build_sql_rejects_a_csv_with_no_marker_columns(tmp_path):
    p = tmp_path / "bad.csv"
    with open(p, "w", newline="") as f:
        csv.writer(f).writerow(["Image", "Object ID", "Cell: Area µm^2"])
    with pytest.raises(SystemExit, match="No 'Cell: <marker>: Mean' columns"):
        cp.build_sql(p, tmp_path / "out", None, 20.0)
