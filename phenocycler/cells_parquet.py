#!/usr/bin/env python3
"""
Step 1 — stream the QuPath per-cell export into a tidy, donor-partitioned
Parquet dataset via DuckDB.

Faithful port of ``scripts/senior/build_cells_parquet.py`` (Islet-Explorer-
Senior).  The only changes vs. upstream are that the input CSV, output root and
thread count come from :class:`phenocycler.config.PipelineConfig` (config.ini)
instead of hardcoded ``~/IO60panc2nd`` paths.

Input  : <cells_csv>                     (the full QuPath ``Cellmeasurements.csv``,
                                           ~23.3M cells, ~1,228 cols, ~138 GB)
Output : <data_dir>/cells/donor_id=<key>/*.parquet   (~65 cols: identity, spatial,
         morphology, region, and the whole-cell marker means)

``<key>`` is a **section** key, not a bare donor id — ``6539`` for a donor's pancreas and
``6539pln`` for its lymph node, derived from the qptiff name by
:data:`phenocycler.sections.SECTION_KEY_SQL`. One partition is one image throughout this
pipeline: REDSEA masks a partition with a single GeoJSON, RESTORE fits one threshold set per
partition, and micron coordinates are measured from one slide origin. Keying on the donor
alone put a donor's two images in one partition, which is not a labelling nuisance but a
silent correctness failure — see ``phenocycler/sections.py`` for the full argument.

Region scheme: Parent ``Islet_N`` = core, ``Islet_N_exp20um`` = peri,
``Annotation (Tissue)`` / ``Root object`` = tissue, else other.

Performance: DuckDB is the one already-parallel stage upstream — it streams the
138 GB CSV out-of-core with ``PARTITION_BY donor_id`` + ZSTD.  We keep that and
expose ``threads`` (and an optional ``memory_limit``) via config.

    python -m phenocycler.cells_parquet                    # full run
    python -m phenocycler.cells_parquet --limit 200000 --out /tmp/senior_test
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Optional

from .config import PipelineConfig, load_config
from .sections import SECTION_KEY_SQL


def header_cols(csv_path: Path) -> list[str]:
    with open(csv_path, newline="") as f:
        return next(csv.reader(f))


def marker_means(cols: list[str]) -> dict[str, str]:
    """Map 'Cell: <marker>: Mean' -> <marker>."""
    out = {}
    for c in cols:
        m = re.match(r"^Cell: (.+): Mean$", c)
        if m:
            out[c] = m.group(1)
    return out


def sanitize(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def build_sql(csv_path: Path, out_dir: Path, limit: Optional[int],
              min_cell_area: float) -> tuple[str, dict]:
    cols = header_cols(csv_path)
    mm = marker_means(cols)
    if not mm:
        raise SystemExit("No 'Cell: <marker>: Mean' columns found — check the CSV header.")
    selects = [
        # Section key, not donor id: '6539_Scan1...' -> 6539, '6539pLN_Scan1...' -> 6539pln.
        # The column keeps the name `donor_id` because it is the hive partition key every
        # downstream stage already globs; `phenocycler.sections.parse` takes it apart.
        f"{SECTION_KEY_SQL} AS donor_id",
        '"Image" AS image',
        '"Object ID" AS object_id',
        "TRY_CAST(\"Centroid X µm\" AS DOUBLE) AS X_centroid",
        "TRY_CAST(\"Centroid Y µm\" AS DOUBLE) AS Y_centroid",
        '"Parent" AS parent_raw',
        '"Classification" AS classification',
        ("CASE WHEN \"Parent\" LIKE '%exp20um' THEN 'peri' "
         "WHEN regexp_matches(\"Parent\", '^Islet_[0-9]+$') THEN 'core' "
         "WHEN \"Parent\" = 'Annotation (Tissue)' OR \"Parent\" LIKE 'Root object%' THEN 'tissue' "
         "ELSE 'other' END AS cell_region"),
        "regexp_extract(\"Parent\", 'Islet_([0-9]+)', 1) AS islet_num",
        "TRY_CAST(\"Cell: Area µm^2\" AS DOUBLE) AS cell_area",
        "TRY_CAST(\"Nucleus: Area µm^2\" AS DOUBLE) AS nucleus_area",
        # extra morphology + signed-distance columns (needed by downstream aggregation)
        "TRY_CAST(\"Cell: Length µm\" AS DOUBLE) AS cell_length",
        "TRY_CAST(\"Cell: Circularity\" AS DOUBLE) AS cell_circularity",
        "TRY_CAST(\"Cell: Solidity\" AS DOUBLE) AS cell_solidity",
        "TRY_CAST(\"Cell: Max diameter µm\" AS DOUBLE) AS cell_maxdiam",
        "TRY_CAST(\"Cell: Min diameter µm\" AS DOUBLE) AS cell_mindiam",
        "TRY_CAST(\"Nucleus: Length µm\" AS DOUBLE) AS nucleus_length",
        "TRY_CAST(\"Nucleus: Circularity\" AS DOUBLE) AS nucleus_circularity",
        "TRY_CAST(\"Nucleus: Solidity\" AS DOUBLE) AS nucleus_solidity",
        "TRY_CAST(\"Nucleus: Max diameter µm\" AS DOUBLE) AS nucleus_maxdiam",
        "TRY_CAST(\"Nucleus: Min diameter µm\" AS DOUBLE) AS nucleus_mindiam",
        "TRY_CAST(\"Nucleus/Cell area ratio\" AS DOUBLE) AS nc_area_ratio",
        "TRY_CAST(\"Signed distance to annotation Islet µm\" AS DOUBLE) AS dist_islet",
        "TRY_CAST(\"Signed distance to annotation Tissue µm\" AS DOUBLE) AS dist_tissue",
    ]
    for orig, marker in mm.items():
        selects.append(f'TRY_CAST("{orig}" AS DOUBLE) AS "{sanitize(marker)}"')
    sel = ",\n  ".join(selects)
    src = (f"read_csv('{csv_path}', header=true, all_varchar=true, "
           f"ignore_errors=true)")
    inner = f"SELECT\n  {sel}\nFROM {src}"
    if min_cell_area and min_cell_area > 0:
        # Drop degenerate sub-cell segmentation objects (QuPath emits NaN
        # intensities for these tiny artifacts; real cells have area >~12 µm²).
        area = "TRY_CAST(\"Cell: Area µm^2\" AS DOUBLE)"
        inner += f"\nWHERE isfinite({area}) AND {area} >= {min_cell_area}"
    if limit:
        inner += f"\nLIMIT {limit}"
    cells_dir = out_dir / "cells"
    sql = (f"COPY (\n{inner}\n) TO '{cells_dir}' "
           f"(FORMAT PARQUET, PARTITION_BY (donor_id), "
           f"OVERWRITE_OR_IGNORE, COMPRESSION ZSTD)")
    return sql, mm


def build_cells_parquet(cfg: PipelineConfig, *, limit: Optional[int] = None,
                        memory_limit: Optional[str] = None) -> Path:
    """Run the DuckDB CSV -> partitioned-parquet build. Returns the cells dir."""
    import duckdb

    if not cfg.cells_csv.exists():
        raise SystemExit(f"[err] QuPath cell CSV not found: {cfg.cells_csv}\n"
                         "      set [paths] cells_csv in config.ini")
    sql, mm = build_sql(cfg.cells_csv, cfg.data_dir, limit, cfg.restore_min_cell_area)
    print(f"[info] markers={len(mm)} out={cfg.cells_dir} limit={limit} "
          f"min_cell_area={cfg.restore_min_cell_area} threads={cfg.duckdb_threads}")
    cfg.data_dir.mkdir(parents=True, exist_ok=True)  # PARTITION_BY needs parent to exist
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={cfg.duckdb_threads}")
    if memory_limit:
        con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET enable_progress_bar=false")
    con.execute(sql)
    print("[done] wrote", cfg.cells_dir)
    return cfg.cells_dir


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None, help="path to config.ini")
    ap.add_argument("--csv", type=Path, default=None, help="override QuPath cell CSV")
    ap.add_argument("--out", type=Path, default=None, help="override data_dir")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--min-cell-area", type=float, default=None,
                    help="drop cells below this Cell:Area µm^2 (0 to disable)")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--memory-limit", default=None, help="DuckDB memory_limit, e.g. '64GB'")
    a = ap.parse_args(argv)

    overrides = {}
    if a.csv is not None:
        overrides["cells_csv"] = a.csv
    if a.out is not None:
        overrides["data_dir"] = a.out
    if a.min_cell_area is not None:
        overrides["restore_min_cell_area"] = a.min_cell_area
    if a.threads is not None:
        overrides["duckdb_threads"] = a.threads
    cfg = load_config(a.config, **overrides)
    build_cells_parquet(cfg, limit=a.limit, memory_limit=a.memory_limit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
