#!/usr/bin/env python3
"""Stream contracted QuPath measurements into canonical donor partitions.

Input  : measurement CSVs named by the QuPath cohort manifest. Panels differ
         between batches, so ingest reconciles headers and keeps the marker
         UNION plus availability.
Output : <run>/cells/donor_id=<id>/*.parquet with exact keys, geometry-QC
         inputs, spatial region and one raw Cell/Nucleus source per marker.

Region scheme: Parent ``Islet_N`` = core, ``Islet_N_exp20um`` = peri,
``Annotation (Tissue)`` / ``Root object`` = tissue, else other.

DuckDB streams large CSVs out of core with ``PARTITION_BY donor_id`` and ZSTD.
Production execution is orchestrated by :mod:`phenocycler.pipeline`.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Mapping, Optional, Sequence

import pandas as pd

from .artifacts import ContractError, QuPathCohortManifest, image_id_to_donor_id
from .config import PipelineConfig


REQUIRED_CONTEXT_COLUMNS = (
    "Image",
    "Object ID",
    "Centroid X µm",
    "Centroid Y µm",
    "Parent",
    "Cell: Area µm^2",
    "Cell: Solidity",
)


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


def compartment_means(cols: list[str], compartment: str) -> dict[str, str]:
    """Map one QuPath compartment's marker means to canonical marker names."""

    pattern = re.compile(rf"^{re.escape(compartment)}: (.+): Mean$")
    out = {}
    for column in cols:
        match = pattern.match(column)
        if match:
            out[column] = sanitize(match.group(1))
    return out


def sanitize(name: str) -> str:
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def reconcile_headers(csv_paths: Sequence[Path]) -> tuple[dict[str, str], dict]:
    """Return the UNION of marker means plus explicit per-export availability.

    A batch-specific marker remains scientifically usable for donors acquired with that panel. The
    old intersection silently discarded that information cohort-wide. Missing measurements are now
    represented as null and accompanied by an availability table; downstream donor-local stages must
    emit ``unavailable`` rather than treating null as negative.
    """
    if not csv_paths:
        raise SystemExit(
            "[err] the cohort manifest declares no QuPath measurement CSV"
        )
    missing = [str(p) for p in csv_paths if not Path(p).exists()]
    if missing:
        raise SystemExit("[err] QuPath cell CSV not found:\n       " + "\n       ".join(missing))

    per_file = {}
    for p in csv_paths:
        columns = header_cols(Path(p))
        absent = sorted(set(REQUIRED_CONTEXT_COLUMNS).difference(columns))
        if absent:
            raise SystemExit(
                f"[err] QuPath export {p} is missing required canonical "
                f"columns: {absent}"
            )
        mm = marker_means(columns)
        if not mm:
            raise SystemExit(f"No 'Cell: <marker>: Mean' columns found in {p} — check the header.")
        per_file[Path(p)] = mm

    # Key on the SANITIZED marker name: the same marker must map to the same output column in every
    # file, and it is the sanitized name every downstream reader uses.
    by_sanitized = {p: {sanitize(v): k for k, v in mm.items()} for p, mm in per_file.items()}
    union = set.union(*(set(d) for d in by_sanitized.values()))
    if not union:
        raise SystemExit("[err] the QuPath exports share no marker means — check they are the same panel")

    union_mm: dict[str, str] = {}
    for marker in sorted(union):
        originals = {mapping[marker] for mapping in by_sanitized.values() if marker in mapping}
        if len(originals) != 1:
            raise SystemExit(
                f"[err] marker {marker!r} resolves from inconsistent QuPath columns "
                f"{sorted(originals)}; add an explicit alias before ingest"
            )
        union_mm[originals.pop()] = marker

    report = {
        "files": {str(p): len(mm) for p, mm in per_file.items()},
        "union": len(union),
        "intersection": len(set.intersection(*(set(d) for d in by_sanitized.values()))),
        "availability": {
            str(p): {marker: marker in mapping for marker in sorted(union)}
            for p, mapping in by_sanitized.items()
        },
    }
    return union_mm, report


def reconcile_compartment_headers(
    csv_paths: Sequence[Path],
    compartment: str,
) -> dict[str, str]:
    """Return a union of one raw QuPath compartment across acquisition batches."""

    per_file = {
        Path(path): compartment_means(header_cols(Path(path)), compartment)
        for path in csv_paths
    }
    by_sanitized = {
        path: {marker: original for original, marker in measurements.items()}
        for path, measurements in per_file.items()
    }
    union = set().union(*(set(mapping) for mapping in by_sanitized.values()))
    result: dict[str, str] = {}
    for marker in sorted(union):
        originals = {
            mapping[marker] for mapping in by_sanitized.values() if marker in mapping
        }
        if len(originals) != 1:
            raise SystemExit(
                f"[err] {compartment} marker {marker!r} resolves from inconsistent "
                f"QuPath columns {sorted(originals)}; add an explicit alias before ingest"
            )
        result[originals.pop()] = f"{compartment}__{marker}"
    return result


def _sql_string(value: object) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def donor_case_expression(image_to_donor: Mapping[str, str]) -> str:
    """Exact QuPath Image -> donor SQL expression; no filename inference."""

    mapping = {str(image): str(donor) for image, donor in image_to_donor.items()}
    if not mapping:
        raise ContractError("image_to_donor cannot be empty")
    branches = " ".join(
        f"WHEN {_sql_string(image)} THEN {_sql_string(mapping[image])}"
        for image in sorted(mapping)
    )
    return (
        f'CASE "Image" {branches} ELSE '
        "error('unmapped QuPath Image: ' || coalesce(\"Image\", '<NULL>')) END"
    )


def contracted_image_source(
    csv_paths: Sequence[Path], image_to_donor: Mapping[str, str]
) -> str:
    """Read shared exports but expose only exact manifest-contracted images."""

    contracted_images = sorted(map(str, image_to_donor))
    if not contracted_images:
        raise ContractError("image_to_donor cannot be empty")
    raw_parts = []
    for index, path in enumerate(csv_paths):
        escaped = str(path).replace("'", "''")
        raw_parts.append(
            f"SELECT * FROM read_csv('{escaped}', header=true, all_varchar=true, "
            f"ignore_errors=true, store_rejects=true, "
            f"rejects_table='phenocycler_ingest_rejects_{index}', "
            f"rejects_scan='phenocycler_ingest_reject_scans_{index}', rejects_limit=0)"
        )
    combined = "(" + " UNION ALL BY NAME ".join(raw_parts) + ")"
    allowed = ", ".join(_sql_string(image) for image in contracted_images)
    return (
        f'(SELECT * FROM {combined} AS all_qupath_exports '
        f'WHERE "Image" IN ({allowed})) AS qupath_exports'
    )


def build_sql(
    csv_paths: Sequence[Path],
    out_dir: Path,
    limit: Optional[int],
    min_cell_area: float,
    *,
    image_to_donor: Mapping[str, str],
) -> tuple[str, dict]:
    csv_paths = [Path(p) for p in ([csv_paths] if isinstance(csv_paths, (str, Path)) else csv_paths)]
    mm, _report = reconcile_headers(csv_paths)
    nucleus_mm = reconcile_compartment_headers(csv_paths, "Nucleus")
    selects = [
        f"{donor_case_expression(image_to_donor)} AS donor_id",
        '"Image" AS image',
        '"Object ID" AS object_id',
        "TRY_CAST(\"Centroid X µm\" AS DOUBLE) AS X_centroid",
        "TRY_CAST(\"Centroid Y µm\" AS DOUBLE) AS Y_centroid",
        '"Parent" AS parent_raw',
        ("CASE WHEN \"Parent\" LIKE '%exp20um' THEN 'peri' "
         "WHEN regexp_matches(\"Parent\", '^Islet_[0-9]+$') THEN 'core' "
         "WHEN \"Parent\" = 'Annotation (Tissue)' OR \"Parent\" LIKE 'Root object%' THEN 'tissue' "
         "ELSE 'other' END AS cell_region"),
        "regexp_extract(\"Parent\", 'Islet_([0-9]+)', 1) AS islet_num",
        "TRY_CAST(\"Cell: Area µm^2\" AS DOUBLE) AS cell_area",
        "TRY_CAST(\"Cell: Solidity\" AS DOUBLE) AS cell_solidity",
    ]
    for orig, marker in mm.items():        # marker is already the sanitized name (reconcile_headers)
        selects.append(f'TRY_CAST("{orig}" AS DOUBLE) AS "{marker}"')
    # Nuclear marker means are the selected raw source for nuclear/state markers.
    # Keeping only Cell + Nucleus avoids carrying four competing raw versions of
    # every marker while preserving every source required by marker_registry.
    for orig, output_column in nucleus_mm.items():
        selects.append(f'TRY_CAST("{orig}" AS DOUBLE) AS "{output_column}"')
    sel = ",\n  ".join(selects)
    # Align batched exports by column NAME, then apply the cohort manifest's exact
    # image allow-list. A shared QuPath export may contain deliberately excluded
    # or otherwise out-of-cohort images; those rows must not reach donor mapping.
    # The CSV reader still scans the batched files because CSV has no image index.
    src = contracted_image_source(csv_paths, image_to_donor)
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
    """Run the DuckDB CSV -> partitioned-parquet build. Returns the cells dir.

    Scans each unique measurement export declared by ``cfg.qupath_manifest``
    once, preserving the panel union while writing only contracted images."""
    import duckdb

    cohort = QuPathCohortManifest.read_json(cfg.qupath_manifest)
    cohort.validate_current(mode="fast")
    contracted_csvs = {image.measurement_csv.path for image in cohort.images}
    csv_paths = tuple(Path(path) for path in sorted(contracted_csvs))
    donor_mapping = image_id_to_donor_id(cohort)
    mm, report = reconcile_headers(csv_paths)
    for path, n in report["files"].items():
        print(f"[info] {Path(path).name}: {n} marker means")
    print(
        f"[info] panel union={report['union']} intersection={report['intersection']}; "
        "batch-specific markers are retained with explicit null availability"
    )
    report["cohort_manifest_content_id"] = cohort.content_id
    report["images"] = {
        image.image_id: {
            "donor_id": image.donor_id,
            "panel_id": image.panel_id,
            "markers": [channel.marker for channel in image.channel_map],
        }
        for image in cohort.images
    }
    sql, mm = build_sql(
        csv_paths,
        cfg.data_dir,
        limit,
        cfg.cells_min_cell_area,
        image_to_donor=donor_mapping,
    )
    print(f"[info] files={len(csv_paths)} markers={len(mm)} out={cfg.cells_dir} limit={limit} "
          f"min_cell_area={cfg.cells_min_cell_area} threads={cfg.duckdb_threads}")
    cfg.data_dir.mkdir(parents=True, exist_ok=True)  # PARTITION_BY needs parent to exist
    con = duckdb.connect()
    con.execute(f"PRAGMA threads={cfg.duckdb_threads}")
    if memory_limit:
        con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET enable_progress_bar=false")
    if cfg.cells_dir.exists() and any(cfg.cells_dir.iterdir()):
        raise FileExistsError(
            f"{cfg.cells_dir} is not empty; ingest is immutable and will not merge with a prior run"
        )
    con.execute(sql)
    reject_frames = [
        con.execute(f"SELECT * FROM phenocycler_ingest_rejects_{index}").fetchdf()
        for index in range(len(csv_paths))
    ]
    rejects = pd.concat(reject_frames, ignore_index=True) if reject_frames else pd.DataFrame()
    rejects.to_parquet(cfg.data_dir / "ingest_rejects.parquet", index=False)
    (cfg.data_dir / "panel_availability.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    if len(rejects):
        raise ValueError(
            f"QuPath CSV parsing rejected {len(rejects):,} row(s); inspect "
            f"{cfg.data_dir / 'ingest_rejects.parquet'}"
        )
    print("[done] wrote", cfg.cells_dir)
    return cfg.cells_dir


def main(argv=None) -> int:
    """Run only the manifest-controlled ingest stage."""

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args(argv)
    command = ["run"]
    if args.config is not None:
        command.extend(["--config", str(args.config)])
    command.extend(["--only", "ingest", "--no-export"])
    from .pipeline import main as pipeline_main

    return pipeline_main(command)


if __name__ == "__main__":
    raise SystemExit(main())
