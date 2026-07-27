#!/usr/bin/env python3
"""Export the cell-QC verdicts as per-image CSVs for the QuPath overlay.

QuPath/Groovy cannot read parquet, so this writes one small CSV per donor keyed by the QuPath
detection UUID (``object_id``) — the same join key REDSEA, RESTORE and the lineage import all use,
so no centroid rounding is involved and the match is exact.

By default only EXCLUDED cells are written, which keeps the CSVs small (a few hundred KB) and makes
the overlay show exactly what the QC removed. ``--include-retained`` writes every cell instead, so
retained cells can be shown in a neutral class and toggled off as background.

    python -m phenocycler.cell_qc_export --all
    python -m phenocycler.cell_qc_export --donor 6450 --donor 6523 --include-retained

Then run ``scripts/groovy/import_cell_qc.groovy`` in QuPath with the image open.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
from pathlib import Path

import pandas as pd

from .cohort import ensure_eligible_donors
from .config import PipelineConfig, load_config

# Kept in sync with import_cell_qc.groovy. Codes are stable so the Groovy colour map does not drift.
REASON_CODE = {
    "no_nucleus": 1,
    "area_below_min": 2,
    "area_above_max": 3,
    "low_solidity": 4,
    "raster_degenerate": 5,
    "raster_empty": 6,
    "no_cytoplasm": 7,   # Tier 2 — retained and still called, but barred from setting thresholds
    "retained": 0,
}


def _donor_image(cfg: PipelineConfig, donor: str) -> str:
    """The QuPath image name for this donor, read from the canonical cells parquet."""
    files = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
    files = [f for f in files if not f.endswith(".preqc.parquet")]
    if not files:
        raise FileNotFoundError(f"donor {donor}: no canonical cells parquet")
    img = pd.read_parquet(files[0], columns=["image"]).image
    names = img.dropna().unique()
    if len(names) != 1:
        raise ValueError(f"donor {donor}: expected one image, found {list(names)}")
    return str(names[0])


def export_donor(
    cfg: PipelineConfig, donor: str, out_dir: Path, *, include_retained: bool = False
) -> Path:
    qc_files = sorted(glob.glob(str(cfg.data_dir / "cell_qc" / f"donor_id={donor}" / "*.parquet")))
    if not qc_files:
        raise FileNotFoundError(
            f"donor {donor}: no QC table — run `python -m phenocycler.cell_qc --all` first"
        )
    qc = pd.read_parquet(qc_files[0], columns=[
        "object_id", "qc_pass", "qc_exclude_reason", "qc_no_cytoplasm",
        "cell_area", "nucleus_area", "cell_solidity", "cell_area_px", "raster_area_ratio",
    ])

    reason = qc.qc_exclude_reason.where(qc.qc_exclude_reason != "", None)
    reason = reason.fillna(qc.qc_no_cytoplasm.map({True: "no_cytoplasm", False: "retained"}))
    qc = qc.assign(qc_reason=reason, qc_code=reason.map(REASON_CODE).astype(int))

    if not include_retained:
        qc = qc[qc.qc_code != REASON_CODE["retained"]]

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"cell_qc_{donor}.csv"
    cols = ["object_id", "qc_code", "qc_reason", "cell_area", "nucleus_area",
            "cell_solidity", "cell_area_px", "raster_area_ratio"]
    out = qc[cols].copy()
    out.insert(1, "image", _donor_image(cfg, donor))
    out.round({"cell_area": 2, "nucleus_area": 2, "cell_solidity": 4,
               "cell_area_px": 1, "raster_area_ratio": 4}).to_csv(path, index=False)

    counts = qc.qc_reason.value_counts().to_dict()
    print(f"[{donor}] {len(out):,} rows -> {path.name} | "
          + ", ".join(f"{k}={v:,}" for k, v in sorted(counts.items())), flush=True)
    return path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Export cell-QC verdicts as QuPath-importable CSVs.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--donor", action="append", default=None)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--include-retained", action="store_true",
                    help="also write passing cells (code 0) so they can be toggled as background")
    ap.add_argument("--output", type=Path, default=None)
    a = ap.parse_args(argv)

    cfg = load_config(a.config)
    donors = list(a.donor) if a.donor else (cfg.discover_donors() if a.all else [])
    if not donors:
        ap.error("pass --donor <id> (repeatable) or --all")
    ensure_eligible_donors(donors)

    out_dir = a.output or (cfg.data_dir / "cell_qc" / "qupath_review")
    written = [export_donor(cfg, d, out_dir, include_retained=a.include_retained) for d in donors]

    manifest = {
        "reason_codes": REASON_CODE,
        "include_retained": bool(a.include_retained),
        "files": [p.name for p in written],
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    # Render the importer with this directory baked in, mirroring how the RESTORE review bundle
    # ships a ready-to-run copy rather than a template the reviewer has to edit.
    template = Path(__file__).resolve().parents[1] / "scripts" / "groovy" / "import_cell_qc.groovy"
    if template.exists():
        rendered = template.read_text().replace("__CELL_QC_CSV_DIR__", str(out_dir))
        (out_dir / "import_cell_qc.groovy").write_text(rendered)
        print(f"\nrunnable importer -> {out_dir / 'import_cell_qc.groovy'}")
    print(f"open the donor's image in QuPath, then run that script.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
