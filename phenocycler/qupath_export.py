#!/usr/bin/env python3
"""
Step 6 (optional QC) — export the broad lineage as per-image CSVs keyed by
QuPath's detection UUID, for the companion ``scripts/groovy/import_broad_lineage.groovy``
(UUID match -> set PathClass directly).

Faithful port of ``scripts/senior/export_broad_class_for_qupath.py``.  QuPath's
``PathObject.getID()`` is exactly our ``object_id``, so matching is exact (no
centroid rounding) and the class is set in one step.

    python -m phenocycler.qupath_export
    -> <phenotype_dir>/qupath_class/pheno_class_<donor>.csv   (object_id, broad_lineage, image)
"""

from __future__ import annotations

import argparse
import glob
from pathlib import Path

import pandas as pd

from .config import PipelineConfig, load_config


def export_qupath_classes(cfg: PipelineConfig) -> Path:
    """Write per-donor pheno_class CSVs for the QuPath importer. Returns the output dir."""
    out = cfg.qupath_class_dir
    out.mkdir(parents=True, exist_ok=True)
    total = 0
    for bf in sorted(glob.glob(str(cfg.broad_dir / "donor_id=*" / "*.parquet"))):
        donor = bf.split("donor_id=")[1].split("/")[0]
        b = pd.read_parquet(bf, columns=["object_id", "broad_lineage"])
        cells = sorted(glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet")))
        if not cells:
            print(f"[{donor}] WARN no cells parquet; skipped")
            continue
        img = pd.read_parquet(cells[0], columns=["object_id", "image"])
        merged = b.merge(img.astype({"object_id": str}), on="object_id", how="inner")
        merged = merged[["object_id", "broad_lineage", "image"]]
        csv = out / f"pheno_class_{donor}.csv"
        merged.to_csv(csv, index=False)
        total += len(merged)
        print(f"[{donor}] {len(merged):,} cells -> {csv.name}", flush=True)
    print(f"\n[done] {total:,} cells across "
          f"{len(glob.glob(str(out/'pheno_class_*.csv')))} images -> {out}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, default=None)
    a = ap.parse_args(argv)
    export_qupath_classes(load_config(a.config))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
