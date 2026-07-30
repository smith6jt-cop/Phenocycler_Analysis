#!/usr/bin/env python3
"""Export hierarchical labels and their uncertainty to QuPath-keyed CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

EXPORT_COLUMNS = (
    "object_id",
    "image",
    "broad_type",
    "specific_type",
    "assignment_status",
    "confidence",
    "best_broad_type",
    "best_specific_type",
    "broad_assignment_status",
    "specific_assignment_status",
    "ranked_broad_probabilities",
    "ranked_specific_probabilities",
    "ranked_probabilities",
)


def export_assignment_frame(
    assignments: pd.DataFrame,
    cell_context: pd.DataFrame,
    out_dir: Path,
) -> list[Path]:
    """Write the current hierarchical assignment schema, one CSV per image.

    The full uncertainty contract is exported instead of collapsing ambiguous or unavailable cells
    into a named QuPath class. ``compartment`` and ``cell_type`` compatibility columns are included
    for the existing importer, but always mirror the authoritative broad/specific fields.
    """
    required = {
        "object_id",
        "donor_id",
        "broad_type",
        "specific_type",
        "assignment_status",
        "confidence",
        "best_broad_type",
        "best_specific_type",
        "broad_assignment_status",
        "specific_assignment_status",
        "ranked_broad_probabilities",
        "ranked_specific_probabilities",
        "ranked_probabilities",
    }
    missing = required - set(assignments)
    if missing:
        raise KeyError(f"assignments missing required columns: {sorted(missing)}")
    if not {"object_id", "image"} <= set(cell_context):
        raise KeyError("cell_context must contain object_id and image")
    left = assignments.copy()
    right = cell_context[["object_id", "image"]].copy()
    left["object_id"] = left["object_id"].astype(str)
    right["object_id"] = right["object_id"].astype(str)
    if left["object_id"].duplicated().any() or right["object_id"].duplicated().any():
        raise ValueError("object_id must be unique in assignments and cell_context")
    merged = left.merge(right, on="object_id", how="left", validate="one_to_one")
    if merged["image"].isna().any():
        raise ValueError(f"{int(merged['image'].isna().sum())} assignments have no image context")
    merged["compartment"] = merged["broad_type"]
    merged["cell_type"] = merged["specific_type"]

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for (donor, image), frame in merged.groupby(
        ["donor_id", "image"], sort=True, dropna=False
    ):
        path = out_dir / f"cell_assignments_{donor}.csv"
        if path.exists():
            raise FileExistsError(f"immutable QuPath export already exists: {path}")
        columns = list(EXPORT_COLUMNS) + ["compartment", "cell_type"]
        frame[columns].to_csv(path, index=False)
        paths.append(path)
    return paths


def main(argv=None):
    """Run the manifest-controlled QuPath export."""

    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args(argv)
    command = ["export"]
    if args.config is not None:
        command.extend(["--config", str(args.config)])
    from .pipeline import main as pipeline_main

    return pipeline_main(command)


if __name__ == "__main__":
    raise SystemExit(main())
