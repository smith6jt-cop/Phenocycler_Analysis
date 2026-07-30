"""Create the immutable cohort contract consumed by the workflow.

The preferred input is a small JSON specification.  Paths are resolved relative
to that file, qptiff channel names may be derived directly from image metadata,
and repeated measurement exports are fingerprinted once then reused across
their donor records.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping

from .artifacts import (
    ChannelMapping,
    FileFingerprint,
    QuPathCohortManifest,
    QuPathImageManifest,
)
from .redsea import SANITIZE, qptiff_channels


SPEC_VERSION = 1


def _path(base: Path, value: str | Path) -> Path:
    result = Path(value).expanduser()
    return result.resolve() if result.is_absolute() else (base / result).resolve()


def _measurement_path(base: Path, row: Mapping[str, Any]) -> Path:
    """Resolve one image's source from a single path or shared path list."""

    value = row["measurement_csv"]
    if isinstance(value, (str, Path)):
        if row.get("measurement_csv_index") is not None:
            raise ValueError(
                "measurement_csv_index is only valid when measurement_csv is a list"
            )
        return _path(base, value)
    if not isinstance(value, (list, tuple)) or not value:
        raise ValueError(
            "measurement_csv must be a path or a non-empty list of paths"
        )
    paths = tuple(_path(base, item) for item in value)
    if len(set(paths)) != len(paths):
        raise ValueError("measurement_csv contains duplicate paths")
    index = row.get("measurement_csv_index")
    if isinstance(index, bool) or not isinstance(index, int):
        raise ValueError(
            "images using multiple measurement_csv paths require an integer "
            "measurement_csv_index"
        )
    if not 0 <= index < len(paths):
        raise ValueError(
            f"measurement_csv_index {index} is outside 0..{len(paths) - 1}"
        )
    return paths[index]


def _channels_from_qptiff(path: Path) -> tuple[ChannelMapping, ...]:
    tf, _pages, names, _shape, _pixel_size = qptiff_channels(path, 0)
    tf.close()
    return tuple(
        ChannelMapping(
            marker=SANITIZE(name),
            channel_name=str(name),
            channel_index=index,
        )
        for index, name in enumerate(names)
    )


def _channel_mapping(
    row: Mapping[str, Any], qptiff: Path
) -> tuple[ChannelMapping, ...]:
    observed = _channels_from_qptiff(qptiff)
    raw = row.get("channel_map")
    if raw is None:
        return observed
    declared = tuple(ChannelMapping.from_dict(item) for item in raw)
    if len(declared) != len(observed):
        raise ValueError(
            f"{qptiff}: declared channel map has {len(declared)} entries; "
            f"qptiff has {len(observed)}"
        )
    by_index = {channel.channel_index: channel for channel in declared}
    if set(by_index) != set(range(len(observed))):
        raise ValueError(
            f"{qptiff}: declared channel indices must exactly cover "
            f"0..{len(observed) - 1}"
        )
    for actual in observed:
        declared_channel = by_index[actual.channel_index]
        if declared_channel.channel_name != actual.channel_name:
            raise ValueError(
                f"{qptiff}: channel {actual.channel_index} is "
                f"{actual.channel_name!r}, not declared "
                f"{declared_channel.channel_name!r}"
            )
    return tuple(by_index[index] for index in range(len(observed)))


def create_cohort_manifest(
    specification: str | Path,
    output: str | Path,
    *,
    measurement_hash_mode: str = "full",
    qptiff_hash_mode: str = "sampled",
) -> QuPathCohortManifest:
    """Create and write a cohort manifest from a versioned JSON specification."""

    spec_path = Path(specification).expanduser().resolve()
    payload = json.loads(spec_path.read_text(encoding="utf-8"))
    if int(payload.get("spec_version", 0)) != SPEC_VERSION:
        raise ValueError(f"manifest specification must have spec_version={SPEC_VERSION}")
    rows = payload.get("images")
    if not isinstance(rows, list) or not rows:
        raise ValueError("manifest specification requires a non-empty images list")
    base = spec_path.parent

    measurement_cache: dict[Path, FileFingerprint] = {}
    qptiff_cache: dict[Path, FileFingerprint] = {}
    images: list[QuPathImageManifest] = []
    for raw_row in rows:
        row = {**payload.get("defaults", {}), **dict(raw_row)}
        measurement = _measurement_path(base, row)
        qptiff = _path(base, row["qptiff"])
        cell_geojson = _path(base, row["cell_geojson"])
        nucleus_value = row.get("nucleus_geojson")
        nucleus_geojson = (
            None if nucleus_value in (None, "") else _path(base, nucleus_value)
        )
        if measurement not in measurement_cache:
            measurement_cache[measurement] = FileFingerprint.capture(
                measurement, hash_mode=measurement_hash_mode
            )
        if qptiff not in qptiff_cache:
            qptiff_cache[qptiff] = FileFingerprint.capture(
                qptiff, hash_mode=qptiff_hash_mode
            )
        measurement_fp = measurement_cache[measurement]
        qptiff_fp = qptiff_cache[qptiff]
        images.append(
            QuPathImageManifest.create(
                donor_id=str(row["donor_id"]),
                image_id=str(row["image_id"]),
                measurement_csv=measurement,
                qptiff=qptiff,
                cell_geojson=cell_geojson,
                nucleus_geojson=nucleus_geojson,
                pixel_size_um_x=float(row["pixel_size_um_x"]),
                pixel_size_um_y=float(row["pixel_size_um_y"]),
                panel_id=str(row["panel_id"]),
                channel_map=_channel_mapping(row, qptiff),
                segmentation_version=str(row["segmentation_version"]),
                qptiff_hash_mode=qptiff_hash_mode,
                measurement_fingerprint=measurement_fp,
                qptiff_fingerprint=qptiff_fp,
            )
        )

    expected = payload.get(
        "expected_donors", [image.donor_id for image in images]
    )
    cohort = QuPathCohortManifest.create(
        cohort_name=str(payload["cohort_name"]),
        expected_donors=expected,
        images=images,
        validate_sources=True,
    )
    cohort.write_json(output)
    return cohort


def example_specification() -> dict[str, Any]:
    """Return a concise, self-documenting specification template."""

    return {
        "spec_version": SPEC_VERSION,
        "cohort_name": "islet-phenocycler",
        "expected_donors": ["DONOR_ID"],
        "defaults": {
            "measurement_csv": [
                "/absolute/path/CellMeasurementsBatch1.csv",
                "/absolute/path/CellMeasurementsBatch2.csv",
            ],
            "segmentation_version": "qupath-project-commit-or-model-version",
        },
        "images": [
            {
                "donor_id": "DONOR_ID",
                "image_id": "exact QuPath Image column value",
                "measurement_csv_index": 0,
                "qptiff": "/absolute/path/image.qptiff",
                "cell_geojson": "/absolute/path/cells__image.geojson",
                "nucleus_geojson": None,
                "pixel_size_um_x": 0.325,
                "pixel_size_um_y": 0.325,
                "panel_id": "panel-v1",
                "channel_map": None,
            }
        ],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    create = sub.add_parser("create", help="create a cohort manifest from JSON")
    create.add_argument("--spec", type=Path, required=True)
    create.add_argument("--out", type=Path, required=True)
    create.add_argument(
        "--measurement-hash-mode", choices=("full", "sampled"), default="full"
    )
    create.add_argument(
        "--qptiff-hash-mode", choices=("full", "sampled"), default="sampled"
    )
    template = sub.add_parser("template", help="print a JSON specification template")
    template.add_argument("--out", type=Path)
    args = parser.parse_args(argv)

    if args.command == "create":
        cohort = create_cohort_manifest(
            args.spec,
            args.out,
            measurement_hash_mode=args.measurement_hash_mode,
            qptiff_hash_mode=args.qptiff_hash_mode,
        )
        print(
            f"wrote {args.out} ({len(cohort.images)} images, "
            f"content_id={cohort.content_id})"
        )
        return 0
    text = json.dumps(example_specification(), indent=2) + "\n"
    if args.out is None:
        print(text, end="")
    else:
        if args.out.exists():
            raise FileExistsError(args.out)
        args.out.write_text(text, encoding="utf-8")
        print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
