"""Compact configuration for the content-addressed workflow.

Biological marker meaning lives in ``marker_registry.json`` and
``typing_rules.json``.  Raw input identity lives in the QuPath cohort manifest.
This file therefore contains only output roots, fixed image/QC parameters and
operational compute settings; it does not duplicate marker pairs or cell-type
rules.
"""

from __future__ import annotations

import os
from configparser import ConfigParser
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Optional

from .cohort import filter_eligible_donors


_REPO_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_INI = _REPO_ROOT / "config.ini"
_PACKAGE_DIR = Path(__file__).resolve().parent


def _boolean(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


@dataclass
class PipelineConfig:
    """Resolved workflow paths and scientific/compute parameters."""

    # Source and output contracts.
    data_dir: Path = _REPO_ROOT / "data"
    qupath_manifest: Path = _REPO_ROOT / "data" / "qupath_manifest.json"
    marker_registry: Path = _PACKAGE_DIR / "marker_registry.json"
    typing_rules: Path = _PACKAGE_DIR / "typing_rules.json"

    # Ingest retains tiny fragments only above this structural floor. Geometry
    # QC is a later immutable flag and never deletes rows.
    cells_min_cell_area: float = 5.0

    # Geometry QC.
    cell_qc_min_cell_area: float = 20.0
    cell_qc_max_area_median_multiple: float = 4.0
    cell_qc_min_cell_solidity: float = 0.70
    cell_qc_min_raster_area_ratio: float = 0.50
    cell_qc_max_raster_area_ratio: float = 2.00
    cell_qc_no_cytoplasm_nc_ratio: float = 0.999
    cell_qc_duplicate_radius_um: float = 1.0
    cell_qc_duplicate_area_tol: float = 0.15

    # REDSEA. Compartment output and fail-fast alignment are mandatory in code.
    redsea_downsample: float = 1.0
    redsea_edge_radius: int = 0
    redsea_comp_mode: int = 1
    redsea_alpha: float = 1.0
    redsea_gap_bridge: int = 1
    redsea_norm_form: str = "donor"
    redsea_exclude_no_spillover: bool = True

    # Operational settings.
    n_jobs: int = 1
    duckdb_threads: int = 8
    use_gpu: bool = False
    gpu_device: int = 0

    _config_path: Optional[Path] = field(default=None, repr=False, compare=False)

    @property
    def runs_dir(self) -> Path:
        return self.data_dir / "runs"

    @property
    def cells_dir(self) -> Path:
        return self.data_dir / "cells"

    @property
    def geometry_qc_dir(self) -> Path:
        return self.data_dir / "geometry_qc"

    @property
    def redsea_dir(self) -> Path:
        return self.data_dir / "redsea"

    @property
    def selected_expression_dir(self) -> Path:
        return self.data_dir / "expression"

    @property
    def reference_controls_dir(self) -> Path:
        return self.data_dir / "reference_controls"

    @property
    def marker_evidence_dir(self) -> Path:
        return self.data_dir / "marker_evidence"

    @property
    def assignments_dir(self) -> Path:
        return self.data_dir / "assignments"

    @property
    def state_dir(self) -> Path:
        return self.data_dir / "cell_states"

    @property
    def qupath_class_dir(self) -> Path:
        return self.data_dir / "qupath_export"

    @property
    def audit_dir(self) -> Path:
        return self.data_dir / "audit"

    @property
    def manifests_dir(self) -> Path:
        return self.data_dir / "manifests"

    @property
    def redsea_scratch(self) -> Path:
        return self.data_dir / "scratch" / "redsea"

    @property
    def mask_dir(self) -> Path:
        return self.redsea_scratch / "masks"

    @property
    def inter_dir(self) -> Path:
        return self.redsea_scratch / "intermediates"

    def discover_donors(self, from_dir: Optional[Path] = None) -> list[str]:
        base = Path(from_dir) if from_dir is not None else self.cells_dir
        discovered = sorted(
            path.name.split("=", 1)[1] for path in base.glob("donor_id=*")
        )
        return filter_eligible_donors(discovered)


_INI_SCHEMA = {
    "paths": {
        "data_dir": ("data_dir", Path),
        "qupath_manifest": ("qupath_manifest", Path),
        "marker_registry": ("marker_registry", Path),
        "typing_rules": ("typing_rules", Path),
    },
    "cells": {
        "min_cell_area": ("cells_min_cell_area", float),
    },
    "geometry_qc": {
        "min_cell_area": ("cell_qc_min_cell_area", float),
        "max_area_median_multiple": (
            "cell_qc_max_area_median_multiple",
            float,
        ),
        "min_cell_solidity": ("cell_qc_min_cell_solidity", float),
        "min_raster_area_ratio": ("cell_qc_min_raster_area_ratio", float),
        "max_raster_area_ratio": ("cell_qc_max_raster_area_ratio", float),
        "no_cytoplasm_nc_ratio": ("cell_qc_no_cytoplasm_nc_ratio", float),
        "duplicate_radius_um": ("cell_qc_duplicate_radius_um", float),
        "duplicate_area_tol": ("cell_qc_duplicate_area_tol", float),
    },
    "redsea": {
        "downsample": ("redsea_downsample", float),
        "edge_radius": ("redsea_edge_radius", int),
        "comp_mode": ("redsea_comp_mode", int),
        "alpha": ("redsea_alpha", float),
        "gap_bridge": ("redsea_gap_bridge", int),
        "norm_form": ("redsea_norm_form", str),
        "exclude_no_spillover": ("redsea_exclude_no_spillover", _boolean),
    },
    "compute": {
        "n_jobs": ("n_jobs", int),
        "duckdb_threads": ("duckdb_threads", int),
        "use_gpu": ("use_gpu", _boolean),
        "gpu_device": ("gpu_device", int),
    },
}

_ENV_OVERRIDES = {
    "data_dir": ("PHENOCYCLER_DATA_DIR", Path),
    "qupath_manifest": ("PHENOCYCLER_QUPATH_MANIFEST", Path),
    "n_jobs": ("PHENOCYCLER_JOBS", int),
    "use_gpu": ("PHENOCYCLER_USE_GPU", _boolean),
}


def load_config(
    config_path: Optional[os.PathLike | str] = None,
    **overrides,
) -> PipelineConfig:
    """Load config.ini, environment variables and explicit overrides."""

    cfg = PipelineConfig()
    ini_path = (
        Path(config_path).expanduser()
        if config_path is not None
        else _DEFAULT_INI
    )
    if ini_path.exists():
        ini_path = ini_path.resolve()
        cfg._config_path = ini_path
        parser = ConfigParser(inline_comment_prefixes=("#", ";"))
        parser.read(ini_path)
        for section, mapping in _INI_SCHEMA.items():
            if not parser.has_section(section):
                continue
            for key, (attribute, caster) in mapping.items():
                if parser.has_option(section, key):
                    raw = parser.get(section, key).strip()
                    if raw:
                        setattr(cfg, attribute, caster(raw))

    for attribute, (environment, caster) in _ENV_OVERRIDES.items():
        raw = os.environ.get(environment)
        if raw:
            setattr(cfg, attribute, caster(raw))

    valid = {item.name for item in fields(cfg)}
    for key, value in overrides.items():
        if key not in valid:
            raise TypeError(f"unknown PipelineConfig override: {key!r}")
        setattr(cfg, key, value)

    base = ini_path.parent if ini_path.exists() else _REPO_ROOT
    for name in (
        "data_dir",
        "qupath_manifest",
        "marker_registry",
        "typing_rules",
    ):
        path = Path(getattr(cfg, name)).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
        setattr(cfg, name, path)

    if cfg.redsea_norm_form != "donor":
        raise ValueError("production REDSEA requires mass-conserving norm_form='donor'")
    return cfg


__all__ = ["PipelineConfig", "load_config"]
