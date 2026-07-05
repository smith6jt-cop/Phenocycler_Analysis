"""
Central configuration for the phenocycler pipeline.

Everything that was a hardcoded absolute path or a magic constant in the
upstream ``scripts/senior/*.py`` scripts lives here instead, resolved from
``config.ini`` (see the repo root) with environment-variable and keyword
overrides.  The scientific defaults (REDSEA subtract-only / α=1 / 1-px band,
RESTORE SSC model + robust guard, the 6 broad lineages) are preserved exactly;
this module only makes the *paths* and *compute knobs* configurable.

Resolution order for every value: explicit keyword override  >  environment
variable  >  ``config.ini``  >  built-in default.
"""

from __future__ import annotations

import os
from configparser import ConfigParser
from dataclasses import dataclass, field, fields
from pathlib import Path
from typing import Optional

# --------------------------------------------------------------------------- #
# Scientific constants (faithful to Islet-Explorer-Senior; not usually tuned)
# --------------------------------------------------------------------------- #

# Broad-lineage gating markers -> the six mutually-exclusive lineages
# (scripts/senior/assign_broad_lineage.py).
LINEAGES: dict[str, list[str]] = {
    "Epithelial": ["Pan_Cytokeratin"],
    "Fibroblast": ["Vimentin"],
    "Immune": ["CD3e", "CD20", "CD163"],
    "Endocrine": ["INS", "GCG", "SST"],
    "Endothelial": ["CD31"],
    "Muscle": ["SMA"],
}
# Structural background is resolved by argmax of these three `_norm` scores.
STRUCT_LINEAGES: list[str] = ["Epithelial", "Fibroblast", "Muscle"]

LINEAGE_COLORS: dict[str, str] = {
    "Epithelial": "#4477AA", "Fibroblast": "#EE6677", "Immune": "#228833",
    "Endocrine": "#CCBB44", "Endothelial": "#66CCEE", "Muscle": "#AA3377",
}
STATUS_ORDER: list[str] = ["ND", "AAB", "T1D"]

# RESTORE mutually-exclusive [target, reference] pairs
# (scripts/senior/restore_normalize.py DEFAULT_MARKER_PAIRS).  Pan_Cytokeratin
# is the universal negative control; CD3e<-CD163 is the documented exception.
DEFAULT_MARKER_PAIRS: list[list[str]] = [
    ["Pan_Cytokeratin", "Vimentin"],
    ["Vimentin", "Pan_Cytokeratin"],
    ["INS", "Pan_Cytokeratin"],
    ["GCG", "Pan_Cytokeratin"],
    ["SST", "Pan_Cytokeratin"],
    ["CD31", "Pan_Cytokeratin"],
    ["SMA", "Pan_Cytokeratin"],
    ["CD3e", "CD163"],
    ["CD20", "Pan_Cytokeratin"],
    ["CD163", "Pan_Cytokeratin"],
]

_REPO_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_INI = _REPO_ROOT / "config.ini"


@dataclass
class PipelineConfig:
    """Resolved paths + compute knobs for a pipeline run."""

    # -- roots (absolute; usually machine-specific -> set in config.ini) ------
    data_dir: Path = _REPO_ROOT / "data"
    images_dir: Path = _REPO_ROOT / "data" / "raw" / "images"
    cells_csv: Path = _REPO_ROOT / "data" / "raw" / "Cellmeasurements.csv"
    donor_metadata: Path = _REPO_ROOT / "data" / "raw" / "donor_metadata_panc.xlsx"
    restore_vendor: Path = _REPO_ROOT / "external" / "RESTORE" / "python_code"

    # -- metadata schema -----------------------------------------------------
    metadata_donor_col: str = "donor_id"
    metadata_status_col: str = "disease.status"

    # -- REDSEA (scripts/senior/redsea_full.py defaults) ---------------------
    redsea_downsample: float = 1.0
    redsea_edge_radius: int = 0      # 0 == the 1-px cell rim
    redsea_comp_mode: int = 0        # 0 == subtract-only
    redsea_alpha: float = 1.0
    redsea_gap_bridge: int = 1

    # -- RESTORE (scripts/senior/restore_normalize.py defaults) --------------
    restore_model: str = "SSC"       # SSC | GMM | KMeans
    restore_subsample: int = 15000
    restore_robust: bool = True
    restore_robust_factor: float = 3.0
    restore_min_cell_area: float = 5.0
    restore_seed: int = 0

    # -- compute -------------------------------------------------------------
    n_jobs: int = 1                  # per-donor process pool size (1 == serial)
    duckdb_threads: int = 8
    use_gpu: bool = False            # opt-in CuPy backend for REDSEA
    gpu_device: int = 0

    # cached derived-path store (not a config field)
    _config_path: Optional[Path] = field(default=None, repr=False, compare=False)

    # ---- derived pipeline directories (hive-partitioned donor_id=* layout) --
    @property
    def cells_dir(self) -> Path:
        return self.data_dir / "cells"

    @property
    def cells_redsea_dir(self) -> Path:
        return self.data_dir / "cells_redsea"

    @property
    def redsea_scratch(self) -> Path:
        return self.data_dir / "redsea_scratch"

    @property
    def geojson_dir(self) -> Path:
        return self.redsea_scratch / "geojson"

    @property
    def mask_dir(self) -> Path:
        return self.redsea_scratch / "masks"

    @property
    def inter_dir(self) -> Path:
        return self.redsea_scratch / "intermediates"

    @property
    def restore_dir(self) -> Path:
        return self.data_dir / "restore_redsea"

    @property
    def restore_gated_dir(self) -> Path:
        return self.data_dir / "restore_gated_redsea"

    @property
    def restore_thresholds_csv(self) -> Path:
        return self.data_dir / "restore_thresholds_redsea.csv"

    @property
    def phenotype_dir(self) -> Path:
        return self.data_dir / "phenotype"

    @property
    def broad_dir(self) -> Path:
        return self.phenotype_dir / "broad"

    @property
    def qupath_class_dir(self) -> Path:
        return self.phenotype_dir / "qupath_class"

    # ---- helpers -----------------------------------------------------------
    def discover_donors(self, from_dir: Optional[Path] = None) -> list[str]:
        """Donor ids by globbing ``<dir>/donor_id=*`` (defaults to cells_dir)."""
        base = Path(from_dir) if from_dir is not None else self.cells_dir
        return sorted(p.name.split("=", 1)[1] for p in base.glob("donor_id=*"))

    def ensure_dirs(self) -> None:
        for d in (self.data_dir, self.geojson_dir, self.mask_dir, self.inter_dir):
            d.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# Loading
# --------------------------------------------------------------------------- #

# config.ini [section] -> {ini_key: (attr_name, caster)}
_INI_SCHEMA = {
    "paths": {
        "data_dir": ("data_dir", Path),
        "images_dir": ("images_dir", Path),
        "cells_csv": ("cells_csv", Path),
        "donor_metadata": ("donor_metadata", Path),
        "restore_vendor": ("restore_vendor", Path),
    },
    "metadata": {
        "donor_col": ("metadata_donor_col", str),
        "status_col": ("metadata_status_col", str),
    },
    "redsea": {
        "downsample": ("redsea_downsample", float),
        "edge_radius": ("redsea_edge_radius", int),
        "comp_mode": ("redsea_comp_mode", int),
        "alpha": ("redsea_alpha", float),
        "gap_bridge": ("redsea_gap_bridge", int),
    },
    "restore": {
        "model": ("restore_model", str),
        "subsample": ("restore_subsample", int),
        "robust": ("restore_robust", lambda s: str(s).lower() in ("1", "true", "yes", "on")),
        "robust_factor": ("restore_robust_factor", float),
        "min_cell_area": ("restore_min_cell_area", float),
        "seed": ("restore_seed", int),
    },
    "compute": {
        "n_jobs": ("n_jobs", int),
        "duckdb_threads": ("duckdb_threads", int),
        "use_gpu": ("use_gpu", lambda s: str(s).lower() in ("1", "true", "yes", "on")),
        "gpu_device": ("gpu_device", int),
    },
}

# attr_name -> environment variable that overrides it
_ENV_OVERRIDES = {
    "data_dir": ("PHENOCYCLER_DATA_DIR", Path),
    "images_dir": ("PHENOCYCLER_IMAGES_DIR", Path),
    "cells_csv": ("PHENOCYCLER_CELLS_CSV", Path),
    "donor_metadata": ("PHENOCYCLER_DONOR_METADATA", Path),
    "restore_vendor": ("PHENOCYCLER_RESTORE_VENDOR", Path),
    "n_jobs": ("PHENOCYCLER_JOBS", int),
    "use_gpu": ("PHENOCYCLER_USE_GPU", lambda s: str(s).lower() in ("1", "true", "yes", "on")),
}


def load_config(config_path: Optional[os.PathLike | str] = None, **overrides) -> PipelineConfig:
    """Build a :class:`PipelineConfig`.

    Parameters
    ----------
    config_path
        Path to a ``config.ini``.  Defaults to ``<repo>/config.ini`` if it
        exists; missing file is fine (built-in defaults are used).
    **overrides
        Keyword overrides that win over everything (values are used verbatim;
        pass ``Path`` objects for path fields).
    """
    cfg = PipelineConfig()
    ini_path = Path(config_path) if config_path is not None else _DEFAULT_INI

    if ini_path.exists():
        cfg._config_path = ini_path
        parser = ConfigParser(inline_comment_prefixes=("#", ";"))
        parser.read(ini_path)
        for section, mapping in _INI_SCHEMA.items():
            if not parser.has_section(section):
                continue
            for key, (attr, caster) in mapping.items():
                if parser.has_option(section, key):
                    raw = parser.get(section, key).strip()
                    if raw != "":
                        setattr(cfg, attr, caster(raw))

    for attr, (env, caster) in _ENV_OVERRIDES.items():
        raw = os.environ.get(env)
        if raw:
            setattr(cfg, attr, caster(raw))

    valid = {f.name for f in fields(cfg)}
    for key, val in overrides.items():
        if key not in valid:
            raise TypeError(f"unknown PipelineConfig override: {key!r}")
        setattr(cfg, key, val)

    # Normalize path fields to absolute Paths. Relative paths resolve against the
    # config file's directory (the repo root) — not the cwd — so the pipeline works
    # from notebooks/, SLURM jobs, or anywhere the package is imported.
    base = ini_path.parent if ini_path.exists() else _REPO_ROOT
    for name in ("data_dir", "images_dir", "cells_csv", "donor_metadata", "restore_vendor"):
        p = Path(getattr(cfg, name)).expanduser()
        if not p.is_absolute():
            p = (base / p).resolve()
        setattr(cfg, name, p)

    return cfg
