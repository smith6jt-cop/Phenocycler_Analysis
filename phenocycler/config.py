"""
Central configuration for the phenocycler pipeline.

Everything that was a hardcoded absolute path or a magic constant in the
upstream ``scripts/senior/*.py`` scripts lives here instead, resolved from
``config.ini`` (see the repo root) with environment-variable and keyword
overrides.  The scientific defaults (REDSEA subtract-only / α=1 / 1-px band,
RESTORE SSC model + robust guard, the 8 broad lineages) are preserved exactly;
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

# Broad-lineage gating markers -> the eight mutually-exclusive lineages
# (scripts/senior/assign_broad_lineage.py).  CD99 (bright-only) is an Endocrine
# gate; Neural (B3TUBB) and Neutrophil (MPO) are carved out of the structural /
# immune background.  CD99/B3TUBB/MPO are gated in a SEPARATE RESTORE pass
# (EXTRA_MARKER_PAIRS) and merged by object_id at assignment time.
LINEAGES: dict[str, list[str]] = {
    "Epithelial": ["Pan_Cytokeratin"],
    "Fibroblast": ["Vimentin"],
    "Immune": ["CD3e", "CD20", "CD163"],
    "Endocrine": ["INS", "GCG", "SST", "CD99"],
    "Endothelial": ["CD31"],
    "Muscle": ["SMA"],
    "Neural": ["B3TUBB"],
    "Neutrophil": ["MPO"],
}
# Structural background is resolved by argmax of these three `_norm` scores.
STRUCT_LINEAGES: list[str] = ["Epithelial", "Fibroblast", "Muscle"]

LINEAGE_COLORS: dict[str, str] = {
    "Epithelial": "#4477AA", "Fibroblast": "#EE6677", "Immune": "#228833",
    "Endocrine": "#CCBB44", "Endothelial": "#66CCEE", "Muscle": "#AA3377",
    "Neural": "#B5838D", "Neutrophil": "#E69F00",   # dusty rose / Okabe-Ito orange (CVD-safe)
}
# unique 3-char codes for the terse per-donor progress line (Endocrine/Endothelial
# and Neural/Neutrophil would otherwise both collide under a naive name[:3]).
LINEAGE_ABBR: dict[str, str] = {
    "Epithelial": "Epi", "Fibroblast": "Fib", "Immune": "Imm", "Endocrine": "Enc",
    "Endothelial": "Eth", "Muscle": "Mus", "Neural": "Nrl", "Neutrophil": "Neu",
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

# Extra RESTORE pass (run separately with --no-robust --no-ref-qc so the validated
# 10-pair gates in restore_gated_redsea stay byte-identical).  These three markers
# (Endocrine-CD99 / Neural-B3TUBB / Neutrophil-MPO) are gated against Pan_Cytokeratin
# into restore_gated_redsea_extra and merged into the lineage call by object_id.
EXTRA_MARKER_PAIRS: list[list[str]] = [
    ["CD99", "Pan_Cytokeratin"],
    ["B3TUBB", "Pan_Cytokeratin"],
    ["MPO", "Pan_Cytokeratin"],
]
EXTRA_MARKERS: list[str] = [p[0] for p in EXTRA_MARKER_PAIRS]   # CD99, B3TUBB, MPO

# {INS,GCG,SST}_pos are floored to _norm >= hormone_min_norm by the hormone_floor
# stage (scripts/senior/apply_hormone_floor.py) before lineage assignment.
HORMONE_MARKERS: list[str] = ["INS", "GCG", "SST"]

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

    # -- lineage (assign_broad_lineage.py + apply_hormone_floor.py + marker_taxonomy.py) --
    hormone_min_norm: float = 5.0    # K: {INS,GCG,SST}_pos := _norm >= K (false-endocrine floor)
    cd99_bright: float = 3.0         # CD99_pos := CD99_norm >= this (bright-only Endocrine gate)

    # -- compute -------------------------------------------------------------
    n_jobs: int = 1                  # per-donor process pool size (1 == serial)
    duckdb_threads: int = 8
    use_gpu: bool = False            # opt-in CuPy backend for REDSEA
    gpu_device: int = 0

    # -- integration (PhenoCycler <-> Xenium; phenocycler/integration/) -------
    # Nothing in steps 1-7 reads these; the core pipeline is unaffected.
    integration_mode: str = "sequential"          # sequential | same_slide
    xenium_paths_csv: Path = _REPO_ROOT / "data" / "integration" / "xenium_paths.csv"
    donor_overrides_csv: Path = _REPO_ROOT / "data" / "integration" / "donor_overrides.csv"
    xenium_root: str = ""                         # optional prefix rewrite for bundle paths
    panel_explorer: Path = _REPO_ROOT / "external" / "XeniumPanelExplorer"
    tissue: str = "pancreas"                      # selects the XeniumPanelExplorer tissue dir
    fixed_modality: str = "phenocycler"           # which frame registration targets
    reg_pixel_um: float = 2.0                     # raster resolution for registration
    reg_nonrigid: bool = True
    reg_max_disp_um: float = 200.0                # displacement cap (prevents folding)
    islet_eps_um: float = 50.0                    # = insulitis_analysis.EPS_UM
    islet_min_samples: int = 10                   # = insulitis_analysis.MIN_SAMPLES
    niche_k: int = 50                             # = nb03 K_COMP
    niche_n: int = 12                             # = nb03 n_niches
    niche_smooth_k: int = 30                      # = nb03 K_SMOOTH
    grid_um: float = 100.0
    match_max_dist_um: float = 200.0
    match_area_ratio: float = 2.5
    crossmodal_min_anchors: int = 8               # refuse pseudo-cell linking below this
    qc_tissue_dice_min: float = 0.80
    qc_islet_rmse_max_um: float = 150.0

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
    def restore_gated_prefloor_dir(self) -> Path:
        # pre-hormone-floor backup; the hormone_floor stage reads this and (re)writes restore_gated_dir
        return self.data_dir / "restore_gated_redsea.pre_hormonefloor"

    @property
    def restore_redsea_extra_dir(self) -> Path:
        return self.data_dir / "restore_redsea_extra"

    @property
    def restore_gated_extra_dir(self) -> Path:
        return self.data_dir / "restore_gated_redsea_extra"

    @property
    def restore_thresholds_csv(self) -> Path:
        return self.data_dir / "restore_thresholds_redsea.csv"

    @property
    def restore_thresholds_extra_csv(self) -> Path:
        return self.data_dir / "restore_thresholds_extra.csv"

    @property
    def redsea_reassess_dir(self) -> Path:
        return self.data_dir / "redsea_reassess"

    @property
    def phenotype_dir(self) -> Path:
        return self.data_dir / "phenotype"

    @property
    def broad_dir(self) -> Path:
        return self.phenotype_dir / "broad"

    @property
    def qupath_class_dir(self) -> Path:
        return self.phenotype_dir / "qupath_class"

    # ---- derived integration directories (data/integration/*) --------------
    @property
    def integration_dir(self) -> Path:
        return self.data_dir / "integration"

    @property
    def manifest_csv(self) -> Path:
        return self.integration_dir / "manifest.csv"

    @property
    def vocab_crosswalk_csv(self) -> Path:
        return self.integration_dir / "vocab_crosswalk.csv"

    @property
    def cells_pheno_dir(self) -> Path:
        return self.integration_dir / "cells_pheno"

    @property
    def cells_xen_dir(self) -> Path:
        return self.integration_dir / "cells_xen"

    @property
    def structures_dir(self) -> Path:
        return self.integration_dir / "structures"

    @property
    def registration_dir(self) -> Path:
        return self.integration_dir / "registration"

    @property
    def paired_dir(self) -> Path:
        return self.integration_dir / "paired"

    @property
    def niches_dir(self) -> Path:
        return self.integration_dir / "niches"

    @property
    def crossmodal_dir(self) -> Path:
        return self.integration_dir / "crossmodal"

    @property
    def integration_qc_dir(self) -> Path:
        return self.integration_dir / "qc"

    @property
    def integration_figures_dir(self) -> Path:
        return self.integration_dir / "figures"

    @property
    def integration_export_dir(self) -> Path:
        return self.integration_dir / "export"

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

def _as_bool(s) -> bool:
    return str(s).lower() in ("1", "true", "yes", "on")


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
    "lineage": {
        "hormone_min_norm": ("hormone_min_norm", float),
        "cd99_bright": ("cd99_bright", float),
    },
    "compute": {
        "n_jobs": ("n_jobs", int),
        "duckdb_threads": ("duckdb_threads", int),
        "use_gpu": ("use_gpu", lambda s: str(s).lower() in ("1", "true", "yes", "on")),
        "gpu_device": ("gpu_device", int),
    },
    "integration": {
        "mode": ("integration_mode", str),
        "xenium_paths_csv": ("xenium_paths_csv", Path),
        "donor_overrides_csv": ("donor_overrides_csv", Path),
        "xenium_root": ("xenium_root", str),
        "panel_explorer": ("panel_explorer", Path),
        "tissue": ("tissue", str),
        "fixed_modality": ("fixed_modality", str),
        "reg_pixel_um": ("reg_pixel_um", float),
        "reg_nonrigid": ("reg_nonrigid", _as_bool),
        "reg_max_disp_um": ("reg_max_disp_um", float),
        "islet_eps_um": ("islet_eps_um", float),
        "islet_min_samples": ("islet_min_samples", int),
        "niche_k": ("niche_k", int),
        "niche_n": ("niche_n", int),
        "niche_smooth_k": ("niche_smooth_k", int),
        "grid_um": ("grid_um", float),
        "match_max_dist_um": ("match_max_dist_um", float),
        "match_area_ratio": ("match_area_ratio", float),
        "crossmodal_min_anchors": ("crossmodal_min_anchors", int),
        "qc_tissue_dice_min": ("qc_tissue_dice_min", float),
        "qc_islet_rmse_max_um": ("qc_islet_rmse_max_um", float),
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
    "use_gpu": ("PHENOCYCLER_USE_GPU", _as_bool),
    "hormone_min_norm": ("PHENOCYCLER_HORMONE_MIN_NORM", float),
    # integration
    "integration_mode": ("PHENOCYCLER_INTEGRATION_MODE", str),
    "xenium_paths_csv": ("PHENOCYCLER_XENIUM_PATHS_CSV", Path),
    "xenium_root": ("PHENOCYCLER_XENIUM_ROOT", str),
    "panel_explorer": ("PHENOCYCLER_PANEL_EXPLORER", Path),
    "tissue": ("PHENOCYCLER_TISSUE", str),
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
    for name in ("data_dir", "images_dir", "cells_csv", "donor_metadata", "restore_vendor",
                 "xenium_paths_csv", "donor_overrides_csv", "panel_explorer"):
        p = Path(getattr(cfg, name)).expanduser()
        if not p.is_absolute():
            p = (base / p).resolve()
        setattr(cfg, name, p)

    if cfg.integration_mode not in ("sequential", "same_slide"):
        raise ValueError(
            f"[integration] mode must be 'sequential' or 'same_slide', got {cfg.integration_mode!r}"
        )
    if cfg.fixed_modality not in ("phenocycler", "xenium"):
        raise ValueError(
            f"[integration] fixed_modality must be 'phenocycler' or 'xenium', "
            f"got {cfg.fixed_modality!r}"
        )

    return cfg
