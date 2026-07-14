"""
Central configuration for the phenocycler pipeline.

Everything that was a hardcoded absolute path or a magic constant in the
upstream ``scripts/senior/*.py`` scripts lives here instead, resolved from
``config.ini`` (see the repo root) with environment-variable and keyword
overrides.  The scientific defaults (REDSEA subtract-only / α=1 / 1-px band,
RESTORE SSC model + robust guard, the 7 broad lineages) are preserved exactly;
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

# Broad phenotyping — 5-compartment ORDERED RESIDUAL partition (+ explicit "Other").
# Each cell is typed by the first matching gate in COMPARTMENT_ORDER (each gate runs on the
# residual of the prior); cells failing every gate -> "Other" (real panel gaps: mast/Schwann/
# adipocyte/quiescent-stellate). Endocrine/exocrine are a SUB-branch of Epithelial (hormone+ vs
# hormone-). Priority: Epithelial > Endothelial > Neural > Immune > Mesenchymal > Other.
# See phenocycler/lineage.py for the tree + sub-splits.
COMPARTMENT_ORDER: list[str] = ["Epithelial", "Endothelial", "Neural", "Immune", "Mesenchymal"]
OTHER_LABEL: str = "Other"

# Anchor gate marker(s) per compartment — a cell enters the compartment if ANY is `_pos`.
# Mesenchymal is ORDERED: SMA (muscle/pericyte) claimed first, then Vimentin (fibroblast/stellate)
# in the SMA-negative residual (Vimentin is the most promiscuous marker -> gated last).
COMPARTMENT_GATES: dict[str, list[str]] = {
    "Epithelial":  ["E_cadherin"],          # only epithelial marker that stays + on endocrine
    "Endothelial": ["CD31"],
    "Neural":      ["B3TUBB"],               # residual (after epithelial) so islet TUBB3 is gone
    "Immune":      ["CD3e", "CD20", "CD79a", "CD68", "CD163", "CD206", "Iba1", "CD11b", "CD11c", "MPO"],
    "Mesenchymal": ["SMA", "Vimentin"],
}

# Endocrine sub-branch of Epithelial: hormone+ = endocrine (beta INS, alpha GCG, delta SST);
# hormone- epithelial = exocrine. (IAPP removed 2026-07-10 — failed marker.)
ENDOCRINE_MARKERS: list[str] = ["INS", "GCG", "SST"]

COMPARTMENT_COLORS: dict[str, str] = {
    "Epithelial": "#4477AA", "Endothelial": "#66CCEE", "Neural": "#B5838D",
    "Immune": "#228833", "Mesenchymal": "#AA3377", "Other": "#BBBBBB",
}
# unique 3-char codes for the terse per-donor progress line.
COMPARTMENT_ABBR: dict[str, str] = {
    "Epithelial": "Epi", "Endothelial": "Eth", "Neural": "Nrl",
    "Immune": "Imm", "Mesenchymal": "Mes", "Other": "Oth",
}
STATUS_ORDER: list[str] = ["ND", "AAB", "T1D"]

# RESTORE mutually-exclusive pairs, DIRECTIONAL: [target, counterpart] — [0] is thresholded, [1] is the
# mutually-exclusive counterpart that defines the negative population (RESTORE keys threshs on the TARGET
# only, so reciprocals like INS<->GCG / CD3e<->CD20 are separate entries). ONE RESTORE pass over all
# markers (the old separate "extra" pass is folded in). Counterparts are biologically curated (user).
MARKER_PAIRS: list[list[str]] = [
    # Epithelial (<- CD31: spatially-separated exclusive negative, avoids intraepithelial-lymphocyte leak)
    ["E_cadherin", "CD31"], ["Pan_Cytokeratin", "CD31"], ["Ker8_18", "CD31"], ["EpCAM", "CD31"],
    ["Keratin_5", "CD31"], ["TP63", "CD31"],
    # Endothelial / stroma (<- E_cadherin: epithelium = dominant negative mass)
    ["CD31", "E_cadherin"], ["CD34", "E_cadherin"], ["Podoplanin", "E_cadherin"], ["SMA", "E_cadherin"],
    ["Vimentin", "E_cadherin"],                 # Vimentin only ever a TARGET, never a counterpart
    # Neural (<- CD3e: TUBB3-negative & non-epithelial; endocrine co-expresses TUBB3)
    ["B3TUBB", "CD3e"],
    # Endocrine, intra-islet reciprocal (negative sits in the target's microenvironment)
    ["INS", "GCG"], ["GCG", "INS"], ["SST", "INS"],   # IAPP removed 2026-07-10 — failed marker
    # Immune (co-localized immune counterparts = shared autofluorescence context)
    ["CD3e", "CD20"], ["CD20", "CD3e"], ["CD79a", "CD3e"],
    ["CD4", "CD20"], ["CD8", "CD20"], ["FOXP3", "CD20"],
    ["CD68", "CD3e"], ["CD163", "CD3e"], ["CD206", "CD3e"], ["Iba1", "CD3e"], ["CD11b", "CD3e"],
    ["CD11c", "CD3e"], ["CD209", "CD3e"], ["MPO", "CD3e"],
    ["Granzyme_B", "CD20"], ["CD57", "CD20"], ["CD56", "CD20"],   # CD56/CD57 for NK; never a counterpart
]

# Optional functional/state markers (activation/exhaustion/IFN-driven) — RESTORE's reliably-negative-
# lineage assumption is weaker here; run separately and validate the cutoffs. These are STATE
# annotations, not compartments. NB: CD38 is batch-1 only (b_Catenin1 replaces it in batch-2).
FUNCTIONAL_PAIRS: list[list[str]] = [
    ["PD_1", "E_cadherin"], ["LAG3", "E_cadherin"], ["TOX", "E_cadherin"], ["TCF_1", "CD68"],
    ["ICOS", "E_cadherin"], ["VISTA", "E_cadherin"], ["CD39", "E_cadherin"], ["CD38", "E_cadherin"],
    ["PD_L1", "CD20"], ["IDO1", "CD20"], ["iNOS", "CD20"], ["HLA_DR", "CD3e"],
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
    # idx_select OFF-marker near-zero ceiling: per-(image,marker) quantile for the corrected
    # MUTUALLY-EXCLUSIVE single-positive selection (REF+/TGT≈0 negative pop | TGT+/REF≈0 signal),
    # replacing RESTORE's co-positive corner. DYNAMIC knob (default median q0.5). 0/negative ->
    # revert to the vendored co-positive >50 prefilter (comparison escape hatch only).
    restore_idx_floor_q: float = 0.5
    # -- 3 opt-in evaluation knobs (default == current behavior; activate per-donor from figures) --
    # Item 1: negative-cluster threshold statistic. "mean3sd" == vendored mu+3sigma (production);
    #   "mad"/"logmad"/"pctile" are robust variants (median+k*1.4826*MAD / coverage-matched pctile).
    restore_neg_stat: str = "mean3sd"
    # Item 2: per-image signal-presence guard. If >0, an (image,marker) whose idx_select cloud shows
    #   cluster-mean separation < this many negative-sigma is called ABSENT (thresh=+inf -> 0 positives).
    restore_presence_min_sep: float = 0.0
    # Item 3: REDSEA zero-clip guard for the idx_select OFF-marker ceilings. If >0, floor c0/c1 at the
    #   nonzero_q quantile of the marker's NONZERO cells (prevents ceiling collapse to 0 when >=50%
    #   of cells are REDSEA-clamped to zero). 0 == disabled (exact current math).
    restore_idx_nonzero_q: float = 0.0

    # -- RESTORE efficacy diagnostic (scripts/R/restore_mxnorm_diagnostics.R via mxnorm) ---
    restore_diag_subsample: int = 20000   # cells/donor for the mxnorm diagnostic (stratified by cell_region)
    restore_diag_seed: int = 0

    # (former [lineage] _norm-floor + cd99_bright knobs removed — floors dropped in the curated redesign)

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
    def restore_mxnorm_dir(self) -> Path:
        """mxnorm before/after RESTORE-efficacy diagnostic outputs (efficacy CSV + figures)."""
        return self.restore_dir / "mxnorm"

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
        "idx_floor_q": ("restore_idx_floor_q", float),
        "neg_stat": ("restore_neg_stat", str),
        "presence_min_sep": ("restore_presence_min_sep", float),
        "idx_nonzero_q": ("restore_idx_nonzero_q", float),
    },
    "restore_diag": {
        "subsample": ("restore_diag_subsample", int),
        "seed": ("restore_diag_seed", int),
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
