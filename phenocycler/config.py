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

from .cohort import filter_eligible_donors

# --------------------------------------------------------------------------- #
# Scientific constants (faithful to Islet-Explorer-Senior; not usually tuned)
# --------------------------------------------------------------------------- #

# Broad phenotyping — 5-compartment ORDERED RESIDUAL partition over the RESTORE tri-state calls
# (+ explicit "Other" and "Unresolved"). Each cell is typed by the first POSITIVE gate in
# COMPARTMENT_ORDER (each gate on the residual of the prior). The first gate that is UNAVAILABLE
# (indeterminate — a required marker has no valid call for that cell) BLOCKS assignment -> "Unresolved":
# we cannot rule the cell into or out of that (higher-priority) compartment, so no lower gate may claim
# it. A cell that is valid-and-NEGATIVE for every gate -> "Other" (a real panel gap, NOT force-assigned).
# Priority high->low: Immune > Epithelial > Endothelial > Neural > Mesenchymal > Other.
# Gate 4 (2026-07-24): Immune moved to FIRST, and the tree runs on tri-state positive|negative|unavailable
# calls emitted by restore_apply from the frozen ACCEPTED_PAIRS. Endocrine/exocrine + immune/mesenchymal
# SUBTYPES are DEFERRED to a validated level-2 pass (docs/level2_hierarchy.md); this first deliverable is
# the 5 broad compartments only. See phenocycler/lineage.py for the tree.
# Endothelial moved AHEAD of Epithelial, 2026-07-27, when the Epithelial gate gained the hormones.
# Measured cause: inside islets 67.8% of CD31-POSITIVE cells were being claimed by Epithelial, and 78%
# of islet-core CD31+ cells were also hormone-positive -- against 1.8% in tissue. Islet capillaries sit
# inside a dense mass of hormone-bright neighbours and pick up their signal; outside islets CD31 and
# hormone are essentially exclusive. That asymmetry is spatial contamination, not biology, and inside an
# islet CD31 (reference_fold ~30x, the panel's most specific structural marker) is the more trustworthy
# evidence. Under the 7-marker gates endothelium was correctly ENRICHED in islets (8-13% vs 2-6% in
# tissue); adding the hormones inverted that to 2.1% vs 10.2%.
#
# Both constraints from docs/BroadLineageRevised.md still hold: Immune stays FIRST (CD31 is expressed on
# monocytes, neutrophils, NK and naive T cells, so endothelium must not precede immune), and Neural
# still follows Epithelial (endocrine cells are TUBB3+).
COMPARTMENT_ORDER: list[str] = ["Immune", "Endothelial", "Epithelial", "Neural", "Mesenchymal"]
OTHER_LABEL: str = "Other"
UNRESOLVED_LABEL: str = "Unresolved"   # blocked by an unavailable higher-priority gate (NOT a negative Other)

# Anchor gate marker(s) per compartment — a cell enters the compartment if ANY gate marker is a POSITIVE
# RESTORE call, and the gate set is EXACTLY the divisor-backed ACCEPTED_PAIRS targets.
#
# Gate markers, expanded 2026-07-27 (maintainer sign-off) from 7 to 24. Every entry is a divisor-backed
# ACCEPTED_PAIRS target, so a tri-state "unavailable" always means a real per-donor QC gap and never
# "never screened" (enforced by tests/test_config.py::test_gate_markers_are_the_accepted_restore_targets).
#
# Why the expansion: `Other` is the ordered residual's terminal bucket, so a cell whose only lineage
# marker was ungated could not be typed. With 7 markers, B cells, plasma cells, DC, neutrophils,
# macrophage subsets, ENDOCRINE cells, lymphatics and smooth muscle all landed in `Other` by
# construction. Measured consequence of the old Epithelial gate: only 1-6% of islet-core cells --
# endocrine by QuPath annotation, hormone-validated -- were called Epithelial, and 82-94% `Other`
# (docs/ANATOMY_CHECK_2026-07-27.md).
COMPARTMENT_GATES: dict[str, list[str]] = {
    # T / B / macrophage / DC / granulocyte. CD163 screened and rejected (arm separation below the 2x
    # floor on every disease group); CD68/CD11b/Iba1/CD206 already cover macrophages.
    "Immune":      ["CD3e", "CD68", "CD11b", "CD79a", "CD11c", "Iba1", "CD206", "MPO", "CD4", "CD8"],
    # Exocrine keratins + EpCAM + the hormones. No single marker spans this compartment: endocrine
    # sits BELOW acinar on E-cadherin and ABOVE it on EpCAM, and the keratins are near-zero in islets.
    "Epithelial":  ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM", "Keratin_5",
                    "INS", "GCG", "SST"],
    "Endothelial": ["CD31", "CD34", "Podoplanin"],   # vascular + lymphatic
    "Neural":      ["B3TUBB"],                       # residual (after epithelial) so islet TUBB3 is gone
    "Mesenchymal": ["Vimentin", "SMA"],              # fibroblast/stellate + mural
}

# Endocrine sub-branch of Epithelial: hormone+ = endocrine (beta INS, alpha GCG, delta SST). DEFERRED to
# the level-2 pass (Gate 4 delivers the 5 broad compartments only); retained here for the level-2 code and
# the config drift guard, NOT consumed by the Gate-4 lineage tree. (IAPP removed 2026-07-10 — failed marker.)
ENDOCRINE_MARKERS: list[str] = ["INS", "GCG", "SST"]

COMPARTMENT_COLORS: dict[str, str] = {
    "Immune": "#228833", "Epithelial": "#4477AA", "Endothelial": "#66CCEE", "Neural": "#B5838D",
    "Mesenchymal": "#AA3377", "Other": "#BBBBBB", "Unresolved": "#666666",
}
# unique 3-char codes for the terse per-donor progress line.
COMPARTMENT_ABBR: dict[str, str] = {
    "Immune": "Imm", "Epithelial": "Epi", "Endothelial": "Eth", "Neural": "Nrl",
    "Mesenchymal": "Mes", "Other": "Oth", "Unresolved": "Unr",
}
STATUS_ORDER: list[str] = ["ND", "AAB", "T1D"]

# RETIRED 2026-07-23 — MARKER_PAIRS and FUNCTIONAL_PAIRS (the curated web fed to the vendored-SSC
# `restore.run_restore` / `run_restore_functional`) are deleted along with those entry points, rather
# than repaired. Two reasons:
#   1. Seven MARKER_PAIRS and three FUNCTIONAL_PAIRS used CD20 as the reference (CD3e, CD4, CD8, FOXP3,
#      Granzyme_B, CD57, CD56, PD_L1, IDO1, iNOS). A reference's positive cells ARE the target's
#      negative control, and B cells are sparse in pancreas even under inflammation or disease, so CD20
#      cannot define one. There is no mechanical substitute either: CD3e is invalid for the T-subset
#      targets (CD4/CD8/FOXP3 cells are CD3e+), and epithelium expresses PD-L1 and IDO1.
#   2. The pair web is being rebuilt from the manuscript in `restore_validation.py` (NNMF separator,
#      maximum-negative-control divisor, per-pair acceptance state), and Gate 4 of
#      docs/restore_faithful_rebuild_plan.md replaces the production wiring outright. Maintaining a
#      second, known-wrong web until then is a hazard, not a fallback.
# The reusable RESTORE machinery in restore.py (idx_select, neg_stat, fit_clusters, threshold LUT and
# apply helpers, the proliferation pass) is untouched.

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
    redsea_comp_mode: int = 1        # 1 == subtract + reinforce (Bai et al. default); 0 == subtract-only
    redsea_alpha: float = 1.0
    redsea_gap_bridge: int = 1
    # "donor" == published mass-conserving operator F[i,j]=contact[i,j]/B[j]; "recipient" == pre-v10.
    redsea_norm_form: str = "donor"
    # Skip compensation for nuclear/ECM channels (marker_taxonomy.NO_SPILLOVER).
    redsea_exclude_no_spillover: bool = True

    # -- cell QC (phenocycler.cell_qc) — applied BEFORE any parameter estimation ---------------
    # Tier 1, hard exclusion (the object is not a cell). See cell_qc.py for the measured evidence
    # behind each value; all are donor-local or scale-free so none normalises across donors.
    cell_qc_min_cell_area: float = 20.0            # um^2; a 5 um disc, still below a lymphocyte
    cell_qc_require_nucleus: bool = True           # DNA gate: no nucleus -> DAPI 0.26-0.41x normal
    cell_qc_max_area_median_multiple: float = 4.0  # x the donor's own median retained cell area
    cell_qc_min_cell_solidity: float = 0.70        # area/convex-hull; concave/multi-lobed objects
    cell_qc_min_raster_area_ratio: float = 0.5     # rasterized vs reported area (geometry defect)
    # Tier 2, excluded from ESTIMATION only — these are real cells and still get a call, but their
    # whole-cell mean is a nuclear mean so they must not set a floor, arm, separator or divisor.
    cell_qc_no_cytoplasm_nc_ratio: float = 0.999
    # Duplicate-detection removal. Detection was run with selectAnnotations(), which selects the
    # Tissue annotation AND every nested Islet_N; InstanSeg segments each selected parent
    # independently, so islet nuclei were segmented twice. 8/20 donors, 100% of duplicates in islets.
    cell_qc_deduplicate: bool = True
    cell_qc_duplicate_radius_um: float = 1.0
    cell_qc_duplicate_area_tol: float = 0.15

    # -- RESTORE pair method policy (Gate 4 production apply) ------------------
    # These two are the METHOD_VERSION v10 knobs on restore_validation.PairValidationConfig. They were
    # previously reachable ONLY from the `restore_validation evaluate-locked` CLI, so a production
    # `restore_apply` run silently used the dataclass defaults regardless of which arm was chosen --
    # i.e. the policy under review and the policy in production could not be made to agree. Threading
    # them here is what makes an accepted sweep arm reproducible by the pipeline.
    #   control_definition: "reference_only" | "reference_and_target"
    #   divisor_statistic:  "max" | "p999" | "p99" | "p95"
    # Defaults deliberately mirror PairValidationConfig so this indirection cannot change behaviour on
    # its own; restore_validation remains the single source of truth for the allowed values.
    restore_control_definition: str = "reference_only"
    restore_divisor_statistic: str = "max"
    # QuPath compartments feeding the NNMF separator, comma-separated. The locked value is
    # "Membrane,Cell", but 14 of the 24 accepted targets post-date the compartment sample and have only
    # whole-cell values, so a production run over the full accepted set REQUIRES "Cell". Measured
    # 2026-07-27 on the 7 original pairs x 3 donors: states identical 21/21, divisor ratio 0.96-1.05,
    # call Jaccard median 0.997, call-rate |delta| median 0.04 pp.
    restore_feature_compartments: str = "Membrane,Cell"

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
    # Partner-low thresholding (ENDOCRINE targets only — INS/GCG/SST): if >0, set the cutoff from a 2-GMM
    #   crossover among ISLET cells that are LOW for the mutually-exclusive reference (reference <= this
    #   quantile), where the target is cleanly bimodal (removes the double-positive + non-islet-background
    #   confounds); +inf (absent) if not bimodal (e.g. beta-loss). Non-endocrine markers fall through to the
    #   negative-cluster statistic. 0 == disabled. (Does NOT generalize beyond islet-resident markers.)
    restore_partner_low_q: float = 0.0
    # Proliferation (Ki67/PCNA) binarization: these are GRADED (not bimodal), so a robust per-donor upper-
    #   tail cutoff = mean + k*std of log10(nonzero+1) (adaptive; no fixed-fraction assumption). k here:
    restore_proliferation_k: float = 3.0

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
    def restore_gated_proliferation_dir(self) -> Path:   # standalone-bimodal Ki67/PCNA positivity
        return self.data_dir / "restore_gated_proliferation"

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

    # ---- RESTORE production-apply inputs (Gate 4) --------------------------
    @property
    def restore_pair_validation_dir(self) -> Path:
        """Root of the gitignored RESTORE pair-validation artifacts (screens, pilots, reference dossier)."""
        return self.data_dir / "restore_pair_validation"

    @property
    def restore_allcells_sample_dir(self) -> Path:
        """All-cells (modulus-1, unsampled) QuPath compartment sample consumed by restore_apply."""
        return self.restore_pair_validation_dir / "qupath_compartment_full_20donor"

    @property
    def restore_pair_reviews_csv(self) -> Path:
        """Maintainer-signed ACCEPTED target/reference rules (Gate-2 sign-off); required by restore_apply."""
        return self.restore_pair_validation_dir / "pair_reviews_accepted.csv"

    @property
    def restore_expression_reviews_csv(self) -> Path:
        """Optional per-(donor,target) image reviews (e.g. confirmed absence); passed through to the screen."""
        return self.restore_pair_validation_dir / "expression_reviews.csv"

    @property
    def restore_apply_manifest_csv(self) -> Path:
        """Cohort donor x accepted-pair QC disposition + divisor manifest written by restore_apply."""
        return self.restore_gated_dir / "apply_manifest.csv"

    # ---- helpers -----------------------------------------------------------
    def discover_donors(self, from_dir: Optional[Path] = None) -> list[str]:
        """Discover donor partitions while enforcing central cohort exclusions."""
        base = Path(from_dir) if from_dir is not None else self.cells_dir
        discovered = sorted(
            p.name.split("=", 1)[1] for p in base.glob("donor_id=*")
        )
        return filter_eligible_donors(discovered)

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
        "norm_form": ("redsea_norm_form", str),
        "exclude_no_spillover": (
            "redsea_exclude_no_spillover",
            lambda s: str(s).lower() in ("1", "true", "yes", "on"),
        ),
    },
    "cell_qc": {
        "min_cell_area": ("cell_qc_min_cell_area", float),
        "require_nucleus": (
            "cell_qc_require_nucleus",
            lambda s: str(s).lower() in ("1", "true", "yes", "on"),
        ),
        "max_area_median_multiple": ("cell_qc_max_area_median_multiple", float),
        "min_cell_solidity": ("cell_qc_min_cell_solidity", float),
        "min_raster_area_ratio": ("cell_qc_min_raster_area_ratio", float),
        "no_cytoplasm_nc_ratio": ("cell_qc_no_cytoplasm_nc_ratio", float),
        "deduplicate": (
            "cell_qc_deduplicate",
            lambda s: str(s).lower() in ("1", "true", "yes", "on"),
        ),
        "duplicate_radius_um": ("cell_qc_duplicate_radius_um", float),
        "duplicate_area_tol": ("cell_qc_duplicate_area_tol", float),
    },
    "restore": {
        # METHOD_VERSION v10 pair-method policy (consumed by restore_apply -> evaluate_locked_pair).
        "control_definition": ("restore_control_definition", str),
        "divisor_statistic": ("restore_divisor_statistic", str),
        "feature_compartments": ("restore_feature_compartments", str),
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
        "partner_low_q": ("restore_partner_low_q", float),
        "proliferation_k": ("restore_proliferation_k", float),
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
