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

from .cohort import filter_eligible_donors

_REPO_ROOT = Path(__file__).resolve().parents[1]
_DEFAULT_INI = _REPO_ROOT / "config.ini"
_PACKAGE_DIR = Path(__file__).resolve().parent


def _boolean(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


# --------------------------------------------------------------------------- #
# Legacy broad-lineage vocabulary — integration layer only.
#
# The core typing pipeline on this branch does NOT use these: it types from
# marker_registry.json + typing_rules.json into its own compartments. They
# remain because phenocycler/integration/ imports them at module scope
# (vocab.py derives the PhenoCycler side of the PhenoCycler<->Xenium crosswalk
# from LINEAGES) and would not import at all without them.
#
# They therefore describe the vocabulary the integration crosswalk was BUILT
# against, not the vocabulary this branch's typing emits. Reconciling the two
# is the open question tracked in docs/INTEGRATION.md; until that lands, do not
# read these as a statement of the current taxonomy.
# --------------------------------------------------------------------------- #
LINEAGES: dict[str, list[str]] = {
    "Epithelial": ["Pan_Cytokeratin"],
    "Fibroblast": ["Vimentin"],
    "Immune": ["CD3e", "CD20", "CD163", "MPO"],
    "Endocrine": ["INS", "GCG", "SST", "CD99"],
    "Endothelial": ["CD31"],
    "Muscle": ["SMA"],
    "Neural": ["B3TUBB"],
}
LINEAGE_COLORS: dict[str, str] = {
    "Epithelial": "#4477AA", "Fibroblast": "#EE6677", "Immune": "#228833",
    "Endocrine": "#CCBB44", "Endothelial": "#66CCEE", "Muscle": "#AA3377",
    "Neural": "#B5838D",
}
STATUS_ORDER: list[str] = ["ND", "AAB", "T1D"]
HORMONE_MARKERS: list[str] = ["INS", "GCG", "SST"]
EXTRA_MARKERS: list[str] = ["CD99", "B3TUBB", "MPO"]


@dataclass
class PipelineConfig:
    """Resolved workflow paths and scientific/compute parameters."""

    # Source and output contracts.
    data_dir: Path = _REPO_ROOT / "data"
    qupath_manifest: Path = _REPO_ROOT / "data" / "qupath_manifest.json"
    marker_registry: Path = _PACKAGE_DIR / "marker_registry.json"
    typing_rules: Path = _PACKAGE_DIR / "typing_rules.json"
    # Optional, immutable probabilistic rescue model.  Keeping this unset is
    # the production rules-only default. Production requires its separately
    # configured immutable promotion manifest.
    typing_model_bundle: Path | None = None
    typing_model_promotion_manifest: Path | None = None

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
    cell_qc_duplicate_policy: str = "fatal"
    # InstanSeg fragments without a nucleus remain auditable source objects but
    # cannot enter analysis, model estimation, or physical spillover context.
    cell_qc_missing_geometry_policy: str = "exclude"

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
    pipeline_geometry_workers: int = 1
    pipeline_redsea_workers: int = 1
    pipeline_downstream_workers: int = 2

    # ---------------------------------------------------------------------- #
    # Integration-layer surface (phenocycler/integration/).
    #
    # Consumed ONLY by the PhenoCycler <-> Xenium layer. The eight-stage core
    # pipeline never reads any of it, and none of it enters _identity_config(),
    # so it cannot change what a run records as its scientific configuration.
    #
    # It lives here because the integration layer arrived on main after this
    # branch diverged: it is 22 modules of shipped, tested code whose only
    # coupling to the core is this config surface. Dropping the fields would
    # delete that work by omission rather than by decision.
    # ---------------------------------------------------------------------- #
    images_dir: Path = _REPO_ROOT / "data" / "raw" / "images"
    donor_metadata: Path = _REPO_ROOT / "data" / "raw" / "donor_metadata_panc.xlsx"
    metadata_donor_col: str = "donor_id"
    metadata_status_col: str = "disease.status"
    hormone_min_norm: float = 5.0    # K: {INS,GCG,SST}_pos := _norm >= K (false-endocrine floor)
    cd99_bright: float = 3.0         # CD99_pos := CD99_norm >= this (bright-only Endocrine gate)
    integration_mode: str = "sequential"          # sequential | same_slide
    xenium_paths_csv: Path = _REPO_ROOT / "data" / "integration" / "xenium_paths.csv"
    donor_overrides_csv: Path = _REPO_ROOT / "data" / "integration" / "donor_overrides.csv"
    xenium_root: str = ""                         # optional prefix rewrite for bundle paths
    panel_explorer: Path = _REPO_ROOT / "external" / "XeniumPanelExplorer"
    tissue: str = "pancreas"                      # default for datasets that do not declare one
    tissues: str = "pancreas,pancreatic_lymph_node"   # comma-separated; the cohort's tissues
    pheno_tissue_dirs: str = ""
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
    xenium_hormone_if: bool = False
    cellmatch_max_dist_um: float = 3.0            # centroid cap for a candidate pair
    cellmatch_max_residual_um: float = 2.0        # skip an islet above this local residual
    cellmatch_min_signal_over_null: float = 3.0   # skip if not clearly above the permuted null
    qc_tissue_dice_min: float = 0.80
    qc_islet_rmse_max_um: float = 150.0

    _config_path: Path | None = field(default=None, repr=False, compare=False)

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

    def discover_donors(self, from_dir: Path | None = None) -> list[str]:
        """Unique, eligible **donor ids** behind the discovered sections.

        Distinct from :meth:`discover_sections` on purpose. Stages that process images want
        sections; the pairing manifest and the donor workbook join want donors. Conflating
        them is how ``6539pln`` becomes a phantom donor with no metadata and no Xenium
        counterpart, which then shows up in the manifest as a real pairing gap.

        Composes the two definitions this branch and ``main`` arrived at independently: the
        section-to-donor parse, and the fail-closed cohort exclusion. Both are load-bearing
        and neither subsumes the other — on a one-image-per-donor cohort the parse is the
        identity, so keeping it costs nothing and stops a multi-section donor from being
        counted twice the moment the lymph-node sections land.
        """
        from .sections import parse

        out = set()
        for key in self.discover_sections(from_dir):
            try:
                out.add(parse(key).donor_id)
            except Exception:  # noqa: BLE001 - an unparseable key is surfaced by the stage
                out.add(key)
        return filter_eligible_donors(sorted(out))

    # ---------------------------------------------------------------------- #
    # Integration-layer derived paths and tissue/ROI helpers. Same rationale as
    # the field block above. `broad_dir` / `restore_gated_dir` still name the
    # legacy layout because integration/export_pheno.py has not yet been
    # repointed at this branch's assignments/ + marker_evidence/ outputs; see
    # docs/INTEGRATION.md for the status of that follow-up.
    # ---------------------------------------------------------------------- #
    @property
    def phenotype_dir(self) -> Path:
        return self.data_dir / "phenotype"

    @property
    def broad_dir(self) -> Path:
        return self.phenotype_dir / "broad"

    @property
    def restore_gated_dir(self) -> Path:
        return self.data_dir / "restore_gated_redsea"

    @property
    def restore_gated_extra_dir(self) -> Path:
        return self.data_dir / "restore_gated_redsea_extra"

    @property
    def geojson_dir(self) -> Path:
        return self.redsea_scratch / "geojson"

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
    def cells_pheno_xif_dir(self) -> Path:
        """Contract tables for PhenoCycler images OF THE XENIUM SECTION (post-Xenium re-stain).

        Kept apart from `cells_pheno_dir` because they share a (donor, roi) but are a
        different piece of tissue: writing them into the same partition would silently mix
        two sections. Their counterpart is `cells_xen_dir` — same cells, same slide.
        """
        return self.integration_dir / "cells_pheno_xif"

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
    def tissue_list(self) -> list[str]:
        """The cohort's tissues, in declaration order."""
        return [t.strip() for t in str(self.tissues).split(",") if t.strip()]

    def tissue_for_roi(self, roi: str) -> str:
        """Which tissue an ROI belongs to.

        ROIs are the Xenium naming (``panc``, ``pln_1`` .. ``pln_3``); tissue is the
        biological class they map onto. Kept as a function rather than a literal dict so a new
        ROI (``pln_4``) needs no config change, and unknown ROIs fall back to the declared
        default rather than silently becoming their own tissue.
        """
        r = str(roi).strip().lower()
        if r.startswith("pln") or r.startswith("ln"):
            return "pancreatic_lymph_node"
        if r.startswith("pan"):
            return "pancreas"
        return self.tissue

    def rois_for_tissue(self, tissue: str) -> list[str]:
        """Known ROI names for a tissue. ``pln_1..3`` because one donor can have several."""
        if tissue == "pancreatic_lymph_node":
            return ["pln_1", "pln_2", "pln_3"]
        if tissue == "pancreas":
            return ["panc"]
        return []

    def pheno_dir_for_tissue(self, tissue: str) -> Path:
        """The PhenoCycler core-pipeline data dir holding this tissue's donors.

        Falls back to ``data_dir`` when no per-tissue mapping is configured. Callers must
        treat two tissues resolving to the same directory as a configuration gap, not as
        permission to export the same cells twice — see ``pheno_tissue_dirs``.
        """
        for entry in str(self.pheno_tissue_dirs).split(","):
            name, _, path = entry.partition("=")
            if name.strip() and name.strip() == tissue and path.strip():
                return Path(path.strip()).expanduser()
        return self.data_dir

    def resolve_rois(self, roi=None, tissues=None) -> list[str]:
        """Turn a ``--roi`` / ``--tissue`` selection into the ROI list to iterate.

        The three ways a run is scoped, in precedence order:

        * ``--roi pln_2``  -> exactly that ROI (the escape hatch for one section);
        * ``--tissue pancreatic_lymph_node`` -> every ROI of that tissue (a *separate* run);
        * neither -> every ROI of every configured tissue (a *combined* run).

        Defaulting to combined is deliberate. The cohort is 20 donors x 2 platforms x 2
        tissues, so a bare ``python -m phenocycler.integration.pipeline`` should process the
        whole dataset rather than silently doing pancreas and leaving the lymph nodes out.
        """
        if roi:
            return [str(roi)]
        names = list(tissues) if tissues else self.tissue_list
        out: list[str] = []
        for t in names:
            for r in self.rois_for_tissue(str(t).strip()):
                if r not in out:
                    out.append(r)
        return out

    def discover_sections(self, from_dir: Path | None = None) -> list[str]:
        """Section keys by globbing ``<dir>/donor_id=*`` (defaults to cells_dir).

        This is the pipeline's **work unit**: one key is one image, one GeoJSON, one
        RESTORE threshold set, one coordinate frame. A donor with a pancreas and a lymph
        node contributes two — ``6539`` and ``6539pln``.
        """
        base = Path(from_dir) if from_dir is not None else self.cells_dir
        return sorted(p.name.split("=", 1)[1] for p in base.glob("donor_id=*"))

    def discover_section_objects(self, from_dir: Path | None = None) -> list:
        """:class:`phenocycler.sections.Section` for every discovered partition."""
        from .sections import parse

        return [parse(k) for k in self.discover_sections(from_dir)]

_INI_SCHEMA = {
    "paths": {
        "data_dir": ("data_dir", Path),
        "qupath_manifest": ("qupath_manifest", Path),
        "marker_registry": ("marker_registry", Path),
        "typing_rules": ("typing_rules", Path),
    },
    "typing": {
        "model_bundle": ("typing_model_bundle", Path),
        "promotion_manifest": ("typing_model_promotion_manifest", Path),
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
        "duplicate_policy": ("cell_qc_duplicate_policy", str),
        "missing_geometry_policy": ("cell_qc_missing_geometry_policy", str),
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
        "geometry_workers": ("pipeline_geometry_workers", int),
        "redsea_workers": ("pipeline_redsea_workers", int),
        "downstream_workers": ("pipeline_downstream_workers", int),
    },
    "metadata": {
        "donor_col": ("metadata_donor_col", str),
        "status_col": ("metadata_status_col", str),
    },
    "lineage": {
        "hormone_min_norm": ("hormone_min_norm", float),
        "immune_min_norm": ("immune_min_norm", float),
        "cd99_bright": ("cd99_bright", float),
    },
    "integration": {
        "mode": ("integration_mode", str),
        "xenium_paths_csv": ("xenium_paths_csv", Path),
        "donor_overrides_csv": ("donor_overrides_csv", Path),
        "xenium_root": ("xenium_root", str),
        "panel_explorer": ("panel_explorer", Path),
        "tissue": ("tissue", str),
        "tissues": ("tissues", str),
        "pheno_tissue_dirs": ("pheno_tissue_dirs", str),
        "fixed_modality": ("fixed_modality", str),
        "reg_pixel_um": ("reg_pixel_um", float),
        "reg_nonrigid": ("reg_nonrigid", _boolean),
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
        "xenium_hormone_if": ("xenium_hormone_if", _boolean),
        "cellmatch_max_dist_um": ("cellmatch_max_dist_um", float),
        "cellmatch_max_residual_um": ("cellmatch_max_residual_um", float),
        "cellmatch_min_signal_over_null": ("cellmatch_min_signal_over_null", float),
        "qc_tissue_dice_min": ("qc_tissue_dice_min", float),
        "qc_islet_rmse_max_um": ("qc_islet_rmse_max_um", float),
    },
}

_ENV_OVERRIDES = {
    "data_dir": ("PHENOCYCLER_DATA_DIR", Path),
    "qupath_manifest": ("PHENOCYCLER_QUPATH_MANIFEST", Path),
    "n_jobs": ("PHENOCYCLER_JOBS", int),
    "use_gpu": ("PHENOCYCLER_USE_GPU", _boolean),
    # Integration layer (see the quarantined block on PipelineConfig).
    "images_dir": ("PHENOCYCLER_IMAGES_DIR", Path),
    "donor_metadata": ("PHENOCYCLER_DONOR_METADATA", Path),
    "hormone_min_norm": ("PHENOCYCLER_HORMONE_MIN_NORM", float),
    "integration_mode": ("PHENOCYCLER_INTEGRATION_MODE", str),
    "xenium_paths_csv": ("PHENOCYCLER_XENIUM_PATHS_CSV", Path),
    "xenium_root": ("PHENOCYCLER_XENIUM_ROOT", str),
    "panel_explorer": ("PHENOCYCLER_PANEL_EXPLORER", Path),
    "tissue": ("PHENOCYCLER_TISSUE", str),
}


def load_config(
    config_path: os.PathLike | str | None = None,
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

    for name in (
        "typing_model_bundle",
        "typing_model_promotion_manifest",
    ):
        value = getattr(cfg, name)
        if value is None:
            continue
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
        setattr(cfg, name, path)
    if (cfg.typing_model_bundle is None) != (
        cfg.typing_model_promotion_manifest is None
    ):
        raise ValueError(
            "[typing] model_bundle and promotion_manifest must be configured together"
        )

    for name in ("images_dir", "donor_metadata", "xenium_paths_csv",
                 "donor_overrides_csv", "panel_explorer"):
        path = Path(getattr(cfg, name)).expanduser()
        if not path.is_absolute():
            path = (base / path).resolve()
        setattr(cfg, name, path)

    if cfg.integration_mode not in ("sequential", "same_slide"):
        raise ValueError(
            f"[integration] mode must be 'sequential' or 'same_slide', "
            f"got {cfg.integration_mode!r}")

    if cfg.redsea_norm_form != "donor":
        raise ValueError("production REDSEA requires mass-conserving norm_form='donor'")
    if cfg.cell_qc_missing_geometry_policy not in {"fatal", "exclude"}:
        raise ValueError(
            "geometry_qc missing_geometry_policy must be 'fatal' or 'exclude'"
        )
    if cfg.cell_qc_duplicate_policy not in {"fatal", "deterministic_keep"}:
        raise ValueError(
            "geometry_qc duplicate_policy must be 'fatal' or "
            "'deterministic_keep'"
        )
    worker_limits = {
        "geometry_workers": cfg.pipeline_geometry_workers,
        "redsea_workers": cfg.pipeline_redsea_workers,
        "downstream_workers": cfg.pipeline_downstream_workers,
    }
    invalid_workers = [
        name for name, value in worker_limits.items() if value < 1
    ]
    if invalid_workers:
        raise ValueError(
            "pipeline worker limits must be >= 1: "
            f"{sorted(invalid_workers)}"
        )
    return cfg


__all__ = [
    "PipelineConfig", "load_config",
    # integration-layer vocabulary (see the block above)
    "LINEAGES", "LINEAGE_COLORS", "STATUS_ORDER", "HORMONE_MARKERS", "EXTRA_MARKERS",
]
