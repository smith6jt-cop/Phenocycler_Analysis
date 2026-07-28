"""
Idempotent orchestrator for the integration stages.

Mirrors ``phenocycler/pipeline.py`` deliberately: the same ``STAGES`` / ``ORDER`` shape, the
same "does the output already exist" skip logic, the same ``--status`` table. Someone who
knows how to drive the core pipeline already knows how to drive this one.

    python -m phenocycler.integration.pipeline --status
    python -m phenocycler.integration.pipeline                    # everything for the mode
    python -m phenocycler.integration.pipeline --only register --force
    python -m phenocycler.integration.pipeline --donor 6539

Stage order depends on the mode, and the difference is not cosmetic:

    sequential   manifest -> export_pheno -> import_xenium -> vocab -> structures
                 -> register -> match -> grid -> donor -> crossmodal -> qc -> figures
    same_slide   manifest -> export_pheno -> import_xenium -> vocab -> structures
                 -> register -> sameslide -> grid -> donor -> qc -> figures

``match`` and ``crossmodal`` are sequential-only because they approximate a pairing that
same-slide mode measures outright; ``sameslide`` is same-slide-only because that measurement
does not exist across serial sections. Each module enforces its own guard, so running a stage
directly cannot bypass the mode check either.

The registration-free stages — ``donor`` and the niche half of ``grid`` — run regardless of
whether registration succeeded. That is the point of ordering them after ``register`` but not
gating them on it: a cohort where half the donors fail to align still produces a complete
donor-level and niche-level result for all of them.

Tissue is a scope, not a mode. The cohort is 20 donors x {PhenoCycler, Xenium} x {pancreas,
pancreatic lymph node}, so a bare invocation runs *both* tissues — every ROI of every tissue
in ``[integration] tissues`` — and ``--tissue pancreas`` narrows it to one. The two tissues
differ in exactly one stage: ``structures`` (and therefore ``match``) needs a structural unit,
and a lymph node has none, so those stages announce the skip and move on rather than failing.
Everything else — export, import, vocab, register, grid, donor, crossmodal, qc — applies to
both. See ``phenocycler.integration.tissues``.
"""

from __future__ import annotations

import argparse
import glob
from typing import Optional

from ..config import PipelineConfig, load_config
from .tissues import DISPLAY_NAME, has_structures

#: (stage, output glob under data_dir, human label)
STAGES: list[tuple[str, str, str]] = [
    ("manifest", "integration/manifest.csv", "Donor <-> Xenium pairing manifest"),
    ("export_pheno", "integration/cells_pheno/donor_id=*", "PhenoCycler -> contract"),
    ("import_xenium", "integration/cells_xen/donor_id=*", "Xenium -> contract"),
    ("vocab", "integration/vocab_crosswalk.csv", "Lineage + marker crosswalk"),
    ("structures", "integration/structures/*/donor_id=*", "Islets / ducts / vessels"),
    ("register", "integration/registration/donor_id=*", "Serial-section registration"),
    ("match", "integration/paired/islets.parquet", "Islet matching [sequential]"),
    ("sameslide", "integration/paired/cells.parquet", "Cell pairing [same_slide]"),
    ("grid", "integration/niches/niche_profiles.csv", "Niches + registered grid"),
    ("donor", "integration/qc/composition.csv", "Donor-level concordance"),
    ("crossmodal", "integration/crossmodal/pseudocell_links.parquet",
     "Pseudo-cell links [sequential, inference]"),
    ("qc", "integration/qc/integration_report.txt", "QC gates + report"),
    ("figures", "integration/figures/composition.png", "Figures"),
]

ORDER_SEQUENTIAL = ["manifest", "export_pheno", "import_xenium", "vocab", "structures",
                    "register", "match", "grid", "donor", "crossmodal", "qc", "figures"]
ORDER_SAME_SLIDE = ["manifest", "export_pheno", "import_xenium", "vocab", "structures",
                    "register", "sameslide", "grid", "donor", "qc", "figures"]

#: Stages that need no registration, so they are never skipped because registration failed.
REGISTRATION_FREE = {"manifest", "export_pheno", "import_xenium", "vocab", "structures",
                     "donor", "qc"}

#: Stages that need a structural unit (islets). Not applicable to a tissue that has none —
#: see ``tissues.STRUCTURE_KINDS``. They are skipped per-tissue, never per-run: a combined
#: pancreas + lymph-node run still matches islets in the pancreas half.
NEEDS_STRUCTURES = {"structures", "match"}

#: Stages that run once for the whole selection rather than once per ROI.
COHORT_WIDE = {"manifest", "vocab", "structures", "donor", "qc"}


def order_for(mode: str) -> list[str]:
    return ORDER_SAME_SLIDE if mode == "same_slide" else ORDER_SEQUENTIAL


def applicable_tissues(stage: str, tissues: list[str]) -> list[str]:
    """Which of the selected tissues this stage actually runs on."""
    if stage not in NEEDS_STRUCTURES:
        return list(tissues)
    return [t for t in tissues if has_structures(t)]


def _has_outputs(cfg: PipelineConfig, pattern: str) -> int:
    return len(glob.glob(str(cfg.data_dir / pattern)))


def status(cfg: PipelineConfig, tissues: Optional[list[str]] = None) -> None:
    """Print the stage table. Needs no data and never writes anything."""
    order = order_for(cfg.integration_mode)
    selected = list(tissues) if tissues else cfg.tissue_list
    print(f"integration mode : {cfg.integration_mode}")
    print(f"data dir         : {cfg.data_dir}")
    print(f"fixed frame      : {cfg.fixed_modality}")
    print(f"tissues          : {', '.join(DISPLAY_NAME.get(t, t) for t in selected)}")
    print(f"rois             : {', '.join(cfg.resolve_rois(None, selected))}")
    print()
    print(f"  {'stage':<15} {'done':>5}  {'':<3} description")
    print("  " + "-" * 74)
    for name, pattern, label in STAGES:
        n = _has_outputs(cfg, pattern)
        mark = "ok " if n else "   "
        scope = ""
        if name not in order:
            mark = "n/a"
        else:
            runs_on = applicable_tissues(name, selected)
            if not runs_on:
                mark = "n/a"
                scope = "  [no selected tissue has a structural unit]"
            elif len(runs_on) < len(selected):
                scope = f"  [{', '.join(DISPLAY_NAME.get(t, t) for t in runs_on)} only]"
        print(f"  {name:<15} {n:>5}  {mark} {label}{scope}")
    print()
    print(f"stage order for {cfg.integration_mode}: {' -> '.join(order)}")


def run_pipeline(
    cfg: PipelineConfig,
    *,
    only: Optional[list[str]] = None,
    force: bool = False,
    donors: Optional[list[str]] = None,
    roi: Optional[str] = None,
    tissues: Optional[list[str]] = None,
    use_images: bool = True,
) -> None:
    """Run the stages for the configured mode, skipping anything already complete.

    ``roi`` selects one section; ``tissues`` selects whole tissues; neither runs the whole
    cohort. Per-ROI stages iterate the resolved ROI list; cohort-wide stages see the
    selection and scope themselves.
    """
    selected = list(tissues) if tissues else cfg.tissue_list
    rois = cfg.resolve_rois(roi, tissues)
    if not rois:
        raise SystemExit(f"no ROIs for tissue selection {selected} — check [integration] "
                         f"tissues in config.ini")
    print(f"[pipeline] tissues: {', '.join(selected)} | rois: {', '.join(rois)}", flush=True)

    order = order_for(cfg.integration_mode)
    if only:
        unknown = [s for s in only if s not in {n for n, _, _ in STAGES}]
        if unknown:
            raise SystemExit(f"unknown stage(s): {unknown}")
        wrong_mode = [s for s in only if s not in order]
        if wrong_mode:
            raise SystemExit(
                f"stage(s) {wrong_mode} do not apply in {cfg.integration_mode!r} mode "
                f"(applicable: {order})")
        order = [s for s in order if s in only]

    patterns = {name: pattern for name, pattern, _ in STAGES}
    # One tissue selected -> pass it down so each stage scopes itself; several -> pass None,
    # which means "the whole selection" and lets a stage span both tissues in one pass.
    one_tissue = selected[0] if len(selected) == 1 else None

    def _run(name: str, fn) -> None:
        """Run one stage over the whole selection.

        Stages take the selection rather than a single ROI and write once at the end. That is
        not a stylistic choice: the output paths are stage-level (``niche_profiles.csv``,
        ``islets.parquet``, ``crossmodal_qc.csv``), so calling a stage once per ROI would have
        each ROI's write silently clobber the previous one's.
        """
        if not force and _has_outputs(cfg, patterns[name]):
            print(f"[SKIP] {name}: outputs already present at "
                  f"{cfg.data_dir / patterns[name]}", flush=True)
            return
        print(f"\n=== {name} ===", flush=True)
        fn()

    for name in order:
        runs_on = applicable_tissues(name, selected)
        if not runs_on:
            print(f"[n/a] {name}: no selected tissue has a structural unit "
                  f"({', '.join(DISPLAY_NAME.get(t, t) for t in selected)})", flush=True)
            continue
        # A stage that only some selected tissues support narrows to those tissues; with one
        # supported tissue that is a plain scope, with several the stage skips per-row itself.
        stage_tissue = runs_on[0] if len(runs_on) == 1 else one_tissue

        if name == "manifest":
            from .manifest import run_manifest
            _run(name, lambda: run_manifest(cfg))
        elif name == "export_pheno":
            from .export_pheno import run_export_pheno
            _run(name, lambda: run_export_pheno(cfg, donors, roi=roi, tissues=runs_on))
        elif name == "import_xenium":
            from .import_xenium import run_import_xenium
            _run(name, lambda: run_import_xenium(cfg, donors=donors, roi=roi,
                                                 tissue=stage_tissue))
        elif name == "vocab":
            from .vocab import run_vocab
            _run(name, lambda: [run_vocab(cfg, _panel_genes(cfg, t), tissue=t)
                                for t in runs_on])
        elif name == "structures":
            from .structures import run_structures
            _run(name, lambda: run_structures(cfg, donors=donors, tissues=runs_on))
        elif name == "register":
            from .register import run_register
            _run(name, lambda: run_register(cfg, donors, roi=roi, tissue=stage_tissue,
                                            use_images=use_images))
        elif name == "match":
            from .match import run_match
            _run(name, lambda: run_match(cfg, donors, roi=roi, tissue=stage_tissue))
        elif name == "sameslide":
            from .sameslide import run_sameslide
            _run(name, lambda: run_sameslide(cfg, donors, roi=roi, tissue=stage_tissue))
        elif name == "grid":
            from .grid import run_grid
            _run(name, lambda: run_grid(cfg, donors, roi=roi, tissue=stage_tissue))
        elif name == "donor":
            from .donor import run_donor
            _run(name, lambda: run_donor(cfg, donors, roi=roi, tissue=stage_tissue))
        elif name == "crossmodal":
            from .crossmodal import run_crossmodal
            _run(name, lambda: run_crossmodal(cfg, donors, roi=roi, tissue=stage_tissue))
        elif name == "qc":
            from .qc import run_qc
            # Always re-run: the report summarises whatever the other stages just produced,
            # so a cached one is stale by construction.
            print(f"\n=== {name} ===", flush=True)
            run_qc(cfg, roi=roi, tissue=stage_tissue)
        elif name == "figures":
            from .figures import run_figures
            _run(name, lambda: run_figures(cfg, donors, roi=roi, tissue=stage_tissue))


def _panel_genes(cfg: PipelineConfig, tissue: str) -> Optional[list[str]]:
    """The run's real gene list for a tissue, when a Xenium table has been imported.

    Matters because the marker crosswalk must be validated against the genes actually
    present, not the nominal panel — several textbook markers live only in the custom add-on,
    and INS/GCG/SST are absent from both. Scoped to the tissue's own ROIs so a pancreas
    partition cannot supply the panel a lymph-node crosswalk is checked against; the panel
    is normally identical across ROIs, but reading it from the wrong tissue would hide a
    genuine per-tissue panel difference rather than surface it.
    """
    from .contract import discover_partitions, feature_name, read_cell_table

    wanted = set(cfg.rois_for_tissue(tissue))
    parts = [p for p in discover_partitions(cfg.cells_xen_dir) if not wanted or p[1] in wanted]
    if not parts:
        return None
    donor, part_roi = parts[0]
    try:
        t = read_cell_table(cfg.cells_xen_dir, donor, part_roi, modality="xenium")
    except FileNotFoundError:
        return None
    genes = [feature_name(c) for c in t.features]
    return [g for g in genes if not g.startswith("score_")] or None


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="PhenoCycler <-> Xenium integration pipeline",
        epilog="Stages are idempotent: re-running skips anything already complete.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--status", action="store_true", help="show the stage table and exit")
    ap.add_argument("--only", action="append", dest="only", help="run only this stage")
    ap.add_argument("--force", action="store_true", help="re-run even if outputs exist")
    ap.add_argument("--donor", action="append", dest="donors")
    ap.add_argument("--roi", default=None,
                    help="run one ROI only (e.g. pln_2); overrides --tissue")
    ap.add_argument("--tissue", action="append", dest="tissues",
                    help="run one tissue's ROIs (repeatable); default is every configured "
                         "tissue, i.e. pancreas and pancreatic lymph node together")
    ap.add_argument("--mode", choices=["sequential", "same_slide"], default=None,
                    help="override [integration] mode")
    ap.add_argument("--no-images", action="store_true",
                    help="registration point-set track only (skip morphology reads)")
    args = ap.parse_args(argv)

    overrides = {}
    if args.mode:
        overrides["integration_mode"] = args.mode
    cfg = load_config(args.config, **overrides)

    if args.tissues:
        unknown = [t for t in args.tissues if t not in cfg.tissue_list]
        if unknown:
            raise SystemExit(f"unknown tissue(s) {unknown}; configured: {cfg.tissue_list}")

    if args.status:
        status(cfg, args.tissues)
        return 0

    run_pipeline(cfg, only=args.only, force=args.force, donors=args.donors,
                 roi=args.roi, tissues=args.tissues, use_images=not args.no_images)
    print("\nintegration pipeline complete.")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
