"""
S0 — the donor <-> Xenium-run pairing manifest.

This is the module that resolves the crosswalk problem, and it is worth being explicit
about what that problem is, because it is not obvious from either repo in isolation:

1. ``donor_id`` is a perfect join key *on paper*. Both repos already share
   ``donor_metadata_panc.xlsx`` verbatim, with the same ``donor_id`` and ``disease.status``
   columns and the same 26 nPOD donors. PhenoCycler derives ``donor_id`` by regex from the
   qptiff image name and hive-partitions everything by it.

2. But the join key never reaches the Xenium data. No Xenium h5ad carries ``obs['donor_id']``
   or ``obs['roi']`` — Xenium artifacts are keyed purely by the XETG slide serial
   (``0041323``), which encodes no donor, no region, and no disease status.

3. And ``xenium_paths.csv``, the one file that *does* map donors to Xenium runs, is read by
   nothing in the Xenium repo. It is hand-maintained and unenforced. Worse, the samples that
   actually have processed h5ads (``0041323``, ``0041326``) do not appear in it at all — the
   manifest cohort and the processed cohort are disjoint.

4. Finally, the slide serial is not a section identifier. One serial can carry both a
   pancreas and a pancreatic-lymph-node region (donor 6539 -> serial ``0059865`` with both a
   ``__Panc__`` and a ``__pLN__`` bundle), so keying anything by serial alone silently
   collides two different tissues.

So this module makes the manifest load-bearing rather than decorative: it normalises donor
ids to strings on both sides, keys Xenium sections by ``{serial}__{roi}`` instead of the
bare serial, classifies every row with a ``pair_status`` so unresolved cases are *visible*
rather than silently dropped, and provides a ``--stamp`` path that writes the join key into
the Xenium cell tables where it belongs.

Rows whose donor cannot be resolved (the processed-sample case) are not an error and do not
halt anything. They are recorded as ``donor_unknown`` and excluded from paired analyses
until someone supplies the mapping in ``donor_overrides.csv`` — a two-column file that is
the right place for knowledge that lives in a lab notebook rather than in the data.

CLI::

    python -m phenocycler.integration.manifest --report
    python -m phenocycler.integration.manifest --build
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from ..config import PipelineConfig, load_config
from .contract import normalize_donor_id

#: Columns of ``data/integration/manifest.csv``, in order.
MANIFEST_COLUMNS = [
    "donor_id",
    "roi",
    "disease_status",
    "pair_status",
    "pheno_image",
    "pheno_geojson",
    "xenium_sample_key",
    "xenium_serial",
    "xenium_batch",
    "xenium_bundle",
    "xenium_zarr",
    "xenium_h5ad",
    "section_gap_um",
    "block_id",
]

#: ``pair_status`` vocabulary.
PAIRED = "paired"                  # both modalities present for this donor
PHENO_ONLY = "pheno_only"          # PhenoCycler has it, Xenium does not
XENIUM_ONLY = "xenium_only"        # Xenium has it, PhenoCycler does not
DONOR_UNKNOWN = "donor_unknown"    # a Xenium run whose donor is not recorded anywhere

#: PhenoCycler is a pancreas-only pipeline; pLN Xenium regions are carried for provenance
#: but never paired against it.
PAIRABLE_ROI = "panc"

#: ``output-XETG00298__0059865__Panc__20250328__202338`` -> serial ``0059865``, label ``Panc``.
_BUNDLE_RE = re.compile(r"output-[A-Za-z0-9]+__(?P<serial>\d+)__(?P<label>.*?)__\d{8}__\d+/?$")


class ManifestError(RuntimeError):
    """The manifest cannot be built."""


@dataclass
class ManifestSummary:
    """Counts and diagnostics from a manifest build — what ``--report`` prints."""

    n_rows: int = 0
    by_status: dict[str, int] = field(default_factory=dict)
    paired_donors: list[str] = field(default_factory=list)
    pheno_only_donors: list[str] = field(default_factory=list)
    xenium_only_donors: list[str] = field(default_factory=list)
    unknown_keys: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    def render(self) -> str:
        lines = [
            "PhenoCycler <-> Xenium pairing manifest",
            "=" * 60,
            f"rows: {self.n_rows}",
        ]
        for status in (PAIRED, PHENO_ONLY, XENIUM_ONLY, DONOR_UNKNOWN):
            lines.append(f"  {status:<14} {self.by_status.get(status, 0)}")
        lines.append("")
        lines.append(f"paired donors ({len(self.paired_donors)}): "
                     f"{', '.join(self.paired_donors) or '-'}")
        if self.pheno_only_donors:
            lines.append(f"PhenoCycler only ({len(self.pheno_only_donors)}): "
                         f"{', '.join(self.pheno_only_donors)}")
        if self.xenium_only_donors:
            lines.append(f"Xenium only ({len(self.xenium_only_donors)}): "
                         f"{', '.join(self.xenium_only_donors)}")
        if self.unknown_keys:
            lines.append("")
            lines.append(f"UNRESOLVED Xenium runs ({len(self.unknown_keys)}) — no donor_id "
                         "recorded; add rows to donor_overrides.csv to pair them:")
            for key in self.unknown_keys:
                lines.append(f"    {key}")
        if self.warnings:
            lines.append("")
            lines.append("warnings:")
            for w in self.warnings:
                lines.append(f"    {w}")
        return "\n".join(lines)


def parse_bundle_path(path: str) -> tuple[str, str]:
    """``.../output-XETG00298__0059865__Panc__20250328__202338`` -> ``('0059865', 'Panc')``.

    Returns ``('', '')`` when the path does not look like a Xenium bundle. The serial is the
    only durable identifier in the whole Xenium layout — everything downstream of ingest
    keeps it and discards the rest — so it is worth recovering carefully.
    """
    text = str(path).strip().rstrip("/")
    if not text or text.upper() == "NA":
        return "", ""
    m = _BUNDLE_RE.search(text + "/")
    if not m:
        return "", ""
    return m.group("serial"), m.group("label")


def sample_key(serial: str, roi: str) -> str:
    """The Xenium section key: ``{serial}__{roi}``.

    Never key by serial alone. A single slide serial can carry a pancreas region and a
    lymph-node region, so ``data/raw/{serial}.zarr`` collides between two different tissues.
    """
    serial = str(serial).strip()
    roi = str(roi).strip()
    if not serial:
        return ""
    return f"{serial}__{roi}" if roi else serial


def normalize_roi(label: str) -> str:
    """Map a free-text Xenium bundle label to the ``roi`` vocabulary.

    The label field is entirely unnormalised upstream — ``Panc``, ``Pan``, ``Pancreas``,
    ``PanT``, ``6380_Pan``, ``pLN``, ``pLN1``, ``6457_pLN1``, ``Region_1`` all occur. Only
    the pancreas rows are pairable against PhenoCycler, so getting this mapping right is
    what decides whether a donor is integrable.
    """
    text = str(label).strip().lower()
    if not text:
        return ""
    if "pln" in text or "_ln" in text:
        m = re.search(r"pln[_]?(\d+)", text)
        return f"pln_{m.group(1)}" if m else "pln_1"
    if "pan" in text:
        return PAIRABLE_ROI
    return text


def read_xenium_paths(path: Path) -> pd.DataFrame:
    """Read ``xenium_paths.csv`` -> normalised rows.

    The vendored file has blank spacer rows between batch pairs and an ``NA`` source path
    for at least one run; both are dropped here rather than propagated.
    """
    cols = ["donor_id", "roi", "batch", "source_path", "xenium_serial", "xenium_sample_key"]
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=cols)

    try:
        raw = pd.read_csv(path, dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        # A zero-byte or header-only file is a legitimate state (nothing vendored yet), not
        # an error — build_manifest reports it as a warning and carries on.
        return pd.DataFrame(columns=cols)

    raw.columns = [c.strip() for c in raw.columns]
    if "donor_id" not in raw.columns:
        raise ManifestError(
            f"{path} has no 'donor_id' column (found {list(raw.columns)}); expected "
            f"donor_id,roi,batch,source_path")
    raw = raw[raw["donor_id"].astype(str).str.strip() != ""].copy()
    if raw.empty:
        return pd.DataFrame(columns=cols)

    raw["donor_id"] = raw["donor_id"].map(normalize_donor_id)
    raw["roi"] = raw["roi"].astype(str).str.strip()
    raw["batch"] = raw.get("batch", "").astype(str).str.strip()
    raw["source_path"] = raw["source_path"].astype(str).str.strip()
    raw = raw[~raw["source_path"].str.upper().isin({"", "NA"})].copy()

    parsed = raw["source_path"].map(parse_bundle_path)
    raw["xenium_serial"] = [p[0] for p in parsed]
    # Trust the manifest's own roi column; fall back to the bundle label when it is blank.
    labels = [p[1] for p in parsed]
    raw["roi"] = [r if r else normalize_roi(lab) for r, lab in zip(raw["roi"], labels)]
    raw["xenium_sample_key"] = [sample_key(s, r)
                                for s, r in zip(raw["xenium_serial"], raw["roi"])]
    return raw.reset_index(drop=True)


def read_donor_overrides(path: Path) -> pd.DataFrame:
    """Read ``donor_overrides.csv``: ``xenium_sample_key, donor_id, roi[, section_gap_um, block_id]``.

    This is where knowledge that exists only outside the data goes — most immediately, which
    nPOD donors ``0041323`` and ``0041326`` are. Overrides win over ``xenium_paths.csv``.
    """
    path = Path(path)
    cols = ["xenium_sample_key", "donor_id", "roi", "section_gap_um", "block_id"]
    if not path.exists():
        return pd.DataFrame(columns=cols)
    try:
        df = pd.read_csv(path, dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=cols)
    df.columns = [c.strip() for c in df.columns]
    for c in cols:
        if c not in df.columns:
            df[c] = ""
    df["donor_id"] = df["donor_id"].map(normalize_donor_id)
    df["xenium_sample_key"] = df["xenium_sample_key"].astype(str).str.strip()
    df["roi"] = df["roi"].astype(str).str.strip()
    return df.loc[:, cols].reset_index(drop=True)


def read_donor_metadata(cfg: PipelineConfig) -> pd.DataFrame:
    """``donor_id -> disease_status`` from the shared workbook.

    Both repos point at the same file with the same column names, so this is the one piece
    of donor-level metadata that needs no reconciliation — only a dtype fix, because the
    workbook stores donor ids numerically.
    """
    path = Path(cfg.donor_metadata)
    if not path.exists():
        return pd.DataFrame(columns=["donor_id", "disease_status"])
    try:
        df = pd.read_excel(path)
    except Exception as exc:  # noqa: BLE001 - missing openpyxl, corrupt file, ...
        print(f"[manifest] WARNING: could not read {path}: {exc}", flush=True)
        return pd.DataFrame(columns=["donor_id", "disease_status"])

    donor_col = cfg.metadata_donor_col
    status_col = cfg.metadata_status_col
    if donor_col not in df.columns:
        raise ManifestError(f"{path} has no {donor_col!r} column (got {list(df.columns)})")
    out = pd.DataFrame({
        "donor_id": df[donor_col].map(normalize_donor_id),
        "disease_status": (df[status_col].astype(str).str.strip()
                           if status_col in df.columns else ""),
    })
    return out.drop_duplicates(subset="donor_id").reset_index(drop=True)


def _pheno_paths(cfg: PipelineConfig, donor_id: str) -> tuple[str, str]:
    """Best-effort ``(qptiff, geojson)`` for a PhenoCycler donor.

    Both are looked up rather than constructed: REDSEA rebuilds the GeoJSON tag from the
    donor's ``image`` column and the upstream data README warns not to rename those files,
    so a glob is safer than a guess. Empty strings when the raw inputs are not mounted —
    the tabular stages do not need them and must not fail without them.
    """
    qptiff = ""
    images_dir = Path(cfg.images_dir)
    if images_dir.exists():
        hits = sorted(images_dir.glob(f"{donor_id}*.qptiff"))
        if hits:
            qptiff = str(hits[0])

    geojson = ""
    gdir = Path(cfg.geojson_dir)
    if gdir.exists():
        hits = sorted(gdir.glob(f"cells__{donor_id}*.geojson"))
        if hits:
            geojson = str(hits[0])
    return qptiff, geojson


def _rewrite_root(path: str, xenium_root: str) -> str:
    """Re-root a bundle path recorded on one filesystem onto another.

    ``xenium_root`` replaces everything up to and including ``.../xenium/``, so a manifest
    written against ``/orange/brusko/xenium/...`` works on a machine that mounts the same
    tree elsewhere without editing the vendored CSV.
    """
    if not xenium_root or not path:
        return path
    marker = "/xenium/"
    idx = path.find(marker)
    if idx == -1:
        return path
    return str(Path(xenium_root) / path[idx + len(marker):])


def build_manifest(
    cfg: PipelineConfig,
    *,
    pheno_donors: Optional[list[str]] = None,
    xenium_processed_root: Optional[Path] = None,
) -> tuple[pd.DataFrame, ManifestSummary]:
    """Build the pairing manifest.

    Parameters
    ----------
    cfg
        Resolved pipeline config (supplies both CSV paths, the donor workbook and
        ``xenium_root``).
    pheno_donors
        Override the PhenoCycler donor list; defaults to ``cfg.discover_donors()``, i.e.
        whatever ``data/cells/donor_id=*`` contains.
    xenium_processed_root
        Where ``{sample}/{sample}_phenotyped.h5ad`` and ``{sample}.zarr`` live, if mounted.
        Used only to fill the ``xenium_h5ad`` / ``xenium_zarr`` columns.

    Returns
    -------
    (DataFrame, ManifestSummary)
    """
    summary = ManifestSummary()

    donors = (list(pheno_donors) if pheno_donors is not None
              else cfg.discover_donors())
    pheno_set = {normalize_donor_id(d) for d in donors if normalize_donor_id(d)}

    xen = read_xenium_paths(cfg.xenium_paths_csv)
    overrides = read_donor_overrides(cfg.donor_overrides_csv)
    meta = read_donor_metadata(cfg)
    status_by_donor = dict(zip(meta["donor_id"], meta["disease_status"]))

    if xen.empty and overrides.empty:
        summary.warnings.append(
            f"no Xenium runs found — {cfg.xenium_paths_csv} is missing or empty, and "
            f"{cfg.donor_overrides_csv} has no rows"
        )

    rows: list[dict] = []
    seen_keys: set[str] = set()

    override_by_key = {r["xenium_sample_key"]: r for _, r in overrides.iterrows()
                       if r["xenium_sample_key"]}

    for _, r in xen.iterrows():
        key = r["xenium_sample_key"]
        ov = override_by_key.get(key, {})
        donor = normalize_donor_id(ov.get("donor_id") or r["donor_id"])
        roi = str(ov.get("roi") or r["roi"])
        bundle = _rewrite_root(r["source_path"], cfg.xenium_root)

        rows.append({
            "donor_id": donor,
            "roi": roi,
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": "",                       # filled below
            "pheno_image": "",
            "pheno_geojson": "",
            "xenium_sample_key": key,
            "xenium_serial": r["xenium_serial"],
            "xenium_batch": r["batch"],
            "xenium_bundle": bundle,
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": str(ov.get("section_gap_um", "") or ""),
            "block_id": str(ov.get("block_id", "") or ""),
        })
        seen_keys.add(key)

    # Override rows that name a Xenium sample absent from xenium_paths.csv. This is exactly
    # the 0041323 / 0041326 case: they have processed h5ads and no manifest row at all.
    for _, r in overrides.iterrows():
        key = r["xenium_sample_key"]
        if not key or key in seen_keys:
            continue
        serial = key.split("__", 1)[0]
        donor = normalize_donor_id(r["donor_id"])
        rows.append({
            "donor_id": donor,
            "roi": r["roi"] or (key.split("__", 1)[1] if "__" in key else ""),
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": "",
            "pheno_image": "",
            "pheno_geojson": "",
            "xenium_sample_key": key,
            "xenium_serial": serial,
            "xenium_batch": "",
            "xenium_bundle": "",
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": str(r.get("section_gap_um", "") or ""),
            "block_id": str(r.get("block_id", "") or ""),
        })
        seen_keys.add(key)

    df = pd.DataFrame(rows, columns=MANIFEST_COLUMNS) if rows else pd.DataFrame(
        columns=MANIFEST_COLUMNS)

    # Fill processed-artifact paths where the Xenium outputs are mounted.
    if xenium_processed_root is not None and len(df):
        root = Path(xenium_processed_root)
        h5ads, zarrs = [], []
        for serial in df["xenium_serial"]:
            h = root / serial / f"{serial}_phenotyped.h5ad"
            z = root.parent / "raw" / f"{serial}.zarr"
            h5ads.append(str(h) if h.exists() else "")
            zarrs.append(str(z) if z.exists() else "")
        df["xenium_h5ad"] = h5ads
        df["xenium_zarr"] = zarrs

    # PhenoCycler-side paths + pairing status.
    xen_donors = {d for d in df["donor_id"] if d}
    if len(df):
        statuses, images, geojsons = [], [], []
        for _, r in df.iterrows():
            donor, roi = r["donor_id"], r["roi"]
            if not donor:
                statuses.append(DONOR_UNKNOWN)
                images.append("")
                geojsons.append("")
                continue
            if donor in pheno_set and roi == PAIRABLE_ROI:
                q, g = _pheno_paths(cfg, donor)
                statuses.append(PAIRED)
                images.append(q)
                geojsons.append(g)
            else:
                statuses.append(XENIUM_ONLY)
                images.append("")
                geojsons.append("")
        df["pair_status"] = statuses
        df["pheno_image"] = images
        df["pheno_geojson"] = geojsons

    # PhenoCycler donors with no Xenium run at all get their own rows, so the manifest is a
    # complete picture of the cohort rather than a Xenium-shaped view of it.
    for donor in sorted(pheno_set - xen_donors):
        q, g = _pheno_paths(cfg, donor)
        df = pd.concat([df, pd.DataFrame([{
            "donor_id": donor,
            "roi": PAIRABLE_ROI,
            "disease_status": status_by_donor.get(donor, ""),
            "pair_status": PHENO_ONLY,
            "pheno_image": q,
            "pheno_geojson": g,
            "xenium_sample_key": "",
            "xenium_serial": "",
            "xenium_batch": "",
            "xenium_bundle": "",
            "xenium_zarr": "",
            "xenium_h5ad": "",
            "section_gap_um": "",
            "block_id": "",
        }])], ignore_index=True)

    if len(df):
        df = df.sort_values(["donor_id", "roi", "xenium_sample_key"],
                            kind="stable").reset_index(drop=True)

    # Collision guard. If this ever fires, two different sections are sharing an identity and
    # every downstream artifact for them would overwrite the other.
    if len(df):
        keys = df.loc[df["xenium_sample_key"] != "", "xenium_sample_key"]
        dupes = sorted(keys[keys.duplicated()].unique())
        if dupes:
            raise ManifestError(
                f"duplicate xenium_sample_key values {dupes} — two Xenium runs share a "
                "(serial, roi) identity, which would silently overwrite one another"
            )

    summary.n_rows = len(df)
    summary.by_status = df["pair_status"].value_counts().to_dict() if len(df) else {}
    summary.paired_donors = sorted(set(df.loc[df["pair_status"] == PAIRED, "donor_id"]))
    summary.pheno_only_donors = sorted(set(df.loc[df["pair_status"] == PHENO_ONLY, "donor_id"]))
    summary.xenium_only_donors = sorted(
        set(df.loc[df["pair_status"] == XENIUM_ONLY, "donor_id"]) - set(summary.paired_donors)
    )
    summary.unknown_keys = sorted(
        df.loc[df["pair_status"] == DONOR_UNKNOWN, "xenium_sample_key"])

    if summary.unknown_keys:
        summary.warnings.append(
            f"{len(summary.unknown_keys)} Xenium run(s) have no donor_id. They are excluded "
            f"from paired analyses. Add them to {cfg.donor_overrides_csv} to include them."
        )
    if not summary.paired_donors:
        summary.warnings.append(
            "no paired donors — donor-level, niche and PhenoCycler-only outputs still build, "
            "but registration, islet matching and pseudo-cell linking have nothing to run on"
        )

    return df, summary


def write_manifest(cfg: PipelineConfig, df: pd.DataFrame) -> Path:
    out = cfg.manifest_csv
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    return out


def load_manifest(cfg: PipelineConfig) -> pd.DataFrame:
    """Read the built manifest, rebuilding it on demand if absent."""
    if cfg.manifest_csv.exists():
        df = pd.read_csv(cfg.manifest_csv, dtype=str).fillna("")
        df["donor_id"] = df["donor_id"].map(normalize_donor_id)
        return df
    df, _ = build_manifest(cfg)
    return df


def paired_rows(df: pd.DataFrame, roi: str = PAIRABLE_ROI) -> pd.DataFrame:
    """The subset that can actually be integrated: both modalities, pairable ROI."""
    sub = df[df["pair_status"] == PAIRED]
    if roi:
        sub = sub[sub["roi"] == roi]
    return sub.reset_index(drop=True)


def run_manifest(cfg: PipelineConfig, *, write: bool = True,
                 xenium_processed_root: Optional[Path] = None) -> pd.DataFrame:
    """Stage entry point used by ``pipeline.py``."""
    df, summary = build_manifest(cfg, xenium_processed_root=xenium_processed_root)
    print(summary.render(), flush=True)
    if write:
        path = write_manifest(cfg, df)
        print(f"\n[manifest] wrote {path} ({len(df)} rows)", flush=True)
    return df


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Build/report the PhenoCycler<->Xenium manifest")
    ap.add_argument("--config", default=None, help="path to config.ini")
    ap.add_argument("--report", action="store_true",
                    help="print the pairing table without writing it")
    ap.add_argument("--build", action="store_true", help="write data/integration/manifest.csv")
    ap.add_argument("--xenium-processed", default=None,
                    help="Xenium data/processed root, to resolve h5ad/zarr paths")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    root = Path(args.xenium_processed) if args.xenium_processed else None
    df, summary = build_manifest(cfg, xenium_processed_root=root)
    print(summary.render())
    if args.build or not args.report:
        path = write_manifest(cfg, df)
        print(f"\nwrote {path} ({len(df)} rows)")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
