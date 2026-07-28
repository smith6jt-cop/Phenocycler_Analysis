"""
The modality-agnostic cell-table contract.

PhenoCycler and Xenium arrive in completely different shapes: hive-partitioned parquet
keyed by a QuPath UUID on one side, an AnnData/SpatialData object keyed by a Xenium cell id
on the other. Rather than teach every downstream module about both, S1a/S1b normalise each
modality into ONE flat parquet schema and everything after that reads only this.

That is what lets ``structures``, ``register``, ``match``, ``grid`` and ``qc`` be written
once instead of twice, and it is what makes them testable without any imaging data — a
synthetic frame with the right columns is a valid input.

Schema::

    cell_id          str      modality-native id. PhenoCycler: QuPath object_id (UUID).
                              Xenium: cell_id from cells.parquet / adata.obs_names.
    donor_id         str      the cross-modality join key. ALWAYS str (the donor metadata
                              workbook stores it numeric; PhenoCycler derives it by regex).
    roi              str      'panc' | 'pln_1' | ... PhenoCycler is pancreas-only.
    modality         str      'phenocycler' | 'xenium'
    x_um, y_um       float32  centroid in MICRONS in the modality's own frame.
    x_um_reg,        float32  centroid in microns in the FIXED frame, written by S5.
    y_um_reg                  NaN until registration runs.
    area_um2         float32  segmented cell area.
    lineage_native   str      the modality's own broad label, unmodified.
    lineage_common   str      the harmonised label from vocab.py.
    celltype_fine    str      finest available subtype.
    confidence       float32  assign_margin (PhenoCycler) | broad_confidence (Xenium).
    evidence_positive bool    False for cells assigned by a default/fallback rule rather
                              than by positive evidence. See below — this one matters.
    region           str      'core' | 'peri' | 'tissue' | 'other' where known.
    structure_id     str      islet id once structures.py has run, else ''.
    feat__<name>     float32  per-cell features. PhenoCycler: '<marker>_norm'.
                              Xenium: identity-core lineage scores.

``evidence_positive`` is not decoration. PhenoCycler's ``Epithelial`` is a *default sink* —
``lineage.py`` assigns it to every cell that fails all the positive gates (flagged there as
``epi_default``) — while Xenium's ``Exocrine`` is a positive call from an identity-core
score, and Xenium separately flags ``broad_lowsignal`` cells whose lineage came from an
argmax fallback rather than the classifier. Comparing raw composition fractions across that
asymmetry inflates PhenoCycler epithelium against Xenium exocrine for reasons that have
nothing to do with biology. Carrying one bool through the contract lets ``donor.py`` report
``frac_all`` and ``frac_evidence_positive`` side by side and headline the latter.
"""

from __future__ import annotations

import glob
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

MODALITIES = ("phenocycler", "xenium")

FEATURE_PREFIX = "feat__"

#: (column, numpy dtype) for every non-feature column, in canonical order.
CELL_TABLE_COLUMNS: dict[str, str] = {
    "cell_id": "object",
    "donor_id": "object",
    "roi": "object",
    "modality": "object",
    "x_um": "float32",
    "y_um": "float32",
    "x_um_reg": "float32",
    "y_um_reg": "float32",
    "area_um2": "float32",
    "lineage_native": "object",
    "lineage_common": "object",
    "celltype_fine": "object",
    "confidence": "float32",
    "evidence_positive": "bool",
    "region": "object",
    "structure_id": "object",
}

#: Columns a producer must supply; the rest are filled with typed nulls.
REQUIRED_COLUMNS = ("cell_id", "donor_id", "modality", "x_um", "y_um", "lineage_native")


class ContractError(ValueError):
    """A cell table violates the contract."""


@dataclass(frozen=True)
class CellTable:
    """A validated cell table plus the identity of what it describes.

    Kept as a thin wrapper over a DataFrame rather than a subclass: pandas subclassing is
    fragile across operations, and every consumer here wants the plain frame.
    """

    df: pd.DataFrame
    modality: str
    donor_id: str
    roi: str = ""

    @property
    def n_cells(self) -> int:
        return len(self.df)

    @property
    def features(self) -> list[str]:
        return feature_columns(self.df)

    def xy(self, registered: bool = False) -> np.ndarray:
        """(n, 2) float64 centroids, in microns.

        ``registered=True`` returns the fixed-frame coordinates and raises if S5 has not
        run — silently falling back to unregistered coordinates would produce a
        plausible-looking but meaningless overlay, which is exactly the failure mode of the
        prototype this replaces.
        """
        if registered:
            if "x_um_reg" not in self.df.columns or self.df["x_um_reg"].isna().all():
                raise ContractError(
                    f"{self.modality}/{self.donor_id}: registered coordinates requested but "
                    "x_um_reg is empty — run the register + transform stages first"
                )
            cols = ("x_um_reg", "y_um_reg")
        else:
            cols = ("x_um", "y_um")
        return self.df.loc[:, list(cols)].to_numpy(dtype=np.float64)


def feature_columns(df: pd.DataFrame) -> list[str]:
    """The ``feat__*`` columns, sorted, so feature order never depends on write order."""
    return sorted(c for c in df.columns if c.startswith(FEATURE_PREFIX))


def feature_name(column: str) -> str:
    """``feat__CD31`` -> ``CD31``."""
    return column[len(FEATURE_PREFIX):] if column.startswith(FEATURE_PREFIX) else column


def as_feature(name: str) -> str:
    """``CD31`` -> ``feat__CD31``."""
    return name if name.startswith(FEATURE_PREFIX) else f"{FEATURE_PREFIX}{name}"


def normalize_donor_id(value) -> str:
    """Coerce a donor id to the canonical string form.

    The donor metadata workbook stores ``6539`` as a number, ``xenium_paths.csv`` as text,
    and PhenoCycler derives it with a regex over the image name. Left alone, ``6539`` and
    ``6539.0`` and ``' 6539'`` are three different join keys and the merge silently produces
    zero paired donors. Everything that touches a donor id goes through here.
    """
    if value is None:
        return ""
    if isinstance(value, float):
        if np.isnan(value):
            return ""
        if float(value).is_integer():
            return str(int(value))
        return str(value)
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    text = str(value).strip()
    if text.endswith(".0") and text[:-2].isdigit():
        text = text[:-2]
    return text


def empty_cell_table() -> pd.DataFrame:
    """A zero-row frame with the full contract schema and correct dtypes."""
    return pd.DataFrame({c: pd.Series(dtype=d) for c, d in CELL_TABLE_COLUMNS.items()})


def coerce_cell_table(
    df: pd.DataFrame,
    *,
    modality: str,
    donor_id: str,
    roi: str = "",
) -> pd.DataFrame:
    """Fill missing optional columns, cast dtypes, and order columns canonically.

    Producers supply what they have; this fills the rest with typed nulls so every consumer
    can assume the full schema exists. Feature columns are preserved and sorted after the
    fixed columns.
    """
    if modality not in MODALITIES:
        raise ContractError(f"unknown modality {modality!r}; expected one of {MODALITIES}")

    out = df.copy()
    missing = [c for c in REQUIRED_COLUMNS if c not in out.columns and c not in ("modality",)]
    # donor_id/modality are supplied as arguments, so only complain about genuinely absent data
    missing = [c for c in missing if c not in ("donor_id",)]
    if missing:
        raise ContractError(f"cell table for {modality}/{donor_id} is missing {missing}")

    out["modality"] = modality
    out["donor_id"] = normalize_donor_id(donor_id)
    if "roi" not in out.columns or roi:
        out["roi"] = roi

    for col, dtype in CELL_TABLE_COLUMNS.items():
        if col not in out.columns:
            if dtype == "bool":
                out[col] = True          # absent evidence flag == positively assigned
            elif dtype.startswith("float"):
                out[col] = np.nan
            else:
                out[col] = ""
        if dtype.startswith("float"):
            out[col] = pd.to_numeric(out[col], errors="coerce").astype(dtype)
        elif dtype == "bool":
            out[col] = out[col].fillna(True).astype(bool)
        else:
            out[col] = out[col].fillna("").astype(str)

    feats = feature_columns(out)
    for c in feats:
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float32")

    # Anything a producer added beyond the contract is kept, not dropped. The contract is a
    # floor that every consumer can rely on, not a ceiling: `immune_subclass`,
    # `endocrine_subtype`, `xenium_sample_key` and `spatial_niche` are all modality-specific
    # provenance that would be expensive to recompute and is cheap to carry.
    extras = [c for c in out.columns if c not in CELL_TABLE_COLUMNS and c not in feats]
    return out.loc[:, list(CELL_TABLE_COLUMNS) + feats + sorted(extras)].reset_index(drop=True)


def validate_cell_table(df: pd.DataFrame, *, modality: Optional[str] = None) -> None:
    """Raise :class:`ContractError` if ``df`` is not a valid cell table.

    Checks the things that actually go wrong in practice: a missing column, coordinates
    that are secretly pixels, a mixed-modality frame, duplicate ids.
    """
    missing = [c for c in CELL_TABLE_COLUMNS if c not in df.columns]
    if missing:
        raise ContractError(f"cell table missing contract columns: {missing}")

    if len(df) == 0:
        return

    mods = set(df["modality"].unique())
    if not mods <= set(MODALITIES):
        raise ContractError(f"cell table has unknown modality values: {mods - set(MODALITIES)}")
    if len(mods) > 1:
        raise ContractError(f"cell table mixes modalities {mods}; keep one per file")
    if modality is not None and mods != {modality}:
        raise ContractError(f"expected modality {modality!r}, got {mods}")

    if df["cell_id"].duplicated().any():
        n = int(df["cell_id"].duplicated().sum())
        raise ContractError(f"cell table has {n} duplicate cell_id values")

    if df["donor_id"].eq("").any():
        raise ContractError("cell table has empty donor_id values")

    for col in ("x_um", "y_um"):
        if not np.isfinite(df[col].to_numpy()).any():
            raise ContractError(f"cell table column {col} has no finite values")

    # Both modalities are natively in microns — PhenoCycler from QuPath's 'Centroid X um',
    # Xenium from obsm['spatial']. A tissue section is at most ~4 cm, so a coordinate span
    # far beyond that means someone handed us pixels. Catching it here is cheap; catching it
    # after registration means debugging a silently wrong overlay.
    span = float(np.nanmax(df["x_um"].to_numpy()) - np.nanmin(df["x_um"].to_numpy()))
    if span > 100_000:
        raise ContractError(
            f"x_um spans {span:,.0f} um (>10 cm) — coordinates look like pixels, not microns"
        )


def write_cell_table(df: pd.DataFrame, path: Path) -> Path:
    """Write one partition, creating parents. Validates first — never persist a bad table."""
    validate_cell_table(df)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, index=False)
    return path


def partition_path(root: Path, donor_id: str, roi: str = "") -> Path:
    """Hive path for a donor (+roi), matching the core pipeline's ``donor_id=<id>`` layout.

    The Xenium side adds ``roi=<roi>`` because one donor can have a pancreas and up to three
    pancreatic lymph node regions, and (worse) a single Xenium slide serial can carry both —
    so donor alone does not identify a section.
    """
    base = Path(root) / f"donor_id={normalize_donor_id(donor_id)}"
    if roi:
        base = base / f"roi={roi}"
    return base / "data_0.parquet"


def read_cell_table(
    root: Path,
    donor_id: str,
    roi: str = "",
    *,
    modality: Optional[str] = None,
) -> CellTable:
    """Read one donor (+roi) partition and return a validated :class:`CellTable`."""
    root = Path(root)
    direct = partition_path(root, donor_id, roi)
    if direct.exists():
        files = [direct]
    else:
        pattern = str(root / f"donor_id={normalize_donor_id(donor_id)}" / "**" / "*.parquet")
        files = sorted(glob.glob(pattern, recursive=True))
        if roi:
            files = [f for f in files if f"roi={roi}" in f]
    if not files:
        raise FileNotFoundError(
            f"no cell table for donor {donor_id!r}"
            + (f" roi {roi!r}" if roi else "")
            + f" under {root}"
        )

    df = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    validate_cell_table(df, modality=modality)
    mod = str(df["modality"].iloc[0])
    got_roi = roi or (str(df["roi"].iloc[0]) if "roi" in df.columns and len(df) else "")
    return CellTable(df=df, modality=mod, donor_id=normalize_donor_id(donor_id), roi=got_roi)


def discover_partitions(root: Path) -> list[tuple[str, str]]:
    """List ``(donor_id, roi)`` pairs present under ``root``. ``roi`` is '' when unpartitioned."""
    root = Path(root)
    if not root.exists():
        return []
    out: set[tuple[str, str]] = set()
    for donor_dir in sorted(root.glob("donor_id=*")):
        donor = donor_dir.name.split("=", 1)[1]
        roi_dirs = sorted(donor_dir.glob("roi=*"))
        if roi_dirs:
            out.update((donor, d.name.split("=", 1)[1]) for d in roi_dirs)
        elif any(donor_dir.glob("*.parquet")):
            out.add((donor, ""))
    return sorted(out)


def concat_cell_tables(tables: Iterable[CellTable]) -> pd.DataFrame:
    """Concatenate cell tables, keeping the union of feature columns."""
    frames = [t.df for t in tables]
    if not frames:
        return empty_cell_table()
    out = pd.concat(frames, ignore_index=True, sort=False)
    for c in feature_columns(out):
        out[c] = pd.to_numeric(out[c], errors="coerce").astype("float32")
    return out
