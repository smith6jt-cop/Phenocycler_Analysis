"""Manifest: the donor/roi crosswalk, and the two blockers it exists to resolve.

1. A single XETG slide serial can carry a pancreas region AND a lymph-node region, so keying
   Xenium artifacts by serial alone collides two different tissues.
2. The Xenium samples that actually have processed h5ads carry no donor_id anywhere, so the
   manifest has to surface them as unresolved rather than dropping them.
"""

from __future__ import annotations

import pandas as pd
import pytest

from phenocycler import load_config
from phenocycler.integration.manifest import (
    DONOR_UNKNOWN, PAIRED, PHENO_ONLY, XENIUM_ONLY, ManifestError, build_manifest,
    normalize_roi, parse_bundle_path, paired_rows, read_xenium_paths, sample_key,
)

BUNDLE = ("/orange/brusko/xenium/HuPrime5kPanTissue/Y9QRJ6/MM044_Xenium_nPODBatch10/"
          "output-XETG00298__0059865__{label}__20250328__202338")


def _paths_csv(tmp_path, rows):
    p = tmp_path / "xenium_paths.csv"
    pd.DataFrame(rows).to_csv(p, index=False)
    return p


def _overrides_csv(tmp_path, rows):
    p = tmp_path / "donor_overrides.csv"
    cols = ["xenium_sample_key", "donor_id", "roi", "section_gap_um", "block_id"]
    pd.DataFrame(rows, columns=cols).to_csv(p, index=False)
    return p


def _cfg(tmp_path, paths_rows, override_rows=()):
    return load_config(
        data_dir=tmp_path / "data",
        xenium_paths_csv=_paths_csv(tmp_path, paths_rows),
        donor_overrides_csv=_overrides_csv(tmp_path, list(override_rows)),
        donor_metadata=tmp_path / "missing.xlsx",
    )


@pytest.mark.parametrize("label,expected_serial", [
    ("Panc", "0059865"), ("pLN", "0059865"), ("6380_Pan", "0059865"),
])
def test_parse_bundle_path(label, expected_serial):
    serial, got_label = parse_bundle_path(BUNDLE.format(label=label))
    assert serial == expected_serial
    assert got_label == label


def test_parse_bundle_path_rejects_junk():
    assert parse_bundle_path("NA") == ("", "")
    assert parse_bundle_path("") == ("", "")


@pytest.mark.parametrize("label,expected", [
    ("Panc", "panc"), ("Pan", "panc"), ("Pancreas", "panc"), ("PanT", "panc"),
    ("6380_Pan", "panc"), ("pLN", "pln_1"), ("pLN1", "pln_1"), ("pLN2", "pln_2"),
    ("6457_pLN1", "pln_1"),
])
def test_roi_normalisation(label, expected):
    """The Xenium bundle label field is entirely unnormalised upstream."""
    assert normalize_roi(label) == expected


def test_sample_key_disambiguates_one_serial_two_regions(tmp_path):
    """The collision that makes `data/raw/{serial}.zarr` unsafe."""
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "Y9QRJ6",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6539", "roi": "pln_1", "batch": "Y9QRJ6",
         "source_path": BUNDLE.format(label="pLN")},
    ]
    df = read_xenium_paths(_paths_csv(tmp_path, rows))
    assert df["xenium_serial"].nunique() == 1          # same slide
    assert df["xenium_sample_key"].nunique() == 2      # different sections
    assert set(df["xenium_sample_key"]) == {"0059865__panc", "0059865__pln_1"}


def test_sample_key_helper():
    assert sample_key("0059865", "panc") == "0059865__panc"
    assert sample_key("0059865", "") == "0059865"
    assert sample_key("", "panc") == ""


def test_pair_status_classification(tmp_path):
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6539", "roi": "pln_1", "batch": "b",
         "source_path": BUNDLE.format(label="pLN")},
        {"donor_id": "9999", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc").replace("0059865", "0011111")},
    ]
    cfg = _cfg(tmp_path, rows)
    df, summary = build_manifest(cfg, pheno_donors=["6539", "6414"])

    def status(donor, roi):
        sub = df[(df["donor_id"] == donor) & (df["roi"] == roi)]
        return sub["pair_status"].iloc[0]

    assert status("6539", "panc") == PAIRED
    # The cohort is 20 donors x {PhenoCycler, Xenium} x {pancreas, pLN}, so a lymph-node
    # section pairs exactly as a pancreas section does. An earlier revision hard-coded
    # PhenoCycler as pancreas-only and would have dropped every pLN pairing on the floor.
    assert status("6539", "pln_1") == PAIRED
    assert status("9999", "panc") == XENIUM_ONLY
    assert status("6414", "panc") == PHENO_ONLY

    assert summary.paired_donors == ["6539"]
    assert summary.pheno_only_donors == ["6414"]
    # Unscoped means the whole cohort — both of 6539's sections, not just the pancreas.
    assert paired_rows(df)["donor_id"].tolist() == ["6539", "6539"]
    assert sorted(paired_rows(df)["roi"]) == ["panc", "pln_1"]


def test_pairing_is_tissue_gated_not_roi_gated(tmp_path):
    """An ROI that maps to no known tissue is not pairable, however the donor pairs.

    `pair_status` asks "is this ROI's tissue one we run?", not "is this ROI literally
    `panc`?" — otherwise adding the lymph node means editing pairing logic. A stray
    `Region_1` label still has no PhenoCycler counterpart and must stay `xenium_only`.
    """
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6539", "roi": "pln_2", "batch": "b",
         "source_path": BUNDLE.format(label="pLN2").replace("0059865", "0011111")},
        {"donor_id": "6539", "roi": "region_1", "batch": "b",
         "source_path": BUNDLE.format(label="Region_1").replace("0059865", "0022222")},
    ]
    cfg = _cfg(tmp_path, rows)
    df, _ = build_manifest(cfg, pheno_donors=["6539"])
    by_roi = dict(zip(df["roi"], df["pair_status"]))
    assert by_roi["panc"] == PAIRED
    assert by_roi["pln_2"] == PAIRED
    assert by_roi["region_1"] == XENIUM_ONLY

    tissues = dict(zip(df["roi"], df["tissue"]))
    assert tissues["panc"] == "pancreas"
    assert tissues["pln_2"] == "pancreatic_lymph_node"
    assert tissues["region_1"] == ""


def test_paired_rows_selects_by_tissue_or_roi(tmp_path):
    """A per-tissue run filters by tissue; a combined run passes `roi=None` for everything."""
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6539", "roi": "pln_1", "batch": "b",
         "source_path": BUNDLE.format(label="pLN")},
        {"donor_id": "6539", "roi": "pln_2", "batch": "b",
         "source_path": BUNDLE.format(label="pLN2").replace("0059865", "0011111")},
    ]
    cfg = _cfg(tmp_path, rows)
    df, _ = build_manifest(cfg, pheno_donors=["6539"])

    # Unscoped is the whole cohort; a tissue narrows to its ROIs; an roi narrows to one.
    assert sorted(paired_rows(df)["roi"]) == ["panc", "pln_1", "pln_2"]
    assert paired_rows(df, tissue="pancreas")["roi"].tolist() == ["panc"]
    assert paired_rows(df, tissue="pancreatic_lymph_node")["roi"].tolist() == ["pln_1", "pln_2"]
    assert paired_rows(df, roi="pln_2")["roi"].tolist() == ["pln_2"]
    # roi beats tissue, matching `resolve_rois` and the CLI.
    assert paired_rows(df, roi="panc", tissue="pancreatic_lymph_node")["roi"].tolist() == ["panc"]


def test_pairing_is_exact_when_pheno_sections_are_known(tmp_path):
    """Knowing which PhenoCycler sections exist turns pairing from an assumption into a fact.

    Donor 6539 has both scans; 6414 has only a pancreas. Pairing on the donor alone would
    report 6414's lymph node as paired on the strength of its pancreas — a pairing that does
    not exist, which registration would then fail on for no discoverable reason.
    """
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6539", "roi": "pln_1", "batch": "b",
         "source_path": BUNDLE.format(label="pLN").replace("0059865", "0011111")},
        {"donor_id": "6414", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc").replace("0059865", "0022222")},
        {"donor_id": "6414", "roi": "pln_1", "batch": "b",
         "source_path": BUNDLE.format(label="pLN").replace("0059865", "0033333")},
    ]
    cfg = _cfg(tmp_path, rows)
    df, _ = build_manifest(cfg, pheno_sections=[("6539", "panc"), ("6539", "pln_1"),
                                                ("6414", "panc")])
    status = {(r["donor_id"], r["roi"]): r["pair_status"] for _, r in df.iterrows()}
    assert status[("6539", "panc")] == PAIRED
    assert status[("6539", "pln_1")] == PAIRED
    assert status[("6414", "panc")] == PAIRED
    assert status[("6414", "pln_1")] == XENIUM_ONLY, "6414's lymph node was never scanned"


def test_pheno_section_without_a_xenium_run_is_pheno_only_per_section(tmp_path):
    """A donor can pair on one tissue and be pheno-only on the other. Recording that per
    donor would lose the distinction entirely."""
    rows = [{"donor_id": "6539", "roi": "panc", "batch": "b",
             "source_path": BUNDLE.format(label="Panc")}]
    cfg = _cfg(tmp_path, rows)
    df, summary = build_manifest(cfg, pheno_sections=[("6539", "panc"), ("6539", "pln_1")])
    status = {(r["donor_id"], r["roi"]): r["pair_status"] for _, r in df.iterrows()}
    assert status[("6539", "panc")] == PAIRED
    assert status[("6539", "pln_1")] == PHENO_ONLY
    assert summary.paired_donors == ["6539"]


def test_donor_list_without_sections_still_pairs_on_tissue(tmp_path):
    """The fallback path: given only a donor list, every ROI of a known tissue is pairable.

    That is the most that can be said without section information, and it stays available so
    a caller with no exported partitions (or a test) is not forced to invent them.
    """
    rows = [{"donor_id": "6539", "roi": "pln_1", "batch": "b",
             "source_path": BUNDLE.format(label="pLN")}]
    cfg = _cfg(tmp_path, rows)
    df, _ = build_manifest(cfg, pheno_donors=["6539"])
    assert df.iloc[0]["pair_status"] == PAIRED


def test_unmapped_xenium_sample_is_surfaced_not_dropped(tmp_path):
    """0041323 / 0041326 have processed h5ads but appear in no manifest and carry no donor.

    They must show up as `donor_unknown` so the gap is visible, rather than being silently
    absent from the pairing table.
    """
    cfg = _cfg(tmp_path, [], [{"xenium_sample_key": "0041323__panc", "donor_id": "",
                               "roi": "panc", "section_gap_um": "", "block_id": ""}])
    df, summary = build_manifest(cfg, pheno_donors=["6539"])
    assert DONOR_UNKNOWN in set(df["pair_status"])
    assert summary.unknown_keys == ["0041323__panc"]
    assert any("donor_overrides" in w for w in summary.warnings)


def test_override_resolves_unknown_donor(tmp_path):
    """Supplying the mapping turns the unresolved row into a real pairing."""
    cfg = _cfg(tmp_path, [], [{"xenium_sample_key": "0041323__panc", "donor_id": "6476",
                               "roi": "panc", "section_gap_um": "5", "block_id": "B1"}])
    df, summary = build_manifest(cfg, pheno_donors=["6476"])
    assert summary.unknown_keys == []
    assert summary.paired_donors == ["6476"]
    row = df[df["xenium_sample_key"] == "0041323__panc"].iloc[0]
    assert row["pair_status"] == PAIRED
    assert row["section_gap_um"] == "5"


def test_override_wins_over_manifest(tmp_path):
    rows = [{"donor_id": "1111", "roi": "panc", "batch": "b",
             "source_path": BUNDLE.format(label="Panc")}]
    cfg = _cfg(tmp_path, rows, [{"xenium_sample_key": "0059865__panc", "donor_id": "6539",
                                 "roi": "panc", "section_gap_um": "", "block_id": ""}])
    df, _ = build_manifest(cfg, pheno_donors=["6539"])
    assert df[df["xenium_sample_key"] == "0059865__panc"]["donor_id"].iloc[0] == "6539"


def test_duplicate_sample_key_raises(tmp_path):
    """Two runs sharing a (serial, roi) identity would overwrite each other's outputs."""
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6593", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Pan")},      # same serial, same roi
    ]
    cfg = _cfg(tmp_path, rows)
    with pytest.raises(ManifestError, match="duplicate xenium_sample_key"):
        build_manifest(cfg, pheno_donors=["6539", "6593"])


def test_donor_id_dtype_is_normalised_across_sources(tmp_path):
    """xenium_paths.csv stores text, the workbook stores numbers; both must join."""
    rows = [{"donor_id": "6539.0", "roi": "panc", "batch": "b",
             "source_path": BUNDLE.format(label="Panc")}]
    cfg = _cfg(tmp_path, rows)
    df, summary = build_manifest(cfg, pheno_donors=[6539])
    assert summary.paired_donors == ["6539"]


def test_blank_and_na_rows_are_dropped(tmp_path):
    rows = [
        {"donor_id": "6539", "roi": "panc", "batch": "b",
         "source_path": BUNDLE.format(label="Panc")},
        {"donor_id": "6591", "roi": "pln_1", "batch": "b", "source_path": "NA"},
        {"donor_id": "", "roi": "", "batch": "", "source_path": ""},
    ]
    df = read_xenium_paths(_paths_csv(tmp_path, rows))
    assert len(df) == 1


def test_empty_inputs_warn_rather_than_crash(tmp_path):
    cfg = _cfg(tmp_path, [])
    df, summary = build_manifest(cfg, pheno_donors=[])
    assert len(df) == 0
    assert any("no Xenium runs" in w for w in summary.warnings)
