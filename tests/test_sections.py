"""Section keys: the partition identity that keeps two of a donor's images apart.

The bug these guard against is silent and severe. Before section keys, ``6539_Scan1.er.qptiff``
and ``6539pLN_Scan1.er.qptiff`` both reduced to ``donor_id=6539``, so one partition held two
physical sections — REDSEA would mask both with one section's polygons, RESTORE would fit one
threshold across two images, and two slide-origin coordinate frames would be interleaved into
one overlapping point cloud. Nothing would raise.
"""

from __future__ import annotations

import pytest

from phenocycler import sections as S

#: Real filename shapes from the cohort, and the partition each must land in.
IMAGES = [
    ("6539_Scan1.er.qptiff - resolution #1", "6539", "6539", "panc"),
    ("6539pLN_Scan1.er.qptiff - resolution #1", "6539pln", "6539", "pln_1"),
    ("6539pLN2_Scan1.er.qptiff - resolution #1", "6539pln2", "6539", "pln_2"),
    ("112_Scan1.er.qptiff", "112", "112", "panc"),
    ("6374pln_Scan1.er.qptiff", "6374pln", "6374", "pln_1"),
]


@pytest.mark.parametrize("image,key,donor,roi", IMAGES)
def test_image_resolves_to_its_section(image, key, donor, roi):
    assert S.key_from_image(image) == key
    sec = S.parse(key)
    assert (sec.donor_id, sec.roi) == (donor, roi)


def test_a_donors_two_images_never_share_a_partition():
    """The regression, stated directly."""
    panc = S.key_from_image("6539_Scan1.er.qptiff - resolution #1")
    pln = S.key_from_image("6539pLN_Scan1.er.qptiff - resolution #1")
    assert panc != pln
    assert S.parse(panc).donor_id == S.parse(pln).donor_id == "6539"


def test_sql_and_python_agree_on_every_shape():
    """The DuckDB expression builds the partitions; the Python parser takes them apart.

    They are written in different languages against the same rule, so they can drift — and a
    drift would mean the integration layer looking for partitions the core pipeline never
    wrote. Checked directly rather than trusted.
    """
    duckdb = pytest.importorskip("duckdb")

    rows = ", ".join(f"('{img}')" for img, _k, _d, _r in IMAGES)
    q = (f'SELECT "Image", {S.SECTION_KEY_SQL} AS key '
         f'FROM (VALUES {rows}) AS t("Image")')
    for image, sql_key in duckdb.sql(q).fetchall():
        assert sql_key == S.key_from_image(image), f"SQL/Python disagree on {image!r}"


def test_case_variants_collapse_to_one_partition():
    """`6539PLN` and `6539pLN` are one section. Two partitions would halve its cell count
    with nothing to indicate the split had happened."""
    assert S.key_from_image("6539PLN_Scan1.er.qptiff") == "6539pln"
    assert S.key_from_image("6539pln_Scan1.er.qptiff") == "6539pln"
    assert S.parse("6539PLN") == S.parse("6539pln")


@pytest.mark.parametrize("donor,roi", [
    ("6539", "panc"), ("6539", "pln_1"), ("6539", "pln_2"), ("112", "panc"),
])
def test_key_round_trip(donor, roi):
    sec = S.parse(S.key_for(donor, roi))
    assert (sec.donor_id, sec.roi) == (donor, roi)


def test_unknown_region_is_refused_not_folded_into_pancreas():
    """An unrecognised token must fail loudly. Defaulting it to pancreas would mix a third
    tissue into pancreatic composition, disease trends and islet statistics silently."""
    with pytest.raises(S.SectionError, match="unknown region token"):
        S.parse("6539spleen")
    with pytest.raises(S.SectionError, match="not <digits>"):
        S.parse("abc")
    with pytest.raises(S.SectionError, match="empty"):
        S.parse("")


def test_group_by_donor_orders_pancreas_first():
    grouped = S.group_by_donor(["6579pln2", "6539", "6579pln", "6539pln", "6414"])
    assert sorted(grouped) == ["6414", "6539", "6579"]
    assert [s.roi for s in grouped["6539"]] == ["panc", "pln_1"]
    assert [s.roi for s in grouped["6579"]] == ["pln_1", "pln_2"]     # no pancreas scanned
    assert [s.roi for s in grouped["6414"]] == ["panc"]


def test_describe_lists_every_section():
    text = S.describe(["6539", "6539pln"])
    assert "2 section(s) across 1 donor(s)" in text
    assert "panc (6539)" in text and "pln_1 (6539pln)" in text


def test_config_separates_sections_from_donors(tmp_path):
    """`discover_sections` is the work unit; `discover_donors` is the cohort.

    Conflating them makes `6539pln` a phantom donor with no workbook metadata and no Xenium
    counterpart, which then surfaces in the manifest as a pairing gap that does not exist.
    """
    from phenocycler import load_config

    cfg = load_config()
    cfg.data_dir = tmp_path
    for key in ("6539", "6539pln", "6414"):
        (cfg.cells_dir / f"donor_id={key}").mkdir(parents=True, exist_ok=True)

    assert cfg.discover_sections() == ["6414", "6539", "6539pln"]
    assert cfg.discover_donors() == ["6414", "6539"]
    assert {(s.donor_id, s.roi) for s in cfg.discover_section_objects()} == {
        ("6414", "panc"), ("6539", "panc"), ("6539", "pln_1")}
