"""
Section keys — what one PhenoCycler partition actually is.

The core pipeline has always been *section*-oriented, not donor-oriented: one
``data/cells/donor_id=*`` partition is one qptiff, one QuPath GeoJSON, one RESTORE
threshold set, one coordinate frame. ``redsea.py`` makes this explicit by reading
``image.iloc[0]`` from a partition and using that single image's GeoJSON to mask every cell
in it, and ``restore.py`` labels its per-image thresholds by the partition name.

That was fine while the cohort was one image per donor. It is not fine now: the cohort is
20 donors x {pancreas, pancreatic lymph node}, and both images of a donor are named for the
same donor::

    6539_Scan1.er.qptiff        pancreas
    6539pLN_Scan1.er.qptiff     pancreatic lymph node

The old key was ``regexp_extract("Image", '[0-9]+')`` — the first digit-run — so both
images landed in ``donor_id=6539``. Two sections in one partition is not a labelling
inconvenience; it is silently wrong. REDSEA would mask both sections' cells with one
section's polygons, RESTORE would compute one threshold across two images, and the two
sections' micron coordinates — each measured from its own slide origin — would be
interleaved into a single overlapping point cloud.

So the partition key becomes the **section key**: the donor's digits plus whatever region
token the filename carries, lowercased.

    6539_Scan1.er.qptiff      -> 6539       -> donor 6539, roi panc
    6539pLN_Scan1.er.qptiff   -> 6539pln    -> donor 6539, roi pln_1
    6539pLN2_Scan1.er.qptiff  -> 6539pln2   -> donor 6539, roi pln_2

This module is the single definition of that key, shared by the DuckDB expression that
builds the partitions (``cells_parquet``) and the integration layer that takes them apart
(``integration/export_pheno``). Keeping both sides on one definition is the point: they are
in different languages and would otherwise drift.

An unrecognised region token is an error, not a pancreas. ``6539spleen`` gets its own
partition and :func:`parse` refuses it, rather than being silently folded in with the
pancreas section and averaged into its statistics.
"""

from __future__ import annotations

import re
from typing import NamedTuple

#: The DuckDB expression that derives the partition key from QuPath's ``Image`` column.
#: Digits, then any immediately-following alpha(+digit) token, lowercased.
#:
#: ``_`` is not in the character class, which is what makes ``6539_Scan1`` stop at the
#: donor id while ``6539pLN_Scan1`` carries its region through. Lowercasing means a stray
#: ``6539PLN`` cannot become a second partition for the same section — a split that would
#: halve a section's cell count with nothing to show for it.
SECTION_KEY_SQL = "lower(regexp_extract(\"Image\", '^[0-9]+[A-Za-z]*[0-9]*'))"

#: ``6539pln2`` -> donor ``6539``, region ``pln``, index ``2``.
_SECTION_RE = re.compile(r"^(?P<donor>\d+)(?:(?P<region>[a-z]+)(?P<index>\d*))?$")

#: The Python twin of :data:`SECTION_KEY_SQL`. Kept adjacent so the two cannot drift: the
#: DuckDB expression builds the partitions and this resolves an image *file* to the same key.
_IMAGE_KEY_RE = re.compile(r"^[0-9]+[A-Za-z]*[0-9]*")


def key_from_image(image: str) -> str:
    """``'6539pLN_Scan1.er.qptiff - resolution #1'`` -> ``'6539pln'``.

    Must agree with :data:`SECTION_KEY_SQL` for every input, which is why the two patterns
    sit next to each other. They are checked against each other in the test suite.
    """
    m = _IMAGE_KEY_RE.match(str(image).strip())
    return m.group(0).lower() if m else ""

#: Which physical section an image is of. Both values are PhenoCycler *images* — the
#: distinction is the slide under the objective, not the instrument.
SLIDE_PHENO = "pheno"      # the PhenoCycler section of the block
SLIDE_XENIUM = "xenium"    # the Xenium section, re-stained on the PhenoCycler afterwards

#: Region token -> ``(roi, slide)``.
#:
#: **roi.** The Xenium manifest names lymph-node regions ``pln_1..pln_3``; PhenoCycler images
#: carry a bare ``pLN``, so a token with no index resolves to ``pln_1``. That is the region
#: every row in ``xenium_paths.csv`` has (26 of the 26 donors listed there, three with a
#: second and one a third) — but that CSV is the vendored upstream manifest, not the run
#: cohort, which is 20 donors and is pinned nowhere. So it is an inference from a superset,
#: and wrong for any donor whose single PhenoCycler lymph-node section is really that donor's
#: ``pln_2``. It would surface as a registration that will not align; fix it with a per-donor
#: entry here rather than by loosening the mapping.
#:
#: **slide.** Post-Xenium INS/GCG staining is done on the PhenoCycler and produces an ordinary
#: qptiff, so it flows through the whole core pipeline unchanged — but it images the *Xenium*
#: section, which makes its cells the same cells as the Xenium transcriptome. That is a
#: genuine same-slide relationship sitting inside a serial-section cohort, and it is the only
#: place in this dataset where protein and RNA are measured on one cell. Keeping the slide in
#: the section key is what stops those cells being written into the PhenoCycler section's
#: partition, where they would be silently averaged with a different piece of tissue.
#:
#: An unlisted token is REFUSED by :func:`parse`, not guessed at. If the post-Xenium images in
#: your cohort use a token that is not here, the error names the file and adding one line
#: fixes it — that is deliberately louder than defaulting.
REGION_TO_ROI: dict[str, tuple[str, str]] = {
    "pln": ("pln_1", SLIDE_PHENO),
    "ln": ("pln_1", SLIDE_PHENO),
    # Post-Xenium re-stain of the Xenium section. `xen` = pancreas, `xenpln` = lymph node.
    "xen": ("panc", SLIDE_XENIUM),
    "xpanc": ("panc", SLIDE_XENIUM),
    "xenpln": ("pln_1", SLIDE_XENIUM),
    "xpln": ("pln_1", SLIDE_XENIUM),
}

#: The ROI of a section whose filename carries no region token at all.
DEFAULT_ROI = "panc"


class SectionError(ValueError):
    """A partition key does not name a section this cohort understands."""


class Section(NamedTuple):
    """One physical section: the partition it lives in, and what it is."""

    key: str          # 6539pln    — the data/cells/donor_id=* partition
    donor_id: str     # 6539       — the join key against the donor workbook
    roi: str          # pln_1      — the ROI vocabulary the Xenium manifest uses
    slide: str = SLIDE_PHENO   # which physical section this image is of

    @property
    def is_default_region(self) -> bool:
        return self.roi == DEFAULT_ROI

    @property
    def on_xenium_slide(self) -> bool:
        """True for a post-Xenium re-stain — the same cells as the Xenium transcriptome.

        Consumers use this to decide whether a pairing against Xenium is a *measurement*
        (same slide) or an approximation (serial sections).
        """
        return self.slide == SLIDE_XENIUM


def parse(key: str) -> Section:
    """Partition key -> :class:`Section`. Raises on anything unrecognised.

    Refusing an unknown region token is deliberate. Folding ``6539spleen`` into the pancreas
    would put a third tissue's cells into pancreatic composition, disease trends and islet
    statistics with no error and no warning — the failure mode this whole key exists to
    prevent, reintroduced one level up.
    """
    text = str(key).strip().lower()
    if not text:
        raise SectionError("empty section key")
    m = _SECTION_RE.match(text)
    if not m:
        raise SectionError(
            f"section key {key!r} is not <digits>[<region><index>] — expected e.g. '6539' "
            f"(pancreas) or '6539pln' (lymph node), from image names like "
            f"'6539_Scan1.er.qptiff' / '6539pLN_Scan1.er.qptiff'")

    donor = m.group("donor")
    region = m.group("region") or ""
    index = m.group("index") or ""
    if not region:
        return Section(key=text, donor_id=donor, roi=DEFAULT_ROI, slide=SLIDE_PHENO)

    entry = REGION_TO_ROI.get(region)
    base = entry[0] if entry else None
    if base is None:
        raise SectionError(
            f"section key {key!r} carries an unknown region token {region!r}; known tokens "
            f"are {sorted(REGION_TO_ROI)} (plus no token at all for {DEFAULT_ROI}). If this "
            f"is a real region — including a post-Xenium re-stain under a token not listed "
            f"there — add one line to phenocycler.sections.REGION_TO_ROI. Do not let it fall "
            f"through to {DEFAULT_ROI}: that would mix two sections in one partition.")
    roi = f"{base.rsplit('_', 1)[0]}_{index}" if index else base
    return Section(key=text, donor_id=donor, roi=roi, slide=entry[1])


def key_for(donor_id: str, roi: str = DEFAULT_ROI, slide: str = SLIDE_PHENO) -> str:
    """The inverse of :func:`parse`: ``('6539', 'pln_2')`` -> ``'6539pln2'``."""
    donor = str(donor_id).strip()
    roi = str(roi).strip().lower() or DEFAULT_ROI
    if roi == DEFAULT_ROI and slide == SLIDE_PHENO:
        return donor
    base = roi.partition("_")[0]
    index = roi.partition("_")[2]
    token = next((tok for tok, (r, sl) in REGION_TO_ROI.items()
                  if sl == slide and r.partition("_")[0] == base), None)
    if token is None:
        raise SectionError(f"no region token maps to roi {roi!r} on the {slide} slide")
    return f"{donor}{token}{index if index not in ('', '1') else ''}"


def donor_id(key: str) -> str:
    """Just the donor. Convenience for callers that do not need the ROI."""
    return parse(key).donor_id


def group_by_donor(keys) -> dict[str, list[Section]]:
    """``['6539', '6539pln', '6414'] -> {'6539': [panc, pln_1], '6414': [panc]}``.

    Sections are ordered ``panc`` first, then lymph nodes by index, so a report iterating
    this reads in the same order for every donor.
    """
    out: dict[str, list[Section]] = {}
    for k in keys:
        s = parse(k)
        out.setdefault(s.donor_id, []).append(s)
    for sections in out.values():
        sections.sort(key=lambda s: (s.roi != DEFAULT_ROI, s.roi))
    return out


def describe(keys) -> str:
    """One-line-per-donor summary of which sections exist — used in stage logs."""
    grouped = group_by_donor(keys)
    lines = [f"{len(keys)} section(s) across {len(grouped)} donor(s):"]
    for donor in sorted(grouped):
        lines.append(f"    {donor}: " + ", ".join(f"{s.roi} ({s.key})" for s in grouped[donor]))
    return "\n".join(lines)
