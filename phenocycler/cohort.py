"""Cohort eligibility rules shared by every phenocycler analysis stage."""

from __future__ import annotations

from collections.abc import Iterable


DONOR_EXCLUSIONS: dict[str, str] = {
    "6457": (
        "Excluded until further notice: poor DAPI staining prevented accurate "
        "Phenocycler channel registration."
    ),
    "6579": (
        "Excluded until further notice by maintainer decision: pancreatitis outlier "
        "must not be included in analysis."
    ),
}


def ensure_eligible_donors(
    donors: Iterable[str],
    *,
    context: str = "analysis",
) -> tuple[str, ...]:
    """Return normalized donor IDs or reject any explicitly excluded donor."""
    normalized = tuple(str(donor) for donor in donors)
    blocked = sorted(set(normalized).intersection(DONOR_EXCLUSIONS))
    if blocked:
        details = "; ".join(
            f"{donor}: {DONOR_EXCLUSIONS[donor]}" for donor in blocked
        )
        raise ValueError(f"{context} includes excluded donor(s): {details}")
    return normalized


def filter_eligible_donors(donors: Iterable[str]) -> list[str]:
    """Remove centrally excluded donors from automatic cohort discovery."""
    return [
        str(donor)
        for donor in donors
        if str(donor) not in DONOR_EXCLUSIONS
    ]
