"""Orthogonal cell-state annotations.

Identity and state are deliberately separate. Proliferation markers do not have a biologically
mutually-exclusive reference, so they are never passed through marker calibration or used as cell
type predictors.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


PROLIFERATION_MARKERS = ("Ki67", "PCNA")


@dataclass(frozen=True)
class StateThreshold:
    donor_id: str
    marker: str
    threshold: float
    n_finite: int
    status: str
    method: str = "log10-mean-plus-k-sd-v1"


def proliferation_threshold(values, *, k: float = 3.0, min_n: int = 200) -> float:
    """Return a donor-local threshold for a graded proliferation marker.

    Proliferation remains an independent state annotation because it has no
    biologically mutually-exclusive reference. This threshold is not presented
    as a frequency-neutral identity call and never participates in typing.
    """
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) < min_n:
        return float("inf")
    transformed = np.log10(values + 1.0)
    return float(10 ** (transformed.mean() + k * transformed.std()) - 1.0)


def annotate_proliferation(
    cells: pd.DataFrame,
    *,
    donor_id: str,
    markers=PROLIFERATION_MARKERS,
    k: float = 3.0,
    min_n: int = 200,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Annotate available proliferation markers without changing cell identity."""
    required = {
        "object_id",
        "qc_analysis_eligible",
        "qc_estimation_eligible",
    }
    missing = sorted(required.difference(cells.columns))
    if missing:
        raise KeyError(f"cells are missing required columns: {missing}")
    analysis_eligible = cells["qc_analysis_eligible"]
    estimation_eligible = cells["qc_estimation_eligible"]
    for name, values in (
        ("qc_analysis_eligible", analysis_eligible),
        ("qc_estimation_eligible", estimation_eligible),
    ):
        if values.isna().any() or not values.isin([True, False, 0, 1]).all():
            raise ValueError(f"{name} must contain explicit booleans")
    analysis_mask = analysis_eligible.astype(bool).to_numpy()
    estimation_mask = estimation_eligible.astype(bool).to_numpy()
    out = pd.DataFrame(
        {
            "object_id": cells["object_id"].astype(str),
            "donor_id": str(donor_id),
            "qc_analysis_eligible": analysis_mask,
        }
    )
    models: list[StateThreshold] = []
    for marker in markers:
        if marker not in cells:
            out[f"{marker}__state"] = "unavailable"
            models.append(StateThreshold(str(donor_id), marker, np.nan, 0, "marker_unavailable"))
            continue
        values = pd.to_numeric(cells[marker], errors="coerce").to_numpy(dtype=float)
        finite = np.isfinite(values)
        n_finite_estimation = int((finite & estimation_mask).sum())
        threshold = proliferation_threshold(
            values[estimation_mask], k=k, min_n=min_n
        )
        available = np.isfinite(threshold)
        state = np.full(len(values), "unavailable", dtype=object)
        if available:
            callable_cells = finite & analysis_mask
            state[callable_cells & (values < threshold)] = "negative"
            state[callable_cells & (values >= threshold)] = "positive"
        out[f"{marker}__value"] = values
        out[f"{marker}__state"] = state
        models.append(
            StateThreshold(
                str(donor_id),
                marker,
                threshold,
                n_finite_estimation,
                (
                    "valid"
                    if available
                    else (
                        "marker_unavailable"
                        if not finite.any()
                        else "insufficient_measurements"
                    )
                ),
            )
        )
    return out, pd.DataFrame([m.__dict__ for m in models])
