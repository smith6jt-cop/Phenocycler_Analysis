#!/usr/bin/env python3
"""
RESTORE production apply — the Gate-4 adapter from the frozen validation method to per-cell calls.

The manuscript-faithful RESTORE method is FROZEN and validated in ``restore_validation.py`` (screen +
Gate-3 pilot on all cells). That module deliberately emits only screen/pilot tables — no normalized values
and no per-cell calls. ``lineage.py``, in turn, consumes a per-donor "gated parquet" of per-marker
calls. This module is the production adapter between them:

  * run the frozen ``evaluate_locked_pair`` on ALL cells (the modulus-1, unsampled compartment sample) for each of
    the 7 ``ACCEPTED_PAIRS`` (one predeclared reference per target, maintainer sign-off 2026-07-24), per
    eligible donor;
  * turn each pair's evaluation into a TRI-STATE per-cell call — ``positive`` / ``negative`` /
    ``unavailable`` — via the pure, unit-testable :func:`cell_calls_from_evaluation`;
  * write ONE gated parquet per donor (``{target}_pos`` / ``{target}_norm`` / ``{target}_state`` for the
    7 accepted targets), plus a per-donor + a cohort QC-disposition manifest.

Call semantics (faithful to the frozen method, no re-implementation of its math):

  * VALID pair  ->  cell is ``positive`` iff its REDSEA target value exceeds the canonical divisor
    (``full_lineage_positive``); ``negative`` if measured but at/below the divisor; ``unavailable`` if it
    lacks a finite nonnegative required measurement for that pair (off ``valid_idx``).
  * NO_TARGET_EXPRESSION_CONFIRMED (a maintainer-confirmed biological absence) -> every cell ``negative``
    (a confident absence, NOT ``unavailable``): the marker cannot pull cells into its compartment, and a
    cell falls through to a lower gate / Other normally.
  * Any other (technical) non-VALID state on an ACCEPTED pair -> HALT LOUDLY. The frozen pairs are
    cohort-validated and their anchor markers are universally expressed, so a technical failure is a
    method/data problem to surface, never a biological negative to absorb. No automatic fallback
    (docs/restore_faithful_rebuild_plan.md, resolved decision 4).

Inputs  : <restore_allcells_sample_dir>/donor_id=*/*.parquet   (modulus-1 QuPath compartment means)
          <cells_dir>/donor_id=*/*.parquet                     (canonical raw whole-cell means)
          <cells_redsea_dir>/donor_id=*/*.parquet              (REDSEA whole-cell target values)
          <restore_pair_reviews_csv>                           (ACCEPTED target/reference rules)
          <restore_expression_reviews_csv>                     (optional per-(donor,target) reviews)
Outputs : <restore_gated_dir>/donor_id=*/data_0.parquet        ({m}_pos, {m}_norm, {m}_state)
          <restore_gated_dir>/donor_id=*/apply_manifest.json    (per-donor divisor/state disposition)
          <restore_gated_dir>/apply_manifest.csv                (cohort donor x pair QC disposition)

    python -m phenocycler.restore_apply --donors 6539 6476 --jobs 4
"""

from __future__ import annotations

import argparse
import functools
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from .config import PipelineConfig, load_config
from .cohort import ensure_eligible_donors
from .parallel import map_donors
from .restore_validation import (
    ACCEPTED_PAIRS,
    CONTROL_DEFINITIONS,
    DIVISOR_STATISTICS,
    METHOD_VERSION,
    ExpressionReview,
    PairReview,
    PairState,
    PairValidationConfig,
    evaluate_locked_pair,
    load_expression_reviews,
    load_pair_reviews,
    load_validation_sample,
    locked_pair_metrics,
)

# tri-state per-cell call encoding (widest label is 11 chars -> "<U11")
CELL_POSITIVE = "positive"
CELL_NEGATIVE = "negative"
CELL_UNAVAILABLE = "unavailable"
_CELL_STATE_DTYPE = "<U11"

# The two pair states that emit a per-cell call. VALID emits divisor-gated positive/negative;
# a maintainer-confirmed absence emits a confident negative for every cell (E-1 decision 2026-07-24).
CALLABLE_STATES = frozenset({PairState.VALID, PairState.NO_TARGET_EXPRESSION_CONFIRMED})
# Every other (technical) state on an ACCEPTED pair halts the run — no automatic fallback.
HARD_FAIL_STATES = frozenset(PairState) - CALLABLE_STATES
# Guard against a newly added PairState silently being treated as callable.
assert CALLABLE_STATES | HARD_FAIL_STATES == frozenset(PairState)

# the 7 target markers emitted into the gated parquet (one per accepted pair), in ACCEPTED_PAIRS order
REQUIRED_TARGETS: tuple[str, ...] = tuple(target for target, _ in ACCEPTED_PAIRS)


def config_from_pipeline(cfg: PipelineConfig) -> PairValidationConfig:
    """Build the frozen-method config from the pipeline config's ``[restore]`` policy keys.

    Only the two METHOD_VERSION v10 policy knobs are configurable; every other
    :class:`PairValidationConfig` field is part of the locked specification and is deliberately not
    exposed to ``config.ini``. Validating here means a typo in the ini fails before any donor is
    evaluated, with the allowed values named, rather than surfacing as a mid-run ``KeyError`` inside
    ``negative_control_statistics``.
    """
    control = str(cfg.restore_control_definition)
    statistic = str(cfg.restore_divisor_statistic)
    if control not in CONTROL_DEFINITIONS:
        raise SystemExit(
            f"[restore] config.ini [restore] control_definition={control!r} is not one of "
            f"{sorted(CONTROL_DEFINITIONS)}"
        )
    if statistic not in DIVISOR_STATISTICS:
        raise SystemExit(
            f"[restore] config.ini [restore] divisor_statistic={statistic!r} is not one of "
            f"{sorted(DIVISOR_STATISTICS)}"
        )
    return PairValidationConfig(control_definition=control, divisor_statistic=statistic)


def cell_calls_from_evaluation(ev: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Map one :func:`evaluate_locked_pair` result to aligned per-cell ``(pos, norm, state)`` arrays.

    Pure and total over every :class:`PairState` (so it is directly unit-testable). The arrays are
    positional to ``donor_df`` rows (length ``ev['n_input']``) and never fabricate a call for a
    non-callable state. The driver decides whether a state is allowed to reach here (see
    :data:`HARD_FAIL_STATES`); this function only encodes the mapping.
    """
    n = int(ev["n_input"])
    state = PairState(ev["state"])
    pos = np.zeros(n, dtype=bool)
    norm = np.full(n, np.nan, dtype=float)
    cell_state = np.empty(n, dtype=_CELL_STATE_DTYPE)

    if state is PairState.VALID:
        available = np.zeros(n, dtype=bool)
        available[ev["valid_idx"]] = True
        pos = np.asarray(ev["full_lineage_positive"], dtype=bool)   # value > divisor, only where available
        norm = np.asarray(ev["full_normalized"], dtype=float)       # value / divisor at valid_idx, NaN off it
        cell_state[~available] = CELL_UNAVAILABLE                    # VALID but cell missing a measurement
        cell_state[available & pos] = CELL_POSITIVE
        cell_state[available & ~pos] = CELL_NEGATIVE
    elif state is PairState.NO_TARGET_EXPRESSION_CONFIRMED:
        cell_state[:] = CELL_NEGATIVE                               # confirmed absence -> confident negative gate
    else:
        cell_state[:] = CELL_UNAVAILABLE                            # technical failure -> no call (driver halts)
    return pos, norm, cell_state


def _native(value):
    """Coerce numpy scalars/NaN to JSON-safe native Python (for the per-donor manifest)."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _atomic_write(dst: Path, write_tmp) -> None:
    """Write via a per-process temp file then ``os.replace`` (atomic on the same filesystem)."""
    dst = Path(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(f"{dst.name}.tmp.{os.getpid()}")
    try:
        write_tmp(tmp)
        os.replace(tmp, dst)
    finally:
        if tmp.exists():
            tmp.unlink()


def _apply_one_donor(
    donor: str,
    *,
    sample_root: str,
    raw_cells: str,
    redsea_cells: str,
    gated_dir: str,
    pair_review_map: dict,
    expr_review_map: dict,
    config: PairValidationConfig,
    require_all_cells: bool,
) -> dict:
    """Type one donor's cells for the 7 accepted pairs and write its tri-state gated parquet.

    Top-level + picklable (crosses the process boundary via :func:`map_donors`). Fails loudly on a
    missing/invalid donor, a coverage mismatch, an unaccepted pair, or a technical non-VALID state.
    """
    donor_df, audit = load_validation_sample(
        Path(sample_root), Path(raw_cells), Path(redsea_cells), donor
    )
    if require_all_cells and int(audit["n"]) != int(audit["canonical_n"]):
        raise ValueError(
            f"donor {donor}: all-cells sample covers {audit['n']:,}/{audit['canonical_n']:,} canonical "
            "cells; a production apply must run on the modulus-1 (exhaustive) sample"
        )

    cols: dict[str, object] = {
        "object_id": donor_df["object_id"].astype(str).to_numpy(),
        "donor_id": str(donor),
    }
    metrics_rows: list[dict] = []
    pair_summaries: list[dict] = []
    for target, reference in ACCEPTED_PAIRS:
        review = pair_review_map.get((target, reference), PairReview.UNREVIEWED)
        if review is not PairReview.ACCEPTED:
            raise ValueError(
                f"donor {donor} {target} <- {reference}: pair-review manifest must mark this accepted pair "
                f"ACCEPTED (got {getattr(review, 'value', review)})"
            )
        ev = evaluate_locked_pair(
            donor_df,
            donor,
            target,
            reference,
            config=config,
            expression_review=expr_review_map.get((donor, target), ExpressionReview.UNREVIEWED),
            pair_review=review,
        )
        state = PairState(ev["state"])
        if state in HARD_FAIL_STATES:
            raise ValueError(
                f"donor {donor} {target} <- {reference}: ACCEPTED pair returned {state.value} "
                f"({ev.get('state_reason', '')}); halting — no automatic fallback (a technical failure on "
                "a validated pair is a method/data problem to surface, not a biological negative)"
            )
        if int(ev["n_input"]) != len(donor_df):
            raise AssertionError(
                f"donor {donor} {target}: evaluation length {ev['n_input']} != {len(donor_df)} cells"
            )
        pos, norm, cell_state = cell_calls_from_evaluation(ev)
        cols[f"{target}_pos"] = pos
        cols[f"{target}_norm"] = norm
        cols[f"{target}_state"] = cell_state

        row = dict(locked_pair_metrics(ev))
        row["image"] = audit["image"]
        row["called_positive_n"] = int((cell_state == CELL_POSITIVE).sum())
        row["called_negative_n"] = int((cell_state == CELL_NEGATIVE).sum())
        row["called_unavailable_n"] = int((cell_state == CELL_UNAVAILABLE).sum())
        metrics_rows.append(row)
        pair_summaries.append(
            {
                "target": target,
                "reference": reference,
                "state": state.value,
                "canonical_divisor": _native(row.get("canonical_divisor")),
                "reference_to_target_ratio": _native(row.get("reference_to_target_ratio")),
                "divisor_reproducibility": _native(row.get("divisor_reproducibility")),
                "target_low_sbr": _native(row.get("target_low_sbr")),
                "called_positive_n": row["called_positive_n"],
                "called_negative_n": row["called_negative_n"],
                "called_unavailable_n": row["called_unavailable_n"],
            }
        )

    gated = pd.DataFrame(cols)
    dst = Path(gated_dir) / f"donor_id={donor}"
    _atomic_write(dst / "data_0.parquet", lambda p: gated.to_parquet(p, index=False))
    manifest = {
        "method_version": METHOD_VERSION,
        "donor": str(donor),
        "image": str(audit["image"]),
        "n_cells": int(len(gated)),
        "canonical_n": int(audit["canonical_n"]),
        "accepted_pairs": [f"{t} <- {r}" for t, r in ACCEPTED_PAIRS],
        "pairs": pair_summaries,
    }
    _atomic_write(
        dst / "apply_manifest.json",
        lambda p: p.write_text(json.dumps(manifest, indent=2)),
    )
    n = len(gated)
    called = {t: int((cols[f"{t}_state"] == CELL_POSITIVE).sum()) for t in REQUIRED_TARGETS}
    print(
        f"[{donor}] {n:,} cells | "
        + " ".join(f"{t}={100 * called[t] / max(n, 1):.1f}%" for t in REQUIRED_TARGETS),
        flush=True,
    )
    return {"donor": str(donor), "n": n, "metrics": metrics_rows}


def run_apply(
    cfg: PipelineConfig,
    *,
    donors=None,
    n_jobs=None,
    config: PairValidationConfig | None = None,
    require_all_cells: bool = True,
) -> pd.DataFrame:
    """Apply the 7 accepted RESTORE pairs on all cells per donor; write gated parquets + a cohort manifest.

    Returns the cohort donor x pair disposition/divisor table (also written to
    ``cfg.restore_apply_manifest_csv``).

    When ``config`` is omitted the pair-method policy comes from ``cfg`` (``[restore]``
    ``control_definition`` / ``divisor_statistic``), so a sweep arm accepted through
    ``restore_validation evaluate-locked`` can be reproduced by the pipeline. Passing ``config``
    explicitly overrides the file, for tests and one-off pilots.
    """
    config = config if config is not None else config_from_pipeline(cfg)

    if not cfg.restore_allcells_sample_dir.exists():
        raise SystemExit(
            f"[restore] all-cells (modulus-1) compartment sample not found: {cfg.restore_allcells_sample_dir}\n"
            "          (the modulus-1 QuPath sample from Gate 3 is the production input)"
        )
    if not cfg.restore_pair_reviews_csv.exists():
        raise SystemExit(
            f"[restore] accepted-pair manifest not found: {cfg.restore_pair_reviews_csv}"
        )

    pair_map, _ = load_pair_reviews(cfg.restore_pair_reviews_csv)
    expr_path = cfg.restore_expression_reviews_csv
    expr_map, _ = load_expression_reviews(expr_path if expr_path.exists() else None)
    missing = [
        f"{t} <- {r}"
        for t, r in ACCEPTED_PAIRS
        if pair_map.get((t, r)) is not PairReview.ACCEPTED
    ]
    if missing:
        raise SystemExit(
            f"[restore] pair-review manifest {cfg.restore_pair_reviews_csv} is missing ACCEPTED rows "
            f"for: {missing}"
        )

    donor_ids = (
        list(ensure_eligible_donors(donors, context="RESTORE apply"))
        if donors
        else cfg.discover_donors(cfg.restore_allcells_sample_dir)
    )
    if not donor_ids:
        raise SystemExit(f"[restore] no eligible donors under {cfg.restore_allcells_sample_dir}")

    # Echo the pair-method policy: it is the one thing about this run that is NOT determined by the
    # frozen specification, and it changes every divisor.
    print(
        f"[restore] {METHOD_VERSION} | control_definition={config.control_definition} "
        f"divisor_statistic={config.divisor_statistic} | {len(donor_ids)} donors "
        f"x {len(ACCEPTED_PAIRS)} accepted pairs"
    )

    cfg.restore_gated_dir.mkdir(parents=True, exist_ok=True)
    fn = functools.partial(
        _apply_one_donor,
        sample_root=str(cfg.restore_allcells_sample_dir),
        raw_cells=str(cfg.cells_dir),
        redsea_cells=str(cfg.cells_redsea_dir),
        gated_dir=str(cfg.restore_gated_dir),
        pair_review_map=pair_map,
        expr_review_map=expr_map,
        config=config,
        require_all_cells=require_all_cells,
    )
    n_jobs = cfg.n_jobs if n_jobs is None else n_jobs
    # on_error='raise' — a per-donor failure must stop the whole run, never be swallowed as a None row.
    results = map_donors(fn, donor_ids, n_jobs=n_jobs, ordered=True, on_error="raise")

    rows = [row for res in results if res for row in res["metrics"]]
    manifest = pd.DataFrame(rows)
    _atomic_write(cfg.restore_apply_manifest_csv, lambda p: manifest.to_csv(p, index=False))
    print(
        f"\n[restore] wrote {len([r for r in results if r])} gated donors -> {cfg.restore_gated_dir}\n"
        f"[restore] cohort manifest ({len(manifest)} donor-pairs) -> {cfg.restore_apply_manifest_csv}"
    )
    return manifest


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--config", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=None, help="per-donor process pool size")
    ap.add_argument("--donors", nargs="*", default=None)
    # Same two knobs as `restore_validation evaluate-locked`, so a sweep arm and its production apply
    # are invoked identically. Default None == take the value from config.ini.
    ap.add_argument(
        "--control-definition",
        choices=sorted(CONTROL_DEFINITIONS),
        default=None,
        help="override config.ini [restore] control_definition",
    )
    ap.add_argument(
        "--divisor-statistic",
        choices=sorted(DIVISOR_STATISTICS),
        default=None,
        help="override config.ini [restore] divisor_statistic",
    )
    a = ap.parse_args(argv)
    cfg = load_config(a.config)
    if a.jobs is not None:
        cfg.n_jobs = a.jobs
    if a.control_definition is not None:
        cfg.restore_control_definition = a.control_definition
    if a.divisor_statistic is not None:
        cfg.restore_divisor_statistic = a.divisor_statistic
    run_apply(cfg, donors=a.donors)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
