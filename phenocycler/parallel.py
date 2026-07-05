"""
Per-donor parallel execution helpers.

Every heavy stage of the Senior pipeline (REDSEA, RESTORE threshold+apply, broad
lineage) processes each donor image completely independently — separate GeoJSON,
qptiff, and output parquet — so the whole cohort is embarrassingly parallel.
Upstream ran them strictly serially; this module fans the per-donor work out
across processes with a bounded pool, falling back to a plain serial loop when
``n_jobs <= 1`` (which keeps logs readable and is what the unit tests use).

The mapped callable and its bound arguments must be picklable (they cross the
process boundary), so pass top-level functions and dataclass/`Path`/str args.
"""

from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, Iterable, Optional


def resolve_jobs(n_jobs: int, n_items: int) -> int:
    """Clamp the requested worker count to something sane.

    ``n_jobs <= 0`` means "use all CPUs".  The result is never larger than the
    number of items or the CPU count.
    """
    cpu = os.cpu_count() or 1
    if n_jobs is None or n_jobs == 1:
        return 1
    if n_jobs <= 0:
        n_jobs = cpu
    return max(1, min(n_jobs, cpu, max(1, n_items)))


def map_donors(
    func: Callable,
    donors: Iterable[str],
    *,
    n_jobs: int = 1,
    ordered: bool = False,
    on_error: str = "log",
) -> list:
    """Apply ``func(donor)`` for each donor, optionally across processes.

    Parameters
    ----------
    func
        Top-level (picklable) callable taking a single donor id.  Bind extra
        arguments with :func:`functools.partial` before calling.
    donors
        Iterable of donor ids.
    n_jobs
        Pool size; ``1`` runs a serial loop, ``<=0`` uses all CPUs.
    ordered
        If True, results are returned in the input donor order; otherwise in
        completion order.
    on_error
        ``"log"`` (default) prints and records ``None`` for a failed donor and
        continues; ``"raise"`` re-raises the first exception.

    Returns
    -------
    list
        Per-donor return values (``None`` for donors that failed under
        ``on_error="log"``).
    """
    donors = list(donors)
    jobs = resolve_jobs(n_jobs, len(donors))

    if jobs == 1:
        results = []
        for d in donors:
            try:
                results.append(func(d))
            except Exception as exc:  # noqa: BLE001
                if on_error == "raise":
                    raise
                print(f"[parallel] donor {d} FAILED: {exc}", flush=True)
                results.append(None)
        return results

    out: dict[str, object] = {}
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        fut_to_donor = {ex.submit(func, d): d for d in donors}
        for fut in as_completed(fut_to_donor):
            d = fut_to_donor[fut]
            try:
                out[d] = fut.result()
            except Exception as exc:  # noqa: BLE001
                if on_error == "raise":
                    raise
                print(f"[parallel] donor {d} FAILED: {exc}", flush=True)
                out[d] = None

    if ordered:
        return [out.get(d) for d in donors]
    return list(out.values())
