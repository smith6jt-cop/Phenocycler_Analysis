"""Unit tests for the per-donor parallel helpers."""

from __future__ import annotations

from phenocycler.parallel import map_donors, resolve_jobs


def _square(x):
    return int(x) ** 2


def _boom(x):
    if x == "3":
        raise ValueError("boom")
    return int(x)


def test_resolve_jobs_clamps():
    assert resolve_jobs(1, 100) == 1
    assert resolve_jobs(4, 2) == 2          # never more than items
    assert resolve_jobs(0, 8) >= 1          # 0 => all cpus, at least 1


def test_map_donors_serial_preserves_order():
    out = map_donors(_square, ["1", "2", "3", "4"], n_jobs=1, ordered=True)
    assert out == [1, 4, 9, 16]


def test_map_donors_parallel_ordered():
    out = map_donors(_square, ["1", "2", "3", "4"], n_jobs=2, ordered=True)
    assert out == [1, 4, 9, 16]


def test_map_donors_error_logged_as_none():
    out = map_donors(_boom, ["1", "2", "3", "4"], n_jobs=1, ordered=True, on_error="log")
    assert out == [1, 2, None, 4]


def test_map_donors_error_raises():
    import pytest
    with pytest.raises(ValueError):
        map_donors(_boom, ["3"], n_jobs=1, on_error="raise")
