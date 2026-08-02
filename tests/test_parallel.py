"""Unit tests for the per-donor parallel helpers."""

from __future__ import annotations

import pytest

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
    with pytest.raises(ValueError):
        map_donors(_boom, ["3"], n_jobs=1, on_error="raise")


def test_parallel_interrupt_terminates_submitted_workers(monkeypatch):
    from phenocycler import parallel

    class Future:
        cancelled = False

        def result(self):
            raise KeyboardInterrupt

        def cancel(self):
            self.cancelled = True

    class Process:
        alive = True
        terminated = False
        killed = False

        def is_alive(self):
            return self.alive

        def terminate(self):
            self.terminated = True

        def join(self, timeout):
            pass

        def kill(self):
            self.killed = True
            self.alive = False

    class Executor:
        def __init__(self):
            self.process = Process()
            self._processes = {1: self.process}
            self.futures = []
            self.shutdown_call = None

        def submit(self, function, donor):
            future = Future()
            self.futures.append(future)
            return future

        def shutdown(self, *, wait, cancel_futures):
            self.shutdown_call = (wait, cancel_futures)

    executor = Executor()
    monkeypatch.setattr(
        parallel,
        "ProcessPoolExecutor",
        lambda max_workers: executor,
    )
    monkeypatch.setattr(
        parallel,
        "as_completed",
        lambda futures: tuple(futures),
    )

    with pytest.raises(KeyboardInterrupt):
        map_donors(_square, ["1", "2"], n_jobs=2)
    assert all(future.cancelled for future in executor.futures)
    assert executor.process.terminated
    assert executor.process.killed
    assert executor.shutdown_call == (False, True)
