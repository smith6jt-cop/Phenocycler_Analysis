"""
Optional GPU backend (CuPy) with a transparent CPU fallback.

REDSEA's two numeric hot spots — the per-channel ``np.bincount`` whole-cell /
boundary sums and the sparse ``F @ edge`` spillover matmul — map cleanly onto
CuPy (``cupy.bincount`` + ``cupyx.scipy.sparse``).  This module detects CuPy at
import time and exposes a tiny ``xp``-style shim so ``redsea.py`` can run the
exact same code path on GPU or CPU.  Modelled on KINTSUGI's ``RAPIDSManager``
(size-gated GPU use + automatic fallback).

Nothing here is required: if CuPy is absent or ``use_gpu=False`` the pure
NumPy/SciPy path is used and results are identical.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

_CUPY = None
_CUPYX_SPARSE = None
try:  # pragma: no cover - environment dependent
    import cupy as _cp
    import cupyx.scipy.sparse as _cpx_sparse
    _CUPY = _cp
    _CUPYX_SPARSE = _cpx_sparse
except Exception:  # ImportError, or CUDA init failure
    _CUPY = None
    _CUPYX_SPARSE = None


def gpu_available() -> bool:
    """True if CuPy imported and at least one CUDA device is visible."""
    if _CUPY is None:
        return False
    try:  # pragma: no cover - environment dependent
        return _CUPY.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


class Backend:
    """Minimal array backend: ``xp`` (numpy or cupy) + helpers used by REDSEA.

    Use :func:`get_backend` to construct.  When ``on_gpu`` is False every method
    is a thin wrapper over NumPy/SciPy, so callers need no branching.
    """

    def __init__(self, on_gpu: bool, device: int = 0):
        self.on_gpu = bool(on_gpu and gpu_available())
        self.device = device
        if self.on_gpu:  # pragma: no cover - environment dependent
            self.xp = _CUPY
            self._spmod = _CUPYX_SPARSE
            _CUPY.cuda.Device(device).use()
        else:
            self.xp = np
            self._spmod = sp

    # -- transfer helpers ----------------------------------------------------
    def asarray(self, a):
        return self.xp.asarray(a)

    def asnumpy(self, a):
        if self.on_gpu:  # pragma: no cover - environment dependent
            return _CUPY.asnumpy(a)
        return np.asarray(a)

    def bincount(self, x, weights=None, minlength=0):
        return self.xp.bincount(x, weights=weights, minlength=minlength)

    def csr_from_scipy(self, m: sp.spmatrix):
        """Move a SciPy CSR to the active backend (identity on CPU)."""
        if self.on_gpu:  # pragma: no cover - environment dependent
            return self._spmod.csr_matrix(m)
        return m.tocsr()

    def diags(self, v):
        return self._spmod.diags(v)

    def clip_(self, a, lo, hi=None):
        self.xp.clip(a, lo, hi, out=a)
        return a


def get_backend(use_gpu: bool, device: int = 0) -> Backend:
    """Return a :class:`Backend`; silently falls back to CPU when GPU is absent."""
    return Backend(on_gpu=use_gpu, device=device)
