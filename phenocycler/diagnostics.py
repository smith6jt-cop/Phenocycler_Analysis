#!/usr/bin/env python3
"""
Interactive RESTORE threshold-tracing diagnostics — the shared COMPUTE behind the step-by-step figures in
``notebooks/03_restore_normalize.ipynb`` (and the composition trace in ``04_broad_lineage.ipynb``).

Traces the per-image RESTORE threshold ONE transformation at a time, validated against each donor's own
data (its within-islet Otsu valley) rather than a single summary metric — a summary metric (whole-tissue
Otsu) once falsely flagged INS, a clean marker, as the worst pair. The per-step PLOTS live inline in the
notebook (short, pedagogical); this module provides the compute so the diagnostic and the pipeline never
drift:

    read_donor    a donor's REDSEA-corrected target/reference MFI + cell_region (islet membership)
    r_otsu        pure Otsu valley (skimage) — the honest within-islet background<->signal boundary
    idx_select    mutually-exclusive single-positive selection             (re-exported from restore)
    neg_stat      threshold statistic on the negative population            (re-exported from restore)
    fit_clusters  RESTORE 2-cluster negative-population fit (SSC/GMM/KMeans) (re-exported from restore)

Also re-exports the production constants/mappings the tracing reuses (so the notebook does not hardcode
them): ``MARKER_PAIRS``, ``COMPARTMENT_ORDER``, ``OTHER_LABEL``, ``STATUS_ORDER`` and ``status_map(cfg)``.

NB: importing this module imports ``phenocycler.restore``, which sets the matplotlib ``Agg`` backend — in a
notebook, re-assert ``%matplotlib inline`` in the cell AFTER importing, so the tracing figures render inline.
"""

from __future__ import annotations

import glob

import numpy as np
import pandas as pd  # noqa: F401  (re-exported convenience for notebook cells)
import pyarrow.parquet as pq

from .config import MARKER_PAIRS, COMPARTMENT_ORDER, OTHER_LABEL, STATUS_ORDER
from .lineage import status_map
from .restore import (idx_select, neg_stat, fit_clusters, r_otsu, partner_low_thresh, bimodal_thresh,
                      logmean_ksigma)

__all__ = [
    "read_donor", "r_otsu",
    "idx_select", "neg_stat", "fit_clusters", "partner_low_thresh", "bimodal_thresh", "logmean_ksigma",
    "MARKER_PAIRS", "COMPARTMENT_ORDER", "OTHER_LABEL", "STATUS_ORDER", "status_map",
]


def read_donor(cfg, donor, markers):
    """REDSEA-corrected raw MFI for ``markers`` (``cfg.cells_redsea_dir``) joined to ``cell_region``
    (``cfg.cells_dir``) on ``object_id``, so the diagnostics can look at the islet core where a marker is
    genuinely bimodal. Returns a DataFrame of ``markers`` + ``cell_region`` (rows with any NaN marker dropped).
    Reads a single parquet file per donor (dodges the hive ``donor_id`` string/int32 inference bug)."""
    donor = str(donor)
    markers = list(markers)
    f = glob.glob(str(cfg.cells_redsea_dir / f"donor_id={donor}" / "*.parquet"))[0]
    df = pq.ParquetFile(f).read(columns=["object_id"] + markers).to_pandas()
    df["object_id"] = df.object_id.astype(str)
    fc = glob.glob(str(cfg.cells_dir / f"donor_id={donor}" / "*.parquet"))[0]
    reg = pq.ParquetFile(fc).read(columns=["object_id", "cell_region"]).to_pandas()
    reg["object_id"] = reg.object_id.astype(str)
    return df.merge(reg, on="object_id", how="left").dropna(subset=markers)


# r_otsu + partner_low_thresh now live in phenocycler.restore (single source; re-exported above) so the
# pipeline and this diagnostics module share one implementation. See restore.py.
