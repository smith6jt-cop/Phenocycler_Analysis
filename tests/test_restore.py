"""Unit tests for the RESTORE robust-threshold guard and threshold flattening."""

from __future__ import annotations

import numpy as np
import pytest
import pandas as pd

from phenocycler import restore


def _thr_df(model="SSC", marker="CD20", values=None):
    values = values or {"img1": 100.0, "img2": 110.0, "img3": 90.0,
                        "img4": 105.0, "img5": 1000.0}
    return pd.DataFrame([{"image": img, "marker": marker, "model": model, "threshold": v}
                         for img, v in values.items()])


def test_build_chosen_lut_imputes_high_outlier():
    """A 10x-cohort-median threshold is replaced by the cohort median; others kept."""
    thr = _thr_df()
    lut, imputed = restore.build_chosen_lut(thr, "SSC", robust=True, factor=3.0)
    # cohort median over non-outliers {100,110,90,105} = 102.5
    assert ("img5", "CD20") in imputed
    assert np.isclose(lut[("img5", "CD20")], 102.5)
    # inliers untouched
    assert np.isclose(lut[("img1", "CD20")], 100.0)
    assert ("img1", "CD20") not in imputed


def test_build_chosen_lut_imputes_low_and_nonfinite():
    thr = _thr_df(values={"img1": 100.0, "img2": 100.0, "img3": 100.0,
                          "img4": 5.0, "img5": np.nan})   # low outlier + nan
    lut, imputed = restore.build_chosen_lut(thr, "SSC", robust=True, factor=3.0)
    assert ("img4", "CD20") in imputed and ("img5", "CD20") in imputed
    assert np.isclose(lut[("img4", "CD20")], 100.0)
    assert np.isclose(lut[("img5", "CD20")], 100.0)


def test_build_chosen_lut_robust_off_is_passthrough():
    thr = _thr_df()
    lut, imputed = restore.build_chosen_lut(thr, "SSC", robust=False)
    assert imputed == {}
    assert np.isclose(lut[("img5", "CD20")], 1000.0)   # outlier NOT imputed


def test_build_chosen_lut_selects_model():
    """Only the chosen model's rows populate the LUT."""
    thr = pd.concat([_thr_df("SSC", "CD20", {"img1": 100.0}),
                     _thr_df("GMM", "CD20", {"img1": 5000.0})], ignore_index=True)
    lut, _ = restore.build_chosen_lut(thr, "SSC", robust=True)
    assert np.isclose(lut[("img1", "CD20")], 100.0)     # SSC, not GMM


def test_flatten_threshs_drops_global_batch():
    threshs = {
        "global": {"6539": {"CD20": {"SSC": 999.0}}},
        0: {"6539": {"CD20": {"SSC": 120.0, "GMM": 800.0}}},
    }
    df = restore.flatten_threshs(threshs)
    assert set(df.model) == {"SSC", "GMM"}
    assert 999.0 not in set(df.threshold)               # 'global' dropped
    row = df[(df.image == "6539") & (df.model == "SSC")]
    assert np.isclose(row.threshold.iloc[0], 120.0)


def test_paired_restore_cli_directs_to_the_replacement():
    """The paired driver is gone; the CLI must say so instead of failing obscurely."""
    with pytest.raises(SystemExit):
        restore.main([])          # no --proliferation -> argparse error
    assert not hasattr(restore, "run_restore")


def test_headless_stubs_importable():
    """The holoviews/selenium/tqdm.notebook stubs install without a display."""
    restore._install_headless_stubs()
    import holoviews as hv
    import selenium.webdriver  # noqa: F401
    from tqdm.notebook import tqdm  # noqa: F401
    # stubbed hv objects support .opts() and the `*` overlay operator
    obj = hv.Curve() * hv.Scatter()
    assert obj.opts().cols() is not None
