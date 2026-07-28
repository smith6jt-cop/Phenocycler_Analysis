"""Xenium import: null handling, label fallback, and the same-slide docstring contract.

Regression guards for the three review findings on PR #5. Two of them fail *silently* — a
null label becomes the string "nan" and a fallback path quietly stops firing — so they need
tests rather than just a fix.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from phenocycler.integration.import_xenium import _obs_str, import_sample

anndata = pytest.importorskip("anndata")


def _h5ad(obs: pd.DataFrame, tmp_path, n_genes: int = 6):
    """Write a minimal phenotyped-style h5ad; `import_sample` reads a path, not an object."""
    rng = np.random.default_rng(0)
    n = len(obs)
    X = rng.poisson(0.5, (n, n_genes)).astype("float32")
    a = anndata.AnnData(X=X, obs=obs,
                        var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))
    a.layers["lognorm"] = np.log1p(X)
    a.obsm["spatial"] = np.column_stack([rng.uniform(0, 2000, n),
                                         rng.uniform(0, 2000, n)]).astype("float32")
    path = Path(tmp_path) / "0059865_phenotyped.h5ad"
    a.write_h5ad(path)
    return path


# --------------------------------------------------------------------------- #
# _obs_str null handling
# --------------------------------------------------------------------------- #

def test_obs_str_preserves_nulls_as_the_default():
    """`astype(str)` before `fillna` turns NaN into the literal string "nan"."""
    obs = pd.DataFrame({"celltype": pd.Series(["Endocrine", None, "Exocrine"], dtype=object)})
    out = _obs_str(obs, "celltype")
    assert list(out) == ["Endocrine", "", "Exocrine"]
    assert "nan" not in set(out)


def test_obs_str_handles_categorical_columns():
    """AnnData obs columns are routinely Categorical, where `fillna` raises outright."""
    obs = pd.DataFrame({"celltype": pd.Categorical(["Beta", None, "Alpha"])})
    out = _obs_str(obs, "celltype")
    assert list(out) == ["Beta", "", "Alpha"]


def test_obs_str_respects_a_custom_default():
    obs = pd.DataFrame({"x": pd.Series([None, "a"], dtype=object)})
    assert list(_obs_str(obs, "x", default="?")) == ["?", "a"]


def test_obs_str_missing_column_returns_default():
    obs = pd.DataFrame({"other": [1, 2]})
    assert list(_obs_str(obs, "celltype", default="-")) == ["-", "-"]


# --------------------------------------------------------------------------- #
# The second-order bug: a silently dead fallback
# --------------------------------------------------------------------------- #

def test_celltype_fine_fallback_fires_when_celltype_is_all_null(tmp_path):
    """`import_sample` falls back to `celltype_fine` when `celltype` is empty.

    With the null bug the column read as all-"nan" rather than all-"", so `.eq("").all()`
    was False and the fallback never ran — an h5ad carrying `celltype_fine` but a null
    `celltype` would have lost its fine labels with no error and no warning.
    """
    # Categorical with unset values, which is how AnnData actually stores a null label
    # (h5py cannot serialise an object column containing None at all).
    n = 20
    obs = pd.DataFrame({
        "celltype_broad": pd.Categorical(["Endocrine"] * n),
        "celltype": pd.Categorical([None] * n, categories=["Beta", "Alpha"]),
        "celltype_fine": pd.Categorical(["Beta"] * 12 + ["Alpha"] * 8),
    }, index=[f"c{i}" for i in range(n)])

    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc")
    assert set(df["celltype_fine"]) == {"Beta", "Alpha"}
    assert "nan" not in set(df["celltype_fine"])


def test_null_labels_do_not_leak_the_string_nan(tmp_path):
    """A null broad label must become '' (and map to Unassigned), never the label "nan"."""
    n = 10
    obs = pd.DataFrame({
        "celltype_broad": pd.Categorical(["Endocrine"] * 6 + [None] * 4),
        "celltype": pd.Categorical(["Beta"] * 6 + [None] * 4),
    }, index=[f"c{i}" for i in range(n)])

    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc")
    assert "nan" not in set(df["lineage_native"])
    assert "nan" not in set(df["celltype_fine"])
    assert (df["lineage_native"] == "").sum() == 4
    assert (df["lineage_common"] == "Unassigned").sum() == 4


def test_import_sample_stamps_the_join_key(tmp_path):
    """donor_id/roi reach the data here and nowhere else in the Xenium world."""
    obs = pd.DataFrame({"celltype_broad": ["Exocrine"] * 5, "celltype": ["Acinar"] * 5},
                       index=[f"c{i}" for i in range(5)])
    df = import_sample(_h5ad(obs, tmp_path), donor_id="6539", roi="panc",
                       sample_key="0059865__panc")
    assert set(df["donor_id"]) == {"6539"}
    assert set(df["roi"]) == {"panc"}
    assert set(df["xenium_sample_key"]) == {"0059865__panc"}


def test_broad_lowsignal_becomes_evidence_negative(tmp_path):
    """Xenium's argmax-fallback cells are the counterpart of PhenoCycler's epi_default."""
    n = 10
    obs = pd.DataFrame({
        "celltype_broad": ["Exocrine"] * n,
        "celltype": ["Acinar"] * n,
        "broad_lowsignal": [True] * 3 + [False] * 7,
    }, index=[f"c{i}" for i in range(n)])
    df = import_sample(_h5ad(obs, tmp_path), donor_id="1", roi="panc")
    assert int((~df["evidence_positive"]).sum()) == 3


# --------------------------------------------------------------------------- #
# same-slide docstring must match the implementation
# --------------------------------------------------------------------------- #

def test_sameslide_docstring_matches_what_pair_donor_does():
    """The docstring claimed IoU was preferred while `pair_donor` only used centroids.

    Asserted rather than eyeballed because the two drift apart silently — the reader has no
    way to tell which criterion actually ran.
    """
    import inspect

    from phenocycler.integration import sameslide

    doc = inspect.getdoc(sameslide) or ""
    src = inspect.getsource(sameslide.pair_donor)

    uses_iou = "match_by_iou(" in src
    uses_centroid = "match_by_centroid(" in src
    assert uses_centroid and not uses_iou, (
        "pair_donor's implementation changed — update the module docstring to match")
    assert "not wired in" in doc or "not yet wired in" in doc, (
        "the docstring must say IoU is available but unused while that remains true")
