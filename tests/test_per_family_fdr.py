"""Per-family BH-FDR matches a hand-computed reference for small p-vectors."""
import numpy as np
import pandas as pd

from utils import bh_fdr, bh_fdr_by_family, bonferroni_by_family


def _bh_reference(p):
    """Plain BH on a list — naive monotone-from-right."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = p[order] * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    return out


def test_bh_matches_reference():
    p = [0.001, 0.01, 0.03, 0.05, 0.20, 0.50]
    np.testing.assert_allclose(bh_fdr(p), _bh_reference(p), rtol=1e-9)


def test_bh_by_family_independent_per_family():
    df = pd.DataFrame({
        "family": ["secondary"] * 4 + ["exploratory"] * 6,
        "p":      [0.001, 0.01, 0.03, 0.05, 0.001, 0.005, 0.01, 0.05, 0.20, 0.50],
    })
    q = bh_fdr_by_family(df, family_col="family", p_col="p")
    # Compare to standalone BH on each family
    expect_sec = _bh_reference(df.loc[df["family"] == "secondary", "p"].values)
    expect_exp = _bh_reference(df.loc[df["family"] == "exploratory", "p"].values)
    np.testing.assert_allclose(q[df["family"] == "secondary"].values, expect_sec, rtol=1e-9)
    np.testing.assert_allclose(q[df["family"] == "exploratory"].values, expect_exp, rtol=1e-9)


def test_bonferroni_by_family_uses_per_family_n():
    df = pd.DataFrame({
        "family": ["A"] * 2 + ["B"] * 4,
        "p":      [0.01, 0.04, 0.01, 0.02, 0.03, 0.05],
    })
    bonf = bonferroni_by_family(df, family_col="family", p_col="p")
    # Family A: n=2 → 0.01*2=0.02 ; 0.04*2=0.08
    # Family B: n=4 → 0.01*4=0.04 ; 0.02*4=0.08 ; 0.03*4=0.12→cap 1 ; 0.05*4=0.20
    expected = pd.Series([0.02, 0.08, 0.04, 0.08, 0.12, 0.20], dtype=float)
    np.testing.assert_allclose(bonf.values, expected.values, rtol=1e-9)
