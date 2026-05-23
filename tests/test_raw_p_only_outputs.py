"""Raw-p-only result annotations."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OLS_PATH = ROOT / "src" / "05a_run_ols.py"
spec = importlib.util.spec_from_file_location("run_ols_05a", OLS_PATH)
run_ols = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(run_ols)


def test_raw_p_annotations_do_not_create_adjusted_p_columns():
    df = pd.DataFrame({
        "family": ["primary_skeleton_tensor", "primary_skeleton_tensor", "controls"],
        "p": [0.02, 0.10, np.nan],
    })

    out = run_ols._add_raw_p_annotations(df, thresholds=[0.05, 0.01])

    assert "raw_p" in out.columns
    assert "minus_log10_raw_p" in out.columns
    assert "raw_p_rank_within_family" in out.columns
    assert "raw_p_lt_0_05" in out.columns
    forbidden = ("fdr", "bonferroni", "q_value", "qval", "adjusted")
    assert not any(any(tok in c.lower() for tok in forbidden) for c in out.columns)


def test_raw_p_summary_counts_nominal_thresholds_only():
    df = pd.DataFrame({
        "family": ["A", "A", "B", "B"],
        "raw_p": [0.001, 0.20, 0.04, np.nan],
    })

    out = run_ols._raw_p_summary(df, thresholds=[0.05, 0.01])
    by_family = out.set_index("family")

    assert by_family.loc["A", "n_tests"] == 2
    assert by_family.loc["A", "n_raw_p_lt_0.05"] == 1
    assert by_family.loc["A", "n_raw_p_lt_0.01"] == 1
    assert by_family.loc["B", "n_tests"] == 1
    assert by_family.loc["B", "n_raw_p_lt_0.05"] == 1
