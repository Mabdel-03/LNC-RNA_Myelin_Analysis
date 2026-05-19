"""Verify composite sign conventions + missing-component handling."""
import numpy as np
import pandas as pd

from synthetic_data import make_synthetic_cohort
from utils import compute_composite_score


def test_demyelination_like_sign_under_injected_signal():
    """We inject FA↓, MD↑, ICVF↓, ISOVF↑ at the dosage in scp_left.

    `demyelination_like` is defined as +RD +MD -FA -ICVF +ISOVF. Each term
    contributes positively under our injection (FA↓ → -(-) > 0; MD↑ > 0;
    ICVF↓ → -(-) > 0; ISOVF↑ > 0). The score should be POSITIVELY correlated
    with dosage.
    """
    wide, manifest = make_synthetic_cohort(n=2000, seed=11,
                                            inject_signal_in_tract="scp_left")
    # Build the per-tract component DataFrame for the injected tract
    tract_cols = {row["metric_guess"]: row["column_name"]
                  for _, row in manifest.iterrows()
                  if row["tract_or_region_guess"] == "scp_left"}
    comp_df = pd.DataFrame({m: wide[c] for m, c in tract_cols.items()})

    weights = {"RD": +1, "MD": +1, "FA": -1, "ICVF": -1, "ISOVF": +1}
    score = compute_composite_score(comp_df, weights, min_components=2,
                                     standardize="zscore", rescale_by_sqrt_n=True)
    r = pd.Series(score).corr(pd.Series(wide["dosage_A"]))
    assert np.isfinite(r), "composite-dosage correlation must be finite"
    assert r > 0.05, f"demyelination_like should correlate POSITIVELY with dosage; got r={r:.3f}"


def test_min_components_rule():
    """If only 1 component is present per eid, score must be NaN."""
    # Construct a tiny synthetic frame with only FA available (others NaN for everyone)
    n = 100
    comp = pd.DataFrame({
        "FA": np.random.default_rng(0).normal(size=n),
        "MD": np.nan,
        "RD": np.nan,
        "ICVF": np.nan,
        "ISOVF": np.nan,
    })
    weights = {"RD": +1, "MD": +1, "FA": -1, "ICVF": -1, "ISOVF": +1}
    score = compute_composite_score(comp, weights, min_components=2,
                                     standardize="zscore")
    # With only FA present (only 1 weighted component), all scores must be NaN
    assert pd.Series(score).notna().sum() == 0


def test_zero_components_returns_all_nan():
    """If none of the weight metrics exist in the frame, score is all NaN."""
    n = 50
    comp = pd.DataFrame({"OTHER_METRIC": np.random.normal(size=n)})
    weights = {"FA": -1, "MD": +1}
    score = compute_composite_score(comp, weights, min_components=2)
    assert pd.Series(score).notna().sum() == 0


def test_axonal_loss_sign_under_fa_injection():
    """axonal_loss_like = -L1 -FA -ICVF. With FA↓ and ICVF↓ (and L1 random),
    the -FA and -ICVF terms become positive → score correlates POSITIVELY
    with dosage."""
    wide, manifest = make_synthetic_cohort(n=2000, seed=12,
                                            inject_signal_in_tract="scp_left")
    tract_cols = {row["metric_guess"]: row["column_name"]
                  for _, row in manifest.iterrows()
                  if row["tract_or_region_guess"] == "scp_left"}
    comp_df = pd.DataFrame({m: wide[c] for m, c in tract_cols.items()})
    weights = {"L1": -1, "FA": -1, "ICVF": -1}
    score = compute_composite_score(comp_df, weights, min_components=2,
                                     standardize="zscore")
    r = pd.Series(score).corr(pd.Series(wide["dosage_A"]))
    assert r > 0.03, f"axonal_loss_like should correlate POSITIVELY; got r={r:.3f}"
