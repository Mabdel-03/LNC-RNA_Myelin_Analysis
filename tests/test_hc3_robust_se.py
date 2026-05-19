"""HC3 robust SE differs from nonrobust under heteroskedasticity."""
import numpy as np
import pandas as pd

from utils import fit_ols_model, make_design_matrix


def test_hc3_differs_under_heteroskedasticity():
    """Heteroskedasticity must depend on the EXPOSURE for HC3 to perturb the
    exposure's SE — variance correlated only with unrelated covariates barely
    moves the dosage SE."""
    rng = np.random.default_rng(2024)
    n = 3000
    dose = rng.choice([0, 1, 2], size=n, p=[0.25, 0.5, 0.25]).astype(float)
    age = rng.normal(55, 7, size=n)
    # Dose-dependent heteroskedasticity: noise variance scales strongly with dose
    sigma = 0.3 + 1.5 * dose
    y = 0.10 * dose + 0.5 * (age - 55) / 7 + rng.normal(0, sigma, size=n)
    df = pd.DataFrame({"age": age})
    X = make_design_matrix(df, cont_cols=["age"], cat_cols=[],
                           add_intercept=True, drop_first=True)
    base = fit_ols_model(pd.Series(y), pd.Series(dose), X,
                          exposure_name="dosage_A", vcov_type="nonrobust")
    robust = fit_ols_model(pd.Series(y), pd.Series(dose), X,
                            exposure_name="dosage_A", vcov_type="HC3")
    assert base["status"] == "ok" and robust["status"] == "ok"
    # Point estimates must match exactly
    np.testing.assert_allclose(base["beta"], robust["beta"], rtol=1e-10)
    # SEs should differ noticeably under dose-dependent heteroskedasticity
    rel_diff = abs(robust["se"] - base["se"]) / base["se"]
    assert rel_diff > 0.05, f"HC3 SE should differ by >5% from nonrobust; got {rel_diff:.3f}"
    # p-values still finite
    assert np.isfinite(robust["p"])


def test_hc3_matches_nonrobust_under_homoskedasticity():
    rng = np.random.default_rng(2025)
    n = 3000
    dose = rng.choice([0, 1, 2], size=n, p=[0.25, 0.5, 0.25]).astype(float)
    age = rng.normal(55, 7, size=n)
    y = 0.10 * dose + 0.5 * (age - 55) / 7 + rng.normal(0, 1.0, size=n)
    df = pd.DataFrame({"age": age})
    X = make_design_matrix(df, cont_cols=["age"], cat_cols=[],
                           add_intercept=True, drop_first=True)
    base = fit_ols_model(pd.Series(y), pd.Series(dose), X,
                          exposure_name="dosage_A", vcov_type="nonrobust")
    robust = fit_ols_model(pd.Series(y), pd.Series(dose), X,
                            exposure_name="dosage_A", vcov_type="HC3")
    # Under homoskedastic noise the two SEs should be close (within ~10%)
    rel_diff = abs(robust["se"] - base["se"]) / base["se"]
    assert rel_diff < 0.15, f"HC3 and nonrobust SE should agree under homoskedasticity; got {rel_diff:.3f}"
