"""OLS recovers a known additive effect on synthetic data."""
import numpy as np
import pandas as pd

from utils import (
    fit_ols_model,
    inverse_rank_normalize,
    make_design_matrix,
)


def test_ols_recovers_known_beta():
    rng = np.random.default_rng(1234)
    n = 5000
    dose = rng.choice([0, 1, 2], size=n, p=[0.25, 0.5, 0.25]).astype(float)
    age = rng.normal(55, 7, size=n)
    sex = rng.integers(0, 2, size=n)
    pcs = rng.normal(0, 0.01, size=(n, 5))
    noise = rng.normal(size=n)
    beta_true = 0.10
    y_raw = beta_true * dose + 0.5 * (age - 55) / 7 + 0.3 * sex + noise
    df = pd.DataFrame({"age": age, "sex": sex,
                        "PC1": pcs[:, 0], "PC2": pcs[:, 1], "PC3": pcs[:, 2],
                        "PC4": pcs[:, 3], "PC5": pcs[:, 4]})
    X = make_design_matrix(df, cont_cols=["age", "PC1", "PC2", "PC3", "PC4", "PC5"],
                           cat_cols=["sex"], add_intercept=True, drop_first=True)
    y_inrt = inverse_rank_normalize(y_raw)
    res = fit_ols_model(y_inrt, pd.Series(dose), X, exposure_name="dosage_A")
    assert res["status"] == "ok"
    # INRT scales y so β is on a roughly comparable but not identical scale.
    # Just verify direction + significance.
    assert res["beta"] > 0, f"expected positive beta; got {res['beta']:.4f}"
    assert res["p"] < 1e-4, f"expected significant p; got {res['p']:.2g}"
