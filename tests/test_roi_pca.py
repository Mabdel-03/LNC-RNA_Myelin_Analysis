"""ROI PCA orientation + reproducibility + PC2 dropping."""
import numpy as np
import pandas as pd

from synthetic_data import make_synthetic_cohort
from utils import compute_roi_pca


def test_pc1_reproducible_under_same_seed():
    wide, manifest = make_synthetic_cohort(n=800, seed=33,
                                            inject_signal_in_tract="scp_left")
    tract_cols = [r["column_name"] for _, r in manifest.iterrows()
                   if r["tract_or_region_guess"] == "scp_left"]
    metric_df = wide[tract_cols].rename(
        columns={c: c.split("_")[-3] for c in tract_cols})  # quick metric label

    res_a = compute_roi_pca(metric_df, random_state=123)
    res_b = compute_roi_pca(metric_df, random_state=123)
    # Absolute loadings should match (sign is arbitrary in PCA)
    la = res_a["loadings"].abs().sort_index()
    lb = res_b["loadings"].abs().sort_index()
    np.testing.assert_allclose(la.values, lb.values, rtol=1e-6)


def test_pc2_dropped_when_below_evr_threshold():
    """Simulate a tract where one direction explains nearly all variance →
    PC2's EVR is below threshold and only PC1 is kept."""
    rng = np.random.default_rng(0)
    n = 500
    x = rng.normal(size=n)
    # Make 4 metrics that are essentially noisy copies of x → only 1 strong direction
    df = pd.DataFrame({
        "m1": x + 0.01 * rng.normal(size=n),
        "m2": x + 0.01 * rng.normal(size=n),
        "m3": x + 0.01 * rng.normal(size=n),
        "m4": x + 0.01 * rng.normal(size=n),
    })
    res = compute_roi_pca(df, pc2_min_evr=0.10, random_state=7)
    assert list(res["scores"].columns) == ["PC1"], \
        f"expected only PC1 kept; got {list(res['scores'].columns)}"
    assert res["evr"][0] > 0.95


def test_pc2_kept_when_meaningful():
    """When two independent signal directions exist, PC2 should be retained."""
    rng = np.random.default_rng(1)
    n = 600
    x = rng.normal(size=n)
    y = rng.normal(size=n)  # independent direction
    df = pd.DataFrame({
        "m1": x + 0.05 * rng.normal(size=n),
        "m2": x + 0.05 * rng.normal(size=n),
        "m3": y + 0.05 * rng.normal(size=n),
        "m4": y + 0.05 * rng.normal(size=n),
    })
    res = compute_roi_pca(df, pc2_min_evr=0.10, random_state=7)
    assert list(res["scores"].columns) == ["PC1", "PC2"]
    assert res["evr"][1] >= 0.10
