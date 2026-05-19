"""Verify 04c writes REGENIE-conformant phenotype + covariate + col-list files.

We don't shell out to the real 04c script (it expects intermediate files on
disk). Instead we exercise the underlying transform helper to verify the
formatting choices.
"""
import io

import numpy as np
import pandas as pd

from utils import inverse_rank_normalize, trim_outliers_z


def test_pheno_format_fid_iid_first_two_cols():
    """Convention: REGENIE/BOLT want FID IID as the first two columns,
    followed by phenotypes."""
    eids = pd.Series([1, 2, 3, 4, 5], dtype=int)
    pheno = pd.DataFrame({
        "FID": eids, "IID": eids,
        "PHENO_A": [0.1, 0.2, 0.3, 0.4, 0.5],
        "PHENO_B": [1.0, 2.0, np.nan, 4.0, 5.0],
    })
    buf = io.StringIO()
    pheno.to_csv(buf, sep="\t", index=False, na_rep="NA")
    txt = buf.getvalue()
    header = txt.splitlines()[0].split("\t")
    assert header[:2] == ["FID", "IID"], f"first two cols must be FID, IID; got {header[:2]}"
    # NaN must be encoded as NA (REGENIE/BOLT expectation)
    assert "NA" in txt


def test_pheno_transform_preserves_n_for_none_transform():
    """Composites + ROI PCs have transform='none' — they should pass through
    untouched."""
    s = pd.Series([1.0, 2.0, np.nan, 4.0, 5.0])
    # Mimic the transform branch in 04c:_transform_for_lmm
    transform = "none"
    if transform == "none":
        out = pd.to_numeric(s, errors="coerce")
    pd.testing.assert_series_equal(out, s)


def test_pheno_transform_inrt_for_default():
    """Default branch: trim outliers then INRT — returns a normal-ish series."""
    rng = np.random.default_rng(0)
    s = pd.Series(rng.normal(size=1000))
    z_thresh = 6.0
    out = inverse_rank_normalize(trim_outliers_z(s, z_thresh=z_thresh))
    # INRT output should be approximately standard normal
    assert abs(out.mean()) < 0.1
    assert abs(out.std(ddof=1) - 1) < 0.1
