"""Genotype-level model variants."""
from __future__ import annotations

import numpy as np
import pandas as pd

from utils import fit_genotype_pairwise_model


def test_pairwise_genotype_model_reports_global_and_raw_pairwise_p_values():
    rng = np.random.default_rng(123)
    n_per_group = 45
    genotype = pd.Series(["GG"] * n_per_group + ["AG"] * n_per_group + ["AA"] * n_per_group)
    age = pd.Series(np.tile(np.linspace(-1, 1, n_per_group), 3))
    covariates = pd.DataFrame({"const": 1.0, "age": age})
    group_effect = genotype.map({"GG": 0.0, "AG": 1.0, "AA": 2.0}).astype(float)
    y = group_effect + 0.2 * age + rng.normal(0, 0.15, len(genotype))

    res = fit_genotype_pairwise_model(y, genotype, covariates, levels=["AA", "AG", "GG"])

    assert res["status"] == "ok"
    assert res["global"]["status"] == "ok"
    assert res["global"]["raw_p"] < 1e-20
    assert len(res["pairwise"]) == 3
    assert all("raw_p" in row for row in res["pairwise"])
    assert not any("fdr" in key or "bonferroni" in key
                   for row in res["pairwise"] for key in row)

    contrasts = {row["contrast"]: row for row in res["pairwise"]}
    assert {"AA_vs_AG", "AA_vs_GG", "AG_vs_GG"} == set(contrasts)
    assert contrasts["AA_vs_GG"]["raw_p"] < 1e-20
    assert contrasts["AA_vs_GG"]["beta"] > 1.8
