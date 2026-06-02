from __future__ import annotations

import numpy as np
import pandas as pd

from disease_subsets import (
    collect_disease_columns,
    derive_disease_flags,
    fit_interaction_model,
)
from utils import make_design_matrix


def test_ms_first_occurrence_uses_g35_not_f31() -> None:
    basket = pd.DataFrame({
        "f.eid": [1, 2, 3, 4],
        "f.130892.0.0": ["2010-01-01", None, None, None],  # F31 bipolar, not MS
        "f.131042.0.0": [None, "2012-02-02", None, None],  # G35 MS
        "f.131043.0.0": [None, "HES", None, None],
        "f.20002.0.0": [None, None, "1261", None],
        "f.41270.0.0": [None, None, None, "G35"],
    })
    groups = collect_disease_columns(list(basket.columns))
    flags, audit = derive_disease_flags(basket, groups)

    assert flags.loc[flags["eid"].eq(1), "ms_case"].item() == 0
    assert flags.loc[flags["eid"].eq(2), "ms_case"].item() == 1
    assert flags.loc[flags["eid"].eq(3), "ms_case"].item() == 1
    assert flags.loc[flags["eid"].eq(4), "ms_case"].item() == 1
    assert audit["n_ms_cases_union"] == 3


def test_broad_neurodegeneration_excludes_ms_only_cases() -> None:
    basket = pd.DataFrame({
        "f.eid": [1, 2, 3],
        "f.41270.0.0": ["G35", "F00", None],
        "f.42032.0.0": [None, None, "2015-05-05"],
        "f.131042.0.0": ["2014-04-04", None, None],
    })
    flags, audit = derive_disease_flags(basket)

    assert flags.loc[flags["eid"].eq(1), "ms_case"].item() == 1
    assert flags.loc[flags["eid"].eq(1), "broad_neurodeg_case"].item() == 0
    assert flags.loc[flags["eid"].eq(2), "broad_neurodeg_case"].item() == 1
    assert flags.loc[flags["eid"].eq(3), "broad_neurodeg_case"].item() == 1
    assert audit["n_broad_neurodeg_cases"] == 2


def test_interaction_model_recovers_case_specific_effect() -> None:
    rng = np.random.default_rng(7)
    dosage = np.tile(np.array([0.0, 1.0, 2.0]), 80)
    genotype = pd.Series(np.tile(np.array(["GG", "AG", "AA"]), 80))
    disease = pd.Series(np.repeat([0.0, 1.0], 120))
    age = pd.Series(rng.normal(60, 5, size=240))
    covars = make_design_matrix(
        pd.DataFrame({"age": age}),
        cont_cols=["age"],
        cat_cols=[],
        add_intercept=True,
    )
    y = pd.Series(
        0.10 * dosage
        + 0.80 * dosage * disease.to_numpy()
        + 0.01 * age.to_numpy()
        + rng.normal(0, 0.05, size=240)
    )

    res = fit_interaction_model(
        y,
        pd.Series(dosage),
        disease,
        covars,
        genotype=genotype,
        min_cases=100,
        min_controls=100,
        min_case_genotype_n=10,
    )

    assert res["status"] == "ok"
    assert res["n_case"] == 120
    assert res["beta_interaction"] > 0.7
    assert res["p_interaction"] < 1e-20


def test_interaction_model_requires_case_genotype_cell_size() -> None:
    rng = np.random.default_rng(11)
    n = 220
    disease = pd.Series(np.r_[np.zeros(110), np.ones(110)])
    dosage = pd.Series(np.r_[np.tile([0.0, 1.0, 2.0], 36), 0.0, 1.0,
                             np.repeat(0.0, 90), np.repeat(1.0, 9), np.repeat(2.0, 11)])
    genotype = dosage.map({0.0: "GG", 1.0: "AG", 2.0: "AA"})
    covars = make_design_matrix(
        pd.DataFrame({"age": rng.normal(60, 5, size=n)}),
        cont_cols=["age"],
        cat_cols=[],
        add_intercept=True,
    )
    y = pd.Series(rng.normal(size=n))

    res = fit_interaction_model(
        y,
        dosage,
        disease,
        covars,
        genotype=genotype,
        min_cases=100,
        min_controls=100,
        min_case_genotype_n=10,
    )

    assert res["status"] == "case_genotype_cell_lt_10"
    assert res["n_case_AG"] == 9
