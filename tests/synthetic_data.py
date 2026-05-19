"""Synthetic UKB-shaped data for unit tests.

We simulate n eids with:
  - dosages at one SNP (HWE-consistent, MAF=0.5)
  - age, sex, head_size, motion_dmri, PC1..PC10 covariates
  - ~30 IDPs spanning 3 tracts × {FA, MD, L1, L2, L3, ICVF, OD, ISOVF, MO}
  - one tract has a composite-shaped signal at the dosage so composite tests
    can verify sign conventions
"""
from __future__ import annotations

import numpy as np
import pandas as pd

METRICS = ["FA", "MD", "L1", "L2", "L3", "ICVF", "OD", "ISOVF", "MO"]
TRACTS = ["genu_corpus_callosum", "internal_capsule_left", "scp_left"]


def make_synthetic_cohort(
    n: int = 1000,
    maf: float = 0.5,
    seed: int = 42,
    inject_signal_in_tract: str = "scp_left",
    signal_beta_fa: float = -0.30,   # FA goes DOWN with dose (myelin loss-like)
    signal_beta_md: float = +0.25,
    signal_beta_icvf: float = -0.20,
    signal_beta_isovf: float = +0.20,
    missing_rate: float = 0.05,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (wide_df, manifest_df) with eid index.

    wide_df has columns: eid, dosage_A, age, age2, age_sex, sex, head_size,
    motion_dmri, PC1..PC10, plus per-tract per-metric IDP columns named
    'dMRI_TBSS_<metric>_<tract>'.

    manifest_df mirrors the real pipeline manifest with column_name, panel,
    family, metric_guess, modality_guess, tract_or_region_guess, source,
    transform, include.
    """
    rng = np.random.default_rng(seed)

    eid = np.arange(1, n + 1, dtype=int)
    # HWE dosages: P(GG)=(1-q)^2, P(AG)=2q(1-q), P(AA)=q^2 with q=maf
    q = float(maf)
    dose_probs = np.array([(1 - q) ** 2, 2 * q * (1 - q), q ** 2])
    dosage = rng.choice([0, 1, 2], size=n, p=dose_probs).astype(float)
    genotype = pd.Series(["GG", "AG", "AA"])[dosage.astype(int)].values

    age = rng.normal(55, 7, size=n)
    sex = rng.integers(0, 2, size=n)
    head_size = rng.normal(1.0, 0.05, size=n)
    motion = np.abs(rng.normal(0, 0.5, size=n))
    pcs = {f"PC{i + 1}": rng.normal(0, 0.01, size=n) for i in range(10)}

    wide = pd.DataFrame({
        "eid": eid,
        "dosage_A": dosage,
        "genotype": genotype,
        "age": age,
        "age2": age ** 2,
        "sex": sex,
        "age_sex": age * sex,
        "head_size": head_size,
        "motion_dmri": motion,
        **pcs,
    })

    manifest_rows = []
    for tract in TRACTS:
        for metric in METRICS:
            col = f"dMRI_TBSS_{metric}_{tract}"
            # Base random metric
            base = rng.normal(0, 1, size=n)
            beta = 0.0
            if tract == inject_signal_in_tract:
                if metric == "FA":
                    beta = signal_beta_fa
                elif metric == "MD":
                    beta = signal_beta_md
                elif metric == "ICVF":
                    beta = signal_beta_icvf
                elif metric == "ISOVF":
                    beta = signal_beta_isovf
            y = base + beta * dosage
            # Inject MCAR missingness
            mask = rng.random(size=n) < missing_rate
            y[mask] = np.nan
            wide[col] = y
            manifest_rows.append({
                "field_id": 25000 + len(manifest_rows),
                "column_name": col,
                "description": f"Mean {metric} in {tract} on FA skeleton",
                "modality_guess": "dMRI_TBSS",
                "metric_guess": metric,
                "tract_or_region_guess": tract,
                "panel": "primary",
                "family": "exploratory",
                "transform": "rank_inverse_normal",
                "include": True,
                "source": "basket",
                "region_priority": 0,
            })
    manifest = pd.DataFrame(manifest_rows)
    return wide, manifest
