"""Clinical disease flags and MRI subset-model helpers.

This module keeps disease derivation testable and separate from the command
line runner. The broad neurodegeneration flag is clinical-only by design:
MRI-derived measures are outcomes in the subset analysis and are not used to
define case status.
"""
from __future__ import annotations

import re
from collections.abc import Sequence

import numpy as np
import pandas as pd
from scipy import stats

from utils import fit_ols_model, standardize_eid_column


ICD_FIELDS = ("41202", "41204", "41270")

MS_SELF_REPORT_FIELD = "20002"
MS_SELF_REPORT_CODE = 1261
MS_FIRST_OCCURRENCE_DATE_FIELDS = ("131042",)
MS_FIRST_OCCURRENCE_SOURCE_FIELDS = ("131043",)

BROAD_ALGORITHMIC_DATE_FIELDS = {
    "broad_alg_dementia": ("42018", "42020", "42022", "42024"),
    "broad_alg_mnd": ("42028",),
    "broad_alg_parkinsonism": ("42030", "42032", "42034", "42036"),
}

BROAD_FIRST_OCCURRENCE_DATE_FIELDS = {
    "broad_fo_dementia_ad": ("130836", "130838", "130840", "130842", "131036"),
    "broad_fo_parkinsonism": ("131022", "131024", "131026", "131028"),
    "broad_fo_other_degeneration": ("131012", "131016", "131038"),
}

BROAD_NEURODEGENERATION_ICD_PREFIXES = (
    "F00", "F01", "F02", "F03",
    "G10", "G12", "G20", "G21", "G22", "G23", "G30", "G31",
)
MS_ICD_PREFIXES = ("G35",)

GENOTYPE_LEVELS = ("AA", "AG", "GG")


def _field_id_from_column(col: str) -> str | None:
    match = re.match(r"^f\.(\d{3,6})\.", col)
    return match.group(1) if match else None


def _columns_for_fields(header: Sequence[str], field_ids: Sequence[str]) -> list[str]:
    wanted = set(str(x) for x in field_ids)
    return [c for c in header if (_field_id_from_column(c) in wanted)]


def collect_disease_columns(header: Sequence[str]) -> dict[str, list[str] | str]:
    """Return grouped basket columns needed for broad neurodegeneration + MS."""
    eid_col = "f.eid" if "f.eid" in header else ("eid" if "eid" in header else None)
    if eid_col is None:
        raise KeyError("basket header has no eid/f.eid column")

    groups: dict[str, list[str] | str] = {"eid": eid_col}
    groups["icd"] = _columns_for_fields(header, ICD_FIELDS)
    groups["ms_self_report"] = _columns_for_fields(header, (MS_SELF_REPORT_FIELD,))
    groups["ms_first_occurrence_date"] = _columns_for_fields(
        header, MS_FIRST_OCCURRENCE_DATE_FIELDS
    )
    groups["ms_first_occurrence_source"] = _columns_for_fields(
        header, MS_FIRST_OCCURRENCE_SOURCE_FIELDS
    )
    for name, fields in BROAD_ALGORITHMIC_DATE_FIELDS.items():
        groups[name] = _columns_for_fields(header, fields)
    for name, fields in BROAD_FIRST_OCCURRENCE_DATE_FIELDS.items():
        groups[name] = _columns_for_fields(header, fields)
    return groups


def disease_use_columns(groups: dict[str, list[str] | str]) -> list[str]:
    """Flatten disease-column groups into a de-duplicated usecols list."""
    out: list[str] = [str(groups["eid"])]
    for key, cols in groups.items():
        if key == "eid":
            continue
        out.extend(cols)
    return list(dict.fromkeys(out))


def _any_notna(df: pd.DataFrame, cols: Sequence[str]) -> pd.Series:
    if not cols:
        return pd.Series(False, index=df.index)
    return df[list(cols)].notna().any(axis=1)


def _any_numeric_code(df: pd.DataFrame, cols: Sequence[str], code: int) -> pd.Series:
    out = pd.Series(False, index=df.index)
    for col in cols:
        out |= pd.to_numeric(df[col], errors="coerce").eq(code)
    return out


def _any_icd_prefix(df: pd.DataFrame,
                    cols: Sequence[str],
                    prefixes: Sequence[str]) -> pd.Series:
    out = pd.Series(False, index=df.index)
    if not cols or not prefixes:
        return out
    pattern = re.compile(r"^(?:" + "|".join(re.escape(p) for p in prefixes) + ")",
                         re.IGNORECASE)
    for col in cols:
        out |= df[col].astype("string").str.match(pattern, na=False)
    return out


def derive_disease_flags(basket: pd.DataFrame,
                         groups: dict[str, list[str] | str] | None = None
                         ) -> tuple[pd.DataFrame, dict[str, int]]:
    """Derive broad neurodegeneration and corrected MS flags from a basket slice."""
    if groups is None:
        groups = collect_disease_columns(list(basket.columns))
    eid_col = str(groups["eid"])
    df = basket.rename(columns={eid_col: "eid"}).copy()
    df = standardize_eid_column(df)

    icd_cols = list(groups.get("icd", []))
    broad_icd = _any_icd_prefix(df, icd_cols, BROAD_NEURODEGENERATION_ICD_PREFIXES)
    ms_icd = _any_icd_prefix(df, icd_cols, MS_ICD_PREFIXES)
    ms_self_report = _any_numeric_code(
        df, list(groups.get("ms_self_report", [])), MS_SELF_REPORT_CODE
    )
    ms_first_occurrence = _any_notna(df, list(groups.get("ms_first_occurrence_date", [])))
    ms_first_occurrence_source = _any_notna(
        df, list(groups.get("ms_first_occurrence_source", []))
    )

    broad_sources: dict[str, pd.Series] = {"broad_icd": broad_icd}
    for name in BROAD_ALGORITHMIC_DATE_FIELDS:
        broad_sources[name] = _any_notna(df, list(groups.get(name, [])))
    for name in BROAD_FIRST_OCCURRENCE_DATE_FIELDS:
        broad_sources[name] = _any_notna(df, list(groups.get(name, [])))

    broad_neurodeg = pd.Series(False, index=df.index)
    for source_flag in broad_sources.values():
        broad_neurodeg |= source_flag
    ms_case = ms_icd | ms_self_report | ms_first_occurrence

    out = pd.DataFrame({
        "eid": df["eid"],
        "broad_neurodeg_case": broad_neurodeg.astype(int),
        "ms_case": ms_case.astype(int),
        "src_broad_icd": broad_icd.astype(int),
        "src_broad_alg_dementia": broad_sources["broad_alg_dementia"].astype(int),
        "src_broad_alg_mnd": broad_sources["broad_alg_mnd"].astype(int),
        "src_broad_alg_parkinsonism": broad_sources["broad_alg_parkinsonism"].astype(int),
        "src_broad_fo_dementia_ad": broad_sources["broad_fo_dementia_ad"].astype(int),
        "src_broad_fo_parkinsonism": broad_sources["broad_fo_parkinsonism"].astype(int),
        "src_broad_fo_other_degeneration": broad_sources[
            "broad_fo_other_degeneration"
        ].astype(int),
        "src_ms_icd_G35": ms_icd.astype(int),
        "src_ms_self_report_1261": ms_self_report.astype(int),
        "src_ms_first_occurrence_G35": ms_first_occurrence.astype(int),
        "src_ms_first_occurrence_source_G35": ms_first_occurrence_source.astype(int),
    })
    audit = {
        "n_basket_rows": int(len(out)),
        "n_broad_neurodeg_cases": int(out["broad_neurodeg_case"].sum()),
        "n_ms_cases_union": int(out["ms_case"].sum()),
        "n_overlap_broad_neurodeg_ms": int(
            ((out["broad_neurodeg_case"] == 1) & (out["ms_case"] == 1)).sum()
        ),
    }
    for col in out.columns:
        if col.startswith("src_"):
            audit[f"n_{col}"] = int(out[col].sum())
    return out, audit


def _case_genotype_counts(genotype: pd.Series,
                          disease_flag: pd.Series,
                          levels: Sequence[str] = GENOTYPE_LEVELS) -> dict[str, int]:
    case_genotype = genotype.loc[disease_flag.astype(int).eq(1)].astype(str)
    counts = case_genotype.value_counts(dropna=False)
    return {level: int(counts.get(level, 0)) for level in levels}


def _gate_subset_model(frame: pd.DataFrame,
                       genotype_col: str | None = "genotype",
                       min_cases: int = 100,
                       min_controls: int = 100,
                       min_case_genotype_n: int = 10) -> dict | None:
    n = int(len(frame))
    if n < 30:
        return {"status": "too_few_obs", "n": n}
    n_case = int(frame["disease_flag"].sum())
    n_control = int(n - n_case)
    if n_case < min_cases or n_control < min_controls:
        return {
            "status": "too_few_cases_or_controls",
            "n": n,
            "n_case": n_case,
            "n_control": n_control,
        }
    if genotype_col and genotype_col in frame.columns:
        counts = _case_genotype_counts(frame[genotype_col], frame["disease_flag"])
        if any(v < min_case_genotype_n for v in counts.values()):
            return {
                "status": f"case_genotype_cell_lt_{min_case_genotype_n}",
                "n": n,
                "n_case": n_case,
                "n_control": n_control,
                **{f"n_case_{k}": v for k, v in counts.items()},
            }
    return None


def fit_interaction_model(y: pd.Series,
                          dosage: pd.Series,
                          disease_flag: pd.Series,
                          covariates: pd.DataFrame,
                          genotype: pd.Series | None = None,
                          vcov_type: str = "nonrobust",
                          min_cases: int = 100,
                          min_controls: int = 100,
                          min_case_genotype_n: int = 10) -> dict:
    """Fit y ~ dosage + disease_flag + dosage:disease_flag + covariates."""
    import statsmodels.api as sm

    parts = [
        y.rename("__y__"),
        dosage.rename("dosage_A"),
        disease_flag.rename("disease_flag").astype(float),
        covariates,
    ]
    if genotype is not None:
        parts.append(genotype.rename("genotype"))
    frame = pd.concat(parts, axis=1).dropna()
    gate = _gate_subset_model(
        frame,
        genotype_col="genotype" if genotype is not None else None,
        min_cases=min_cases,
        min_controls=min_controls,
        min_case_genotype_n=min_case_genotype_n,
    )
    if gate is not None:
        return gate

    frame["dosage_A_x_disease_flag"] = frame["dosage_A"] * frame["disease_flag"]
    drop_cols = ["__y__"] + (["genotype"] if "genotype" in frame.columns else [])
    X = frame.drop(columns=drop_cols)
    Y = frame["__y__"].astype(float)
    try:
        if vcov_type == "nonrobust":
            model = sm.OLS(Y.values, X.values, hasconst=True).fit()
        else:
            model = sm.OLS(Y.values, X.values, hasconst=True).fit(cov_type=vcov_type)
    except Exception as exc:  # noqa: BLE001
        return {"status": f"fit_error:{exc.__class__.__name__}", "n": int(len(frame))}

    cols = list(X.columns)
    i_dose = cols.index("dosage_A")
    i_int = cols.index("dosage_A_x_disease_flag")
    cov = np.asarray(model.cov_params())
    beta_case = float(model.params[i_dose] + model.params[i_int])
    var_case = float(cov[i_dose, i_dose] + cov[i_int, i_int] + 2.0 * cov[i_dose, i_int])
    se_case = float(np.sqrt(var_case)) if var_case >= 0 and np.isfinite(var_case) else np.nan
    t_case = beta_case / se_case if se_case and np.isfinite(se_case) else np.nan
    p_case = float(2.0 * stats.t.sf(abs(t_case), model.df_resid)) if np.isfinite(t_case) else np.nan

    out = {
        "status": "ok",
        "n": int(len(frame)),
        "n_case": int(frame["disease_flag"].sum()),
        "n_control": int((1 - frame["disease_flag"]).sum()),
        "beta_per_A_controls": float(model.params[i_dose]),
        "se_per_A_controls": float(model.bse[i_dose]),
        "p_per_A_controls": float(model.pvalues[i_dose]),
        "beta_interaction": float(model.params[i_int]),
        "se_interaction": float(model.bse[i_int]),
        "p_interaction": float(model.pvalues[i_int]),
        "beta_per_A_cases_from_interaction": beta_case,
        "se_per_A_cases_from_interaction": se_case,
        "p_per_A_cases_from_interaction": p_case,
        "vcov_type": vcov_type,
    }
    if genotype is not None:
        counts = _case_genotype_counts(frame["genotype"], frame["disease_flag"])
        out.update({f"n_case_{k}": v for k, v in counts.items()})
    return out


def fit_stratified_model(y: pd.Series,
                         dosage: pd.Series,
                         disease_flag: pd.Series,
                         covariates: pd.DataFrame,
                         stratum: int,
                         genotype: pd.Series | None = None,
                         vcov_type: str = "nonrobust",
                         min_case_obs: int = 100,
                         min_case_genotype_n: int = 10) -> dict:
    """Fit y ~ dosage + covariates inside one disease stratum."""
    parts = [
        y.rename("__y__"),
        dosage.rename("dosage_A"),
        disease_flag.rename("disease_flag").astype(float),
        covariates,
    ]
    if genotype is not None:
        parts.append(genotype.rename("genotype"))
    frame = pd.concat(parts, axis=1).dropna()
    frame = frame[frame["disease_flag"].astype(int).eq(int(stratum))].copy()
    n = int(len(frame))
    if stratum == 1 and n < min_case_obs:
        return {"status": "too_few_case_obs", "n": n, "stratum": int(stratum)}
    if genotype is not None and stratum == 1:
        counts = frame["genotype"].astype(str).value_counts(dropna=False)
        cell_counts = {level: int(counts.get(level, 0)) for level in GENOTYPE_LEVELS}
        if any(v < min_case_genotype_n for v in cell_counts.values()):
            return {
                "status": f"case_genotype_cell_lt_{min_case_genotype_n}",
                "n": n,
                "stratum": int(stratum),
                **{f"n_case_{k}": v for k, v in cell_counts.items()},
            }
    x_covar = frame.drop(columns=["__y__", "dosage_A", "disease_flag"]
                         + (["genotype"] if "genotype" in frame.columns else []))
    res = fit_ols_model(
        frame["__y__"],
        frame["dosage_A"],
        x_covar,
        exposure_name="dosage_A",
        vcov_type=vcov_type,
    )
    return {
        "status": res.get("status"),
        "n": res.get("n", n),
        "stratum": int(stratum),
        "beta_per_A": res.get("beta"),
        "se_per_A": res.get("se"),
        "p_per_A": res.get("p"),
        "vcov_type": res.get("vcov_type", vcov_type),
    }
