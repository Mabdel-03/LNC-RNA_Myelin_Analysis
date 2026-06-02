#!/usr/bin/env python
"""Step 07 — Disease-risk regression for rs2546890-A on multiple sclerosis.

MS cases are derived from the UKB basket using ICD-10 G35:
  - main hospital diagnoses     (field 41202)
  - secondary hospital diagnoses (field 41204)
  - all-source ICD-10 summary    (field 41270)
  - self-report non-cancer       (field 20002 == 1261)
  - first-occurrence G35 date    (field 131042; source field 131043 retained
                                  as provenance only)

Headline regression (additive dosage, logistic):
    logit(P(ms_case)) = b0 + b_A * dosage_A + covariates
Pairwise: AA-vs-GG, AG-vs-GG, AA-vs-AG (each as a 2-group Logit with the
listed reference removed). Dominant / recessive sensitivities follow the
same pattern as `05a_run_ols.py`.

Reads:
  results/intermediate/variant_dosage.tsv  (eid, dosage_A, genotype)
  results/intermediate/covariates.tsv      (eid + age + sex + array + site + PCs)
  results/intermediate/keep_ols.txt        (FID IID — OLS keep set)
  inputs.basket_tab                        (UKB basket TSV for ICD-10 / self-report)

Writes:
  results/ms_phenotype.tsv          (per-eid ms_case + source provenance; gitignored)
  results/association_ms_risk.csv   (one row per model variant: additive, dominant,
                                     recessive, pairwise contrasts)
  results/ms_phenotype_audit.csv    (case-count breakdown per source field)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    make_design_matrix,
    resolve_under_repo,
    safe_read_table,
    setup_logger,
    set_seed,
)

# G35 ICD-10 = multiple sclerosis; self-report code 1261 = MS.
# Note: UKB fields 130892/130893 are F31 bipolar affective disorder, not MS.
ICD_FIELDS = ("41202", "41204", "41270")
SELF_REPORT_FIELD = "20002"
SELF_REPORT_MS_CODE = 1261
FIRST_OCCURRENCE_FIELDS = ("131042",)
FIRST_OCCURRENCE_SOURCE_FIELDS = ("131043",)
G35_RE = re.compile(r"^G35", re.IGNORECASE)


def _collect_ms_columns(header: list[str]) -> dict[str, list[str]]:
    """Group basket columns by their G35-source role."""
    out: dict[str, list[str]] = {
        "icd": [],
        "self_report": [],
        "first_occurrence": [],
        "first_occurrence_source": [],
    }
    for c in header:
        m = re.match(r"^f\.(\d{3,6})\.", c)
        if not m:
            continue
        fid = m.group(1)
        if fid in ICD_FIELDS:
            out["icd"].append(c)
        elif fid == SELF_REPORT_FIELD:
            out["self_report"].append(c)
        elif fid in FIRST_OCCURRENCE_FIELDS:
            out["first_occurrence"].append(c)
        elif fid in FIRST_OCCURRENCE_SOURCE_FIELDS:
            out["first_occurrence_source"].append(c)
    return out


def _derive_ms_case(basket: pd.DataFrame,
                     groups: dict[str, list[str]],
                     log) -> tuple[pd.DataFrame, dict[str, int]]:
    """Per-eid binary MS case + per-source case counts."""
    n = len(basket)
    flag_icd = pd.Series(False, index=basket.index)
    flag_sr = pd.Series(False, index=basket.index)
    flag_fo = pd.Series(False, index=basket.index)
    flag_fo_source = pd.Series(False, index=basket.index)

    for c in groups["icd"]:
        s = basket[c].astype(str).str.upper()
        flag_icd |= s.str.match(r"^G35", na=False)

    for c in groups["self_report"]:
        v = pd.to_numeric(basket[c], errors="coerce")
        flag_sr |= (v == SELF_REPORT_MS_CODE)

    for c in groups["first_occurrence"]:
        flag_fo |= basket[c].notna()

    for c in groups.get("first_occurrence_source", []):
        flag_fo_source |= basket[c].notna()

    any_flag = flag_icd | flag_sr | flag_fo

    audit = {
        "n_basket_rows": int(n),
        "n_icd_G35": int(flag_icd.sum()),
        "n_self_report_1261": int(flag_sr.sum()),
        "n_first_occurrence": int(flag_fo.sum()),
        "n_first_occurrence_source": int(flag_fo_source.sum()),
        "n_ms_cases_union": int(any_flag.sum()),
    }
    log.info(
        f"MS case derivation: ICD={audit['n_icd_G35']:,}  "
        f"self-report={audit['n_self_report_1261']:,}  "
        f"first-occ={audit['n_first_occurrence']:,}  "
        f"union={audit['n_ms_cases_union']:,} of {n:,} basket rows"
    )

    out = pd.DataFrame({
        "eid": basket["eid"],
        "ms_case": any_flag.astype(int),
        "src_icd_G35": flag_icd.astype(int),
        "src_self_report_1261": flag_sr.astype(int),
        "src_first_occurrence": flag_fo.astype(int),
        "src_first_occurrence_source": flag_fo_source.astype(int),
    })
    return out, audit


def _covariate_design(merged: pd.DataFrame, cfg: dict,
                      max_missing_frac: float = 0.30,
                      log=None) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Mirror of _covariate_design in 05a_run_ols.py — but without imaging
    confounds (which only exist for the imaging subsample). MS-risk is
    estimated on the full WB-unrelated sample, not the imaging subset.

    Continuous covars with > max_missing_frac NaN in the analysis sample are
    dropped (catches imaging-only fields like bmi/head_size/motion that
    bleed into covariates.tsv).
    """
    cont = []
    cat = []
    base_cont = ("age", "age2", "age_sex", "bmi")
    for c in base_cont:
        if c in merged.columns:
            cont.append(c)
    pc_cols = [c for c in merged.columns
               if c.startswith("PC") and c[2:].isdigit()
               and int(c[2:]) <= int(cfg["covariates"].get("genetic_pc_count", 20))]
    cont.extend(sorted(pc_cols, key=lambda c: int(c[2:])))
    cont = list(dict.fromkeys(cont))
    # Drop continuous covars with high missingness in this analysis sample.
    keep_cont = []
    for c in cont:
        frac_nan = float(merged[c].isna().mean()) if c in merged.columns else 1.0
        if frac_nan > max_missing_frac:
            if log:
                log.info(f"  dropping continuous covar '{c}' "
                         f"(missing fraction {frac_nan:.1%} > {max_missing_frac:.0%})")
            continue
        keep_cont.append(c)
    cont = keep_cont
    # sex + array are sufficient. `site` is the imaging-center field (only
    # populated for the ~40K MRI participants) — including it would restrict
    # the MS-risk regression to the imaging subsample and drop ~85% of cases.
    # `genotype_batch` (~95 levels) drives perfect separation in this
    # low-prevalence binary outcome.
    for c in ("sex", "array"):
        if c in merged.columns:
            cat.append(c)
    X = make_design_matrix(merged, cont_cols=cont, cat_cols=cat,
                            add_intercept=True, drop_first=True)
    return X, cont, cat


def _fit_logit(y: pd.Series,
               exposure: pd.Series | pd.DataFrame,
               covariates: pd.DataFrame,
               exposure_names: list[str]) -> dict:
    """Logit fit with NaN-aware joint dropna. Returns per-exposure OR/CI/raw_p."""
    import statsmodels.api as sm

    if isinstance(exposure, pd.Series):
        exp_df = exposure.to_frame(exposure_names[0])
    else:
        exp_df = exposure.copy()
        exp_df.columns = exposure_names
    frame = pd.concat([y.rename("__y__"), exp_df, covariates], axis=1).dropna()
    n = len(frame)
    if n < 30:
        return {"status": "too_few_obs", "n": n}
    y_ = frame["__y__"].astype(int).values
    n_case = int(y_.sum())
    n_ctrl = int((1 - y_).sum())
    if n_case < 20 or n_ctrl < 20:
        return {"status": "too_few_cases_or_controls", "n": n,
                "n_case": n_case, "n_ctrl": n_ctrl}
    X = frame.drop(columns="__y__")
    try:
        model = sm.Logit(y_, X.values).fit(disp=0, maxiter=200)
    except Exception as exc:  # noqa: BLE001
        return {"status": f"fit_error:{exc.__class__.__name__}", "n": n,
                "n_case": n_case, "n_ctrl": n_ctrl}
    cols = list(X.columns)
    out: dict = {
        "status": "ok",
        "n": int(n),
        "n_case": n_case,
        "n_ctrl": n_ctrl,
        "converged": bool(model.mle_retvals.get("converged", True)),
    }
    ci = model.conf_int()
    for ename in exposure_names:
        if ename not in cols:
            out[f"beta_{ename}"] = np.nan
            out[f"se_{ename}"] = np.nan
            out[f"or_{ename}"] = np.nan
            out[f"or_lo_{ename}"] = np.nan
            out[f"or_hi_{ename}"] = np.nan
            out[f"p_{ename}"] = np.nan
            continue
        i = cols.index(ename)
        beta = float(model.params[i])
        se = float(model.bse[i])
        out[f"beta_{ename}"] = beta
        out[f"se_{ename}"] = se
        out[f"or_{ename}"] = float(np.exp(beta))
        out[f"or_lo_{ename}"] = float(np.exp(ci.iloc[i, 0] if hasattr(ci, "iloc") else ci[i, 0]))
        out[f"or_hi_{ename}"] = float(np.exp(ci.iloc[i, 1] if hasattr(ci, "iloc") else ci[i, 1]))
        out[f"p_{ename}"] = float(model.pvalues[i])
    return out


def _row(model: str, contrast: str, res: dict, exposure_name: str, cfg: dict) -> dict:
    ea = cfg["variant"]["effect_allele"]
    oa = cfg["variant"]["other_allele"]
    return {
        "phenotype": "MS_G35",
        "variant": cfg["variant"]["rsid"],
        "effect_allele": ea,
        "other_allele": oa,
        "engine": "OLS-equivalent_logistic",
        "model": model,
        "contrast": contrast,
        "exposure": exposure_name,
        "n": res.get("n"),
        "n_case": res.get("n_case"),
        "n_ctrl": res.get("n_ctrl"),
        "converged": res.get("converged"),
        "beta": res.get(f"beta_{exposure_name}"),
        "se": res.get(f"se_{exposure_name}"),
        "OR": res.get(f"or_{exposure_name}"),
        "OR_lo95": res.get(f"or_lo_{exposure_name}"),
        "OR_hi95": res.get(f"or_hi_{exposure_name}"),
        "raw_p": res.get(f"p_{exposure_name}"),
        "status": res.get("status"),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--keep", default="ols", choices=["ols", "lmm"],
                        help="Which keep list to filter to (default: ols / unrelated)")
    parser.add_argument("--reuse-phenotype", action="store_true",
                        help="Skip basket re-read when results/ms_phenotype.tsv exists")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    set_seed(int(cfg["project"].get("random_seed", 42)))
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", False))

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    log = setup_logger("ms_risk", logs_dir / "run_log.txt")

    basket = cfg["inputs"]["basket_tab"]
    dosage_path = inter_dir / "variant_dosage.tsv"
    covars_path = inter_dir / "covariates.tsv"
    keep_path = inter_dir / (f"keep_{args.keep}.txt")

    out_ms_pheno = results_dir / "ms_phenotype.tsv"
    out_assoc = results_dir / "association_ms_risk.csv"
    out_audit = results_dir / "ms_phenotype_audit.csv"

    for p in (Path(basket), dosage_path, covars_path, keep_path):
        if not Path(p).is_file():
            log.error(f"required input missing: {p}")
            return 1

    if dry:
        log.info("DRY RUN — would derive MS cases and fit logistic regression")
        return 0

    if args.reuse_phenotype and out_ms_pheno.is_file():
        log.info(f"--reuse-phenotype set; loading existing {out_ms_pheno} "
                 "(skipping basket re-read)")
        ms = safe_read_table(out_ms_pheno)
        ms["eid"] = pd.to_numeric(ms["eid"], errors="coerce").astype("Int64")
        audit = {
            "n_basket_rows": int(len(ms)),
            "n_ms_cases_union": int(ms["ms_case"].sum()),
            "source": "reused_from_disk",
        }
        log.info(f"reused {len(ms):,} rows, {audit['n_ms_cases_union']:,} MS cases")
    else:
        # 1. Discover G35-relevant columns in basket header
        from utils import read_header_only
        header = read_header_only(basket)
        eid_col = "f.eid" if "f.eid" in header else ("eid" if "eid" in header else None)
        if eid_col is None:
            log.error("basket has no eid / f.eid column")
            return 2
        groups = _collect_ms_columns(header)
        use_cols = (
            [eid_col]
            + groups["icd"]
            + groups["self_report"]
            + groups["first_occurrence"]
            + groups["first_occurrence_source"]
        )
        log.info(f"basket subset: {len(use_cols):,} cols "
                 f"(icd={len(groups['icd'])}, sr={len(groups['self_report'])}, "
                 f"fo={len(groups['first_occurrence'])}, "
                 f"fo_source={len(groups['first_occurrence_source'])})")

        # 2. Read basket subset (object dtype for ICD strings, numeric for sr / fo dates)
        log.info(f"reading {len(use_cols):,} basket columns; this may take ~30-60s")
        basket_df = safe_read_table(basket, usecols=use_cols, dtype=str)
        basket_df = basket_df.rename(columns={eid_col: "eid"})
        basket_df["eid"] = pd.to_numeric(basket_df["eid"], errors="coerce").astype("Int64")
        log.info(f"basket subset loaded: {basket_df.shape}")

        # 3. Derive MS case flag
        ms, audit = _derive_ms_case(basket_df, groups, log)
        ms.to_csv(out_ms_pheno, sep="\t", index=False)
        pd.DataFrame([audit]).to_csv(out_audit, index=False)
        log.info(f"wrote {out_ms_pheno}")
        log.info(f"wrote {out_audit}")

    # 4. Read dosage + covariates and join on eid
    dose = safe_read_table(dosage_path)
    dose["eid"] = pd.to_numeric(dose["eid"], errors="coerce").astype("Int64")
    covars = safe_read_table(covars_path)
    covars["eid"] = pd.to_numeric(covars["eid"], errors="coerce").astype("Int64")
    # age2 if missing
    if "age" in covars.columns and "age2" not in covars.columns:
        covars["age2"] = covars["age"].astype(float) ** 2
    if "age" in covars.columns and "sex" in covars.columns and "age_sex" not in covars.columns:
        try:
            covars["age_sex"] = covars["age"].astype(float) * covars["sex"].astype(float)
        except Exception:
            pass

    merged = (ms[["eid", "ms_case"]]
              .merge(dose, on="eid", how="inner")
              .merge(covars, on="eid", how="inner"))
    log.info(f"merged (ms + dose + covars): {merged.shape}")
    # drop rows with sentinel negative eids from variant extract
    merged = merged[merged["eid"] >= 0].copy()

    # 5. Filter to keep set
    keep = pd.read_csv(keep_path, sep=r"\s+", header=None, names=["FID", "IID"])
    keep_set = set(int(x) for x in keep["IID"])
    before = len(merged)
    merged = merged[merged["eid"].astype("Int64").isin(keep_set)].copy()
    log.info(f"keep_{args.keep} filter: {len(merged):,} of {before:,} rows retained")
    n_cases_keep = int(merged["ms_case"].sum())
    log.info(f"after keep filter: {n_cases_keep:,} MS cases / "
             f"{len(merged)-n_cases_keep:,} controls")

    # 6. Build covariate design
    X_covar, cont_cols, cat_cols = _covariate_design(merged, cfg, log=log)
    log.info(f"covariate design: {len(cont_cols)} cont + {len(cat_cols)} cat "
             f"-> ncols={X_covar.shape[1]}")

    ea = cfg["variant"]["effect_allele"]
    oa = cfg["variant"]["other_allele"]
    aa = ea + ea
    ag = "".join(sorted([ea, oa]))
    gg = oa + oa

    y = merged["ms_case"].astype(float)
    rows: list[dict] = []

    def _fmt(x, spec=".3g"):
        return format(x, spec) if isinstance(x, (int, float)) and not pd.isna(x) else "NA"

    # A. Additive dosage
    res = _fit_logit(y, merged["dosage_A"], X_covar, ["dosage_A"])
    rows.append(_row("additive_dosage", f"per-{ea}-allele OR", res, "dosage_A", cfg))
    log.info(f"additive: status={res.get('status')}  OR={_fmt(res.get('or_dosage_A'))} "
             f"[{_fmt(res.get('or_lo_dosage_A'))}-{_fmt(res.get('or_hi_dosage_A'))}] "
             f"p={_fmt(res.get('p_dosage_A'))}  n_case={res.get('n_case')} "
             f"n_ctrl={res.get('n_ctrl')}")

    # B. Dominant: AA|AG vs GG
    g = merged.get("genotype")
    if g is not None:
        dom = g.map({gg: 0.0, ag: 1.0, aa: 1.0}).rename("dominant_A")
        res = _fit_logit(y, dom, X_covar, ["dominant_A"])
        rows.append(_row("dominant_A", f"{aa}|{ag}_vs_{gg}", res, "dominant_A", cfg))
        log.info(f"dominant: status={res.get('status')} OR={_fmt(res.get('or_dominant_A'))} "
                 f"p={_fmt(res.get('p_dominant_A'))}")

        # C. Recessive: AA vs AG|GG
        rec = g.map({gg: 0.0, ag: 0.0, aa: 1.0}).rename("recessive_A")
        res = _fit_logit(y, rec, X_covar, ["recessive_A"])
        rows.append(_row("recessive_A", f"{aa}_vs_{ag}|{gg}", res, "recessive_A", cfg))
        log.info(f"recessive: status={res.get('status')} OR={_fmt(res.get('or_recessive_A'))} "
                 f"p={_fmt(res.get('p_recessive_A'))}")

        # D. Pairwise: each pair (AA-GG, AG-GG, AA-AG) restricted to the 2 genotypes
        for a, b in (("AA", "GG"), ("AG", "GG"), ("AA", "AG")):
            sub_mask = g.isin([a, b])
            if sub_mask.sum() < 30:
                rows.append(_row(f"pairwise_{a}_vs_{b}", f"{a}_vs_{b}",
                                  {"status": "too_few_obs", "n": int(sub_mask.sum())},
                                  "geno_indicator", cfg))
                continue
            sub_y = y[sub_mask]
            sub_x = X_covar.loc[sub_mask]
            indicator = (g[sub_mask] == a).astype(float).rename("geno_indicator")
            res = _fit_logit(sub_y, indicator, sub_x, ["geno_indicator"])
            rows.append(_row(f"pairwise_{a}_vs_{b}", f"{a}_vs_{b}",
                              res, "geno_indicator", cfg))
            log.info(f"pairwise {a} vs {b}: status={res.get('status')} "
                     f"OR={_fmt(res.get('or_geno_indicator'))} "
                     f"p={_fmt(res.get('p_geno_indicator'))}")

    # 7. Write
    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_assoc, index=False)
    log.info(f"wrote {out_assoc}  ({len(out_df)} rows)")

    print("\n=== MS-risk summary (rs2546890-A) ===")
    print(out_df[["model", "contrast", "n", "n_case", "n_ctrl",
                  "OR", "OR_lo95", "OR_hi95", "raw_p", "status"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
