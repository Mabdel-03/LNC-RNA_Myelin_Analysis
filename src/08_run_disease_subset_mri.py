#!/usr/bin/env python
"""Step 08 — clinical neurodegeneration/MS subset MRI interaction models.

Derives clinical-only disease flags, then tests whether rs2546890-A MRI effects
differ in disease-enriched subsets. The primary model is an interaction:

    IDP ~ dosage_A + disease_flag + dosage_A:disease_flag + covariates

The runner also reports within-case and within-noncase additive estimates.
Outputs use raw p values only, matching the rest of the project.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from disease_subsets import (  # noqa: E402
    collect_disease_columns,
    derive_disease_flags,
    disease_use_columns,
    fit_interaction_model,
    fit_stratified_model,
)
from utils import (  # noqa: E402
    inverse_rank_normalize,
    load_config,
    make_design_matrix,
    read_header_only,
    resolve_under_repo,
    safe_read_table,
    setup_logger,
    standardize_eid_column,
    trim_outliers_z,
)


DEFAULT_FLAGS = ("broad_neurodeg_case", "ms_case")
FLAG_LABELS = {
    "broad_neurodeg_case": "Broad clinical neurodegeneration, excluding MS definition",
    "ms_case": "Multiple sclerosis current union, corrected G35-first-occurrence",
}
DEFAULT_EXTRA_SCOPE_FAMILIES = ("secondary_composite", "secondary_roi_pca")


def _transform_phenotype(series: pd.Series, transform: str, z_thresh: float) -> pd.Series:
    if transform == "none":
        return pd.to_numeric(series, errors="coerce")
    out = trim_outliers_z(pd.to_numeric(series, errors="coerce"), z_thresh=z_thresh)
    if transform.startswith("log1p"):
        out = pd.Series(np.where(out >= 0, np.log1p(out), np.nan), index=out.index)
    return inverse_rank_normalize(out)


def _add_raw_p(row: dict) -> dict:
    if row.get("analysis") == "interaction":
        row["raw_p"] = row.get("p_interaction")
    else:
        row["raw_p"] = row.get("p_per_A")
    raw_p = row.get("raw_p")
    row["minus_log10_raw_p"] = (
        -np.log10(max(float(raw_p), np.finfo(float).tiny))
        if raw_p is not None and pd.notna(raw_p)
        else np.nan
    )
    return row


def _metadata(row: pd.Series, disease_flag: str, cfg: dict) -> dict:
    return {
        "disease_flag": disease_flag,
        "disease_label": FLAG_LABELS.get(disease_flag, disease_flag),
        "phenotype": row.get("column_name", ""),
        "variant": cfg["variant"]["rsid"],
        "effect_allele": cfg["variant"]["effect_allele"],
        "other_allele": cfg["variant"]["other_allele"],
        "family": row.get("family", ""),
        "analysis_tier": row.get("analysis_tier", ""),
        "confirmatory_role": row.get("confirmatory_role", ""),
        "panel": row.get("panel", ""),
        "metric": row.get("metric_guess", row.get("metric", "")),
        "region": row.get("tract_or_region_guess", row.get("region", "")),
        "source": row.get("source", ""),
    }


def _primary_covariate_design(merged: pd.DataFrame,
                              cfg: dict,
                              mode: str = "primary") -> tuple[pd.DataFrame, list[str], list[str]]:
    n_pcs = int(cfg.get("covariates", {}).get("genetic_pc_count", 20))
    cont = []
    base_cont = ("age", "age2", "age_sex", "head_size", "motion_dmri", "t1_motion")
    if mode == "full":
        base_cont = base_cont + ("bmi",)
    for col in base_cont:
        if col in merged.columns:
            cont.append(col)

    primary_prefixes = ("scanner_pos_", "dmri_qc_", "t1_qc_", "modality_discrepancy_")
    full_prefixes = primary_prefixes + ("wmh_",)
    prefixes = full_prefixes if mode == "full" else primary_prefixes
    cont.extend([c for c in merged.columns if c.startswith(prefixes)])
    pc_cols = [
        c for c in merged.columns
        if c.startswith("PC") and c[2:].isdigit() and int(c[2:]) <= n_pcs
    ]
    cont.extend(sorted(pc_cols, key=lambda c: int(c[2:])))
    cont = list(dict.fromkeys(cont))

    cat = []
    base_cat = ("sex", "site", "array")
    if mode == "full":
        base_cat = base_cat + (
            "smoking_status",
            "diabetes",
            "hypertension_6150",
            "hypertension_6177",
            "hypertension_self_report",
        )
    for col in base_cat:
        if col in merged.columns:
            cat.append(col)
    cat.extend([c for c in merged.columns if c.startswith("protocol_")])
    cat = list(dict.fromkeys(cat))
    return make_design_matrix(
        merged, cont_cols=cont, cat_cols=cat, add_intercept=True, drop_first=True
    ), cont, cat


def _filter_keep_ols(df: pd.DataFrame, keep_path: Path, log) -> pd.DataFrame:
    if not keep_path.is_file():
        log.warning(f"keep_ols file not found: {keep_path}; using all rows")
        return df
    keep = pd.read_csv(keep_path, sep=r"\s+", header=None, names=["FID", "IID"])
    keep_set = set(int(x) for x in keep["IID"])
    before = len(df)
    out = df[df["eid"].astype("Int64").isin(keep_set)].copy()
    log.info(f"OLS keep filter: retained {len(out):,} of {before:,} rows")
    return out


def _select_manifest(manifest: pd.DataFrame,
                     available_cols: set[str],
                     cfg: dict,
                     max_phenotypes: int | None = None) -> pd.DataFrame:
    out = manifest.copy()
    if "include" in out.columns:
        out = out[out["include"].astype(str).str.lower().isin(("true", "1"))].copy()
    subset_cfg = cfg.get("disease_subset", {})
    scope = list(subset_cfg.get("scope_families", []))
    if not scope:
        scope = list(cfg.get("lmm", {}).get("scope_families", []))
        scope.extend(DEFAULT_EXTRA_SCOPE_FAMILIES)
    if scope and "family" in out.columns:
        out = out[out["family"].astype(str).isin(scope)].copy()
    out = out[out["column_name"].astype(str).isin(available_cols)].copy()
    if cfg.get("project", {}).get("test_mode", False):
        cap = int(cfg["project"].get("max_test_phenotypes", 3))
        out = out.head(cap)
    if max_phenotypes is not None and max_phenotypes > 0:
        out = out.head(max_phenotypes)
    return out


def _load_or_derive_flags(cfg: dict,
                          results_dir: Path,
                          log,
                          reuse: bool = False) -> tuple[pd.DataFrame, dict[str, int]]:
    pheno_path = results_dir / "disease_subset_phenotypes.tsv"
    audit_path = results_dir / "disease_subset_audit.csv"
    if reuse and pheno_path.is_file():
        pheno = safe_read_table(pheno_path)
        pheno = standardize_eid_column(pheno)
        audit = {
            "n_basket_rows": int(len(pheno)),
            "n_broad_neurodeg_cases": int(pheno["broad_neurodeg_case"].sum()),
            "n_ms_cases_union": int(pheno["ms_case"].sum()),
            "source": "reused_from_disk",
        }
        return pheno, audit

    basket = cfg["inputs"]["basket_tab"]
    header = read_header_only(basket)
    groups = collect_disease_columns(header)
    use_cols = disease_use_columns(groups)
    log.info(
        f"reading {len(use_cols):,} disease basket columns "
        f"(icd={len(groups['icd'])}, ms_sr={len(groups['ms_self_report'])})"
    )
    basket_df = safe_read_table(basket, usecols=use_cols, dtype=str)
    pheno, audit = derive_disease_flags(basket_df, groups)
    pheno.to_csv(pheno_path, sep="\t", index=False)
    pd.DataFrame([audit]).to_csv(audit_path, index=False)
    log.info(f"wrote {pheno_path}")
    log.info(f"wrote {audit_path}")
    return pheno, audit


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--reuse-phenotypes", action="store_true",
                        help="Reuse results/disease_subset_phenotypes.tsv when present")
    parser.add_argument("--derive-only", action="store_true",
                        help="Only derive disease flags and audit counts")
    parser.add_argument("--flags", nargs="+", default=list(DEFAULT_FLAGS),
                        help="Disease flag columns to analyze")
    parser.add_argument("--covariate-mode", choices=["primary", "full"], default="primary",
                        help="primary excludes pathology/vascular covariates; full includes them")
    parser.add_argument("--max-phenotypes", type=int, default=0,
                        help="Optional cap for smoke tests")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", False))
    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    log = setup_logger("disease_subset_mri", logs_dir / "run_log.txt")

    if dry:
        log.info("DRY RUN — would derive disease flags and fit subset MRI models")
        return 0

    disease_pheno, audit = _load_or_derive_flags(
        cfg, results_dir, log, reuse=args.reuse_phenotypes
    )
    log.info(
        f"disease flags: broad={audit.get('n_broad_neurodeg_cases'):,}, "
        f"MS={audit.get('n_ms_cases_union'):,}"
    )
    if args.derive_only:
        return 0

    manifest_path = inter_dir / "phenotype_manifest_extended.csv"
    if not manifest_path.is_file():
        manifest_path = results_dir / "phenotype_manifest.csv"
    analysis_path = inter_dir / "phenotypes_extended.tsv"
    if not analysis_path.is_file():
        analysis_path = inter_dir / "analysis_ready.tsv"
    keep_path = inter_dir / "keep_ols.txt"
    for path in (manifest_path, analysis_path, keep_path):
        if not path.is_file():
            log.error(f"required input missing: {path}")
            return 1

    manifest = pd.read_csv(manifest_path)
    merged = pd.read_csv(analysis_path, sep="\t")
    merged = standardize_eid_column(merged)
    disease_pheno = standardize_eid_column(disease_pheno)
    merged = merged.merge(disease_pheno, on="eid", how="inner")
    merged = _filter_keep_ols(merged, keep_path, log)

    flag_cols = [flag for flag in args.flags if flag in merged.columns]
    missing_flags = sorted(set(args.flags) - set(flag_cols))
    if missing_flags:
        log.error(f"requested disease flags missing from merged table: {missing_flags}")
        return 2
    selected = _select_manifest(
        manifest,
        available_cols=set(merged.columns),
        cfg=cfg,
        max_phenotypes=args.max_phenotypes if args.max_phenotypes > 0 else None,
    )
    if selected.empty:
        log.error("no phenotypes selected for disease-subset analysis")
        return 3
    log.info(f"selected {len(selected):,} phenotypes for disease-subset models")

    x_covar, cont_cols, cat_cols = _primary_covariate_design(
        merged, cfg, mode=args.covariate_mode
    )
    log.info(
        f"covariate design ({args.covariate_mode}): {len(cont_cols)} cont + "
        f"{len(cat_cols)} cat -> ncols={x_covar.shape[1]}"
    )

    subset_cfg = cfg.get("disease_subset", {})
    vcov_type = str(subset_cfg.get("vcov_type", cfg.get("models", {}).get("ols_robust_vcov", "nonrobust")))
    min_cases = int(subset_cfg.get("min_cases_per_phenotype", 100))
    min_controls = int(subset_cfg.get("min_controls_per_phenotype", 100))
    min_cell = int(subset_cfg.get("min_case_genotype_n", 10))
    z_thresh = float(cfg.get("models", {}).get("outlier_z_threshold", 6))

    rows: list[dict] = []
    excluded: list[dict] = []
    for n_done, (_, row) in enumerate(selected.iterrows(), start=1):
        phenotype = row["column_name"]
        y = _transform_phenotype(
            merged[phenotype],
            str(row.get("transform", "rank_inverse_normal")),
            z_thresh,
        )
        for flag in flag_cols:
            meta = _metadata(row, flag, cfg)
            interaction = fit_interaction_model(
                y,
                merged["dosage_A"],
                merged[flag],
                x_covar,
                genotype=merged["genotype"] if "genotype" in merged.columns else None,
                vcov_type=vcov_type,
                min_cases=min_cases,
                min_controls=min_controls,
                min_case_genotype_n=min_cell,
            )
            if interaction.get("status") == "ok":
                rows.append(_add_raw_p({
                    **meta,
                    "analysis": "interaction",
                    "model": "dosage_A_x_disease_flag",
                    **interaction,
                }))
            else:
                excluded.append({**meta, "analysis": "interaction", **interaction})

            for stratum, label in ((1, "stratified_case"), (0, "stratified_noncase")):
                res = fit_stratified_model(
                    y,
                    merged["dosage_A"],
                    merged[flag],
                    x_covar,
                    stratum=stratum,
                    genotype=merged["genotype"] if "genotype" in merged.columns else None,
                    vcov_type=vcov_type,
                    min_case_obs=min_cases,
                    min_case_genotype_n=min_cell,
                )
                if res.get("status") == "ok":
                    rows.append(_add_raw_p({
                        **meta,
                        "analysis": label,
                        "model": "additive_dosage_within_stratum",
                        **res,
                    }))
                else:
                    excluded.append({**meta, "analysis": label, **res})
        if n_done % 25 == 0 or n_done == len(selected):
            log.info(f"processed {n_done}/{len(selected)} manifest rows")

    out_assoc = results_dir / "association_disease_subset_mri.csv"
    out_excluded = results_dir / "disease_subset_excluded_models.csv"
    rows_df = pd.DataFrame(rows)
    if not rows_df.empty:
        rows_df = rows_df.sort_values(["disease_flag", "analysis", "raw_p"],
                                      na_position="last")
    rows_df.to_csv(out_assoc, index=False)
    pd.DataFrame(excluded).to_csv(out_excluded, index=False)
    log.info(f"wrote {out_assoc} ({len(rows):,} rows)")
    log.info(f"wrote {out_excluded} ({len(excluded):,} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
