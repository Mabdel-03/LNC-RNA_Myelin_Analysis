#!/usr/bin/env python
"""Step 04c — write REGENIE/BOLT-formatted phenotype + covariate + keep files.

Reads:
  results/intermediate/phenotypes_extended.tsv (analysis_ready + composites + ROI PCs)
  results/intermediate/phenotype_manifest_extended.csv
  results/intermediate/keep_lmm.txt
Phenotypes are pre-transformed the same way OLS does so REGENIE/BOLT effects
are directly comparable: |z|>6 trimmed, WMH log1p'd, then inverse-rank-normal.
Composite + ROI-PC phenotypes are already standardized — no extra transform.

Writes:
  results/lmm_inputs/phenotypes.tsv         FID IID + every analysis phenotype
  results/lmm_inputs/covariates.tsv         FID IID + quant + categorical covars
  results/lmm_inputs/keep_lmm.txt           (copied from intermediate/)
  results/lmm_inputs/pheno_col_list.txt     one phenotype per line
  results/lmm_inputs/covar_col_quant.txt    quantitative covariate names
  results/lmm_inputs/covar_col_cat.txt      categorical covariate names
  results/lmm_inputs/excluded_phenotypes.csv  phenotypes dropped (with reason)
  results/lmm_inputs/target_variants.txt    REGENIE step-2 --extract content
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    inverse_rank_normalize,
    load_config,
    resolve_under_repo,
    setup_logger,
    standardize_eid_column,
    trim_outliers_z,
)


def _transform_for_lmm(series: pd.Series, transform: str, z_thresh: float) -> pd.Series:
    """Apply the same preprocessing the OLS step would apply, so the LMM β
    is on the same scale as the OLS β.

    - "none"               → return as-is (composites + ROI PCs already standardized)
    - "log1p_then_rank_inverse" → trim outliers, log1p, INRT
    - default              → trim outliers, INRT
    """
    if transform == "none":
        return pd.to_numeric(series, errors="coerce")
    s = pd.to_numeric(series, errors="coerce")
    s = trim_outliers_z(s, z_thresh=z_thresh)
    if transform == "log1p_then_rank_inverse":
        s = np.log1p(s.where(s >= 0))
    return inverse_rank_normalize(s)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    lmm_dir = results_dir / "lmm_inputs"
    lmm_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logger("prepare_lmm", logs_dir / "run_log.txt")

    ext_path = inter_dir / "phenotypes_extended.tsv"
    manifest_path = inter_dir / "phenotype_manifest_extended.csv"
    keep_src = inter_dir / "keep_lmm.txt"

    pheno_out = lmm_dir / "phenotypes.tsv"
    covar_out = lmm_dir / "covariates.tsv"
    keep_out = lmm_dir / "keep_lmm.txt"
    pheno_list_out = lmm_dir / "pheno_col_list.txt"
    quant_list_out = lmm_dir / "covar_col_quant.txt"
    cat_list_out = lmm_dir / "covar_col_cat.txt"
    excluded_out = lmm_dir / "excluded_phenotypes.csv"
    targets_out = lmm_dir / "target_variants.txt"

    if dry:
        log.info("DRY RUN — checking expected inputs")
        for p in (ext_path, manifest_path, keep_src):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    for p in (ext_path, manifest_path, keep_src):
        if not p.is_file():
            log.error(f"required input missing: {p}")
            return 1

    log.info(f"loading {ext_path}")
    wide = pd.read_csv(ext_path, sep="\t")
    wide = standardize_eid_column(wide)
    manifest = pd.read_csv(manifest_path)

    # Filter to LMM keep
    log.info(f"applying LMM keep list: {keep_src}")
    keep_df = pd.read_csv(keep_src, sep=r"\s+", header=None, names=["FID", "IID"])
    keep_eids = set(int(x) for x in keep_df["IID"])
    wide = wide[wide["eid"].astype("Int64").isin(keep_eids)].copy()
    log.info(f"  kept {len(wide):,} rows (LMM superset post-keep)")

    # ---- phenotype side
    included = manifest[manifest["include"].astype(bool)].copy() \
        if "include" in manifest.columns else manifest.copy()
    if "include_in_lmm" in included.columns:
        before_lmm_flag = len(included)
        included = included[included["include_in_lmm"].astype(str).str.lower().isin(("true", "1"))].copy()
        log.info(f"LMM include_in_lmm filter: kept {len(included)} of {before_lmm_flag} phenotypes")
    # Apply LMM scope: which families go to REGENIE/BOLT. Exploratory (single
    # IDPs) is opt-in via lmm.scope_include_exploratory because it adds 600+
    # phenotypes and dominates the wall-clock.
    scope_families = list(cfg.get("lmm", {}).get("scope_families", []))
    include_explor = bool(cfg.get("lmm", {}).get("scope_include_exploratory", False))
    if include_explor and "exploratory" not in scope_families:
        scope_families.append("exploratory")
    if scope_families and "family" in included.columns:
        before = len(included)
        included = included[included["family"].astype(str).isin(scope_families)].copy()
        log.info(f"LMM scope filter: kept {len(included)} of {before} phenotypes "
                 f"(families={scope_families})")
    max_pheno = int(cfg.get("lmm", {}).get("max_phenotypes", 0) or 0)
    if max_pheno > 0 and len(included) > max_pheno:
        included = included.head(max_pheno)
        log.info(f"LMM max_phenotypes cap: limited to first {max_pheno} phenotypes")
    pheno_cols = [c for c in included["column_name"] if c in wide.columns]
    transforms = dict(zip(included["column_name"], included.get("transform", "rank_inverse_normal")))
    z_thresh = float(cfg.get("models", {}).get("outlier_z_threshold", 6))
    min_n = int(cfg.get("lmm", {}).get("min_n_per_phenotype", 500))

    pheno_df = pd.DataFrame({"FID": wide["eid"].astype(int).values,
                              "IID": wide["eid"].astype(int).values})
    kept_pheno: list[str] = []
    excluded_rows: list[dict] = []
    for col in pheno_cols:
        trf = transforms.get(col, "rank_inverse_normal")
        s = _transform_for_lmm(wide[col], trf, z_thresh)
        n_nonna = int(s.notna().sum())
        if n_nonna < min_n:
            excluded_rows.append({"column_name": col, "n_nonna": n_nonna,
                                   "reason": f"n_nonna<{min_n}"})
            continue
        pheno_df[col] = s.values
        kept_pheno.append(col)
    log.info(f"  kept {len(kept_pheno)} phenotypes (dropped {len(excluded_rows)} below n={min_n})")

    pheno_df.to_csv(pheno_out, sep="\t", index=False, na_rep="NA")
    log.info(f"wrote {pheno_out} ({len(pheno_df):,} rows × {len(pheno_df.columns)} cols)")

    pheno_list_out.write_text("\n".join(kept_pheno) + "\n")
    excluded_out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(excluded_rows).to_csv(excluded_out, index=False)
    log.info(f"wrote {pheno_list_out} and {excluded_out}")

    # ---- covariate side
    # Quantitative: lean covariates plus focused imaging-QC/pathology/vascular fields.
    n_pcs = int(cfg.get("covariates", {}).get("genetic_pc_count", 20))
    quant_candidates = ["age", "age2", "age_sex", "head_size", "motion_dmri",
                        "t1_motion", "bmi"] + [f"PC{i}" for i in range(1, n_pcs + 1)]
    quant_prefixes = ("scanner_pos_", "dmri_qc_", "t1_qc_", "modality_discrepancy_", "wmh_")
    quant_present = [c for c in quant_candidates if c in wide.columns]
    quant_present.extend([c for c in wide.columns if c.startswith(quant_prefixes)])
    quant_present = list(dict.fromkeys(quant_present))
    # NOTE: `genotype_batch` has ~106 unique levels in this UKB extract, which
    # exceeds REGENIE's --maxCatLevels (default 10). Batch effects are captured
    # adequately via the PCs and `array`, so we omit genotype_batch as a
    # categorical here. Same reasoning for any cat covar with >50 levels —
    # filtered out below after population check.
    cat_candidates = ["sex", "array", "site", "smoking_status",
                      "diabetes", "hypertension_6150", "hypertension_6177",
                      "hypertension_self_report"]
    cat_prefixes = ("protocol_",)
    cat_present = [c for c in cat_candidates if c in wide.columns]
    cat_present.extend([c for c in wide.columns if c.startswith(cat_prefixes)])
    cat_present = list(dict.fromkeys(cat_present))
    # Drop cat covars that REGENIE will reject (default --maxCatLevels=10) or
    # that have only one level (REGENIE silently warns and ignores them).
    keep_cat = []
    for c in cat_present:
        n_levels = wide[c].dropna().nunique()
        if n_levels < 2:
            log.info(f"  dropping cat covar '{c}' (only {n_levels} unique level)")
            continue
        if n_levels > 10:
            log.info(f"  dropping cat covar '{c}' ({n_levels} levels — exceeds REGENIE default maxCatLevels)")
            continue
        keep_cat.append(c)
    cat_present = keep_cat
    if not quant_present and not cat_present:
        log.error("no recognized covariates found in extended table — check step 03 output")
        return 2

    covar_df = pd.DataFrame({"FID": wide["eid"].astype(int).values,
                              "IID": wide["eid"].astype(int).values})
    for c in quant_present + cat_present:
        covar_df[c] = wide[c].values
    covar_df.to_csv(covar_out, sep="\t", index=False, na_rep="NA")
    log.info(f"wrote {covar_out} ({len(covar_df):,} rows × {len(covar_df.columns)} cols; "
              f"quant={len(quant_present)}, cat={len(cat_present)})")

    quant_list_out.write_text("\n".join(quant_present) + "\n")
    cat_list_out.write_text("\n".join(cat_present) + "\n")
    log.info(f"wrote {quant_list_out} and {cat_list_out}")

    # ---- keep file
    shutil.copyfile(keep_src, keep_out)
    log.info(f"copied keep list to {keep_out}")

    # ---- target variant for step 2
    target = cfg.get("lmm", {}).get("target_variant_id", "5:158759900:A:G")
    targets_out.write_text(str(target) + "\n")
    log.info(f"wrote {targets_out} (target = {target})")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
