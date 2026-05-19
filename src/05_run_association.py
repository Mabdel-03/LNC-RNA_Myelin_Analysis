#!/usr/bin/env python
"""Step 05 — run per-IDP OLS + (optional) genotypic-model sensitivity.

For each IDP in phenotype_manifest.csv with include=true:
  1. Trim |z|>outlier_z_threshold on raw scale.
  2. log1p for WMH-like rows whose transform == 'log1p_then_rank_inverse'.
  3. Inverse-rank-normalize.
  4. PRIMARY model: y_inrt ~ dosage_A + age + age2 + sex + age*sex
                              + C(site) + head_size + motion_dmri
                              + C(array) + PC1..PCk
  5. GENOTYPIC sensitivity: y_inrt ~ C(genotype, ref=GG) + same covariates
                            (requires genotype column from extract step).

Multiple testing:
  - Bonferroni and BH-FDR computed across each panel separately.

Outputs:
  results/association_results_primary.csv
  results/association_results_genotypic.csv
  logs/model_warnings.txt
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    bh_fdr,
    fit_genotypic_model,
    fit_ols_model,
    inverse_rank_normalize,
    load_config,
    make_design_matrix,
    resolve_under_repo,
    setup_logger,
    set_seed,
    trim_outliers_z,
)


def _covariate_design(merged: pd.DataFrame, cfg: dict) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Decide which columns are continuous vs categorical covariates,
    return (design_df, cont_cols, cat_cols)."""
    cont = []
    cat = []
    if "age" in merged.columns:
        cont.append("age")
    if "age2" in merged.columns:
        cont.append("age2")
    if "age_sex" in merged.columns:
        cont.append("age_sex")
    if "head_size" in merged.columns:
        cont.append("head_size")
    if "motion_dmri" in merged.columns:
        cont.append("motion_dmri")
    pc_cols = [c for c in merged.columns
               if c.startswith("PC") and c[2:].isdigit()
               and int(c[2:]) <= int(cfg["covariates"].get("genetic_pc_count", 20))]
    cont.extend(sorted(pc_cols, key=lambda c: int(c[2:])))
    if "sex" in merged.columns:
        cat.append("sex")
    if "site" in merged.columns:
        cat.append("site")
    if "array" in merged.columns:
        cat.append("array")
    # Build design matrix (without exposure y; constant added)
    X = make_design_matrix(merged, cont_cols=cont, cat_cols=cat,
                            add_intercept=True, drop_first=True)
    return X, cont, cat


def _apply_correction(df: pd.DataFrame, p_col: str = "p") -> pd.DataFrame:
    out = df.copy()
    panels = out["panel"].unique()
    out["bonferroni"] = np.nan
    out["fdr_bh"] = np.nan
    for panel in panels:
        mask = (out["panel"] == panel) & out[p_col].notna()
        if not mask.any():
            continue
        n = int(mask.sum())
        out.loc[mask, "bonferroni"] = np.minimum(out.loc[mask, p_col] * n, 1.0)
        out.loc[mask, "fdr_bh"] = bh_fdr(out.loc[mask, p_col].values)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    set_seed(int(cfg["project"]["random_seed"]))
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    log = setup_logger("assoc", logs_dir / "run_log.txt")

    manifest_path = results_dir / "phenotype_manifest.csv"
    analysis_path = inter_dir / "analysis_ready.tsv"
    out_primary = results_dir / "association_results_primary.csv"
    out_genotypic = results_dir / "association_results_genotypic.csv"
    warn_path = logs_dir / "model_warnings.txt"

    if dry:
        log.info("DRY RUN — would fit OLS + genotypic models per IDP")
        for p in (manifest_path, analysis_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    if not manifest_path.is_file():
        log.error(f"manifest not found: {manifest_path}")
        return 1
    if not analysis_path.is_file():
        log.error(f"analysis-ready table not found: {analysis_path}")
        return 1

    manifest = pd.read_csv(manifest_path)
    merged = pd.read_csv(analysis_path, sep="\t")
    log.info(f"manifest: {len(manifest)} rows, analysis table: {len(merged):,} rows")

    # restrict to include=true (manifest)
    if "include" in manifest.columns:
        manifest = manifest[manifest["include"].astype(str).str.lower().isin(("true", "1"))].copy()
    log.info(f"manifest after include filter: {len(manifest)} rows")

    # test-mode cap
    if cfg.get("project", {}).get("test_mode", False):
        cap = int(cfg["project"].get("max_test_phenotypes", 3))
        manifest = manifest.head(cap)
        log.info(f"TEST MODE — capped to first {cap} phenotypes")

    # build covariate design once
    X_covar, cont_cols, cat_cols = _covariate_design(merged, cfg)
    log.info(f"covariate design: {len(cont_cols)} continuous, "
             f"{len(cat_cols)} categorical, total ncols={X_covar.shape[1]}")

    warnings: list[str] = []
    primary_rows = []
    geno_rows = []
    do_genotypic = bool(cfg["models"].get("run_genotypic_sensitivity", True)) \
        and ("genotype" in merged.columns)
    if not do_genotypic:
        log.warning("genotypic-model sensitivity disabled "
                    "(config flag false OR no hard-call genotype column).")

    z_thresh = float(cfg["models"].get("outlier_z_threshold", 6))

    for _, row in manifest.iterrows():
        idp = row["column_name"]
        if idp not in merged.columns:
            warnings.append(f"{idp}: column not in analysis table; skipped")
            continue
        y_raw = pd.to_numeric(merged[idp], errors="coerce")
        y_trim = trim_outliers_z(y_raw, z_thresh=z_thresh)
        if str(row.get("transform", "")).startswith("log1p"):
            # only valid for non-negative values
            y_trim = pd.Series(np.where(y_trim >= 0, np.log1p(y_trim), np.nan),
                                index=y_trim.index)
        y_inrt = inverse_rank_normalize(y_trim)
        # PRIMARY: additive dosage
        res = fit_ols_model(y_inrt, merged["dosage_A"], X_covar,
                            exposure_name="dosage_A")
        rec = {
            "column_name": idp,
            "field_id": row.get("field_id", ""),
            "panel": row.get("panel", ""),
            "modality": row.get("modality_guess", ""),
            "metric": row.get("metric_guess", ""),
            "region": row.get("tract_or_region_guess", ""),
            "effect_allele": cfg["variant"]["effect_allele"],
            "other_allele": cfg["variant"]["other_allele"],
            "n": res.get("n"),
            "beta_per_A": res.get("beta"),
            "se": res.get("se"),
            "t": res.get("t"),
            "p": res.get("p"),
            "aa_vs_gg_additive_2beta": (res.get("beta") * 2.0) if res.get("status") == "ok" else None,
            "status": res.get("status"),
        }
        primary_rows.append(rec)
        if res.get("status") != "ok":
            warnings.append(f"{idp}: primary model status={res.get('status')}")

        # GENOTYPIC sensitivity
        if do_genotypic:
            gres = fit_genotypic_model(y_inrt, merged["genotype"], X_covar,
                                       reference=cfg["variant"]["other_allele"] * 2)
            grec = {
                "column_name": idp,
                "panel": row.get("panel", ""),
                "modality": row.get("modality_guess", ""),
                "metric": row.get("metric_guess", ""),
                "reference": gres.get("reference"),
                "n": gres.get("n"),
                "status": gres.get("status"),
            }
            for k in gres:
                if k.startswith("beta_") or k.startswith("se_") or k.startswith("p_") \
                   or k in ("wald_p", "wald_chi2", "wald_df"):
                    grec[k] = gres[k]
            geno_rows.append(grec)
            if gres.get("status") != "ok":
                warnings.append(f"{idp}: genotypic model status={gres.get('status')}")

    primary_df = pd.DataFrame(primary_rows)
    primary_df = _apply_correction(primary_df, p_col="p")
    primary_df = primary_df.sort_values(["panel", "p"], na_position="last")
    out_primary.parent.mkdir(parents=True, exist_ok=True)
    primary_df.to_csv(out_primary, index=False)
    log.info(f"wrote {out_primary} ({len(primary_df)} rows)")

    if geno_rows:
        geno_df = pd.DataFrame(geno_rows)
        if "wald_p" in geno_df.columns:
            geno_df = _apply_correction(geno_df, p_col="wald_p")
        geno_df.to_csv(out_genotypic, index=False)
        log.info(f"wrote {out_genotypic} ({len(geno_df)} rows)")

    warn_path.parent.mkdir(parents=True, exist_ok=True)
    warn_path.write_text("\n".join(warnings) + ("\n" if warnings else ""))
    log.info(f"wrote {warn_path} ({len(warnings)} warnings)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
