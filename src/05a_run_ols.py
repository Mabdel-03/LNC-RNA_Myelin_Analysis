#!/usr/bin/env python
"""Step 05a — per-phenotype OLS additive (+ genotypic) on the OLS keep set.

Reads:
  results/intermediate/phenotype_manifest_extended.csv   (incl. composites + ROI PCs)
  results/intermediate/phenotypes_extended.tsv           (analysis_ready + composites + ROI PCs)
  results/intermediate/keep_ols.txt                       (FID IID — OLS-eligible subset)

For each phenotype with include=true:
  1. Trim |z|>outlier_z_threshold on raw scale (composites/PCs: no trim, transform="none").
  2. Apply log1p (WMH-style rows whose transform == 'log1p_then_rank_inverse').
  3. Inverse-rank-normalize (skipped when transform=="none").
  4. PRIMARY model: y ~ dosage_A + age + age2 + sex + age*sex + C(site)
                       + head_size + motion_dmri + C(array) + PC1..PCk
     SE estimator follows config.models.ols_robust_vcov (default "nonrobust").
  5. GENOTYPIC sensitivity: y ~ C(genotype, ref=GG) + same covariates.
  6. (optional) Per-phenotype diagnostics → model_diagnostics_summary.csv.

Multiple testing:
  - Per-panel Bonferroni + BH-FDR (back-compat — `bonferroni`, `fdr_bh` columns).
  - Per-family Bonferroni + BH-FDR (new — `family_bonferroni`, `family_fdr_bh`).

Outputs (split by phenotype source so the existing PheWAS CSV schema is preserved):
  results/association_results_primary.csv     — single IDPs (flat PheWAS, back-compat)
  results/association_results_composites.csv  — composite phenotypes
  results/association_results_roi_pca.csv     — ROI-PCA phenotypes
  results/association_results_genotypic.csv   — genotypic sensitivity (single IDPs)
  results/multiple_testing_summary.csv        — per-family test counts and hit summary
  results/model_diagnostics_summary.csv       — optional residual diagnostics
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
    bh_fdr_by_family,
    bonferroni_by_family,
    compute_ols_diagnostics,
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
    cont = []
    cat = []
    for c in ("age", "age2", "age_sex", "head_size", "motion_dmri"):
        if c in merged.columns:
            cont.append(c)
    pc_cols = [c for c in merged.columns
               if c.startswith("PC") and c[2:].isdigit()
               and int(c[2:]) <= int(cfg["covariates"].get("genetic_pc_count", 20))]
    cont.extend(sorted(pc_cols, key=lambda c: int(c[2:])))
    for c in ("sex", "site", "array"):
        if c in merged.columns:
            cat.append(c)
    X = make_design_matrix(merged, cont_cols=cont, cat_cols=cat,
                            add_intercept=True, drop_first=True)
    return X, cont, cat


def _apply_per_panel_correction(df: pd.DataFrame, p_col: str = "p") -> pd.DataFrame:
    """Back-compat: per-panel Bonferroni + BH-FDR (writes `bonferroni`, `fdr_bh`)."""
    out = df.copy()
    out["bonferroni"] = np.nan
    out["fdr_bh"] = np.nan
    if "panel" not in out.columns:
        return out
    for panel in out["panel"].unique():
        mask = (out["panel"] == panel) & out[p_col].notna()
        if not mask.any():
            continue
        n = int(mask.sum())
        out.loc[mask, "bonferroni"] = np.minimum(out.loc[mask, p_col] * n, 1.0)
        out.loc[mask, "fdr_bh"] = bh_fdr(out.loc[mask, p_col].values)
    return out


def _apply_per_family_correction(df: pd.DataFrame, p_col: str = "p") -> pd.DataFrame:
    out = df.copy()
    if "family" not in out.columns:
        return out
    out["family_bonferroni"] = bonferroni_by_family(out, family_col="family", p_col=p_col)
    out["family_fdr_bh"] = bh_fdr_by_family(out, family_col="family", p_col=p_col)
    return out


def _filter_to_keep(df: pd.DataFrame, keep_path: Path, log) -> pd.DataFrame:
    if not keep_path.is_file():
        log.warning(f"keep file not found: {keep_path}; using all rows")
        return df
    keep_df = pd.read_csv(keep_path, sep=r"\s+", header=None, names=["FID", "IID"])
    keep_set = set(int(x) for x in keep_df["IID"])
    before = len(df)
    out = df[df["eid"].astype("Int64").isin(keep_set)].copy()
    log.info(f"OLS keep filter: {len(out):,} of {before:,} rows retained")
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
    log = setup_logger("assoc_ols", logs_dir / "run_log.txt")

    # Prefer extended manifest/table if present; fall back to originals for back-compat
    manifest_path = inter_dir / "phenotype_manifest_extended.csv"
    if not manifest_path.is_file():
        manifest_path = results_dir / "phenotype_manifest.csv"
    analysis_path = inter_dir / "phenotypes_extended.tsv"
    if not analysis_path.is_file():
        analysis_path = inter_dir / "analysis_ready.tsv"
    keep_ols_path = inter_dir / "keep_ols.txt"

    out_primary = results_dir / "association_results_primary.csv"
    out_composites = results_dir / "association_results_composites.csv"
    out_roi = results_dir / "association_results_roi_pca.csv"
    out_genotypic = results_dir / "association_results_genotypic.csv"
    out_mt = results_dir / "multiple_testing_summary.csv"
    out_diag = results_dir / "model_diagnostics_summary.csv"
    warn_path = logs_dir / "model_warnings.txt"

    if dry:
        log.info("DRY RUN — would fit OLS + genotypic models per phenotype")
        for p in (manifest_path, analysis_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    if not manifest_path.is_file():
        log.error(f"manifest not found: {manifest_path}")
        return 1
    if not analysis_path.is_file():
        log.error(f"phenotype table not found: {analysis_path}")
        return 1

    manifest = pd.read_csv(manifest_path)
    if "family" not in manifest.columns:
        log.warning("manifest missing `family` column — defaulting all to 'exploratory'")
        manifest["family"] = "exploratory"
    if "source" not in manifest.columns:
        manifest["source"] = "basket"

    merged = pd.read_csv(analysis_path, sep="\t")
    log.info(f"manifest: {len(manifest)} rows, phenotype table: {len(merged):,} rows")

    # Apply OLS keep filter (drop relatives)
    if "eid" not in merged.columns:
        log.error("phenotype table missing 'eid' column")
        return 2
    merged = _filter_to_keep(merged, keep_ols_path, log)

    if "include" in manifest.columns:
        manifest = manifest[manifest["include"].astype(str).str.lower().isin(("true", "1"))].copy()
    log.info(f"manifest after include filter: {len(manifest)} rows")

    if cfg.get("project", {}).get("test_mode", False):
        cap = int(cfg["project"].get("max_test_phenotypes", 3))
        # Prioritize secondary (composites + ROI PCs) and controls in test mode
        # so the demo exercises the new family paths, not just single IDPs.
        if "family" in manifest.columns:
            priority = {"primary": 0, "secondary": 1, "controls": 2, "exploratory": 3}
            manifest = manifest.assign(
                _prio=manifest["family"].map(priority).fillna(9)
            ).sort_values("_prio").drop(columns="_prio")
        manifest = manifest.head(cap)
        log.info(f"TEST MODE — capped to first {cap} phenotypes (family priority order)")

    X_covar, cont_cols, cat_cols = _covariate_design(merged, cfg)
    log.info(f"covariate design: {len(cont_cols)} continuous, "
             f"{len(cat_cols)} categorical, total ncols={X_covar.shape[1]}")

    vcov = str(cfg.get("models", {}).get("ols_robust_vcov", "nonrobust"))
    compute_diag = bool(cfg.get("models", {}).get("compute_diagnostics", True))
    z_thresh = float(cfg.get("models", {}).get("outlier_z_threshold", 6))

    warnings: list[str] = []
    primary_rows: list[dict] = []
    geno_rows: list[dict] = []
    diag_rows: list[dict] = []

    do_genotypic = bool(cfg["models"].get("run_genotypic_sensitivity", True)) \
        and ("genotype" in merged.columns)
    if not do_genotypic:
        log.warning("genotypic-model sensitivity disabled "
                    "(config flag false OR no hard-call genotype column).")

    for _, row in manifest.iterrows():
        idp = row["column_name"]
        if idp not in merged.columns:
            warnings.append(f"{idp}: column not in analysis table; skipped")
            continue
        y_raw = pd.to_numeric(merged[idp], errors="coerce")
        transform = str(row.get("transform", "rank_inverse_normal"))

        if transform == "none":
            # composites + ROI PCs are already standardized; no further preproc
            y_inrt = y_raw
        else:
            y_trim = trim_outliers_z(y_raw, z_thresh=z_thresh)
            if transform.startswith("log1p"):
                y_trim = pd.Series(np.where(y_trim >= 0, np.log1p(y_trim), np.nan),
                                    index=y_trim.index)
            y_inrt = inverse_rank_normalize(y_trim)

        # PRIMARY additive model
        want_diag = compute_diag and (row.get("source") in ("composite", "roi_pca"))
        res = fit_ols_model(y_inrt, merged["dosage_A"], X_covar,
                            exposure_name="dosage_A",
                            vcov_type=vcov,
                            return_model=want_diag)
        rec = {
            "column_name": idp,
            "field_id": row.get("field_id", ""),
            "family": row.get("family", ""),
            "panel": row.get("panel", ""),
            "modality": row.get("modality_guess", ""),
            "metric": row.get("metric_guess", ""),
            "region": row.get("tract_or_region_guess", ""),
            "source": row.get("source", ""),
            "effect_allele": cfg["variant"]["effect_allele"],
            "other_allele": cfg["variant"]["other_allele"],
            "n": res.get("n"),
            "beta_per_A": res.get("beta"),
            "se": res.get("se"),
            "t": res.get("t"),
            "p": res.get("p"),
            "vcov_type": res.get("vcov_type", vcov),
            "aa_vs_gg_additive_2beta": (res.get("beta") * 2.0) if res.get("status") == "ok" else None,
            "status": res.get("status"),
        }
        primary_rows.append(rec)
        if res.get("status") != "ok":
            warnings.append(f"{idp}: primary model status={res.get('status')}")
        if want_diag and res.get("status") == "ok":
            d = compute_ols_diagnostics(res)
            d["column_name"] = idp
            d["family"] = row.get("family", "")
            diag_rows.append(d)

        # GENOTYPIC sensitivity (single IDPs only — meaningless for composites/PCs)
        if do_genotypic and row.get("source", "") not in ("composite", "roi_pca"):
            gres = fit_genotypic_model(y_inrt, merged["genotype"], X_covar,
                                        reference=cfg["variant"]["other_allele"] * 2)
            grec = {
                "column_name": idp,
                "family": row.get("family", ""),
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

    all_df = pd.DataFrame(primary_rows)
    # Per-panel correction (back-compat) within each phenotype source group
    all_df = _apply_per_panel_correction(all_df, p_col="p")
    # Per-family correction across the entire results set (the headline)
    all_df = _apply_per_family_correction(all_df, p_col="p")
    all_df = all_df.sort_values(["family", "p"], na_position="last")

    # ---- write split outputs by phenotype source for back-compat
    out_primary.parent.mkdir(parents=True, exist_ok=True)
    src = all_df["source"].fillna("basket")
    single_df = all_df[~src.isin(["composite", "roi_pca"])].copy()
    comp_df = all_df[src.eq("composite")].copy()
    pca_df = all_df[src.eq("roi_pca")].copy()
    single_df.to_csv(out_primary, index=False)
    log.info(f"wrote {out_primary} ({len(single_df)} rows — single-IDP back-compat)")
    if not comp_df.empty:
        comp_df.to_csv(out_composites, index=False)
        log.info(f"wrote {out_composites} ({len(comp_df)} rows)")
    if not pca_df.empty:
        pca_df.to_csv(out_roi, index=False)
        log.info(f"wrote {out_roi} ({len(pca_df)} rows)")

    if geno_rows:
        geno_df = pd.DataFrame(geno_rows)
        if "wald_p" in geno_df.columns:
            geno_df = _apply_per_panel_correction(geno_df, p_col="wald_p")
            geno_df = _apply_per_family_correction(geno_df, p_col="wald_p")
        geno_df.to_csv(out_genotypic, index=False)
        log.info(f"wrote {out_genotypic} ({len(geno_df)} rows)")

    # ---- multiple_testing_summary.csv
    mt_rows = []
    if "family" in all_df.columns:
        for fam, sub in all_df.groupby("family", dropna=False):
            n = int(sub["p"].notna().sum())
            if n == 0:
                mt_rows.append({"family": fam, "n_tests": 0,
                                "alpha_bonf": np.nan,
                                "n_sig_bonf": 0, "n_sig_fdr_05": 0, "n_sig_fdr_10": 0,
                                "min_p": np.nan, "q_at_top_hit": np.nan})
                continue
            min_p = float(sub["p"].min())
            n_sig_bonf = int((sub["p"] * n <= 0.05).sum())
            n_sig_fdr05 = int((sub.get("family_fdr_bh", pd.Series(np.nan, index=sub.index)) <= 0.05).sum())
            n_sig_fdr10 = int((sub.get("family_fdr_bh", pd.Series(np.nan, index=sub.index)) <= 0.10).sum())
            top_q = float(sub.sort_values("p").iloc[0].get("family_fdr_bh", np.nan)) \
                if "family_fdr_bh" in sub.columns else np.nan
            mt_rows.append({
                "family": fam, "n_tests": n,
                "alpha_bonf": 0.05 / n,
                "n_sig_bonf": n_sig_bonf,
                "n_sig_fdr_05": n_sig_fdr05,
                "n_sig_fdr_10": n_sig_fdr10,
                "min_p": min_p,
                "q_at_top_hit": top_q,
            })
    pd.DataFrame(mt_rows).to_csv(out_mt, index=False)
    log.info(f"wrote {out_mt} ({len(mt_rows)} family rows)")

    # ---- model_diagnostics_summary.csv
    if diag_rows:
        diag_df = pd.DataFrame(diag_rows)
        diag_df.to_csv(out_diag, index=False)
        log.info(f"wrote {out_diag} ({len(diag_df)} rows)")

    warn_path.parent.mkdir(parents=True, exist_ok=True)
    warn_path.write_text("\n".join(warnings) + ("\n" if warnings else ""))
    log.info(f"wrote {warn_path} ({len(warnings)} warnings)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
