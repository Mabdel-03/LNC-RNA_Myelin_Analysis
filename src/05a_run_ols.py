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

Outputs (split by phenotype source so the existing PheWAS CSV schema is preserved):
  results/association_results_primary.csv     — single IDPs (flat PheWAS, back-compat)
  results/association_results_composites.csv  — composite phenotypes
  results/association_results_roi_pca.csv     — ROI-PCA phenotypes
  results/association_results_genotypic.csv   — genotypic sensitivity (single IDPs)
  results/association_results_model_matrix.csv — additive/dominant/recessive model variants
  results/association_results_pairwise_genotype.csv — AA/AG/GG pairwise contrasts
  results/raw_p_testing_summary.csv           — per-family raw-p test counts and thresholds
  results/model_diagnostics_summary.csv       — optional residual diagnostics
  logs/model_warnings.txt
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import multiprocessing as mp
import os
import sys
from pathlib import Path

# Keep each worker's BLAS footprint small. Slurm wrappers set these explicitly;
# the defaults here protect ad hoc worker launches from oversubscribing cores.
for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_thread_var, os.environ.get("OLS_BLAS_THREADS", "1"))

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    compute_ols_diagnostics,
    fit_genotype_pairwise_model,
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
    base_cont = ("age", "age2", "age_sex", "head_size", "motion_dmri",
                 "t1_motion", "bmi")
    for c in base_cont:
        if c in merged.columns:
            cont.append(c)
    cont.extend([c for c in merged.columns if c.startswith((
        "scanner_pos_", "dmri_qc_", "t1_qc_", "modality_discrepancy_", "wmh_"
    ))])
    pc_cols = [c for c in merged.columns
               if c.startswith("PC") and c[2:].isdigit()
               and int(c[2:]) <= int(cfg["covariates"].get("genetic_pc_count", 20))]
    cont.extend(sorted(pc_cols, key=lambda c: int(c[2:])))
    cont = list(dict.fromkeys(cont))
    for c in ("sex", "site", "array", "genotype_batch", "smoking_status",
              "diabetes", "hypertension_6150", "hypertension_6177",
              "hypertension_self_report"):
        if c in merged.columns:
            cat.append(c)
    cat.extend([c for c in merged.columns if c.startswith("protocol_")])
    cat = list(dict.fromkeys(cat))
    X = make_design_matrix(merged, cont_cols=cont, cat_cols=cat,
                            add_intercept=True, drop_first=True)
    return X, cont, cat


def _add_raw_p_annotations(df: pd.DataFrame,
                           p_col: str = "p",
                           family_col: str = "family",
                           thresholds: list[float] | None = None) -> pd.DataFrame:
    out = df.copy()
    if p_col not in out.columns:
        return out
    out["raw_p"] = pd.to_numeric(out[p_col], errors="coerce")
    out["minus_log10_raw_p"] = -np.log10(out["raw_p"].clip(lower=np.finfo(float).tiny))
    if family_col in out.columns:
        out["n_tests_in_family"] = out.groupby(family_col)["raw_p"].transform(lambda s: int(s.notna().sum()))
        out["raw_p_rank_within_family"] = out.groupby(family_col)["raw_p"].rank(method="min", na_option="bottom")
    for thr in thresholds or []:
        label = str(thr).replace(".", "_")
        out[f"raw_p_lt_{label}"] = out["raw_p"] < float(thr)
    return out


def _raw_p_summary(df: pd.DataFrame,
                   family_col: str = "family",
                   p_col: str = "raw_p",
                   thresholds: list[float] | None = None) -> pd.DataFrame:
    rows = []
    if family_col not in df.columns or p_col not in df.columns:
        return pd.DataFrame(rows)
    for fam, sub in df.groupby(family_col, dropna=False):
        p = pd.to_numeric(sub[p_col], errors="coerce")
        n = int(p.notna().sum())
        rec = {"family": fam, "n_tests": n, "min_raw_p": float(p.min()) if n else np.nan}
        for thr in thresholds or []:
            rec[f"n_raw_p_lt_{thr:g}"] = int((p < float(thr)).sum())
        rows.append(rec)
    return pd.DataFrame(rows).sort_values("min_raw_p", na_position="last")


def _metadata_from_manifest_row(row: pd.Series, cfg: dict) -> dict:
    return {
        "column_name": row["column_name"],
        "field_id": row.get("field_id", ""),
        "field_id_base": row.get("field_id_base", ""),
        "family": row.get("family", ""),
        "analysis_tier": row.get("analysis_tier", ""),
        "confirmatory_role": row.get("confirmatory_role", ""),
        "panel": row.get("panel", ""),
        "modality": row.get("modality_guess", ""),
        "metric": row.get("metric_guess", ""),
        "region": row.get("tract_or_region_guess", ""),
        "source": row.get("source", ""),
        "effect_allele": cfg["variant"]["effect_allele"],
        "other_allele": cfg["variant"]["other_allele"],
    }


def _binary_model_rows(y: pd.Series,
                       genotype: pd.Series,
                       x_covar: pd.DataFrame,
                       meta: dict,
                       cfg: dict,
                       vcov_types: list[str]) -> list[dict]:
    rows = []
    ea = cfg["variant"]["effect_allele"]
    oa = cfg["variant"]["other_allele"]
    aa = ea + ea
    ag = "".join(sorted([ea, oa]))
    gg = oa + oa
    g = genotype.astype(object)
    encodings = []
    if cfg["models"].get("run_dominant_sensitivity", False):
        encodings.append(("dominant_A", g.map({gg: 0.0, ag: 1.0, aa: 1.0}), f"{aa}/{ag}_vs_{gg}"))
    if cfg["models"].get("run_recessive_sensitivity", False):
        encodings.append(("recessive_A", g.map({gg: 0.0, ag: 0.0, aa: 1.0}), f"{aa}_vs_{ag}/{gg}"))
    for model_name, exposure, contrast in encodings:
        for vcov in vcov_types:
            res = fit_ols_model(y, exposure, x_covar,
                                exposure_name=model_name,
                                vcov_type=vcov)
            rows.append({
                **meta,
                "model": model_name,
                "contrast": contrast,
                "n": res.get("n"),
                "beta": res.get("beta"),
                "se": res.get("se"),
                "stat": res.get("t"),
                "p": res.get("p"),
                "vcov_type": res.get("vcov_type", vcov),
                "status": res.get("status"),
            })
    return rows


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


_OLS_WORKER: dict = {}
_THREADPOOL_LIMITER = None


def _resolve_parallel_jobs(cfg: dict, n_phenotypes: int, log) -> int:
    model_cfg = cfg.get("models", {})
    raw = model_cfg.get("ols_parallel_jobs", model_cfg.get("ols_n_jobs", "auto"))
    if str(raw).lower() == "auto":
        raw = (
            os.environ.get("OLS_N_JOBS")
            or os.environ.get("SLURM_CPUS_PER_TASK")
            or os.environ.get("SLURM_CPUS_ON_NODE")
            or os.environ.get("SLURM_NTASKS")
            or "1"
        )
    try:
        n_jobs = int(raw)
    except (TypeError, ValueError):
        log.warning(f"invalid models.ols_parallel_jobs={raw!r}; using 1")
        n_jobs = 1
    n_jobs = max(1, min(n_jobs, max(1, n_phenotypes)))
    if n_jobs > 1 and not hasattr(os, "fork"):
        log.warning("phenotype-level OLS parallelism requires fork; using 1 worker")
        return 1
    return n_jobs


def _init_ols_worker(merged: pd.DataFrame,
                     x_covar: pd.DataFrame,
                     cfg: dict,
                     vcov_types: list[str],
                     default_vcov: str,
                     compute_diag: bool,
                     z_thresh: float,
                     do_genotypic: bool) -> None:
    global _OLS_WORKER, _THREADPOOL_LIMITER
    _OLS_WORKER = {
        "merged": merged,
        "x_covar": x_covar,
        "cfg": cfg,
        "vcov_types": vcov_types,
        "default_vcov": default_vcov,
        "compute_diag": compute_diag,
        "z_thresh": z_thresh,
        "do_genotypic": do_genotypic,
    }
    try:
        from threadpoolctl import threadpool_limits

        blas_threads = int(os.environ.get("OLS_BLAS_THREADS", "1"))
        _THREADPOOL_LIMITER = threadpool_limits(limits=max(1, blas_threads))
        _THREADPOOL_LIMITER.__enter__()
    except Exception:
        _THREADPOOL_LIMITER = None


def _analyze_phenotype_record(row: dict) -> dict:
    merged: pd.DataFrame = _OLS_WORKER["merged"]
    x_covar: pd.DataFrame = _OLS_WORKER["x_covar"]
    cfg: dict = _OLS_WORKER["cfg"]
    vcov_types: list[str] = _OLS_WORKER["vcov_types"]
    default_vcov: str = _OLS_WORKER["default_vcov"]
    compute_diag: bool = _OLS_WORKER["compute_diag"]
    z_thresh: float = _OLS_WORKER["z_thresh"]
    do_genotypic: bool = _OLS_WORKER["do_genotypic"]

    warnings: list[str] = []
    primary_rows: list[dict] = []
    matrix_rows: list[dict] = []
    geno_rows: list[dict] = []
    pairwise_rows: list[dict] = []
    diag_rows: list[dict] = []

    idp = row["column_name"]
    if idp not in merged.columns:
        warnings.append(f"{idp}: column not in analysis table; skipped")
        return {
            "primary_rows": primary_rows,
            "matrix_rows": matrix_rows,
            "geno_rows": geno_rows,
            "pairwise_rows": pairwise_rows,
            "diag_rows": diag_rows,
            "warnings": warnings,
        }

    y_raw = pd.to_numeric(merged[idp], errors="coerce")
    transform = str(row.get("transform", "rank_inverse_normal"))

    if transform == "none":
        y_inrt = y_raw
    else:
        y_trim = trim_outliers_z(y_raw, z_thresh=z_thresh)
        if transform.startswith("log1p"):
            y_trim = pd.Series(np.where(y_trim >= 0, np.log1p(y_trim), np.nan),
                               index=y_trim.index)
        y_inrt = inverse_rank_normalize(y_trim)

    meta = _metadata_from_manifest_row(row, cfg)
    want_diag = compute_diag and (row.get("source") in ("composite", "roi_pca"))
    default_rec = None

    for vcov in vcov_types:
        res = fit_ols_model(y_inrt, merged["dosage_A"], x_covar,
                            exposure_name="dosage_A",
                            vcov_type=vcov,
                            return_model=(want_diag and vcov == default_vcov))
        rec = {
            **meta,
            "n": res.get("n"),
            "beta_per_A": res.get("beta"),
            "se": res.get("se"),
            "t": res.get("t"),
            "p": res.get("p"),
            "vcov_type": res.get("vcov_type", vcov),
            "aa_vs_gg_additive_2beta": (res.get("beta") * 2.0) if res.get("status") == "ok" else None,
            "status": res.get("status"),
        }
        if vcov == default_vcov:
            default_rec = rec
        matrix_rows.append({
            **meta,
            "model": "additive_dosage",
            "contrast": "per_A_allele",
            "n": res.get("n"),
            "beta": res.get("beta"),
            "se": res.get("se"),
            "stat": res.get("t"),
            "p": res.get("p"),
            "vcov_type": res.get("vcov_type", vcov),
            "status": res.get("status"),
        })
        if res.get("status") != "ok":
            warnings.append(f"{idp}: additive {vcov} model status={res.get('status')}")
        if want_diag and vcov == default_vcov and res.get("status") == "ok":
            d = compute_ols_diagnostics(res)
            d["column_name"] = idp
            d["family"] = row.get("family", "")
            d["analysis_tier"] = row.get("analysis_tier", "")
            diag_rows.append(d)
    if default_rec is not None:
        primary_rows.append(default_rec)

    if do_genotypic:
        ea = cfg["variant"]["effect_allele"]
        oa = cfg["variant"]["other_allele"]
        gg = oa + oa
        ag = "".join(sorted([ea, oa]))
        aa = ea + ea
        gres = fit_genotypic_model(y_inrt, merged["genotype"], x_covar,
                                   reference=gg)
        grec = {
            **meta,
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

        if cfg["models"].get("run_pairwise_genotype", True):
            pres = fit_genotype_pairwise_model(
                y_inrt, merged["genotype"], x_covar,
                levels=[aa, ag, gg],
                vcov_type=default_vcov,
            )
            glob = pres.get("global", {})
            matrix_rows.append({
                **meta,
                "model": "genotypic_global",
                "contrast": f"{aa}/{ag}/{gg}_2df",
                "n": glob.get("n", pres.get("n")),
                "beta": np.nan,
                "se": np.nan,
                "stat": glob.get("wald_chi2"),
                "p": glob.get("raw_p"),
                "vcov_type": default_vcov,
                "status": glob.get("status", pres.get("status")),
            })
            for prow in pres.get("pairwise", []):
                out = {
                    **meta,
                    "model": "genotype_pairwise",
                    "contrast": prow.get("contrast"),
                    "level_a": prow.get("level_a"),
                    "level_b": prow.get("level_b"),
                    "n": prow.get("n"),
                    "beta": prow.get("beta"),
                    "se": prow.get("se"),
                    "t": prow.get("t"),
                    "p": prow.get("raw_p"),
                    "ci95_low": prow.get("ci95_low"),
                    "ci95_high": prow.get("ci95_high"),
                    "vcov_type": prow.get("vcov_type"),
                    "status": prow.get("status"),
                }
                pairwise_rows.append(out)
                matrix_rows.append({
                    **meta,
                    "model": "genotype_pairwise",
                    "contrast": prow.get("contrast"),
                    "n": prow.get("n"),
                    "beta": prow.get("beta"),
                    "se": prow.get("se"),
                    "stat": prow.get("t"),
                    "p": prow.get("raw_p"),
                    "vcov_type": prow.get("vcov_type"),
                    "status": prow.get("status"),
                })
        matrix_rows.extend(_binary_model_rows(
            y_inrt, merged["genotype"], x_covar, meta, cfg, vcov_types
        ))

    return {
        "primary_rows": primary_rows,
        "matrix_rows": matrix_rows,
        "geno_rows": geno_rows,
        "pairwise_rows": pairwise_rows,
        "diag_rows": diag_rows,
        "warnings": warnings,
    }


def _extend_results(bundle: dict,
                    primary_rows: list[dict],
                    matrix_rows: list[dict],
                    geno_rows: list[dict],
                    pairwise_rows: list[dict],
                    diag_rows: list[dict],
                    warnings: list[str]) -> None:
    primary_rows.extend(bundle.get("primary_rows", []))
    matrix_rows.extend(bundle.get("matrix_rows", []))
    geno_rows.extend(bundle.get("geno_rows", []))
    pairwise_rows.extend(bundle.get("pairwise_rows", []))
    diag_rows.extend(bundle.get("diag_rows", []))
    warnings.extend(bundle.get("warnings", []))


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
    out_model_matrix = results_dir / "association_results_model_matrix.csv"
    out_pairwise = results_dir / "association_results_pairwise_genotype.csv"
    out_raw_summary = results_dir / "raw_p_testing_summary.csv"
    out_mt = results_dir / "multiple_testing_summary.csv"  # compatibility name; raw-p summary only
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

    model_cfg = cfg.get("models", {})
    os.environ["OLS_BLAS_THREADS"] = str(model_cfg.get(
        "ols_blas_threads", os.environ.get("OLS_BLAS_THREADS", "1")
    ))
    vcov_types = list(model_cfg.get("ols_vcov_types", [])) or [str(model_cfg.get("ols_robust_vcov", "nonrobust"))]
    vcov_types = [str(v) for v in vcov_types]
    default_vcov = vcov_types[0]
    compute_diag = bool(cfg.get("models", {}).get("compute_diagnostics", True))
    z_thresh = float(cfg.get("models", {}).get("outlier_z_threshold", 6))
    raw_thresholds = [float(x) for x in model_cfg.get("raw_p_thresholds", [0.05, 0.01, 0.001])]

    warnings: list[str] = []
    primary_rows: list[dict] = []
    matrix_rows: list[dict] = []
    geno_rows: list[dict] = []
    pairwise_rows: list[dict] = []
    diag_rows: list[dict] = []

    do_genotypic = bool(cfg["models"].get("run_genotypic_sensitivity", True)) \
        and ("genotype" in merged.columns)
    if not do_genotypic:
        log.warning("genotypic-model sensitivity disabled "
                    "(config flag false OR no hard-call genotype column).")

    manifest_records = manifest.to_dict("records")
    n_jobs = _resolve_parallel_jobs(cfg, len(manifest_records), log)
    blas_threads = int(os.environ.get("OLS_BLAS_THREADS", "1"))
    log.info(f"OLS phenotype parallelism: n_jobs={n_jobs}, blas_threads_per_worker={blas_threads}")

    if n_jobs == 1:
        _init_ols_worker(merged, X_covar, cfg, vcov_types, default_vcov,
                         compute_diag, z_thresh, do_genotypic)
        for i, record in enumerate(manifest_records, start=1):
            bundle = _analyze_phenotype_record(record)
            _extend_results(bundle, primary_rows, matrix_rows, geno_rows,
                            pairwise_rows, diag_rows, warnings)
            if i % 25 == 0 or i == len(manifest_records):
                log.info(f"OLS phenotypes completed: {i}/{len(manifest_records)}")
    else:
        ctx = mp.get_context("fork")
        with cf.ProcessPoolExecutor(
            max_workers=n_jobs,
            mp_context=ctx,
            initializer=_init_ols_worker,
            initargs=(merged, X_covar, cfg, vcov_types, default_vcov,
                      compute_diag, z_thresh, do_genotypic),
        ) as executor:
            futures = [executor.submit(_analyze_phenotype_record, record)
                       for record in manifest_records]
            for i, fut in enumerate(cf.as_completed(futures), start=1):
                bundle = fut.result()
                _extend_results(bundle, primary_rows, matrix_rows, geno_rows,
                                pairwise_rows, diag_rows, warnings)
                if i % 25 == 0 or i == len(futures):
                    log.info(f"OLS phenotypes completed: {i}/{len(futures)}")

    all_df = pd.DataFrame(primary_rows)
    all_df = _add_raw_p_annotations(all_df, p_col="p", thresholds=raw_thresholds)
    all_df = all_df.sort_values(["family", "raw_p"], na_position="last")

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
            geno_df = _add_raw_p_annotations(geno_df, p_col="wald_p", thresholds=raw_thresholds)
        geno_df.to_csv(out_genotypic, index=False)
        log.info(f"wrote {out_genotypic} ({len(geno_df)} rows)")

    if pairwise_rows:
        pairwise_df = pd.DataFrame(pairwise_rows)
        pairwise_df = _add_raw_p_annotations(pairwise_df, p_col="p", thresholds=raw_thresholds)
        pairwise_df.sort_values(["family", "raw_p"], na_position="last").to_csv(out_pairwise, index=False)
        log.info(f"wrote {out_pairwise} ({len(pairwise_df)} rows)")

    if matrix_rows:
        matrix_df = pd.DataFrame(matrix_rows)
        matrix_df = _add_raw_p_annotations(matrix_df, p_col="p", thresholds=raw_thresholds)
        matrix_df.sort_values(["family", "model", "raw_p"], na_position="last").to_csv(out_model_matrix, index=False)
        log.info(f"wrote {out_model_matrix} ({len(matrix_df)} rows)")

    # ---- raw-p family summary (compat copy at multiple_testing_summary.csv)
    raw_summary = _raw_p_summary(all_df, thresholds=raw_thresholds)
    raw_summary.to_csv(out_raw_summary, index=False)
    raw_summary.to_csv(out_mt, index=False)
    log.info(f"wrote {out_raw_summary} and raw-p compatibility copy {out_mt} "
             f"({len(raw_summary)} family rows)")

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
