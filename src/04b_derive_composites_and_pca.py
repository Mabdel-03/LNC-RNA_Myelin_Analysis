#!/usr/bin/env python
"""Step 04b — build dMRI composite phenotypes and ROI-level PCA scores.

Reads:
  results/intermediate/analysis_ready.tsv   (LMM superset from step 04)
  results/phenotype_manifest.csv            (with `family` column from step 02)

For each (tract, composite) the score is a sqrt(n)-rescaled weighted sum of
z-scored (or INRT) component metrics from `config.composites.definitions`.
Tracts must have ≥ `config.composites.min_components` non-missing metrics for
that composite to produce a score for an eid.

ROI PCA: per `config.roi_pca.tracts` substring, collect all metric columns for
that tract, drop high-missing cols/rows, z-score, mean-impute residual NA,
fit `sklearn.decomposition.PCA(n_components=2)`. PC1 always kept; PC2 kept
when explained-variance-ratio ≥ `config.roi_pca.pc2_min_evr`.

Writes:
  results/composite_phenotypes.csv          (eid + composite__tract cols)
  results/roi_pca_phenotypes.csv            (eid + tract__PC1[/PC2] cols)
  results/roi_pca_loadings.csv              (tract, pc, metric, loading, evr)
  results/phenotype_tiers.csv               (column_name, family, tier_source, ...)
  results/intermediate/phenotypes_extended.tsv (analysis_ready + new cols)
  results/intermediate/phenotype_manifest_extended.csv (manifest + new rows)
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
    compute_composite_score,
    compute_roi_pca,
    load_config,
    resolve_under_repo,
    setup_logger,
    standardize_eid_column,
)


def _slugify(text: str, max_len: int = 60) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_").lower()
    return s[:max_len] if max_len else s


def _tract_matches(tract_query: str, tract_disp: str) -> bool:
    """Case-insensitive substring match between a config query and a manifest tract."""
    if not tract_disp:
        return False
    return tract_query.strip().lower() in str(tract_disp).strip().lower()


def _collect_tract_metrics(manifest: pd.DataFrame,
                            wide: pd.DataFrame,
                            tract_disp: str,
                            modality: str | None = None) -> dict[str, str]:
    """For a single tract display name, return {metric_name: wide_col} for cols
    present in `wide`."""
    sub = manifest[manifest["tract_or_region_guess"].astype(str).str.strip()
                   == str(tract_disp).strip()]
    if modality:
        sub = sub[sub["modality_guess"] == modality]
    out = {}
    for _, row in sub.iterrows():
        metric = str(row["metric_guess"]).upper()
        col = str(row["column_name"])
        if col in wide.columns and metric not in out:
            out[metric] = col
    return out


def _build_composites(wide: pd.DataFrame,
                       manifest: pd.DataFrame,
                       cfg: dict,
                       logger) -> tuple[pd.DataFrame, list[dict]]:
    """Returns (composite_wide_df, list_of_manifest_rows)."""
    comp_cfg = cfg.get("composites", {})
    definitions = comp_cfg.get("definitions", {}) or {}
    standardize = comp_cfg.get("standardize", "zscore")
    min_components = int(comp_cfg.get("min_components", 2))
    rescale = bool(comp_cfg.get("rescale_by_sqrt_n", True))

    if not definitions:
        logger.info("no composite definitions in config; skipping composites")
        return pd.DataFrame({"eid": wide["eid"]}), []

    # Unique tract display names from the single-IDP manifest rows
    tracts = sorted(set(
        str(t).strip() for t in manifest["tract_or_region_guess"].fillna("")
        if t and str(t).strip()
    ))

    out = pd.DataFrame({"eid": wide["eid"].values})
    new_rows: list[dict] = []

    for comp_name, weights in definitions.items():
        weights_norm = {str(k).upper(): float(v) for k, v in weights.items()}
        n_made = 0
        for tract in tracts:
            # Composites are constructed within a single modality so we don't
            # mix TBSS and weighted-mean values. Try TBSS first, then WM.
            for modality in ("dMRI_TBSS", "dMRI_weighted_mean"):
                metric_to_col = _collect_tract_metrics(manifest, wide, tract, modality)
                if not metric_to_col:
                    continue
                # Build component DataFrame indexed by row position (matches wide.index)
                comp_df = pd.DataFrame(index=wide.index)
                for metric, col in metric_to_col.items():
                    comp_df[metric] = pd.to_numeric(wide[col], errors="coerce")
                score = compute_composite_score(
                    comp_df, weights_norm,
                    min_components=min_components,
                    standardize=standardize,
                    rescale_by_sqrt_n=rescale,
                )
                # Only keep if at least one eid has a non-NaN score
                if score.notna().sum() == 0:
                    continue
                modality_slug = "tbss" if modality == "dMRI_TBSS" else "wm"
                col_name = f"{comp_name}__{modality_slug}__{_slugify(tract, 40)}"
                out[col_name] = score.values
                new_rows.append({
                    "field_id": np.nan,
                    "instance": np.nan,
                    "array": np.nan,
                    "basket_col": "",
                    "column_name": col_name,
                    "description": f"Composite {comp_name} (modality={modality}) for tract '{tract}'",
                    "modality_guess": f"dMRI_composite_{modality_slug}",
                    "metric_guess": comp_name,
                    "tract_or_region_guess": tract,
                    "panel": "secondary",
                    "family": "secondary",
                    "transform": "none",  # composites already standardized
                    "include": True,
                    "notes": f"min_components={min_components}, standardize={standardize}",
                    "source": "composite",
                    "region_priority": 0,
                })
                n_made += 1
        logger.info(f"composite {comp_name}: built {n_made} tract×modality scores")

    return out, new_rows


def _build_roi_pca(wide: pd.DataFrame,
                    manifest: pd.DataFrame,
                    cfg: dict,
                    logger) -> tuple[pd.DataFrame, list[dict], pd.DataFrame]:
    """Returns (pc_wide_df, list_of_manifest_rows, loadings_df)."""
    rp = cfg.get("roi_pca", {})
    if not rp.get("enabled", True):
        logger.info("ROI PCA disabled in config")
        return pd.DataFrame({"eid": wide["eid"]}), [], pd.DataFrame()
    queries = list(rp.get("tracts", []))
    if not queries:
        logger.info("no roi_pca.tracts in config; skipping ROI PCA")
        return pd.DataFrame({"eid": wide["eid"]}), [], pd.DataFrame()

    max_col = float(rp.get("max_col_missing_rate", 0.5))
    max_row = float(rp.get("max_row_missing_rate", 0.3))
    pc2_min_evr = float(rp.get("pc2_min_evr", 0.10))
    seed = int(cfg.get("project", {}).get("random_seed", 12345))

    out = pd.DataFrame({"eid": wide["eid"].values})
    new_rows: list[dict] = []
    all_loadings: list[pd.DataFrame] = []

    # All distinct tracts in the manifest
    tracts_in_manifest = sorted(set(
        str(t).strip() for t in manifest["tract_or_region_guess"].fillna("")
        if t and str(t).strip()
    ))

    for q in queries:
        matched_tracts = [t for t in tracts_in_manifest if _tract_matches(q, t)]
        if not matched_tracts:
            logger.info(f"  ROI PCA query '{q}' matched 0 tracts")
            continue
        for tract in matched_tracts:
            # Combine metrics from BOTH TBSS and WM for a richer PCA (the
            # signals are usually correlated and we want to summarize all
            # information about the tract).
            metric_to_col: dict[str, str] = {}
            for modality in ("dMRI_TBSS", "dMRI_weighted_mean"):
                m = _collect_tract_metrics(manifest, wide, tract, modality)
                for k, v in m.items():
                    key = f"{k}_{'tbss' if modality == 'dMRI_TBSS' else 'wm'}"
                    metric_to_col[key] = v
            if len(metric_to_col) < 2:
                continue
            metric_df = pd.DataFrame(index=wide.index)
            for k, col in metric_to_col.items():
                metric_df[k] = pd.to_numeric(wide[col], errors="coerce")
            res = compute_roi_pca(
                metric_df,
                max_col_missing_rate=max_col,
                max_row_missing_rate=max_row,
                n_components=2,
                pc2_min_evr=pc2_min_evr,
                random_state=seed,
            )
            scores = res["scores"]
            if scores.empty or scores.shape[1] == 0:
                logger.info(f"  ROI PCA tract '{tract}': insufficient data after filters")
                continue
            for pc in scores.columns:
                col_name = f"roi__{_slugify(tract, 40)}__{pc}"
                out[col_name] = scores[pc].values
                new_rows.append({
                    "field_id": np.nan,
                    "instance": np.nan,
                    "array": np.nan,
                    "basket_col": "",
                    "column_name": col_name,
                    "description": f"ROI PCA {pc} for tract '{tract}' (combined TBSS+WM)",
                    "modality_guess": "dMRI_roi_pca",
                    "metric_guess": pc,
                    "tract_or_region_guess": tract,
                    "panel": "secondary",
                    "family": "secondary",
                    "transform": "none",
                    "include": True,
                    "notes": f"n_used={res['n_used']}, evr={res['evr'].tolist()}",
                    "source": "roi_pca",
                    "region_priority": 0,
                })
            loadings = res["loadings"].copy()
            loadings["tract"] = tract
            loadings["tract_query"] = q
            loadings = loadings.reset_index().rename(columns={"index": "metric"})
            evr_arr = res["evr"]
            evr_map = {f"PC{i + 1}": float(evr_arr[i]) for i in range(len(evr_arr))}
            for pc in [c for c in loadings.columns if c.startswith("PC")]:
                long = loadings[["tract", "tract_query", "metric", pc]].rename(
                    columns={pc: "loading"})
                long["pc"] = pc
                long["explained_variance_ratio"] = evr_map.get(pc, np.nan)
                all_loadings.append(long)
            logger.info(f"  ROI PCA tract '{tract}': {scores.shape[1]} PC(s) kept "
                        f"(evr={res['evr'].tolist()}, n_used={res['n_used']})")

    loadings_df = (pd.concat(all_loadings, ignore_index=True)
                   if all_loadings else pd.DataFrame())
    return out, new_rows, loadings_df


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
    inter_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logger("derive_composites", logs_dir / "run_log.txt")

    analysis_path = inter_dir / "analysis_ready.tsv"
    manifest_path = results_dir / "phenotype_manifest.csv"
    extended_path = inter_dir / "phenotypes_extended.tsv"
    extended_manifest_path = inter_dir / "phenotype_manifest_extended.csv"
    composite_out = results_dir / "composite_phenotypes.csv"
    roi_out = results_dir / "roi_pca_phenotypes.csv"
    loadings_out = results_dir / "roi_pca_loadings.csv"
    tiers_out = results_dir / "phenotype_tiers.csv"

    if dry:
        log.info("DRY RUN — checking expected inputs exist")
        for p in (analysis_path, manifest_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    for p in (analysis_path, manifest_path):
        if not p.is_file():
            log.error(f"required input missing: {p}")
            return 1

    log.info(f"loading analysis_ready: {analysis_path}")
    wide = pd.read_csv(analysis_path, sep="\t")
    wide = standardize_eid_column(wide)
    log.info(f"  {len(wide):,} rows × {len(wide.columns)} cols")

    log.info(f"loading manifest: {manifest_path}")
    manifest = pd.read_csv(manifest_path)
    if "family" not in manifest.columns:
        log.warning("manifest is missing `family` column — rerun step 02. Defaulting all to exploratory.")
        manifest["family"] = "exploratory"

    # ---- composites
    comp_df, comp_rows = _build_composites(wide, manifest, cfg, log)
    if comp_rows:
        composite_out.parent.mkdir(parents=True, exist_ok=True)
        comp_df.to_csv(composite_out, index=False)
        log.info(f"wrote {composite_out} ({comp_df.shape[1] - 1} composite phenotypes)")
    else:
        # Still write an empty header for downstream consistency
        pd.DataFrame({"eid": wide["eid"]}).to_csv(composite_out, index=False)
        log.info("no composites built; wrote empty composite_phenotypes.csv")

    # ---- ROI PCA
    pc_df, pc_rows, loadings = _build_roi_pca(wide, manifest, cfg, log)
    if pc_rows:
        roi_out.parent.mkdir(parents=True, exist_ok=True)
        pc_df.to_csv(roi_out, index=False)
        log.info(f"wrote {roi_out} ({pc_df.shape[1] - 1} ROI-PC phenotypes)")
        loadings.to_csv(loadings_out, index=False)
        log.info(f"wrote {loadings_out} ({len(loadings)} loadings rows)")
    else:
        pd.DataFrame({"eid": wide["eid"]}).to_csv(roi_out, index=False)
        pd.DataFrame().to_csv(loadings_out, index=False)
        log.info("no ROI-PC phenotypes built; wrote empty CSVs")

    # ---- extended manifest
    new_rows = comp_rows + pc_rows
    if new_rows:
        ext_manifest = pd.concat([manifest, pd.DataFrame(new_rows)], ignore_index=True)
    else:
        ext_manifest = manifest.copy()
    ext_manifest.to_csv(extended_manifest_path, index=False)
    log.info(f"wrote {extended_manifest_path} ({len(ext_manifest)} rows; "
             f"{len(manifest)} original + {len(new_rows)} new)")

    # ---- extended wide table: analysis_ready + composites + ROI PCs
    merged = wide.copy()
    if not comp_df.drop(columns=["eid"]).empty:
        merged = merged.merge(comp_df, on="eid", how="left")
    if not pc_df.drop(columns=["eid"]).empty:
        merged = merged.merge(pc_df, on="eid", how="left")
    extended_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(extended_path, sep="\t", index=False)
    log.info(f"wrote {extended_path} ({len(merged):,} rows × {len(merged.columns)} cols)")

    # ---- phenotype_tiers.csv
    tiers = ext_manifest[["column_name", "family", "panel", "modality_guess",
                           "metric_guess", "tract_or_region_guess", "source", "include"]].copy()
    tiers.columns = ["column_name", "family", "panel", "modality",
                      "metric", "region", "tier_source", "include"]
    tiers.to_csv(tiers_out, index=False)
    log.info(f"wrote {tiers_out}")

    # Brief summary
    by_family = ext_manifest.groupby("family").size().to_dict()
    log.info(f"phenotype family counts: {by_family}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
