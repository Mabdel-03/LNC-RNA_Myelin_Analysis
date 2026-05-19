#!/usr/bin/env python
"""Step 05d — ingest REGENIE / BOLT-LMM results into the OLS-shaped CSV schema.

For each LMM result on disk, extract the rs2546890 row (target variant) and
emit a per-source association table:

  results/association_results_primary_lmm.csv     — single IDPs
  results/association_results_composites_lmm.csv  — composites
  results/association_results_roi_pca_lmm.csv     — ROI-PCA scores
  results/association_results_controls_lmm.csv    — WMH/QC outcomes
  results/multiple_testing_summary_lmm.csv        — per-family Bonferroni/BH

The output schema mirrors the OLS results so the report can render LMM and
OLS side-by-side: column_name, family, panel, modality, metric, region, n,
beta_per_A, se, chisq, p, log10p, engine, info, eaf, status.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    bh_fdr_by_family,
    bonferroni_by_family,
    load_config,
    parse_bolt_stats,
    parse_regenie_output,
    resolve_under_repo,
    setup_logger,
)


def _manifest_lookup(manifest: pd.DataFrame, column_name: str) -> dict:
    row = manifest[manifest["column_name"] == column_name]
    if row.empty:
        return {"family": "unknown", "panel": "unknown", "modality": "",
                "metric": "", "region": "", "source": ""}
    r = row.iloc[0]
    return {
        "family": r.get("family", "unknown"),
        "panel": r.get("panel", "unknown"),
        "modality": r.get("modality_guess", ""),
        "metric": r.get("metric_guess", ""),
        "region": r.get("tract_or_region_guess", ""),
        "source": r.get("source", ""),
    }


def _ingest_regenie(lmm_dir: Path, target_id: str, manifest: pd.DataFrame,
                     log) -> pd.DataFrame:
    """Walk results/lmm/regenie_step2_*.regenie* and pull rs2546890 effect per pheno."""
    rows = []
    files = sorted(list(lmm_dir.glob("regenie_step2_*.regenie.gz")) +
                   list(lmm_dir.glob("regenie_step2_*.regenie")))
    if not files:
        log.info("no REGENIE step-2 output files found")
        return pd.DataFrame()
    for f in files:
        # Filename convention: regenie_step2_<pheno>.regenie[.gz]
        stem = f.name
        if stem.endswith(".gz"):
            stem = stem[:-3]
        if stem.endswith(".regenie"):
            stem = stem[:-len(".regenie")]
        pheno = stem.replace("regenie_step2_", "", 1)
        try:
            df = parse_regenie_output(f, target_variant_id=target_id)
        except Exception as exc:  # noqa: BLE001
            log.warning(f"  parse error for {f.name}: {exc}")
            continue
        if df.empty:
            log.warning(f"  {f.name}: target variant {target_id} not found")
            continue
        r = df.iloc[0]
        meta = _manifest_lookup(manifest, pheno)
        rows.append({
            "column_name": pheno,
            "family": meta["family"],
            "panel": meta["panel"],
            "modality": meta["modality"],
            "metric": meta["metric"],
            "region": meta["region"],
            "source": meta["source"],
            "n": r.get("n", np.nan),
            "beta_per_A": r.get("beta", np.nan),
            "se": r.get("se", np.nan),
            "chisq": r.get("chisq", np.nan),
            "p": r.get("p", np.nan),
            "log10p": r.get("log10p", np.nan),
            "engine": "regenie",
            "info": r.get("info", np.nan),
            "eaf": r.get("eaf", np.nan),
            "status": "ok",
        })
    return pd.DataFrame(rows)


def _ingest_bolt(lmm_dir: Path, target_id: str, target_chrpos: tuple[int, int] | None,
                  target_rsid: str | None, manifest: pd.DataFrame, log) -> pd.DataFrame:
    """Walk results/lmm/bolt_<pheno>.stats and pull rs2546890 effect per pheno."""
    rows = []
    files = sorted(lmm_dir.glob("bolt_*.stats"))
    if not files:
        log.info("no BOLT .stats files found")
        return pd.DataFrame()
    for f in files:
        pheno = f.name[len("bolt_"):-len(".stats")]
        try:
            df = parse_bolt_stats(f, target_variant_id=target_id,
                                    target_rsid=target_rsid,
                                    target_chrpos=target_chrpos)
        except Exception as exc:  # noqa: BLE001
            log.warning(f"  parse error for {f.name}: {exc}")
            continue
        if df.empty:
            log.warning(f"  {f.name}: target variant not found")
            continue
        r = df.iloc[0]
        meta = _manifest_lookup(manifest, pheno)
        rows.append({
            "column_name": pheno,
            "family": meta["family"],
            "panel": meta["panel"],
            "modality": meta["modality"],
            "metric": meta["metric"],
            "region": meta["region"],
            "source": meta["source"],
            "n": np.nan,
            "beta_per_A": r.get("beta", np.nan),
            "se": r.get("se", np.nan),
            "chisq": r.get("chisq", np.nan),
            "p": r.get("p", np.nan),
            "log10p": r.get("log10p", np.nan),
            "engine": "bolt",
            "info": r.get("info", np.nan),
            "eaf": r.get("eaf", np.nan),
            "status": "ok",
        })
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    lmm_cfg = cfg.get("lmm", {})
    target_id = lmm_cfg.get("target_variant_id", "5:158759900:A:G")
    chrpos = None
    parts = target_id.split(":")
    if len(parts) >= 2:
        try:
            chrpos = (int(parts[0]), int(parts[1]))
        except ValueError:
            chrpos = None

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    lmm_dir = results_dir / "lmm"
    log = setup_logger("lmm_ingest", logs_dir / "run_log.txt")

    manifest_path = inter_dir / "phenotype_manifest_extended.csv"
    if not manifest_path.is_file():
        manifest_path = results_dir / "phenotype_manifest.csv"

    if dry:
        log.info("DRY RUN — would parse REGENIE/BOLT outputs from {lmm_dir}")
        log.info(f"  manifest: {manifest_path} {'OK' if manifest_path.is_file() else 'MISSING'}")
        return 0

    if not manifest_path.is_file():
        log.error(f"manifest not found: {manifest_path}")
        return 1
    manifest = pd.read_csv(manifest_path)
    if "family" not in manifest.columns:
        manifest["family"] = "exploratory"
    if "source" not in manifest.columns:
        manifest["source"] = "basket"

    pieces = []
    reg_df = _ingest_regenie(lmm_dir, target_id, manifest, log)
    if not reg_df.empty:
        pieces.append(reg_df)
    bolt_df = _ingest_bolt(lmm_dir, target_id, chrpos,
                            target_rsid=cfg.get("variant", {}).get("rsid"),
                            manifest=manifest, log=log)
    if not bolt_df.empty:
        pieces.append(bolt_df)

    if not pieces:
        log.warning("no LMM results found — nothing to ingest")
        return 0

    all_df = pd.concat(pieces, ignore_index=True)
    all_df["effect_allele"] = cfg["variant"]["effect_allele"]
    all_df["other_allele"] = cfg["variant"]["other_allele"]

    # Per-family multiple testing on LMM p-values
    all_df["family_bonferroni"] = bonferroni_by_family(all_df, family_col="family", p_col="p")
    all_df["family_fdr_bh"] = bh_fdr_by_family(all_df, family_col="family", p_col="p")

    # Split by phenotype source for parallel-to-OLS output layout
    src = all_df["source"].fillna("basket")
    fam = all_df["family"].fillna("exploratory")
    splits = {
        "association_results_primary_lmm.csv":    all_df[~src.isin(["composite", "roi_pca"]) & (fam != "controls")],
        "association_results_composites_lmm.csv": all_df[src.eq("composite")],
        "association_results_roi_pca_lmm.csv":    all_df[src.eq("roi_pca")],
        "association_results_controls_lmm.csv":   all_df[fam.eq("controls")],
    }
    for name, sub in splits.items():
        if sub.empty:
            continue
        path = results_dir / name
        sub.sort_values("p", na_position="last").to_csv(path, index=False)
        log.info(f"wrote {path} ({len(sub)} rows)")

    # multiple_testing_summary_lmm.csv
    mt_rows = []
    for fam_name, sub in all_df.groupby("family", dropna=False):
        n = int(sub["p"].notna().sum())
        if n == 0:
            mt_rows.append({"family": fam_name, "n_tests": 0,
                             "alpha_bonf": np.nan,
                             "n_sig_bonf": 0, "n_sig_fdr_05": 0, "n_sig_fdr_10": 0,
                             "min_p": np.nan, "q_at_top_hit": np.nan})
            continue
        n_sig_bonf = int((sub["family_bonferroni"] <= 0.05).sum())
        n_sig_fdr05 = int((sub["family_fdr_bh"] <= 0.05).sum())
        n_sig_fdr10 = int((sub["family_fdr_bh"] <= 0.10).sum())
        top_q = float(sub.sort_values("p").iloc[0].get("family_fdr_bh", np.nan))
        mt_rows.append({
            "family": fam_name, "n_tests": n,
            "alpha_bonf": 0.05 / n,
            "n_sig_bonf": n_sig_bonf,
            "n_sig_fdr_05": n_sig_fdr05,
            "n_sig_fdr_10": n_sig_fdr10,
            "min_p": float(sub["p"].min()),
            "q_at_top_hit": top_q,
        })
    pd.DataFrame(mt_rows).to_csv(results_dir / "multiple_testing_summary_lmm.csv", index=False)
    log.info("wrote multiple_testing_summary_lmm.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
