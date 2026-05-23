#!/usr/bin/env python
"""Targeted check: rs2546890-A effect on FA + axial diffusivity (L1) + radial
diffusivity (RD), reported per-tract WITHOUT any multiple-testing correction.

Reads (already produced by 05a_run_ols.py):
  results/association_results_primary.csv      (single-IDP OLS)
  results/association_results_genotypic.csv    (genotypic AA-vs-GG sensitivity)

Writes:
  results/targeted_diffusion_check.csv         (focused per-tract effects, raw p)
  figures/targeted_diffusion_forest.png         (forest plot by metric)

Usage:
  python src/targeted_diffusion_check.py --config config.yaml

Notes:
  - "Axial diffusivity" = L1 (principal-axis eigenvalue of the diffusion tensor).
  - RD is the (L2+L3)/2 derived phenotype already in the manifest from step 02.
  - This script does NOT recompute models — it just filters the existing OLS
    output. Refresh by re-running 05a_run_ols.py if the underlying analysis
    changed.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import load_config, resolve_under_repo, setup_logger  # noqa: E402


# Targeted metric set. We treat L1 as axial diffusivity (AD).
TARGET_METRICS = {"FA": "Fractional anisotropy",
                  "L1": "Axial diffusivity (L1)",
                  "RD": "Radial diffusivity (L2+L3)/2"}


def _setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 150
    return plt


def _forest_plot(df: pd.DataFrame, out_path: Path) -> None:
    plt = _setup_mpl()
    metrics = [m for m in ("FA", "L1", "RD") if m in df["metric"].unique()]
    if not metrics:
        return
    fig, axes = plt.subplots(len(metrics), 1,
                              figsize=(8, 1.0 + 0.20 * len(df)),
                              squeeze=False, sharex=True)
    for i, m in enumerate(metrics):
        sub = df[df["metric"] == m].sort_values("beta_per_A").copy()
        ax = axes[i, 0]
        ci = 1.96 * sub["se"]
        ax.errorbar(sub["beta_per_A"], np.arange(len(sub)),
                    xerr=ci, fmt="o", color="#1f77b4", ecolor="grey",
                    capsize=2, markersize=4)
        # Highlight nominally significant (raw p < 0.05) — NO FDR here
        sig = sub["p"] < 0.05
        if sig.any():
            ax.scatter(sub.loc[sig, "beta_per_A"], np.where(sig)[0],
                       color="#d62728", s=30, zorder=3,
                       label=f"raw p < 0.05 ({int(sig.sum())} tracts)")
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels([f"{r['modality']} | {r['region']}" for _, r in sub.iterrows()],
                            fontsize=6)
        ax.set_title(TARGET_METRICS[m], fontsize=10, loc="left")
        if sig.any():
            ax.legend(loc="lower right", fontsize=8)
    axes[-1, 0].set_xlabel("β per A allele (95% CI) — raw p, NO FDR correction")
    fig.suptitle("rs2546890-A → diffusion metrics (FA, axial, radial)",
                  fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    figures_dir = resolve_under_repo(cfg, cfg["project"]["figures_dir"])
    log = setup_logger("targeted_diff", logs_dir / "run_log.txt")

    primary_path = results_dir / "association_results_primary.csv"
    primary_lmm_path = results_dir / "association_results_primary_lmm.csv"
    geno_path = results_dir / "association_results_genotypic.csv"
    out_path = results_dir / "targeted_diffusion_check.csv"
    fig_path = figures_dir / "targeted_diffusion_forest.png"

    if not primary_path.is_file():
        log.error(f"missing {primary_path}; run src/05a_run_ols.py first")
        return 1
    primary = pd.read_csv(primary_path)
    geno = pd.read_csv(geno_path) if geno_path.is_file() else None
    primary_lmm = pd.read_csv(primary_lmm_path) if primary_lmm_path.is_file() else None

    # Filter to the targeted metric set
    df = primary[primary["metric"].astype(str).isin(TARGET_METRICS.keys())].copy()
    if df.empty:
        log.error("no rows match target metrics (FA, L1, RD); is the manifest "
                  "missing the dMRI single-IDP rows?")
        return 2

    # Drop rows where the model didn't fit
    if "status" in df.columns:
        df = df[df["status"].astype(str) == "ok"]

    # Join in genotypic AA-vs-GG sensitivity when available
    ea = cfg["variant"]["effect_allele"]
    oa = cfg["variant"]["other_allele"]
    col_aa = f"beta_{ea + ea}_vs_{oa + oa}"
    col_aa_p = f"p_{ea + ea}_vs_{oa + oa}"
    if geno is not None and col_aa in geno.columns:
        df = df.merge(
            geno[["column_name", col_aa, col_aa_p, "wald_p"]].rename(
                columns={col_aa: "aa_vs_gg_beta_genotypic",
                          col_aa_p: "aa_vs_gg_p_genotypic"}),
            on="column_name", how="left")

    # Join in REGENIE-LMM single-IDP stats when available (engine-independent
    # confirmation under a mixed model). Same per-tract row, raw p only.
    if primary_lmm is not None and "column_name" in primary_lmm.columns:
        lmm_cols = ["column_name", "n", "beta_per_A", "se", "p", "log10p"]
        lmm_cols = [c for c in lmm_cols if c in primary_lmm.columns]
        sub = primary_lmm[lmm_cols].rename(columns={
            "n": "n_lmm",
            "beta_per_A": "beta_per_A_lmm",
            "se": "se_lmm",
            "p": "p_lmm",
            "log10p": "log10p_lmm",
        })
        df = df.merge(sub, on="column_name", how="left")
        log.info(f"joined LMM column for {df['p_lmm'].notna().sum() if 'p_lmm' in df.columns else 0} "
                 f"of {len(df)} rows")

    keep_cols = [c for c in [
        "metric", "modality", "region", "column_name",
        "n", "beta_per_A", "se", "t", "p",
        "n_lmm", "beta_per_A_lmm", "se_lmm", "p_lmm", "log10p_lmm",
        "aa_vs_gg_additive_2beta",
        "aa_vs_gg_beta_genotypic", "aa_vs_gg_p_genotypic", "wald_p",
        "field_id",
    ] if c in df.columns]
    df = df[keep_cols].sort_values(["metric", "p"], na_position="last").reset_index(drop=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    log.info(f"wrote {out_path} ({len(df)} rows: "
              f"{(df['metric'] == 'FA').sum()} FA, "
              f"{(df['metric'] == 'L1').sum()} L1/axial, "
              f"{(df['metric'] == 'RD').sum()} RD)")

    # Forest plot
    figures_dir.mkdir(parents=True, exist_ok=True)
    _forest_plot(df, fig_path)
    log.info(f"wrote {fig_path}")

    # Console summary of nominally significant (raw p < 0.05) hits per metric
    log.info("")
    log.info("=" * 72)
    log.info("Raw p < 0.05 hits (NO multiple-testing correction)")
    log.info("=" * 72)
    has_lmm = "p_lmm" in df.columns
    for m, label in TARGET_METRICS.items():
        sub = df[(df["metric"] == m) & (df["p"] < 0.05)]
        if sub.empty:
            log.info(f"  {label}: 0/{(df['metric'] == m).sum()}")
            continue
        log.info(f"  {label}: {len(sub)}/{(df['metric'] == m).sum()} tracts at OLS p<0.05")
        for _, r in sub.sort_values("p").iterrows():
            lmm_str = ""
            if has_lmm and pd.notna(r.get("p_lmm")):
                lmm_str = (f"   [LMM β={r['beta_per_A_lmm']:+.4f}±"
                           f"{r['se_lmm']:.4f}  p={r['p_lmm']:.2e}]")
            log.info(f"    {r['modality']:24s} {r['region']:50s} "
                      f"OLS β={r['beta_per_A']:+.4f}±{r['se']:.4f}  p={r['p']:.2e}{lmm_str}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
