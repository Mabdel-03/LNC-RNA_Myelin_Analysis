#!/usr/bin/env python
"""Step 06 — plots + Markdown report from association results.

Inputs:
  results/association_results_primary.csv
  results/association_results_genotypic.csv (optional)
  results/genotype_qc_summary.csv
  results/sample_counts.csv
  results/intermediate/analysis_ready.tsv (for optional violin plot)

Outputs:
  figures/qqplot_pvalues.png
  figures/manhattan_like_idp_results.png
  figures/effect_size_forest_top_hits.png
  figures/genotype_violin_top3.png (optional, if genotype + top hits exist)
  results/report.md
"""
from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import load_config, resolve_under_repo, setup_logger  # noqa: E402


def _setup_mpl():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 150
    return plt


def _qq_plot(primary: pd.DataFrame, out_path: Path) -> None:
    plt = _setup_mpl()
    p = primary["p"].dropna().values
    if len(p) == 0:
        return
    n = len(p)
    expected = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    observed = -np.log10(np.sort(p))
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(expected, observed, s=12, alpha=0.7)
    lim = max(expected.max(), observed.max())
    ax.plot([0, lim], [0, lim], color="grey", linestyle="--", lw=1)
    ax.set_xlabel("Expected -log10(p)")
    ax.set_ylabel("Observed -log10(p)")
    ax.set_title(f"QQ — primary OLS (n IDPs = {n})")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _manhattan_by_idp(primary: pd.DataFrame, out_path: Path) -> None:
    plt = _setup_mpl()
    df = primary.dropna(subset=["p"]).copy()
    if df.empty:
        return
    df = df.sort_values(["panel", "modality", "metric", "column_name"]).reset_index(drop=True)
    df["x"] = np.arange(len(df))
    df["nlp"] = -np.log10(df["p"])
    fig, ax = plt.subplots(figsize=(max(6, len(df) * 0.08), 4))
    palette = {"primary": "#1f77b4", "secondary": "#ff7f0e", "exploratory": "#7f7f7f"}
    for panel, sub in df.groupby("panel"):
        ax.scatter(sub["x"], sub["nlp"], s=14, label=panel,
                   color=palette.get(panel, "#444"), alpha=0.85)
    if "bonferroni" in df.columns and df["bonferroni"].notna().any():
        thr = 0.05 / max(1, df["panel"].value_counts().max())
        ax.axhline(-np.log10(thr), color="red", linestyle="--", lw=0.8,
                   label=f"Bonferroni 0.05 (per-panel n={df['panel'].value_counts().max()})")
    ax.set_xlabel("IDP (ordered by panel/modality)")
    ax.set_ylabel("-log10(p)")
    ax.set_title("Per-IDP association p-values for rs2546890 dosage")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _forest_top(primary: pd.DataFrame, out_path: Path, top_n: int = 20) -> None:
    plt = _setup_mpl()
    df = primary.dropna(subset=["beta_per_A", "se", "p"]).copy()
    if df.empty:
        return
    df = df.sort_values("p").head(top_n).iloc[::-1]
    ci = 1.96 * df["se"]
    fig, ax = plt.subplots(figsize=(6, max(3, 0.3 * len(df))))
    ax.errorbar(df["beta_per_A"], np.arange(len(df)),
                xerr=ci, fmt="o", color="#1f77b4", ecolor="grey", capsize=2)
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["column_name"], fontsize=7)
    ax.set_xlabel("beta per A allele (95% CI)")
    ax.set_title(f"Top {len(df)} IDPs by p")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _violin_top(primary: pd.DataFrame, analysis: pd.DataFrame,
                 out_path: Path, top_n: int = 3) -> None:
    if "genotype" not in analysis.columns:
        return
    plt = _setup_mpl()
    df = primary.dropna(subset=["p"]).sort_values("p").head(top_n)
    if df.empty:
        return
    fig, axes = plt.subplots(1, len(df), figsize=(4 * len(df), 4), squeeze=False)
    for i, (_, row) in enumerate(df.iterrows()):
        idp = row["column_name"]
        if idp not in analysis.columns:
            continue
        sub = analysis[["genotype", idp]].dropna()
        if sub.empty:
            continue
        order = sorted(sub["genotype"].unique())
        data = [sub.loc[sub["genotype"] == g, idp].values for g in order]
        ax = axes[0, i]
        ax.violinplot(data, showmeans=True)
        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(order)
        ax.set_title(idp, fontsize=8)
        ax.set_ylabel("raw value")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _make_report(cfg: dict, primary: pd.DataFrame,
                  geno: pd.DataFrame | None,
                  qc: pd.DataFrame | None,
                  counts: pd.DataFrame | None,
                  out_path: Path) -> None:
    var = cfg["variant"]
    lines: list[str] = []
    lines.append(f"# rs2546890 → UKBB MRI IDP association — report")
    lines.append("")
    lines.append(f"_generated: {_dt.datetime.now().isoformat(timespec='seconds')}_")
    lines.append("")
    lines.append("## Variant")
    lines.append(f"- rsID: **{var['rsid']}**")
    lines.append(f"- position: chr{var['chr']}:{var['pos_grch38']} (GRCh38)")
    lines.append(f"- effect allele: **{var['effect_allele']}** (modeled in additive dosage)")
    lines.append(f"- other allele:  **{var['other_allele']}**")
    if qc is not None and not qc.empty:
        q = qc.iloc[0]
        lines.append(f"- INFO (MFI):    {q.get('info_mfi', 'NA')}")
        lines.append(f"- est. EAF:      {q.get('eaf_estimated', 'NA')}")
        lines.append(f"- HWE p:         {q.get('hwe_p', 'NA')}")
        ea = var["effect_allele"]; oa = var["other_allele"]
        for lab in (ea + ea, "".join(sorted([ea, oa])), oa + oa):
            col = f"hardcall_{lab}"
            lines.append(f"- hardcall {lab}: {q.get(col, 'NA')}")
    lines.append("")
    lines.append("## Sample funnel")
    if counts is not None and not counts.empty:
        lines.append("| step | n | n dropped |")
        lines.append("|---|---:|---:|")
        for _, r in counts.iterrows():
            lines.append(f"| {r['step']} | {int(r['n']):,} | {int(r['n_dropped_this_step']):,} |")
    lines.append("")
    lines.append("## Phenotypes tested")
    counts_panel = primary["panel"].value_counts().to_dict()
    for k, v in counts_panel.items():
        lines.append(f"- {k}: {v}")
    lines.append("")
    lines.append("## Top 20 associations (primary additive model)")
    show_cols = ["column_name", "panel", "modality", "metric", "region",
                 "n", "beta_per_A", "se", "p", "bonferroni", "fdr_bh",
                 "aa_vs_gg_additive_2beta"]
    show_cols = [c for c in show_cols if c in primary.columns]
    top = primary.dropna(subset=["p"]).sort_values("p").head(20)[show_cols]
    if not top.empty:
        lines.append(top.to_markdown(index=False, floatfmt=".3g"))
    else:
        lines.append("_(no successful models)_")
    lines.append("")
    sig = primary[(primary.get("fdr_bh", 1) <= 0.1)]
    lines.append(f"## FDR-q ≤ 0.1: {len(sig)} IDP(s)")
    if not sig.empty and geno is not None and not geno.empty:
        lines.append("")
        lines.append("### Sensitivity (additive 2β vs genotypic AA-vs-GG) for FDR hits")
        ea = cfg["variant"]["effect_allele"]
        oa = cfg["variant"]["other_allele"]
        col_aa = f"beta_{ea+ea}_vs_{oa+oa}"
        cmp = sig.merge(geno[["column_name"] + [c for c in geno.columns if c.startswith("beta_") or c == "wald_p"]],
                        on="column_name", how="left")
        keep = ["column_name", "beta_per_A", "aa_vs_gg_additive_2beta",
                col_aa, "p", "wald_p"]
        keep = [c for c in keep if c in cmp.columns]
        if keep:
            lines.append(cmp[keep].to_markdown(index=False, floatfmt=".3g"))
    lines.append("")
    lines.append("## Caveats")
    inputs = cfg["inputs"]
    has_anc = (bool(inputs.get("ancestry_relatedness_file"))
               or bool(inputs.get("genetic_covariates_file")))
    has_unrel = (bool(inputs.get("ancestry_relatedness_file"))
                 or bool(inputs.get("relatedness_file"))
                 or bool(inputs.get("genetic_covariates_file")))
    if not has_anc:
        lines.append("- No ancestry source supplied (neither ancestry_relatedness_file nor "
                     "genetic_covariates_file with SQC); analysis is PC-adjusted only.")
    if not has_unrel:
        lines.append("- No relatedness source supplied; cryptic relatedness not pruned.")
    if not inputs.get("imaging_confounds_file"):
        lines.append("- UKB imaging confounds file (Smith et al. 2020 Resource 1977) not supplied; "
                     "only core imaging covariates (age, sex, site, head size, dMRI outlier slices) included.")
    if not inputs.get("data_dictionary_file"):
        lines.append("- No UKB data dictionary supplied; tract/region labels were inferred "
                     "from field-ID position within known IDP blocks, not authoritative names.")
    lines.append("- OLS additive model is primary; AA-vs-GG read off as 2 × β. "
                 "A 2-df genotypic model provides the sensitivity check.")
    lines.append(f"- Test mode capped phenotypes to {cfg['project'].get('max_test_phenotypes', '?')}; "
                 f"set project.test_mode=false to run all IDPs."
                 if cfg.get("project", {}).get("test_mode", False) else
                 "- Full panel run (test_mode disabled).")
    out_path.write_text("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    figures_dir = resolve_under_repo(cfg, cfg["project"]["figures_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    log = setup_logger("plots_report", logs_dir / "run_log.txt")

    primary_path = results_dir / "association_results_primary.csv"
    geno_path = results_dir / "association_results_genotypic.csv"
    qc_path = results_dir / "genotype_qc_summary.csv"
    counts_path = results_dir / "sample_counts.csv"
    analysis_path = inter_dir / "analysis_ready.tsv"

    if dry:
        log.info("DRY RUN — checking inputs only")
        for p in (primary_path, geno_path, qc_path, counts_path, analysis_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    if not primary_path.is_file():
        log.error(f"primary results missing: {primary_path}")
        return 1
    primary = pd.read_csv(primary_path)
    geno = pd.read_csv(geno_path) if geno_path.is_file() else None
    qc = pd.read_csv(qc_path) if qc_path.is_file() else None
    counts = pd.read_csv(counts_path) if counts_path.is_file() else None

    figures_dir.mkdir(parents=True, exist_ok=True)
    _qq_plot(primary, figures_dir / "qqplot_pvalues.png")
    _manhattan_by_idp(primary, figures_dir / "manhattan_like_idp_results.png")
    _forest_top(primary, figures_dir / "effect_size_forest_top_hits.png")
    if analysis_path.is_file():
        try:
            analysis = pd.read_csv(analysis_path, sep="\t")
            _violin_top(primary, analysis, figures_dir / "genotype_violin_top3.png")
        except Exception as exc:
            log.warning(f"violin plot skipped: {exc}")

    _make_report(cfg, primary, geno, qc, counts, results_dir / "report.md")
    log.info(f"wrote {results_dir / 'report.md'} and figures under {figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
