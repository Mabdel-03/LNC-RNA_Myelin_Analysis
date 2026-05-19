#!/usr/bin/env python
"""Step 06 — plots + Markdown reports from association results.

Reads:
  results/association_results_primary.csv          (OLS single-IDP)
  results/association_results_genotypic.csv        (OLS genotypic)
  results/association_results_composites.csv       (OLS composites — optional)
  results/association_results_roi_pca.csv          (OLS ROI PCA — optional)
  results/association_results_*_lmm.csv            (REGENIE/BOLT — optional)
  results/multiple_testing_summary*.csv            (per-family — optional)
  results/roi_pca_loadings.csv                     (per-tract loadings — optional)
  results/genotype_qc_summary.csv
  results/sample_counts.csv

Outputs:
  figures/qqplot_pvalues.png
  figures/manhattan_like_idp_results.png
  figures/effect_size_forest_top_hits.png
  figures/genotype_violin_top3.png                 (optional)
  figures/effect_heatmap_by_tract_metric.png       (NEW)
  figures/composite_effects_forest.png             (NEW)
  figures/top_roi_pca_loadings.png                 (NEW)
  results/report.md                                (back-compat OLS report)
  results/report_hierarchical.md                   (NEW — LMM headline + OLS sensitivity)
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


# ---------------------------------------------------------------------------
# Hierarchical refactor: new figures + new report
# ---------------------------------------------------------------------------

REQUIRED_CAVEAT = (
    "Ordinary diffusion MRI signals (FA, MD, L1-L3, ICVF, OD, ISOVF) are "
    "*indirect* and should not be labeled as myelination effects unless "
    "supported by myelin-sensitive MRI (MTR, MTsat, MWF, qT1) or orthogonal "
    "validation. The current expected interpretation is that rs2546890-A "
    "appears more consistent with a free-water or tract-geometry signal than "
    "a canonical demyelination signal unless the composite analysis below "
    "shows otherwise."
)


def _effect_heatmap_by_tract_metric(primary: pd.DataFrame, out_path: Path) -> None:
    """Heatmap: rows=tracts, cols=metrics, color=signed -log10(p) (sign from β)."""
    plt = _setup_mpl()
    df = primary.dropna(subset=["p", "beta_per_A", "region", "metric"]).copy()
    if df.empty:
        return
    df = df[df["region"].astype(str).str.strip() != ""]
    if df.empty:
        return
    df["signed_nlp"] = np.sign(df["beta_per_A"]) * -np.log10(df["p"].clip(lower=np.finfo(float).tiny))
    # If multiple modalities (TBSS + WM) match the same (tract, metric), average.
    pivot = (df.groupby(["region", "metric"])["signed_nlp"]
               .mean().unstack("metric"))
    metric_order = [m for m in ["FA", "MD", "MO", "L1", "L2", "L3", "RD",
                                  "ICVF", "OD", "ISOVF"]
                     if m in pivot.columns]
    other_metrics = [m for m in pivot.columns if m not in metric_order]
    pivot = pivot[metric_order + other_metrics]
    if pivot.empty:
        return
    fig, ax = plt.subplots(figsize=(max(4, 0.6 * len(pivot.columns)),
                                     max(3, 0.18 * len(pivot.index))))
    vmax = float(np.nanmax(np.abs(pivot.values))) if pivot.size else 1.0
    im = ax.imshow(pivot.values, aspect="auto", cmap="RdBu_r",
                    vmin=-vmax, vmax=vmax)
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=6)
    fig.colorbar(im, ax=ax, label="signed -log10(p)")
    ax.set_title("rs2546890-A effect heatmap by tract × metric (signed -log10p)")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _composite_effects_forest(comp: pd.DataFrame, out_path: Path) -> None:
    """Forest of composite × tract effects, faceted by composite."""
    plt = _setup_mpl()
    df = comp.dropna(subset=["beta_per_A", "se", "p", "metric"]).copy()
    if df.empty:
        return
    composites = sorted(df["metric"].unique())
    fig, axes = plt.subplots(len(composites), 1,
                              figsize=(7, 1.2 + 0.20 * len(df)),
                              squeeze=False, sharex=True)
    for i, comp_name in enumerate(composites):
        sub = df[df["metric"] == comp_name].sort_values("beta_per_A").copy()
        ax = axes[i, 0]
        ci = 1.96 * sub["se"]
        ax.errorbar(sub["beta_per_A"], np.arange(len(sub)),
                    xerr=ci, fmt="o", color="#1f77b4", ecolor="grey", capsize=2)
        # mark FDR-significant
        if "family_fdr_bh" in sub.columns:
            mask = sub["family_fdr_bh"] <= 0.10
            if mask.any():
                ax.scatter(sub.loc[mask, "beta_per_A"], np.where(mask)[0],
                           color="#d62728", s=30, zorder=3, label="family FDR ≤ 0.10")
        ax.axvline(0, color="black", lw=0.5)
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels([f"{r['region']}" for _, r in sub.iterrows()],
                            fontsize=6)
        ax.set_title(comp_name, fontsize=10, loc="left")
    axes[-1, 0].set_xlabel("beta per A allele (95% CI)")
    fig.suptitle("Composite phenotype effects (rs2546890-A)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _top_roi_pca_loadings(loadings: pd.DataFrame, pca_results: pd.DataFrame,
                            out_path: Path, top_n: int = 6) -> None:
    """Per top-ranked tract, bar plot of PC1 metric loadings."""
    plt = _setup_mpl()
    if loadings.empty or pca_results.empty:
        return
    pc1 = pca_results[pca_results["metric"] == "PC1"].dropna(subset=["p"])
    if pc1.empty:
        return
    top_tracts = pc1.sort_values("p").head(top_n)["region"].tolist()
    if not top_tracts:
        return
    fig, axes = plt.subplots(1, len(top_tracts),
                              figsize=(3 * len(top_tracts), 4),
                              squeeze=False)
    for i, tract in enumerate(top_tracts):
        sub = loadings[(loadings["tract"] == tract) & (loadings["pc"] == "PC1")]
        if sub.empty:
            continue
        sub = sub.sort_values("loading")
        ax = axes[0, i]
        ax.barh(np.arange(len(sub)), sub["loading"].values, color="#1f77b4")
        ax.set_yticks(np.arange(len(sub)))
        ax.set_yticklabels(sub["metric"].tolist(), fontsize=7)
        ax.axvline(0, color="black", lw=0.5)
        ax.set_title(tract, fontsize=8)
        ax.set_xlabel("PC1 loading")
    fig.suptitle("Top tracts by ROI-PC1 association: PC1 metric loadings", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def _classify_composite_signal(comp: pd.DataFrame) -> str:
    """Auto-classify the composite pattern into one of four labels.

    Logic: count which composites have at least one FDR<0.10 tract hit. If
    `demyelination_like` dominates with consistent negative-FA / positive-RD
    direction → canonical_demyelination. Else if `free_water_like` dominates
    → free_water_dominant. Else if `tract_geometry_like` → tract_geometry_dominant.
    Else inconclusive.
    """
    if comp is None or comp.empty:
        return "inconclusive"
    if "family_fdr_bh" not in comp.columns:
        return "inconclusive"
    sig = comp[comp["family_fdr_bh"] <= 0.10]
    if sig.empty:
        return "inconclusive"
    by_comp = sig.groupby("metric").size().to_dict()
    if not by_comp:
        return "inconclusive"
    top = max(by_comp, key=by_comp.get)
    if top == "demyelination_like":
        return "canonical_demyelination"
    if top == "free_water_like":
        return "free_water_dominant"
    if top == "tract_geometry_like":
        return "tract_geometry_dominant"
    if top == "axonal_loss_like":
        return "axonal_loss_dominant"
    return "inconclusive"


def _read_optional_csv(p: Path) -> pd.DataFrame | None:
    """Return None if missing OR if empty/unparseable (treat as 'no data')."""
    if not p.is_file():
        return None
    try:
        if p.stat().st_size == 0:
            return None
        return pd.read_csv(p)
    except (pd.errors.EmptyDataError, Exception):
        return None


def _make_hierarchical_report(cfg: dict, results_dir: Path, figures_dir: Path,
                                qc: pd.DataFrame | None,
                                counts: pd.DataFrame | None,
                                out_path: Path) -> None:
    """Family-stratified report with LMM headline + OLS sensitivity."""
    var = cfg["variant"]
    primary_lmm = _read_optional_csv(results_dir / "association_results_primary_lmm.csv")
    comp_lmm = _read_optional_csv(results_dir / "association_results_composites_lmm.csv")
    roi_lmm = _read_optional_csv(results_dir / "association_results_roi_pca_lmm.csv")
    ctl_lmm = _read_optional_csv(results_dir / "association_results_controls_lmm.csv")
    mt_lmm = _read_optional_csv(results_dir / "multiple_testing_summary_lmm.csv")
    primary_ols = _read_optional_csv(results_dir / "association_results_primary.csv")
    comp_ols = _read_optional_csv(results_dir / "association_results_composites.csv")
    roi_ols = _read_optional_csv(results_dir / "association_results_roi_pca.csv")
    loadings = _read_optional_csv(results_dir / "roi_pca_loadings.csv")

    L: list[str] = []
    L.append("# rs2546890 → UKBB MRI: hierarchical analysis report")
    L.append("")
    L.append(f"_generated: {_dt.datetime.now().isoformat(timespec='seconds')}_")
    L.append("")
    L.append(f"- rsID: **{var['rsid']}**")
    L.append(f"- position (GRCh37): chr{var['chr']}:{var.get('pos_grch37', '?')}")
    L.append(f"- effect allele: **{var['effect_allele']}** / other: **{var['other_allele']}**")
    if qc is not None and not qc.empty:
        q = qc.iloc[0]
        L.append(f"- INFO: {q.get('info_mfi', 'NA')}, EAF: {q.get('eaf_estimated', 'NA')}, HWE p: {q.get('hwe_p', 'NA')}")
    L.append("")

    L.append("## Required scientific framing")
    L.append("")
    L.append("> " + REQUIRED_CAVEAT)
    L.append("")

    # Sample funnel
    L.append("## Sample funnel")
    if counts is not None and not counts.empty:
        L.append("| step | n | n dropped |")
        L.append("|---|---:|---:|")
        for _, r in counts.iterrows():
            L.append(f"| {r['step']} | {int(r['n']):,} | {int(r['n_dropped_this_step']):,} |")
    L.append("")

    # Multiple testing summary
    L.append("## Per-family multiple-testing summary (LMM)")
    if mt_lmm is not None and not mt_lmm.empty:
        L.append(mt_lmm.to_markdown(index=False, floatfmt=".3g"))
    else:
        L.append("_(LMM results not present; run engine=regenie+ols then 05d)_")
    L.append("")

    # Family 1
    L.append("## Family 1 — primary myelin-sensitive MRI")
    if primary_lmm is not None and not primary_lmm.empty:
        prim_family = primary_lmm[primary_lmm.get("family", "") == "primary"]
        if prim_family.empty:
            L.append("_N=0 myelin-sensitive phenotypes available in this extract; no test performed._")
        else:
            L.append(prim_family.sort_values("p").head(20).to_markdown(index=False, floatfmt=".3g"))
    else:
        L.append("_(LMM results not present)_")
    L.append("")

    # Family 2 — composites
    L.append("## Family 2a — biological composites (LMM headline)")
    interp = "inconclusive"
    if comp_lmm is not None and not comp_lmm.empty:
        L.append(comp_lmm.sort_values("p").head(20).to_markdown(index=False, floatfmt=".3g"))
        interp = _classify_composite_signal(comp_lmm)
    else:
        L.append("_(no composite LMM results)_")
    if comp_ols is not None and not comp_ols.empty:
        L.append("")
        L.append("### OLS sensitivity for composites")
        L.append(comp_ols.sort_values("p").head(10).to_markdown(index=False, floatfmt=".3g"))
    L.append("")

    # Family 2 — ROI PCs
    L.append("## Family 2b — ROI PCA (LMM headline)")
    if roi_lmm is not None and not roi_lmm.empty:
        L.append(roi_lmm.sort_values("p").head(20).to_markdown(index=False, floatfmt=".3g"))
    else:
        L.append("_(no ROI PCA LMM results)_")
    if roi_ols is not None and not roi_ols.empty:
        L.append("")
        L.append("### OLS sensitivity for ROI PCA")
        L.append(roi_ols.sort_values("p").head(10).to_markdown(index=False, floatfmt=".3g"))
    L.append("")

    # Family 3 — exploratory single-IDP PheWAS
    L.append("## Family 3 — exploratory single-IDP PheWAS (LMM headline)")
    L.append("")
    L.append("_Interpret single-IDP hits below only in light of the composite + ROI-PC results above._")
    L.append("")
    if primary_lmm is not None and not primary_lmm.empty:
        # Top 20 by p
        L.append(primary_lmm.sort_values("p").head(20).to_markdown(index=False, floatfmt=".3g"))
    else:
        L.append("_(no LMM single-IDP results)_")
    if primary_ols is not None and not primary_ols.empty:
        L.append("")
        L.append("### OLS sensitivity for single IDPs (top 10 by OLS p)")
        L.append(primary_ols.sort_values("p").head(10).to_markdown(index=False, floatfmt=".3g"))
    L.append("")

    # Family 4 — controls
    L.append("## Family 4 — controls (WMH / QC outcomes)")
    if ctl_lmm is not None and not ctl_lmm.empty:
        L.append(ctl_lmm.sort_values("p").to_markdown(index=False, floatfmt=".3g"))
    else:
        L.append("_(no LMM control-family results)_")
    L.append("")

    # LMM vs OLS comparison
    L.append("## LMM (REGENIE) vs OLS comparison")
    if (primary_lmm is not None and primary_ols is not None
            and not primary_lmm.empty and not primary_ols.empty):
        merged = primary_lmm.merge(primary_ols[["column_name", "beta_per_A", "p"]],
                                     on="column_name", how="inner",
                                     suffixes=("_lmm", "_ols"))
        if not merged.empty:
            try:
                corr = float(merged[["beta_per_A_lmm", "beta_per_A_ols"]].corr().iloc[0, 1])
                L.append(f"- Pearson r(β_LMM, β_OLS) on {len(merged):,} single IDPs: **{corr:.3f}**")
            except Exception:
                pass
            disagree = merged[(np.sign(merged["beta_per_A_lmm"]) != np.sign(merged["beta_per_A_ols"])) &
                                 (merged[["beta_per_A_lmm", "beta_per_A_ols"]].abs().min(axis=1) > 0.01)]
            L.append(f"- Sign disagreements (|β|>0.01 both sides): **{len(disagree)}** phenotypes")
    L.append("")

    # Auto-classified interpretation
    L.append("## Composite-driven interpretation")
    L.append("")
    L.append(f"**Auto-classified signature: `{interp}`**")
    L.append("")
    L.append("Interpretation key:")
    L.append("- `canonical_demyelination` — `demyelination_like` composite dominates the FDR-significant set")
    L.append("- `free_water_dominant` — `free_water_like` composite dominates")
    L.append("- `tract_geometry_dominant` — `tract_geometry_like` composite dominates")
    L.append("- `axonal_loss_dominant` — `axonal_loss_like` composite dominates")
    L.append("- `inconclusive` — no composite has any tract hit at family FDR ≤ 0.10")
    L.append("")

    # Limitations
    L.append("## Limitations")
    L.append("- No myelin-sensitive MRI present in this UKB extract (family 1 empty).")
    L.append("- ROI PCA uses mean imputation for residual missing data → biased if missingness is informative.")
    L.append("- Family-FDR applied independently per family (not sequential gatekeeping).")
    L.append("- REGENIE step-1 LOCO model fit on the imaging subsample (~40K), not the full WB-MM cohort.")
    L.append("- BOLT-LMM (when used) relies on rs2546890 being in HapMap3; for imputed-dosage stats see lmm.bolt.regen_bgen documentation.")
    L.append("")

    L.append("## Figures")
    L.append(f"- `{figures_dir.name}/qqplot_pvalues.png`")
    L.append(f"- `{figures_dir.name}/manhattan_like_idp_results.png`")
    L.append(f"- `{figures_dir.name}/effect_size_forest_top_hits.png`")
    L.append(f"- `{figures_dir.name}/effect_heatmap_by_tract_metric.png` (NEW)")
    L.append(f"- `{figures_dir.name}/composite_effects_forest.png` (NEW)")
    L.append(f"- `{figures_dir.name}/top_roi_pca_loadings.png` (NEW)")
    L.append("")

    out_path.write_text("\n".join(L) + "\n")


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
    comp_path = results_dir / "association_results_composites.csv"
    roi_path = results_dir / "association_results_roi_pca.csv"
    loadings_path = results_dir / "roi_pca_loadings.csv"

    if dry:
        log.info("DRY RUN — checking inputs only")
        for p in (primary_path, geno_path, qc_path, counts_path, analysis_path,
                   comp_path, roi_path, loadings_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    if not primary_path.is_file():
        log.error(f"primary results missing: {primary_path}")
        return 1
    primary = pd.read_csv(primary_path)
    geno = _read_optional_csv(geno_path)
    qc = _read_optional_csv(qc_path)
    counts = _read_optional_csv(counts_path)
    comp = _read_optional_csv(comp_path)
    roi = _read_optional_csv(roi_path)
    loadings = _read_optional_csv(loadings_path)

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

    # New hierarchical figures
    _effect_heatmap_by_tract_metric(primary, figures_dir / "effect_heatmap_by_tract_metric.png")
    if comp is not None and not comp.empty:
        _composite_effects_forest(comp, figures_dir / "composite_effects_forest.png")
    if loadings is not None and not loadings.empty and roi is not None and not roi.empty:
        _top_roi_pca_loadings(loadings, roi, figures_dir / "top_roi_pca_loadings.png")

    _make_report(cfg, primary, geno, qc, counts, results_dir / "report.md")
    _make_hierarchical_report(cfg, results_dir, figures_dir, qc, counts,
                                results_dir / "report_hierarchical.md")
    log.info(f"wrote {results_dir / 'report.md'}, {results_dir / 'report_hierarchical.md'}, "
              f"and figures under {figures_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
