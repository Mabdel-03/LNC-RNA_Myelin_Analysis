#!/usr/bin/env python
"""Step 04 — join genotype + phenotype + covariates on eid, apply sample QC,
write the analysis-ready table and the funnel counts.

Funnel order (each step logged):
  1. inner-join eid across genotype + phenotype + covariates
  2. ancestry restriction (if enabled) — gives the LMM-eligible pool
  3. drop rows missing required covariates / dosage — gives the LMM superset
     (this is what `analysis_ready.tsv` contains)
  4. KING relatedness pruning at sample.kinship_threshold — gives the OLS subset
  5. optional MZ-twin filter at lmm.drop_kinship_above — applied to LMM keep list
  6. genotype QC: variant-wide INFO check (halt unless --force)

Outputs:
  results/intermediate/analysis_ready.tsv  (LMM superset; OLS filters to keep_ols)
  results/intermediate/keep_ols.txt        (FID IID — KING-unrelated subset)
  results/intermediate/keep_lmm.txt        (FID IID — LMM keep, may exclude MZ twins)
  results/sample_counts.csv                (per-step funnel counts)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    resolve_under_repo,
    safe_read_table,
    setup_logger,
    standardize_eid_column,
)


REQUIRED_COVARS = ["age", "sex"]  # minimum non-PC covariates


def _funnel_row(name: str, n: int, drop_from: int | None = None) -> dict:
    return {"step": name, "n": int(n),
            "n_dropped_this_step": (int(drop_from - n) if drop_from is not None else 0)}


# Map config.sample.ancestry_label → predicate over SQC columns.
def _ancestry_keep_mask(df: pd.DataFrame, label: str) -> pd.Series | None:
    """Return a boolean mask using SQC columns, or None if no SQC info present."""
    lab = (label or "").strip().lower()
    cols = set(df.columns)
    pop = df["population"].astype(str).str.upper() if "population" in cols else None
    pop_mm = df["population_MM"].astype(str).str.upper() if "population_MM" in cols else None
    sr_wb = df["sr_WB"] if "sr_WB" in cols else None
    if lab in ("white british", "wb", "wb_mm"):
        m = pd.Series(False, index=df.index)
        if pop_mm is not None:
            m |= pop_mm.eq("WB_MM")
        if pop is not None:
            m |= pop.eq("WB")
        if sr_wb is not None:
            m |= pd.to_numeric(sr_wb, errors="coerce").eq(1)
        return m
    if lab in ("eur", "european", "europeans"):
        m = pd.Series(False, index=df.index)
        if pop is not None:
            m |= pop.isin(["WB", "NBW"])
        if pop_mm is not None:
            m |= pop_mm.isin(["WB_MM", "NBW_MM"])
        return m
    if lab in ("afr", "african"):
        m = pd.Series(False, index=df.index)
        if pop is not None:
            m |= pop.eq("AFR")
        if pop_mm is not None:
            m |= pop_mm.eq("AFR_MM")
        return m
    if lab in ("sa", "south asian"):
        m = pd.Series(False, index=df.index)
        if pop is not None:
            m |= pop.eq("SA")
        if pop_mm is not None:
            m |= pop_mm.eq("SA_MM")
        return m
    # unknown label
    return None


def _king_unrelated_keep(merged_eids: pd.Series, rel_path: str,
                          kinship_thresh: float, logger) -> set:
    """Greedy prune of KING relatedness pairs to the unrelated subset.

    The .dat[.gz] is whitespace-separated with cols: ID1 ID2 HetHet IBS0 Kinship.
    We restrict to pairs where both IDs are in `merged_eids` and Kinship >= threshold,
    then remove vertices with highest degree until no edges remain.
    """
    eid_set = set(int(x) for x in merged_eids.dropna().astype(int))
    logger.info(f"  loading relatedness file {rel_path}")
    df = pd.read_csv(rel_path, sep=r"\s+", engine="python")
    df = df.dropna(subset=["ID1", "ID2", "Kinship"])
    df["ID1"] = pd.to_numeric(df["ID1"], errors="coerce").astype("Int64")
    df["ID2"] = pd.to_numeric(df["ID2"], errors="coerce").astype("Int64")
    df["Kinship"] = pd.to_numeric(df["Kinship"], errors="coerce")
    df = df.dropna(subset=["ID1", "ID2", "Kinship"])
    df = df[(df["Kinship"] >= kinship_thresh)
            & df["ID1"].isin(eid_set) & df["ID2"].isin(eid_set)]
    logger.info(f"  {len(df):,} kin pairs above threshold {kinship_thresh}")
    if df.empty:
        return eid_set
    # Build adjacency
    from collections import defaultdict
    adj: dict[int, set] = defaultdict(set)
    for a, b in zip(df["ID1"].astype(int), df["ID2"].astype(int)):
        adj[a].add(b); adj[b].add(a)
    # Greedy: repeatedly remove highest-degree vertex
    removed: set = set()
    while adj:
        v, neighbors = max(adj.items(), key=lambda kv: len(kv[1]))
        if not neighbors:
            del adj[v]
            continue
        removed.add(v)
        for n in list(neighbors):
            adj[n].discard(v)
            if not adj[n]:
                del adj[n]
        del adj[v]
    keep = eid_set - removed
    logger.info(f"  greedy prune removed {len(removed):,}; "
                 f"unrelated keep set = {len(keep):,}")
    return keep


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--force", action="store_true",
                        help="proceed even if variant-wide QC (e.g. INFO) fails")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    inputs = cfg["inputs"]

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    inter_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logger("merge_qc", logs_dir / "run_log.txt")

    geno_path = inter_dir / "variant_dosage.tsv"
    pheno_path = inter_dir / "phenotypes_wide.tsv"
    covar_path = inter_dir / "covariates.tsv"
    qc_path = results_dir / "genotype_qc_summary.csv"
    out_path = inter_dir / "analysis_ready.tsv"
    counts_path = results_dir / "sample_counts.csv"
    keep_ols_path = inter_dir / "keep_ols.txt"
    keep_lmm_path = inter_dir / "keep_lmm.txt"

    if dry:
        log.info("DRY RUN — checking required intermediate files exist")
        for p in (geno_path, pheno_path, covar_path):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        return 0

    for p in (geno_path, pheno_path, covar_path):
        if not p.is_file():
            log.error(f"required intermediate missing: {p}")
            return 1

    log.info("loading intermediates")
    geno = pd.read_csv(geno_path, sep="\t")
    geno = standardize_eid_column(geno)
    pheno = pd.read_csv(pheno_path, sep="\t")
    pheno = standardize_eid_column(pheno)
    covar = pd.read_csv(covar_path, sep="\t")
    covar = standardize_eid_column(covar)

    n_geno = len(geno)
    n_pheno = len(pheno)
    n_covar = len(covar)
    log.info(f"  genotype rows: {n_geno:,}, phenotype rows: {n_pheno:,}, "
             f"covariate rows: {n_covar:,}")

    funnel = [
        _funnel_row("genotype_input", n_geno),
        _funnel_row("phenotype_input", n_pheno),
        _funnel_row("covariate_input", n_covar),
    ]

    # ---- step 1: inner join
    merged = geno.merge(pheno, on="eid", how="inner").merge(covar, on="eid", how="inner")
    funnel.append(_funnel_row("merged_inner_join", len(merged)))
    log.info(f"after inner join: {len(merged):,}")

    # ---- step 2: ancestry restriction
    label = cfg["sample"].get("ancestry_label", "White British")
    anc_path = inputs.get("ancestry_relatedness_file") or ""
    if cfg["sample"].get("restrict_to_ancestry", False):
        if anc_path and Path(anc_path).is_file():
            anc = safe_read_table(anc_path)
            anc = standardize_eid_column(anc)
            if "ancestry_label" in anc.columns:
                keep_ids = anc.loc[anc["ancestry_label"].astype(str) == label, "eid"]
                before = len(merged)
                merged = merged[merged["eid"].isin(keep_ids)].copy()
                funnel.append(_funnel_row(f"ancestry={label}", len(merged), before))
                log.info(f"after ancestry restriction ({label}, via file): {len(merged):,}")
            else:
                log.warning(f"{anc_path} missing 'ancestry_label' column; skipping ancestry filter")
                funnel.append(_funnel_row("ancestry_restriction_skipped", len(merged)))
        else:
            # use SQC columns merged in from step 03
            mask = _ancestry_keep_mask(merged, label)
            if mask is None:
                log.warning(f"no SQC ancestry columns and no ancestry_relatedness_file; "
                            f"skipping ancestry filter for label '{label}'")
                funnel.append(_funnel_row("ancestry_restriction_skipped", len(merged)))
            else:
                before = len(merged)
                merged = merged[mask.values].copy()
                funnel.append(_funnel_row(f"ancestry={label}_via_sqc", len(merged), before))
                log.info(f"after ancestry restriction ({label} via SQC flags): {len(merged):,}")
    else:
        funnel.append(_funnel_row("ancestry_restriction_skipped", len(merged)))
        log.info("ancestry restriction disabled in config")

    # ---- step 3: drop rows missing required covariates / dosage → LMM superset
    missing_required = [c for c in REQUIRED_COVARS if c not in merged.columns]
    if missing_required:
        log.error(f"required covariates missing from table: {missing_required}")
        return 2
    before = len(merged)
    merged = merged.dropna(subset=REQUIRED_COVARS + ["dosage_A"])
    funnel.append(_funnel_row("required_covars_nonmissing", len(merged), before))
    log.info(f"after dropping missing required covars/dosage: {len(merged):,} (LMM superset)")

    # Snapshot: this is the LMM-eligible sample (ancestry-filtered, has dosage + required covars).
    lmm_eids = set(int(x) for x in merged["eid"].dropna().astype(int))

    # ---- step 4: KING relatedness pruning → OLS subset
    rel_path = inputs.get("relatedness_file") or ""
    kin_thresh = float(cfg["sample"].get("kinship_threshold", 0.0884))
    ols_eids = set(lmm_eids)  # default = no pruning
    if cfg["sample"].get("restrict_to_unrelated", False):
        if rel_path and Path(rel_path).is_file():
            keep = _king_unrelated_keep(merged["eid"], rel_path, kin_thresh, log)
            ols_eids = set(int(x) for x in keep)
            funnel.append(_funnel_row(f"ols_unrelated_king≥{kin_thresh}", len(ols_eids),
                                       len(lmm_eids)))
            log.info(f"OLS sample after KING prune: {len(ols_eids):,}")
        elif "used_in_pca_calculation" in merged.columns:
            sub = merged[merged["used_in_pca_calculation"].astype(str).str.upper().isin(("TRUE", "1", "T"))]
            ols_eids = set(int(x) for x in sub["eid"])
            funnel.append(_funnel_row("ols_unrelated_sqc_pca_set", len(ols_eids), len(lmm_eids)))
            log.info(f"OLS sample after SQC PCA-set filter: {len(ols_eids):,}")
        elif anc_path and Path(anc_path).is_file():
            anc2 = safe_read_table(anc_path)
            anc2 = standardize_eid_column(anc2)
            if "kin_keep" in anc2.columns:
                keep_ids = anc2.loc[pd.to_numeric(anc2["kin_keep"], errors="coerce") == 1, "eid"]
                ols_eids = set(int(x) for x in keep_ids if int(x) in lmm_eids)
                funnel.append(_funnel_row("ols_unrelated_via_file", len(ols_eids), len(lmm_eids)))
                log.info(f"OLS sample after relatedness file: {len(ols_eids):,}")
            else:
                log.warning(f"{anc_path} missing 'kin_keep'; OLS sample = LMM superset")
                funnel.append(_funnel_row("ols_unrelated_skipped", len(ols_eids), len(lmm_eids)))
        else:
            log.warning("no relatedness source; OLS sample = LMM superset")
            funnel.append(_funnel_row("ols_unrelated_skipped", len(ols_eids), len(lmm_eids)))
    else:
        funnel.append(_funnel_row("ols_unrelated_skipped", len(ols_eids), len(lmm_eids)))

    # ---- step 5: MZ-twin filter for LMM keep list (drop_kinship_above)
    lmm_keep_eids = set(lmm_eids)
    drop_high = float(cfg.get("lmm", {}).get("drop_kinship_above", 0.354))
    if rel_path and Path(rel_path).is_file() and drop_high > 0:
        # _king_unrelated_keep returns the KEEP set after greedy prune at the
        # given threshold; any eid in lmm_eids not in that set was dropped.
        kept_after_mz_prune = _king_unrelated_keep(merged["eid"], rel_path, drop_high, log)
        dropped = lmm_eids - set(int(x) for x in kept_after_mz_prune)
        if dropped:
            log.info(f"LMM keep list: dropping {len(dropped):,} eids with kinship ≥ {drop_high}")
            lmm_keep_eids = lmm_eids - dropped
        funnel.append(_funnel_row(f"lmm_keep_kinship<{drop_high}", len(lmm_keep_eids), len(lmm_eids)))
    else:
        funnel.append(_funnel_row("lmm_keep_no_mz_filter", len(lmm_keep_eids), len(lmm_eids)))

    # ---- step 6: variant-wide QC
    if qc_path.is_file():
        qc = pd.read_csv(qc_path)
        info = qc["info_mfi"].iloc[0] if "info_mfi" in qc.columns and len(qc) else np.nan
        miss_pct = qc["missing_pct"].iloc[0] if "missing_pct" in qc.columns and len(qc) else np.nan
        log.info(f"variant QC: INFO={info}, missing_pct={miss_pct}")
        if np.isfinite(info) and info < cfg["variant"]["min_info"] and not args.force:
            log.error(f"variant INFO {info} < min_info {cfg['variant']['min_info']}; "
                      f"pass --force to override.")
            return 3

    # ---- write outputs
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False)
    log.info(f"wrote {out_path} (LMM superset: {len(merged):,} rows, {len(merged.columns)} cols)")

    # Two-column FID IID keep lists (PLINK/REGENIE/BOLT format).
    def _write_keep(path: Path, eids: set) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as fh:
            for e in sorted(eids):
                fh.write(f"{e}\t{e}\n")

    _write_keep(keep_ols_path, ols_eids)
    _write_keep(keep_lmm_path, lmm_keep_eids)
    log.info(f"wrote {keep_ols_path} ({len(ols_eids):,} eids)")
    log.info(f"wrote {keep_lmm_path} ({len(lmm_keep_eids):,} eids)")

    counts_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(funnel).to_csv(counts_path, index=False)
    log.info(f"wrote {counts_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
