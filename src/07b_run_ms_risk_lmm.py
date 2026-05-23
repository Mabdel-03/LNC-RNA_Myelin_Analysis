#!/usr/bin/env python
"""Step 07b — REGENIE binary-trait LMM for MS (G35) risk on rs2546890.

Independent of the main QT pipeline (05b). Builds its own minimal
phenotypes / covariates / keep_lmm inputs in `results/lmm_inputs_ms/`,
refits REGENIE step 1 with --bt --firth --approx (Firth correction
recommended at low case fractions), and runs step 2 for rs2546890 only.

Reads:
  results/ms_phenotype.tsv                 (per-eid ms_case from step 07)
  results/intermediate/covariates.tsv      (eid + age + sex + array + site + PCs)
  results/intermediate/keep_lmm.txt        (FID IID — LMM keep set)
  config.lmm.regenie.*                     (paths to HM3 BED, model SNPs, step-2 pgen, etc.)

Writes:
  results/lmm_inputs_ms/{phenotypes,covariates,keep_lmm,pheno_col_list,
                          covar_col_quant,covar_col_cat,target_variants}.*
  results/lmm/ms_step1_*.{loco.gz,pred.list,pheno_sha256}
  results/lmm/ms_step2_ms_case.regenie.gz
  results/association_ms_risk_lmm.csv

The MS step-1 cache hash deliberately differs from the QT one (different
SHA256 over a different pheno_col_list), so the existing QT cache at
results/lmm/regenie_step1_* is untouched. The two runs can co-exist.
"""
from __future__ import annotations

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    parse_regenie_output,
    resolve_under_repo,
    run_regenie_step1,
    run_regenie_step2,
    safe_read_table,
    setup_logger,
)


def _sha256_of_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_list_file(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def _build_inputs(cfg: dict, ms_pheno: Path, covars_tsv: Path,
                  keep_lmm: Path, out_dir: Path, log) -> tuple[Path, Path, Path, Path, Path, Path, Path]:
    """Materialize the 7 REGENIE input files in out_dir for the MS BT run."""
    out_dir.mkdir(parents=True, exist_ok=True)

    ms = safe_read_table(ms_pheno)
    ms["eid"] = pd.to_numeric(ms["eid"], errors="coerce").astype("Int64")
    ms = ms.dropna(subset=["eid"])
    ms = ms[["eid", "ms_case"]].copy()

    covars = safe_read_table(covars_tsv)
    covars["eid"] = pd.to_numeric(covars["eid"], errors="coerce").astype("Int64")
    if "age2" not in covars.columns and "age" in covars.columns:
        covars["age2"] = covars["age"].astype(float) ** 2

    keep_set = set()
    with keep_lmm.open() as fh:
        for ln in fh:
            parts = ln.split()
            if len(parts) >= 2:
                try:
                    keep_set.add(int(parts[1]))
                except ValueError:
                    pass
    log.info(f"keep_lmm set size: {len(keep_set):,}")

    merged = ms.merge(covars, on="eid", how="inner")
    merged = merged[merged["eid"].astype("Int64").isin(keep_set)].copy()
    n_case = int(merged["ms_case"].sum())
    log.info(f"merged + keep_lmm: {len(merged):,} rows, {n_case:,} MS cases "
             f"({100.0 * n_case / max(len(merged), 1):.2f}%)")

    # quantitative + categorical covar lists — same convention as 04c, but
    # for MS-risk on the full WB cohort we drop imaging-only fields and
    # high-cardinality categoricals that REGENIE refuses (genotype_batch ~105
    # levels). `site` is the imaging-center field — only populated for the
    # ~40K MRI participants, would restrict the analysis.
    quant_candidates = ["age", "age2"] + [c for c in merged.columns
                                            if c.startswith("PC") and c[2:].isdigit()
                                            and int(c[2:]) <= int(cfg["covariates"].get("genetic_pc_count", 20))]
    quant_cols = [c for c in quant_candidates if c in merged.columns]
    cat_candidates = ["sex", "array"]
    cat_cols = []
    for c in cat_candidates:
        if c not in merged.columns:
            continue
        # Drop singletons — REGENIE warns and ignores them anyway
        n_levels = merged[c].dropna().nunique()
        if n_levels < 2:
            log.info(f"  dropping '{c}' (only {n_levels} unique level)")
            continue
        cat_cols.append(c)
    # Drop rows missing any required covariate (REGENIE rejects rows with NaN)
    before = len(merged)
    merged = merged.dropna(subset=quant_cols + cat_cols + ["ms_case"]).copy()
    log.info(f"after dropna on required covars: {len(merged):,} of {before:,} rows")
    n_case = int(merged["ms_case"].sum())
    log.info(f"final analysis set: {len(merged):,} rows, {n_case:,} MS cases")

    # Cast categorical covars to int for REGENIE (which expects categorical levels as integers/strings)
    for c in cat_cols:
        merged[c] = pd.to_numeric(merged[c], errors="coerce").astype("Int64").astype(str)

    pheno_file = out_dir / "phenotypes.tsv"
    covar_file = out_dir / "covariates.tsv"
    keep_file = out_dir / "keep_lmm.txt"
    pheno_list_file = out_dir / "pheno_col_list.txt"
    quant_list_file = out_dir / "covar_col_quant.txt"
    cat_list_file = out_dir / "covar_col_cat.txt"
    targets_file = out_dir / "target_variants.txt"

    pheno_out = merged[["eid", "ms_case"]].copy()
    pheno_out.insert(0, "FID", pheno_out["eid"])
    pheno_out = pheno_out.rename(columns={"eid": "IID"})
    pheno_out = pheno_out[["FID", "IID", "ms_case"]]
    pheno_out.to_csv(pheno_file, sep="\t", index=False)

    covar_out = merged[["eid"] + quant_cols + cat_cols].copy()
    covar_out.insert(0, "FID", covar_out["eid"])
    covar_out = covar_out.rename(columns={"eid": "IID"})
    covar_out = covar_out[["FID", "IID"] + quant_cols + cat_cols]
    covar_out.to_csv(covar_file, sep="\t", index=False)

    keep_out = merged[["eid"]].copy()
    keep_out.insert(0, "FID", keep_out["eid"])
    keep_out = keep_out.rename(columns={"eid": "IID"})
    keep_out.to_csv(keep_file, sep=" ", index=False, header=False)

    pheno_list_file.write_text("ms_case\n")
    quant_list_file.write_text("\n".join(quant_cols) + "\n")
    cat_list_file.write_text("\n".join(cat_cols) + "\n")
    targets_file.write_text(str(cfg["lmm"]["target_variant_id"]) + "\n")

    for p in (pheno_file, covar_file, keep_file, pheno_list_file,
              quant_list_file, cat_list_file, targets_file):
        log.info(f"wrote {p}")
    return pheno_file, covar_file, keep_file, pheno_list_file, quant_list_file, cat_list_file, targets_file


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-step1", action="store_true")
    parser.add_argument("--skip-step2", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", False))
    lmm_cfg = cfg.get("lmm", {})
    rcfg = lmm_cfg.get("regenie", {})

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    lmm_in = results_dir / "lmm_inputs_ms"
    lmm_out = results_dir / "lmm"
    lmm_out.mkdir(parents=True, exist_ok=True)

    log = setup_logger("regenie_ms", logs_dir / "run_log.txt")
    step1_log = setup_logger("regenie_step1_ms", logs_dir / "regenie_step1_ms.log")
    step2_log = setup_logger("regenie_step2_ms", logs_dir / "regenie_step2_ms.log")

    ms_pheno_path = results_dir / "ms_phenotype.tsv"
    covars_path = inter_dir / "covariates.tsv"
    keep_lmm_path = inter_dir / "keep_lmm.txt"

    for p in (ms_pheno_path, covars_path, keep_lmm_path):
        if not p.is_file():
            log.error(f"required input missing: {p}")
            log.error("  (07b expects results/ms_phenotype.tsv from step 07; run 07 first)")
            return 1

    if dry:
        log.info("DRY RUN — would build BT REGENIE inputs and refit step 1 for MS")
        return 0

    (pheno_file, covar_file, keep_file, pheno_list_file,
     quant_list_file, cat_list_file, targets_file) = _build_inputs(
        cfg, ms_pheno_path, covars_path, keep_lmm_path, lmm_in, log)

    binary = rcfg.get("binary", "")
    pheno_cols = _read_list_file(pheno_list_file)
    quant_cols = _read_list_file(quant_list_file)
    cat_cols = _read_list_file(cat_list_file)

    step1_prefix = lmm_out / "ms_step1"
    step1_pred = lmm_out / "ms_step1_pred.list"
    step1_hash = lmm_out / "ms_step1.pheno_sha256"
    step2_prefix = lmm_out / "ms_step2"
    tmp_prefix = str(lmm_out / "ms_step1_tmp")

    # Cache check
    skip_step1 = args.skip_step1
    if rcfg.get("reuse_step1", True) and step1_pred.is_file() and step1_hash.is_file():
        cur_hash = _sha256_of_text(pheno_list_file)
        if cur_hash == step1_hash.read_text().strip():
            log.info("step 1 CACHE HIT — reusing ms_step1_pred.list")
            skip_step1 = True

    if not skip_step1:
        log.info("== REGENIE STEP 1 (binary trait, --bt --firth --approx) ==")
        run_regenie_step1(
            binary=binary,
            bed_prefix=rcfg["hm3_bed_prefix"],
            pheno_file=str(pheno_file),
            pheno_col_list=pheno_cols,
            covar_file=str(covar_file),
            covar_col_list=quant_cols,
            cat_covar_list=cat_cols,
            keep_file=str(keep_file),
            extract_file=rcfg["model_snps"],
            out_prefix=str(step1_prefix),
            bsize=int(lmm_cfg.get("bsize_step1", 1000)),
            threads=int(lmm_cfg.get("threads_step1", 16)),
            tmp_prefix=tmp_prefix,
            force_qt=False,
            bt=True,
            firth=True,
            min_mac=int(rcfg.get("min_mac", 20)),
            extra_ld_library_path=str(rcfg.get("extra_ld_library_path", "")),
            logger=step1_log,
        )
        step1_hash.write_text(_sha256_of_text(pheno_list_file))
        log.info(f"step 1 done; cached pheno hash to {step1_hash}")
    else:
        log.info("step 1 skipped")

    if args.skip_step2:
        return 0

    log.info("== REGENIE STEP 2 (BT, rs2546890) ==")
    run_regenie_step2(
        binary=binary,
        pgen_prefix=rcfg["step2_pgen_prefix"],
        pheno_file=str(pheno_file),
        pheno_col_list=pheno_cols,
        covar_file=str(covar_file),
        covar_col_list=quant_cols,
        cat_covar_list=cat_cols,
        keep_file=str(keep_file),
        pred_file=str(step1_pred),
        out_prefix=str(step2_prefix),
        extract_file=str(targets_file),
        bsize=int(lmm_cfg.get("bsize_step2", 400)),
        min_mac=int(rcfg.get("min_mac", 20)),
        min_info=float(rcfg.get("min_info", 0.8)),
        threads=int(lmm_cfg.get("threads_step2", 8)),
        bt=True,
        firth=True,
        extra_ld_library_path=str(rcfg.get("extra_ld_library_path", "")),
        logger=step2_log,
    )

    # Parse the single output file
    produced = sorted(lmm_out.glob("ms_step2_ms_case.regenie*"))
    if not produced:
        log.error("step 2 produced no ms_case output")
        return 3
    res = parse_regenie_output(produced[0],
                               target_variant_id=str(cfg["lmm"]["target_variant_id"]))
    if res.empty:
        log.error(f"target variant {cfg['lmm']['target_variant_id']} not found in step 2 output")
        return 4

    row = res.iloc[0].to_dict()
    # REGENIE reports BETA per-ALLELE1 (the alt allele in the pgen), which may
    # or may not match the user's "effect_allele" from config. Detect and flip
    # the sign so the reported OR is always per-effect-allele.
    beta_raw = float(row.get("beta", np.nan))
    se = float(row.get("se", np.nan))
    p = float(row.get("p", np.nan))
    a0 = str(row.get("allele_ref", "")).upper()  # REGENIE ALLELE0
    a1 = str(row.get("allele_alt", "")).upper()  # REGENIE ALLELE1 (BETA reference)
    ea = str(cfg["variant"]["effect_allele"]).upper()
    if a1 == ea:
        beta = beta_raw            # already per-effect-allele
        orientation = f"REGENIE ALLELE1={a1} matches effect_allele={ea}; BETA used as-is"
    elif a0 == ea:
        beta = -beta_raw           # flip
        orientation = f"REGENIE ALLELE1={a1} != effect_allele={ea}; BETA sign flipped"
    else:
        beta = beta_raw
        orientation = (f"WARNING: neither REGENIE allele ({a0}/{a1}) matches "
                       f"effect_allele={ea}; reporting BETA unflipped")
        log.warning(orientation)
    log.info(orientation)
    odds_ratio = float(np.exp(beta))
    odds_lo = float(np.exp(beta - 1.959964 * se))
    odds_hi = float(np.exp(beta + 1.959964 * se))

    log.info(f"REGENIE-BT MS result: OR={odds_ratio:.3g} "
             f"[{odds_lo:.3g}-{odds_hi:.3g}] p={p:.3g} n={int(row.get('n', 0)):,}")

    out_row = {
        "phenotype": "MS_G35",
        "variant": cfg["variant"]["rsid"],
        "effect_allele": cfg["variant"]["effect_allele"],
        "other_allele": cfg["variant"]["other_allele"],
        "engine": "regenie_BT_firth",
        "model": "additive_dosage",
        "contrast": f"per-{cfg['variant']['effect_allele']}-allele OR",
        "exposure": "dosage_A",
        "n": int(row.get("n", 0)),
        "beta": beta,
        "se": se,
        "OR": odds_ratio,
        "OR_lo95": odds_lo,
        "OR_hi95": odds_hi,
        "raw_p": p,
        "eaf": float(row.get("eaf", np.nan)),
        "info": float(row.get("info", np.nan)),
        "chisq": float(row.get("chisq", np.nan)),
        "log10p": float(row.get("log10p", np.nan)),
        "status": "ok",
    }
    out_path = results_dir / "association_ms_risk_lmm.csv"
    pd.DataFrame([out_row]).to_csv(out_path, index=False)
    log.info(f"wrote {out_path}")
    print(f"\nMS-risk REGENIE-BT: OR={odds_ratio:.3g} "
          f"[{odds_lo:.3g}-{odds_hi:.3g}] raw_p={p:.3g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
