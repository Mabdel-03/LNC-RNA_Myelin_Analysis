#!/usr/bin/env python
"""Step 05c — drive BOLT-LMM as an optional secondary engine.

Reuses the prior SI-loneliness BOLT-LMM infrastructure:
  - --bfile:           config.lmm.bolt.hm3_bed_prefix (defaults to lmm.regenie.hm3_bed_prefix)
  - --modelSnps:       config.lmm.bolt.model_snps      (defaults to lmm.regenie.model_snps)
  - --LDscoresFile:    config.lmm.bolt.ld_scores
  - --geneticMapFile:  config.lmm.bolt.genetic_map

BOLT processes one phenotype per invocation. By default we run only composite
phenotypes + ROI PCs (15 phenotypes ≈ a few hours wall-clock). Single IDPs
are skipped unless config.lmm.bolt.include_single_idps: true.

rs2546890 lookup strategy:
  - Default (Option A): BOLT outputs SNP stats for the HM3 BED variants only.
    rs2546890 IS in HapMap3 (1000G CEU MAF~50%), so we pluck it from the
    .stats file by rsid or by 5:158759900.
  - Option B (lmm.bolt.regen_bgen: true): materialize a chr5 BGEN from the
    pgen via plink2, pass --bgenFile, and read .bgen.stats. Not implemented
    here — emit a clear instructional message instead, since this takes ~1 hr
    of plink2 conversion.

Outputs:
  results/lmm/bolt_<pheno>.stats
  logs/bolt_<pheno>.log
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    resolve_under_repo,
    run_cmd,
    setup_logger,
)


def _read_list_file(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    lmm_cfg = cfg.get("lmm", {})
    bcfg = lmm_cfg.get("bolt", {})

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    lmm_in = results_dir / "lmm_inputs"
    lmm_out = results_dir / "lmm"
    lmm_out.mkdir(parents=True, exist_ok=True)
    log = setup_logger("bolt", logs_dir / "run_log.txt")

    pheno_file = lmm_in / "phenotypes.tsv"
    covar_file = lmm_in / "covariates.tsv"
    keep_file = lmm_in / "keep_lmm.txt"
    pheno_list_file = lmm_in / "pheno_col_list.txt"
    quant_list_file = lmm_in / "covar_col_quant.txt"
    cat_list_file = lmm_in / "covar_col_cat.txt"
    ext_manifest = inter_dir / "phenotype_manifest_extended.csv"

    binary = bcfg.get("binary", "")
    bfile = bcfg.get("hm3_bed_prefix") or lmm_cfg.get("regenie", {}).get("hm3_bed_prefix", "")
    model_snps = bcfg.get("model_snps") or lmm_cfg.get("regenie", {}).get("model_snps", "")
    ld_scores = bcfg.get("ld_scores", "")
    genetic_map = bcfg.get("genetic_map", "")
    threads = int(bcfg.get("threads", 16))

    if dry:
        log.info("DRY RUN — BOLT-LMM invocation planning only")
        for p in (pheno_file, covar_file, keep_file, pheno_list_file,
                  quant_list_file, cat_list_file):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        log.info(f"  binary: {binary} ({'OK' if Path(binary).is_file() else 'MISSING'})")
        log.info(f"  bfile: {bfile}.bed ({'OK' if Path(bfile + '.bed').is_file() else 'MISSING'})")
        return 0

    for p in (pheno_file, covar_file, keep_file, pheno_list_file,
              quant_list_file, cat_list_file):
        if not p.is_file():
            log.error(f"required BOLT input missing: {p}")
            return 1

    if not Path(binary).is_file():
        log.error(f"BOLT binary not found: {binary}")
        return 2

    if bcfg.get("regen_bgen", False):
        log.error("lmm.bolt.regen_bgen=true is templated but not implemented in 05c. "
                  "Use plink2 manually: "
                  f"  plink2 --pfile {lmm_cfg.get('regenie',{}).get('step2_pgen_prefix','')} vzs "
                  "--chr 5 --export bgen-1.2 ref-first --out chr5_bgen "
                  "then pass --bgenFile=chr5_bgen.bgen --sampleFile=chr5_bgen.sample to BOLT.")
        return 3

    pheno_cols_all = _read_list_file(pheno_list_file)
    quant_cols = _read_list_file(quant_list_file)
    cat_cols = _read_list_file(cat_list_file)

    # Scope: composites + ROI PCs by default, optionally single IDPs too
    include_single = bool(bcfg.get("include_single_idps", False))
    if not include_single and ext_manifest.is_file():
        m = pd.read_csv(ext_manifest)
        scope_cols = set(m.loc[m["source"].isin(["composite", "roi_pca"]), "column_name"])
        pheno_cols = [c for c in pheno_cols_all if c in scope_cols]
        log.info(f"BOLT scope: composites + ROI PCs only → {len(pheno_cols)} phenotypes "
                 f"(out of {len(pheno_cols_all)}). Set lmm.bolt.include_single_idps=true to expand.")
    else:
        pheno_cols = pheno_cols_all
        log.info(f"BOLT scope: ALL {len(pheno_cols)} phenotypes (single IDPs included)")

    if not pheno_cols:
        log.warning("no BOLT phenotypes after scope filter; nothing to do")
        return 0

    # Build the --remove file (BED eids minus keep_lmm eids).
    fam_path = Path(bfile + ".fam")
    if not fam_path.is_file():
        log.error(f"BED .fam not found at {fam_path}")
        return 4
    fam = pd.read_csv(fam_path, sep=r"\s+", header=None,
                       names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
    keep_df = pd.read_csv(keep_file, sep=r"\s+", header=None, names=["FID", "IID"])
    keep_set = set(int(x) for x in keep_df["IID"])
    remove_df = fam[~fam["IID"].astype(int).isin(keep_set)][["FID", "IID"]]
    remove_path = lmm_in / "bolt_remove.txt"
    remove_df.to_csv(remove_path, sep="\t", index=False, header=False)
    log.info(f"BOLT --remove list: {len(remove_df):,} eids (BED total {len(fam):,}, "
              f"kept {len(keep_set):,})")

    for pheno in pheno_cols:
        out_stats = lmm_out / f"bolt_{pheno}.stats"
        cmd = [
            binary,
            f"--bfile={bfile}",
            f"--phenoFile={pheno_file}",
            f"--phenoCol={pheno}",
            f"--covarFile={covar_file}",
            "--covarMaxLevels=30",
            "--lmm",
            "--LDscoresMatchBp",
            f"--LDscoresFile={ld_scores}",
            f"--geneticMapFile={genetic_map}",
            f"--modelSnps={model_snps}",
            f"--remove={remove_path}",
            f"--numThreads={threads}",
            f"--statsFile={out_stats}",
        ]
        for q in quant_cols:
            cmd.append(f"--qCovarCol={q}")
        for c in cat_cols:
            cmd.append(f"--covarCol={c}")
        log.info(f"BOLT → {pheno}")
        bolt_log = setup_logger(f"bolt_{pheno}", logs_dir / f"bolt_{pheno}.log")
        try:
            run_cmd(cmd, logger=bolt_log, check=True)
        except RuntimeError as exc:
            log.error(f"BOLT failed for {pheno}: {exc}")
            continue
        log.info(f"  → {out_stats}")
    log.info(f"BOLT-LMM done; wrote {sum(1 for _ in lmm_out.glob('bolt_*.stats'))} stats files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
