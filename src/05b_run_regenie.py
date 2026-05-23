#!/usr/bin/env python
"""Step 05b — drive REGENIE step 1 (LOCO null model) + step 2 (rs2546890).

Reuses the prior SI-loneliness REGENIE infrastructure:
  - step 1 BED:        config.lmm.regenie.hm3_bed_prefix
  - step 1 extract:    config.lmm.regenie.model_snps
  - step 2 pgen:       config.lmm.regenie.step2_pgen_prefix
  - target variant:    config.lmm.target_variant_id (rs2546890 → 5:158759900:A:G)

Step 1 caching: writes the SHA256 of pheno_col_list.txt next to the predictor
file. If the hash matches on a subsequent run AND config.lmm.regenie.reuse_step1
is true, step 1 is skipped.

Outputs:
  results/lmm/regenie_step1_pred.list                 (REGENIE step-1 manifest)
  results/lmm/regenie_step1_<pheno>.loco              (per-pheno LOCO predictions)
  results/lmm/regenie_step2_<pheno>.regenie.gz        (per-pheno association)
  logs/regenie_step1.log, logs/regenie_step2.log
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    resolve_under_repo,
    run_regenie_step1,
    run_regenie_step2,
    setup_logger,
)


def _sha256_of_text(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_list_file(path: Path) -> list[str]:
    if not path.is_file():
        return []
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def _check_regenie(binary: str, logger) -> bool:
    """Return True if the REGENIE binary launches; False on missing libmvec etc."""
    try:
        proc = subprocess.run([binary, "--help"],
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                              text=True, timeout=30)
    except FileNotFoundError:
        logger.error(f"REGENIE binary not found: {binary}")
        return False
    except Exception as exc:  # noqa: BLE001
        logger.error(f"REGENIE failed to launch: {exc}")
        return False
    if proc.returncode != 0:
        stderr = proc.stderr.strip()
        logger.error(f"REGENIE rc={proc.returncode}; stderr: {stderr}")
        if "libmvec" in stderr or "libmvec.so" in stderr:
            logger.error("Set lmm.regenie.extra_ld_library_path to a dir containing libmvec.so.1, "
                          "OR activate the conda env used for the prior SI-loneliness REGENIE work.")
        return False
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-step1", action="store_true",
                        help="force-skip step 1 even without a hash match")
    parser.add_argument("--skip-step2", action="store_true",
                        help="run step 1 only (debug)")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    lmm_cfg = cfg.get("lmm", {})
    rcfg = lmm_cfg.get("regenie", {})

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    lmm_in = results_dir / "lmm_inputs"
    lmm_out = results_dir / "lmm"
    lmm_out.mkdir(parents=True, exist_ok=True)
    log = setup_logger("regenie", logs_dir / "run_log.txt")
    step1_log = setup_logger("regenie_step1", logs_dir / "regenie_step1.log")
    step2_log = setup_logger("regenie_step2", logs_dir / "regenie_step2.log")

    pheno_file = lmm_in / "phenotypes.tsv"
    covar_file = lmm_in / "covariates.tsv"
    keep_file = lmm_in / "keep_lmm.txt"
    pheno_list_file = lmm_in / "pheno_col_list.txt"
    quant_list_file = lmm_in / "covar_col_quant.txt"
    cat_list_file = lmm_in / "covar_col_cat.txt"
    targets_file = lmm_in / "target_variants.txt"

    step1_prefix = lmm_out / "regenie_step1"
    step1_pred = lmm_out / "regenie_step1_pred.list"
    step1_hash = lmm_out / "regenie_step1.pheno_sha256"
    step2_prefix = lmm_out / "regenie_step2"
    tmp_prefix = str(lmm_out / "regenie_step1_tmp")

    binary = rcfg.get("binary", "")
    if dry:
        log.info("DRY RUN — REGENIE invocation planning only")
        for p in (pheno_file, covar_file, keep_file, pheno_list_file,
                  quant_list_file, cat_list_file, targets_file):
            log.info(f"  {p}: {'OK' if p.is_file() else 'MISSING'}")
        log.info(f"  binary: {binary}")
        return 0

    for p in (pheno_file, covar_file, keep_file, pheno_list_file,
              quant_list_file, cat_list_file, targets_file):
        if not p.is_file():
            log.error(f"required REGENIE input missing: {p}")
            return 1

    if not _check_regenie(binary, log):
        return 2

    pheno_cols = _read_list_file(pheno_list_file)
    quant_cols = _read_list_file(quant_list_file)
    cat_cols = _read_list_file(cat_list_file)
    if not pheno_cols:
        log.error(f"empty phenotype list: {pheno_list_file}")
        return 1
    log.info(f"REGENIE will run on {len(pheno_cols)} phenotypes "
              f"(quant covars={len(quant_cols)}, cat={len(cat_cols)})")

    # ----- STEP 1: cached?
    skip_step1 = args.skip_step1 or bool(rcfg.get("skip_step1", False))
    if not skip_step1 and rcfg.get("reuse_step1", True) and step1_pred.is_file() and step1_hash.is_file():
        cur_hash = _sha256_of_text(pheno_list_file)
        cached_hash = step1_hash.read_text().strip()
        if cur_hash == cached_hash:
            log.info(f"step 1 CACHE HIT — reusing {step1_pred}")
            skip_step1 = True
        else:
            log.info("step 1 cache miss — phenotype list changed; re-running")

    if not skip_step1:
        log.info("== REGENIE STEP 1 ==")
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
            force_qt=True,
            min_mac=int(rcfg.get("min_mac", 20)),
            extra_ld_library_path=str(rcfg.get("extra_ld_library_path", "")),
            logger=step1_log,
        )
        # Cache the phenotype-list hash so the next run can reuse step 1
        step1_hash.write_text(_sha256_of_text(pheno_list_file))
        log.info(f"step 1 done; cached pheno hash to {step1_hash}")
    else:
        log.info("step 1 skipped — using existing pred.list")

    if args.skip_step2:
        log.info("--skip-step2 set; exiting after step 1")
        return 0

    # ----- STEP 2
    log.info("== REGENIE STEP 2 ==")
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
        extra_ld_library_path=str(rcfg.get("extra_ld_library_path", "")),
        logger=step2_log,
    )

    # Summarize what landed on disk
    produced = sorted(lmm_out.glob("regenie_step2_*.regenie*"))
    log.info(f"step 2 produced {len(produced)} output files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
