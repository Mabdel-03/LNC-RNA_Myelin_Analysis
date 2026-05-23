#!/usr/bin/env python
"""Step 01 — extract rs2546890 dosage and write a clean per-eid table.

Extraction priority:
  1) inputs.preextracted_variant_dosage_file (validate columns).
  2) inputs.pgen_prefix + plink2 (PRIMARY on this machine).
  3) inputs.chr5_bgen + inputs.sample_file + (plink2|qctool|bgenix).
  4) None available → write templated commands to logs/ and exit 1
     (unless dry_run).

Outputs:
  results/intermediate/variant_dosage.tsv  (cols: eid, dosage_A, [genotype])
  results/genotype_qc_summary.csv
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    dosage_to_hardcall,
    load_config,
    resolve_under_repo,
    run_cmd,
    safe_read_table,
    setup_logger,
    standardize_eid_column,
    tool_available,
)


# ---------------------------------------------------------------------------
# Branch 1: pre-extracted CSV
# ---------------------------------------------------------------------------

def _load_pre_extracted(path: str, effect_allele: str, other_allele: str,
                         logger) -> pd.DataFrame:
    logger.info(f"reading pre-extracted dosage from {path}")
    df = safe_read_table(path)
    df = standardize_eid_column(df)
    if "dosage_A" not in df.columns and "dosage" in df.columns:
        df = df.rename(columns={"dosage": "dosage_A"})
    if "dosage_A" not in df.columns:
        raise ValueError(f"pre-extracted file must contain 'dosage_A' (or 'dosage'); "
                         f"found {list(df.columns)}")
    df["dosage_A"] = pd.to_numeric(df["dosage_A"], errors="coerce")
    keep = ["eid", "dosage_A"]
    if "genotype" in df.columns:
        keep.append("genotype")
    df = df[keep].dropna(subset=["eid"])
    df["source"] = "preextracted"
    return df


# ---------------------------------------------------------------------------
# Branch 2: plink2 on pgen (PRIMARY)
# ---------------------------------------------------------------------------

def _slurm_threads(default: int = 1) -> int:
    raw = (
        os.environ.get("PLINK2_THREADS")
        or os.environ.get("SLURM_CPUS_PER_TASK")
        or os.environ.get("SLURM_CPUS_ON_NODE")
        or os.environ.get("SLURM_NTASKS")
        or str(default)
    )
    try:
        return max(1, int(raw))
    except ValueError:
        return max(1, default)


def _run_plink2_pgen(pgen_prefix: str, rsid: str, chrom: int, pos: int,
                     out_prefix: Path, logger,
                     effect_allele: str = "A", other_allele: str = "G") -> tuple[Path, str]:
    """Try rsID → chr:pos:ref:alt ID → coord range. Returns (raw_path, allele_modeled)."""
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    base_cmd = ["plink2", "--threads", str(_slurm_threads()), "--pfile", pgen_prefix]
    if Path(pgen_prefix + ".pvar.zst").is_file():
        base_cmd += ["vzs"]
    # candidate variant IDs to try with --snp (UKB-style colon format included)
    snp_candidates = [
        rsid,
        f"{chrom}:{pos}:{effect_allele}:{other_allele}",
        f"{chrom}:{pos}:{other_allele}:{effect_allele}",
        f"chr{chrom}:{pos}:{effect_allele}:{other_allele}",
        f"chr{chrom}:{pos}:{other_allele}:{effect_allele}",
    ]
    last_exc = None
    for snp in snp_candidates:
        cmd = base_cmd + ["--snp", snp, "--export", "A", "--out", str(out_prefix)]
        try:
            run_cmd(cmd, logger=logger, check=True)
            logger.info(f"matched variant by ID '{snp}'")
            break
        except RuntimeError as exc:
            last_exc = exc
            logger.info(f"  '{snp}' not found; trying next")
    else:
        # coordinate range fallback
        logger.warning(f"all ID candidates failed ({last_exc}); falling back to chr:bp range")
        cmd = base_cmd + ["--chr", str(chrom), "--from-bp", str(pos),
                          "--to-bp", str(pos), "--export", "A",
                          "--out", str(out_prefix)]
        run_cmd(cmd, logger=logger, check=True)
    raw = Path(str(out_prefix) + ".raw")
    if not raw.is_file():
        raise FileNotFoundError(f"plink2 did not produce {raw}")
    # Parse header to find the dosage column. plink2 --export A names columns
    # like "rsid_<countedAllele>".
    with raw.open("r") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    dosage_cols = [c for c in header if "_" in c and c not in (
        "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    # the first allele-suffix column corresponds to our variant
    allele_modeled = None
    target_col = None
    for c in dosage_cols:
        # last underscore separates allele
        allele = c.rsplit("_", 1)[-1]
        if len(allele) == 1 and allele.upper() in ("A", "C", "G", "T"):
            allele_modeled = allele.upper()
            target_col = c
            break
    if target_col is None:
        raise RuntimeError(f"could not find dosage column in {raw} header: {header}")
    logger.info(f"plink2 .raw dosage column: {target_col} (allele={allele_modeled})")
    return raw, allele_modeled or "?"


def _load_plink2_raw(raw_path: Path, allele_modeled: str,
                     effect_allele: str, other_allele: str,
                     hardcall_threshold: float,
                     logger) -> pd.DataFrame:
    df = pd.read_csv(raw_path, sep="\t")
    df = standardize_eid_column(df, prefer="IID")
    # find dosage column again here (header just re-read)
    dosage_cols = [c for c in df.columns if c.endswith(f"_{allele_modeled}")
                   and c not in ("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
    if not dosage_cols:
        raise RuntimeError(f"no dosage column ending '_{allele_modeled}' in {raw_path}")
    dose_col = dosage_cols[0]
    dose = pd.to_numeric(df[dose_col], errors="coerce")
    if allele_modeled == effect_allele:
        dosage_A = dose
    elif allele_modeled == other_allele:
        logger.warning(f"plink2 counted '{allele_modeled}'; flipping to count '{effect_allele}'")
        dosage_A = 2.0 - dose
    else:
        raise RuntimeError(
            f"plink2 counted allele '{allele_modeled}' is neither effect "
            f"'{effect_allele}' nor other '{other_allele}'. Check variant config.")
    out = pd.DataFrame({"eid": df["eid"], "dosage_A": dosage_A})
    out["genotype"] = dosage_to_hardcall(out["dosage_A"], effect_allele,
                                         other_allele, hardcall_threshold)
    out["source"] = "plink2_pgen"
    return out


# ---------------------------------------------------------------------------
# Branch 3: BGEN fallback (templated)
# ---------------------------------------------------------------------------

def _try_bgen(chr5_bgen: str, sample_file: str, rsid: str,
              chrom: int, pos: int, out_prefix: Path, logger,
              effect_allele: str, other_allele: str,
              hardcall_threshold: float) -> pd.DataFrame:
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    if tool_available("plink2"):
        cmd = ["plink2", "--threads", str(_slurm_threads()),
               "--bgen", chr5_bgen, "ref-first",
               "--sample", sample_file, "--snp", rsid,
               "--export", "A", "--out", str(out_prefix)]
        try:
            run_cmd(cmd, logger=logger, check=True)
        except RuntimeError as exc:
            logger.warning(f"rsID bgen extraction failed ({exc}); retry by coord")
            cmd = ["plink2", "--threads", str(_slurm_threads()),
                   "--bgen", chr5_bgen, "ref-first",
                   "--sample", sample_file,
                   "--chr", str(chrom), "--from-bp", str(pos), "--to-bp", str(pos),
                   "--export", "A", "--out", str(out_prefix)]
            run_cmd(cmd, logger=logger, check=True)
        raw = Path(str(out_prefix) + ".raw")
        # reuse plink2 raw reader
        with raw.open("r") as fh:
            header = fh.readline().rstrip("\n").split("\t")
        dosage_cols = [c for c in header if "_" in c and c not in (
            "FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")]
        target = dosage_cols[0]
        allele_modeled = target.rsplit("_", 1)[-1].upper()
        return _load_plink2_raw(raw, allele_modeled, effect_allele,
                                other_allele, hardcall_threshold, logger)
    raise RuntimeError("BGEN provided but plink2 not on PATH (qctool/bgenix paths "
                       "would require an extra parser; not implemented).")


# ---------------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------------

def _lookup_mfi_info(mfi_path: str | None, chrom: int, pos: int, rsid: str,
                     logger) -> dict:
    """Best-effort: grep the MFI for our variant and return INFO + allele freq.

    Supports two MFI schemas:
      (a) UKB v3 standard (8 cols): chr rsid pos a1 a2 maf minor_allele info
      (b) Tanigawa local (5 cols):  ID(chr:pos:a1:a2)  UKB_VAR_ID  ORIGINAL_VAR_ID  AF_A1  INFO
    Detected by inspecting the first non-comment line.
    """
    out = {"info": np.nan, "mfi_a1": None, "mfi_a2": None, "mfi_maf": np.nan}
    if not mfi_path or not Path(mfi_path).is_file():
        return out
    try:
        opener = f"zstdcat {mfi_path}" if mfi_path.endswith(".zst") else f"cat {mfi_path}"
        # First detect schema by peeking
        peek = subprocess.run(f"{opener} | head -5", shell=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        peek_lines = [l for l in (peek.stdout or "").splitlines() if l and not l.startswith("#")]
        schema = "unknown"
        if peek_lines:
            first = peek_lines[0]
            ncol = len(first.split("\t")) if "\t" in first else len(first.split())
            schema = "tanigawa_5col" if ncol == 5 else ("ukb_v3_8col" if ncol >= 8 else "unknown")
            logger.info(f"MFI schema detected: {schema} (ncols={ncol})")
        # Build awk query
        coord_id = f"{chrom}:{pos}"
        if schema == "tanigawa_5col":
            # match either first col ID prefix (chr:pos:...) or rsID anywhere
            awk = f"$1 ~ /^{chrom}:{pos}:/ || $2==\"{rsid}\" || $3==\"{rsid}\""
        else:
            awk = f"($1==\"{chrom}\" && $3=={pos}) || $2==\"{rsid}\""
        proc = subprocess.run(f"{opener} | awk -F'\\t' '{awk}'",
                              shell=True, text=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        line = (proc.stdout or "").strip().splitlines()
        if not line:
            logger.warning(f"variant not found in MFI {mfi_path} (coord={coord_id}, rsid={rsid})")
            return out
        parts = line[0].split("\t") if "\t" in line[0] else line[0].split()
        if schema == "tanigawa_5col" and len(parts) >= 5:
            # parts[0] = "5:158759900:A:G"
            id_parts = parts[0].split(":")
            if len(id_parts) >= 4:
                out["mfi_a1"] = id_parts[2]
                out["mfi_a2"] = id_parts[3]
            try:
                out["mfi_maf"] = float(parts[3])
            except Exception:
                pass
            try:
                out["info"] = float(parts[4])
            except Exception:
                pass
        elif len(parts) >= 8:
            out["mfi_a1"] = parts[3]
            out["mfi_a2"] = parts[4]
            try:
                out["mfi_maf"] = float(parts[5])
            except Exception:
                pass
            try:
                out["info"] = float(parts[7])
            except Exception:
                pass
    except Exception as exc:
        logger.warning(f"MFI lookup failed: {exc}")
    return out


def _hwe_p(counts_aa: int, counts_ab: int, counts_bb: int) -> float:
    """Standard exact-style chi-square HWE; returns NaN if any cell is 0 total."""
    n = counts_aa + counts_ab + counts_bb
    if n == 0:
        return float("nan")
    p = (2 * counts_aa + counts_ab) / (2 * n)
    q = 1 - p
    exp_aa = n * p * p
    exp_ab = 2 * n * p * q
    exp_bb = n * q * q
    chi2 = 0.0
    for o, e in ((counts_aa, exp_aa), (counts_ab, exp_ab), (counts_bb, exp_bb)):
        if e > 0:
            chi2 += (o - e) ** 2 / e
    from scipy.stats import chi2 as _c2
    return float(1.0 - _c2.cdf(chi2, df=1))


def _write_qc(df: pd.DataFrame, cfg: dict, mfi_info: dict, out_csv: Path,
              extraction_method: str) -> None:
    var = cfg["variant"]
    ea, oa = var["effect_allele"], var["other_allele"]
    aa_lab, ag_lab, gg_lab = ea + ea, "".join(sorted([ea, oa])), oa + oa
    dose = df["dosage_A"].astype(float)
    n_total = len(df)
    n_miss = int(dose.isna().sum())
    counts = df["genotype"].value_counts(dropna=True).to_dict() if "genotype" in df else {}
    n_aa = int(counts.get(aa_lab, 0))
    n_ag = int(counts.get(ag_lab, 0))
    n_gg = int(counts.get(gg_lab, 0))
    hwe = _hwe_p(n_aa, n_ag, n_gg) if (n_aa + n_ag + n_gg) > 0 else float("nan")
    eaf = float(dose.mean(skipna=True) / 2.0) if dose.notna().any() else float("nan")
    row = {
        "rsid": var["rsid"],
        "chr": var["chr"],
        "pos_grch38": var["pos_grch38"],
        "effect_allele": ea,
        "other_allele": oa,
        "extraction_method": extraction_method,
        "dosage_mean": float(dose.mean(skipna=True)) if dose.notna().any() else float("nan"),
        "dosage_sd": float(dose.std(skipna=True, ddof=1)) if dose.notna().sum() > 1 else float("nan"),
        "eaf_estimated": eaf,
        "missing_n": n_miss,
        "missing_pct": (n_miss / n_total * 100.0) if n_total else float("nan"),
        "n_total": n_total,
        f"hardcall_{aa_lab}": n_aa,
        f"hardcall_{ag_lab}": n_ag,
        f"hardcall_{gg_lab}": n_gg,
        "hwe_p": hwe,
        "info_mfi": mfi_info.get("info"),
        "mfi_a1": mfi_info.get("mfi_a1"),
        "mfi_a2": mfi_info.get("mfi_a2"),
        "mfi_maf": mfi_info.get("mfi_maf"),
    }
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([row]).to_csv(out_csv, index=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true",
                        help="check tools/paths only; do not run plink2")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    variant_dir = results_dir / "variant"
    variant_dir.mkdir(parents=True, exist_ok=True)
    out_dose = inter_dir / "variant_dosage.tsv"
    out_qc = results_dir / "genotype_qc_summary.csv"

    log = setup_logger("extract_variant", logs_dir / "run_log.txt")
    inputs = cfg["inputs"]
    var = cfg["variant"]

    # decide branch
    pre = inputs.get("preextracted_variant_dosage_file") or ""
    pgen = inputs.get("pgen_prefix") or ""
    bgen = inputs.get("chr5_bgen") or ""
    sample = inputs.get("sample_file") or ""
    plink2_path = tool_available("plink2")

    method = None
    if pre and Path(pre).is_file():
        method = "preextracted"
    elif pgen and Path(pgen + ".pgen").is_file() and plink2_path:
        method = "plink2_pgen"
    elif bgen and sample and Path(bgen).is_file() and Path(sample).is_file() and (
            plink2_path or tool_available("qctool") or tool_available("bgenix")):
        method = "bgen_fallback"
    else:
        method = "none"

    log.info(f"selected extraction branch: {method}")
    if dry:
        log.info("DRY RUN — not executing extraction; writing template commands only")
        # Emit the would-be commands so the user can run them manually if needed.
        cmd_path = logs_dir / "variant_extraction_commands.sh"
        lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
        if pgen:
            vzs = " vzs" if Path(pgen + ".pvar.zst").is_file() else ""
            lines.append(f"# Primary: plink2 on pgen")
            lines.append(f"plink2 --pfile {pgen}{vzs} \\")
            lines.append(f"       --snp {var['rsid']} \\")
            lines.append(f"       --export A \\")
            lines.append(f"       --out {variant_dir}/{var['rsid']}")
            lines.append("")
            lines.append(f"# Coordinate fallback (if rsID miss)")
            lines.append(f"plink2 --pfile {pgen}{vzs} \\")
            lines.append(f"       --chr {var['chr']} --from-bp {var['pos_grch38']} "
                          f"--to-bp {var['pos_grch38']} \\")
            lines.append(f"       --export A \\")
            lines.append(f"       --out {variant_dir}/{var['rsid']}_bycoord")
            lines.append("")
        if bgen and sample:
            lines.append(f"# BGEN fallback")
            lines.append(f"plink2 --bgen {bgen} ref-first --sample {sample} \\")
            lines.append(f"       --snp {var['rsid']} --export A \\")
            lines.append(f"       --out {variant_dir}/{var['rsid']}_bgen")
            lines.append("")
        cmd_path.write_text("\n".join(lines) + "\n")
        log.info(f"wrote {cmd_path}")
        return 0

    if method == "none":
        log.error("no viable extraction source. Fill one of: "
                  "inputs.preextracted_variant_dosage_file, "
                  "inputs.pgen_prefix (+plink2 on PATH), "
                  "or inputs.chr5_bgen + inputs.sample_file (+plink2|qctool|bgenix).")
        return 1

    # Run the chosen branch
    if method == "preextracted":
        df = _load_pre_extracted(pre, var["effect_allele"], var["other_allele"], log)
        # ensure genotype if absent
        if "genotype" not in df.columns:
            df["genotype"] = dosage_to_hardcall(
                df["dosage_A"], var["effect_allele"], var["other_allele"],
                var["hardcall_probability_threshold"])
    elif method == "plink2_pgen":
        out_prefix = variant_dir / var["rsid"]
        build = str(var.get("pgen_build", "GRCh38")).upper().replace("HG", "GRCH")
        if build in ("GRCH37", "HG19"):
            pos = int(var.get("pos_grch37") or var["pos_grch38"])
            log.info(f"pgen_build={build}; using pos_grch37={pos}")
        else:
            pos = int(var["pos_grch38"])
            log.info(f"pgen_build={build}; using pos_grch38={pos}")
        raw, allele_modeled = _run_plink2_pgen(
            pgen, var["rsid"], var["chr"], pos, out_prefix, log,
            effect_allele=var["effect_allele"], other_allele=var["other_allele"])
        df = _load_plink2_raw(raw, allele_modeled, var["effect_allele"],
                              var["other_allele"], var["hardcall_probability_threshold"], log)
    elif method == "bgen_fallback":
        out_prefix = variant_dir / (var["rsid"] + "_bgen")
        df = _try_bgen(bgen, sample, var["rsid"], var["chr"], var["pos_grch38"],
                       out_prefix, log, var["effect_allele"], var["other_allele"],
                       var["hardcall_probability_threshold"])
    else:
        raise RuntimeError(f"unhandled method {method}")

    # write outputs
    out_dose.parent.mkdir(parents=True, exist_ok=True)
    cols = ["eid", "dosage_A"] + (["genotype"] if "genotype" in df.columns else [])
    df[cols].to_csv(out_dose, sep="\t", index=False)
    log.info(f"wrote {out_dose} ({len(df):,} rows)")

    # MFI is on the same build as the pgen it accompanies; use that build's coord.
    build = str(var.get("pgen_build", "GRCh38")).upper().replace("HG", "GRCH")
    mfi_pos = int(var.get("pos_grch37") or var["pos_grch38"]) if build in ("GRCH37", "HG19") \
              else int(var["pos_grch38"])
    mfi_info = _lookup_mfi_info(inputs.get("variant_metadata_file"),
                                 var["chr"], mfi_pos, var["rsid"], log)
    _write_qc(df, cfg, mfi_info, out_qc, method)
    log.info(f"wrote {out_qc}")

    info = mfi_info.get("info")
    if info is not None and np.isfinite(info) and info < cfg["variant"]["min_info"]:
        log.warning(f"INFO={info:.3f} < min_info={cfg['variant']['min_info']}; "
                    f"04_merge_qc_sample.py will halt unless threshold is relaxed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
