#!/usr/bin/env python
"""Step 00 — validate every input declared in config.yaml.

Writes logs/inspect_report.txt with a per-path status (FOUND / MISSING /
EMPTY), per-tool availability (plink2, qctool, bgenix, regenie), and
Python package versions. Exits 0 unless an input the pipeline cannot run
without is missing (in dry-run mode it always exits 0).
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# ensure repo `src/` on path so plain `python src/00_*.py` works
sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    resolve_under_repo,
    setup_logger,
    tool_available,
    write_version_log,
)


def _status(p: str | None) -> str:
    if not p:
        return "NOT_SET"
    pp = Path(p)
    if not pp.exists():
        return "MISSING"
    if pp.is_file() and pp.stat().st_size == 0:
        return "EMPTY"
    if pp.is_dir():
        try:
            next(iter(pp.iterdir()))
        except StopIteration:
            return "DIR_EMPTY"
        return "DIR_OK"
    return "FOUND"


def _check_pgen_trio(prefix: str | None) -> dict:
    """A pgen 'prefix' implies .pgen + .psam + .pvar (or .pvar.zst)."""
    out = {}
    if not prefix:
        out["pgen_prefix"] = "NOT_SET"
        return out
    for suf in (".pgen", ".psam", ".pvar", ".pvar.zst"):
        out[prefix + suf] = _status(prefix + suf)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="path to config.yaml")
    parser.add_argument("--report", default=None,
                        help="override path for inspect_report.txt")
    parser.add_argument("--dry-run", action="store_true",
                        help="(accepted for orchestrator parity; inspect always runs fully)")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    repo = Path(cfg["_repo_root"])
    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    logs_dir.mkdir(parents=True, exist_ok=True)
    report_path = Path(args.report) if args.report else (logs_dir / "inspect_report.txt")

    log = setup_logger("inspect", logs_dir / "run_log.txt")
    log.info(f"repo root: {repo}")
    log.info(f"config:    {args.config}")
    log.info(f"report:    {report_path}")

    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("rs2546890 → UKBB MRI IDP pipeline — input inspection")
    lines.append("=" * 72)
    lines.append("")
    lines.append(f"repo_root: {repo}")
    lines.append(f"config:    {args.config}")
    lines.append("")

    # --- single-file inputs
    inputs = cfg.get("inputs", {})
    single_file_keys = [
        "phenotype_file", "data_dictionary_file", "basket_tab",
        "extract_phenotype_sh", "variant_metadata_file",
        "chr5_bgen", "sample_file", "preextracted_variant_dosage_file",
        "genetic_covariates_file", "imaging_confounds_file",
        "ancestry_relatedness_file",
    ]
    lines.append("[inputs] single files")
    for k in single_file_keys:
        v = inputs.get(k, "")
        lines.append(f"  {k:40s}  {_status(v):10s}  {v or '(empty)'}")
    lines.append("")

    # --- directory list inputs
    lines.append("[inputs] mri_pheno_tsv_dirs")
    for d in inputs.get("mri_pheno_tsv_dirs", []) or []:
        s = _status(d)
        n_tsv = 0
        if s == "DIR_OK":
            n_tsv = sum(1 for _ in Path(d).rglob("*.tsv"))
        lines.append(f"  {s:10s} ({n_tsv:>5d} tsv)  {d}")
    lines.append("")

    # --- pgen trio
    lines.append("[inputs] pgen trio")
    for k, v in _check_pgen_trio(inputs.get("pgen_prefix", "")).items():
        lines.append(f"  {v:10s}  {k}")
    lines.append("")

    # --- tools
    lines.append("[tools]")
    for tool in ("plink2", "qctool", "bgenix", "regenie", "bcftools", "tabix"):
        loc = tool_available(tool)
        lines.append(f"  {tool:10s}  {loc or 'NOT_ON_PATH'}")
    lines.append("")

    # --- python packages snapshot (delegate to write_version_log into its own file)
    versions_path = logs_dir / "versions.txt"
    write_version_log(versions_path)
    lines.append(f"[versions] written to {versions_path}")
    lines.append("")

    # --- pre-flight decision summary
    lines.append("[pre-flight]")
    pgen_ok = (inputs.get("pgen_prefix") and
               Path(inputs["pgen_prefix"] + ".pgen").is_file() and
               Path(inputs["pgen_prefix"] + ".psam").is_file() and
               (Path(inputs["pgen_prefix"] + ".pvar").is_file() or
                Path(inputs["pgen_prefix"] + ".pvar.zst").is_file()))
    pre_ok = bool(inputs.get("preextracted_variant_dosage_file") and
                  Path(inputs["preextracted_variant_dosage_file"]).is_file())
    bgen_ok = bool(inputs.get("chr5_bgen") and inputs.get("sample_file") and
                   Path(inputs["chr5_bgen"]).is_file() and
                   Path(inputs["sample_file"]).is_file())
    plink2_ok = tool_available("plink2") is not None

    if pre_ok:
        lines.append("  variant extraction → pre-extracted CSV")
    elif pgen_ok and plink2_ok:
        lines.append("  variant extraction → plink2 on pgen (PRIMARY)")
    elif bgen_ok and (plink2_ok or tool_available("qctool") or tool_available("bgenix")):
        lines.append("  variant extraction → BGEN fallback")
    else:
        lines.append("  variant extraction → BLOCKED. provide pgen_prefix + plink2, "
                     "OR chr5_bgen + sample_file + (plink2|qctool|bgenix), "
                     "OR preextracted_variant_dosage_file.")

    pheno_sources_ok = (
        bool(inputs.get("phenotype_file") and Path(inputs["phenotype_file"]).is_file())
        or any(Path(d).is_dir() for d in (inputs.get("mri_pheno_tsv_dirs") or []))
        or bool(inputs.get("basket_tab") and Path(inputs["basket_tab"]).is_file())
    )
    lines.append("  phenotype source available: " + ("YES" if pheno_sources_ok else "NO"))

    covar_source_ok = (
        bool(inputs.get("genetic_covariates_file") and Path(inputs["genetic_covariates_file"]).is_file())
        or bool(inputs.get("basket_tab") and Path(inputs["basket_tab"]).is_file())
    )
    lines.append("  covariate source available:  " + ("YES" if covar_source_ok else "NO"))

    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    lines.append(f"  dry_run effective:            {dry}")
    lines.append("")
    lines.append("=" * 72)

    text = "\n".join(lines) + "\n"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(text)
    log.info(f"wrote {report_path}")
    # Also echo the summary block to stdout so it shows up in run_all.sh log
    print(text)

    # Exit nonzero only if not in dry-run AND a required source is missing
    blocked = (not dry) and (not pheno_sources_ok or not covar_source_ok or
                             not (pre_ok or (pgen_ok and plink2_ok) or
                                  (bgen_ok and (plink2_ok or tool_available("qctool")
                                                or tool_available("bgenix")))))
    return 1 if blocked else 0


if __name__ == "__main__":
    raise SystemExit(main())
