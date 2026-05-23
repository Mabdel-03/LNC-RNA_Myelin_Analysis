#!/usr/bin/env python
"""Step 03 — build the covariate table for downstream OLS.

Pulls (in order of preference):
  - inputs.genetic_covariates_file (a pre-built TSV/CSV).
  - inputs.basket_tab columns identified via config.covariates.field_hints.
  - inputs.imaging_confounds_file is consumed if given (joined as-is).

Outputs:
  results/intermediate/covariates.tsv   (eid + covariate columns)
  results/covariate_missingness.csv     (per-covariate missingness summary)
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from utils import (  # noqa: E402
    load_config,
    read_header_only,
    resolve_under_repo,
    safe_read_table,
    setup_logger,
    standardize_eid_column,
)


def _basket_col_pattern(field_hint: str) -> re.Pattern:
    """Convert '21003-2.0' or '22009-0.' style hints into a basket-col regex.

    Basket cols are f.<id>.<inst>.<arr> (dot-separated).
    """
    if field_hint.endswith("-"):
        # bare field id with trailing '-' → match any instance/array
        fid = field_hint.rstrip("-")
        return re.compile(rf"^f\.{re.escape(fid)}\.\d+\.\d+$")
    # split on '-' and '.'
    parts = re.split(r"[-.]", field_hint)
    if len(parts) == 1:
        return re.compile(rf"^f\.{re.escape(parts[0])}\.\d+\.\d+$")
    if len(parts) == 2:
        # field-instance, any array
        return re.compile(rf"^f\.{re.escape(parts[0])}\.{re.escape(parts[1])}\.\d+$")
    return re.compile(rf"^f\.{re.escape(parts[0])}\.{re.escape(parts[1])}\.{re.escape(parts[2])}$")


def _pc_pattern(pc_prefix: str, count: int) -> list[re.Pattern]:
    """pc_prefix like '22009-0.' → matches f.22009.0.<1..count>."""
    prefix = pc_prefix.rstrip(".").rstrip("-")
    parts = re.split(r"[-.]", prefix)
    if len(parts) >= 2:
        fid, inst = parts[0], parts[1]
    else:
        fid, inst = parts[0], "0"
    return [re.compile(rf"^f\.{re.escape(fid)}\.{re.escape(inst)}\.{i}$")
            for i in range(1, count + 1)]


def _find_basket_cols(header: list[str], hints: list[str]) -> list[str]:
    found = []
    for h in hints:
        pat = _basket_col_pattern(h)
        for col in header:
            if pat.match(col):
                found.append(col)
                break
    return found


def _field_id_from_hint(hint: str) -> str:
    return re.split(r"[-.]", str(hint))[0]


def _load_extra_covariates_from_basket(basket: str, cfg: dict, logger) -> pd.DataFrame:
    """Load focused-analysis imaging QC, WMH, and vascular covariates.

    Missing requested fields are logged but are not fatal because not every UKB
    basket has every processing/QC release field.
    """
    field_hints = cfg["covariates"].get("field_hints", {})
    header = read_header_only(basket)
    wanted: dict[str, str] = {}
    categorical: set[str] = set()

    def add_one(key: str, friendly: str, categorical_flag: bool = False) -> None:
        hints = field_hints.get(key, [])
        cols = _find_basket_cols(header, hints)
        if cols:
            wanted[cols[0]] = friendly
            if categorical_flag:
                categorical.add(friendly)
            logger.info(f"  extra {friendly:28s} ← {cols[0]}")
        else:
            logger.warning(f"  extra {friendly:28s} ← NOT FOUND (tried {hints})")

    def add_many(key: str, prefix: str, categorical_flag: bool = False) -> None:
        for hint in field_hints.get(key, []):
            cols = _find_basket_cols(header, [hint])
            fid = _field_id_from_hint(hint)
            name = f"{prefix}_{fid}"
            if cols:
                wanted[cols[0]] = name
                if categorical_flag:
                    categorical.add(name)
                logger.info(f"  extra {name:28s} ← {cols[0]}")
            else:
                logger.warning(f"  extra {name:28s} ← NOT FOUND (tried {hint})")

    add_one("t1_motion", "t1_motion")
    add_many("scanner_position", "scanner_pos")
    add_many("dmri_qc", "dmri_qc")
    add_many("t1_qc", "t1_qc")
    add_many("modality_discrepancy", "modality_discrepancy")
    add_many("protocol_flags", "protocol", categorical_flag=True)
    add_many("wmh_burden", "wmh")

    vascular_names = {
        "21001": "bmi",
        "20116": "smoking_status",
        "2443": "diabetes",
        "6150": "hypertension_6150",
        "6177": "hypertension_6177",
        "20002": "hypertension_self_report",
    }
    for hint in field_hints.get("vascular_risk", []):
        cols = _find_basket_cols(header, [hint])
        fid = _field_id_from_hint(hint)
        name = vascular_names.get(fid, f"vascular_{fid}")
        if cols:
            wanted[cols[0]] = name
            if fid in {"20116", "2443", "6150", "6177", "20002"}:
                categorical.add(name)
            logger.info(f"  extra {name:28s} ← {cols[0]}")
        else:
            logger.warning(f"  extra {name:28s} ← NOT FOUND (tried {hint})")

    if not wanted:
        return pd.DataFrame({"eid": pd.Series([], dtype="Int64")})
    cols_to_read = ["f.eid"] + list(wanted.keys())
    df = safe_read_table(basket, usecols=cols_to_read)
    df = standardize_eid_column(df, prefer="f.eid")
    df = df.rename(columns=wanted)
    for c in df.columns:
        if c == "eid":
            continue
        if c in categorical:
            df[c] = df[c].astype("category")
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _load_sqc_or_generic(path: str, cfg: dict, logger) -> pd.DataFrame:
    """Detect SQC (`#FID` header, UKB_PC1..PCk cols) vs generic covariate TSV.

    For SQC: rename UKB_PCi → PCi, IID → eid, genotyping_array → array,
    pick up age/sex/population/population_MM/used_in_pca_calculation/pass_QC_filter/sr_WB
    and pass them through so step 04 can filter on the flags.
    """
    p = Path(path)
    # gz / plain detection handled by pandas
    sep = "\t"
    # Peek header to detect format (handle gzip transparently)
    if str(p).endswith(".gz"):
        import gzip
        with gzip.open(p, "rt") as fh:
            head = fh.readline().rstrip("\n")
    else:
        with p.open("r") as fh:
            head = fh.readline().rstrip("\n")
    cols = head.lstrip("#").split(sep)
    is_sqc = ("UKB_PC1" in cols) or ("population_MM" in cols) or ("FID" in cols and "IID" in cols and "population" in cols)
    pc_count = int(cfg["covariates"].get("genetic_pc_count", 20))
    if is_sqc:
        logger.info("  detected UKB SQC format")
        # decide which columns to read
        wanted = ["FID", "IID"]
        for c in ("population", "population_MM", "sex", "age", "genotyping_array", "array",
                  "used_in_pca_calculation", "pass_QC_filter", "het_missing_outliers",
                  "sr_WB", "pop_initial", "genetic_kinship_to_others", "genotype_batch"):
            if c in cols:
                wanted.append(c)
        for i in range(1, pc_count + 1):
            c = f"UKB_PC{i}"
            if c in cols:
                wanted.append(c)
        # pandas read with comment='#' would treat the whole header line as a comment;
        # instead we read with header=None and set names manually.
        df = pd.read_csv(p, sep="\t", comment=None, header=0, low_memory=False,
                         usecols=lambda c: c.lstrip("#") in wanted,
                         na_values=["NA", "", "N/A"])
        df.columns = [c.lstrip("#") for c in df.columns]
        # IID → eid (drop sentinel rows with FID=-1)
        if "FID" in df.columns:
            df = df[pd.to_numeric(df["FID"], errors="coerce") >= 0]
        df = standardize_eid_column(df, prefer="IID")
        # rename PCs
        ren = {f"UKB_PC{i}": f"PC{i}" for i in range(1, pc_count + 1)}
        ren["genotyping_array"] = "array"
        df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
        # Note: SQC's `array` column is also present as 0/1; if both exist after rename
        # we prefer the more descriptive one. Drop the legacy column to avoid collision.
        if "array" in df.columns and df.columns.tolist().count("array") > 1:
            # pandas duplicates → keep first
            df = df.loc[:, ~df.columns.duplicated(keep="first")]
        # sex: SQC's `sex` column is 0/1 (0=female, 1=male per .psam)
        if "sex" in df.columns:
            df["sex"] = pd.to_numeric(df["sex"], errors="coerce").astype("Int64")
        if "age" in df.columns:
            df["age"] = pd.to_numeric(df["age"], errors="coerce")
        logger.info(f"  SQC loaded: {len(df):,} rows, "
                     f"{sum(c.startswith('PC') for c in df.columns)} PCs")
        return df
    # generic path
    logger.info("  treating as generic covariate TSV/CSV")
    df = safe_read_table(p)
    df = standardize_eid_column(df)
    return df


def _load_imaging_covariates_from_basket(basket: str, cfg: dict, logger) -> pd.DataFrame:
    """Pull just the imaging-specific covariates (site, head_size, motion, acq_date)
    from the basket — used to augment an SQC that lacks them."""
    field_hints = cfg["covariates"].get("field_hints", {})
    header = read_header_only(basket)
    wanted: dict[str, str] = {}
    for friendly_meta, friendly in (
        ("imaging_centre", "site"),
        ("head_size_scaling", "head_size"),
        ("head_motion_dmri", "motion_dmri"),
        ("acquisition_date", "acq_date"),
    ):
        hints = field_hints.get(friendly_meta, [])
        cols = _find_basket_cols(header, hints)
        if cols:
            wanted[cols[0]] = friendly
            logger.info(f"  imaging {friendly:14s} ← {cols[0]}")
        else:
            logger.warning(f"  imaging {friendly:14s} ← NOT FOUND (tried {hints})")
    if not wanted:
        return pd.DataFrame({"eid": pd.Series([], dtype="Int64")})
    cols_to_read = ["f.eid"] + list(wanted.keys())
    df = safe_read_table(basket, usecols=cols_to_read)
    df = standardize_eid_column(df, prefer="f.eid")
    df = df.rename(columns=wanted)
    for c in df.columns:
        if c == "eid":
            continue
        if c in ("site",):
            df[c] = df[c].astype("category")
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _load_basket_covariates(basket: str, cfg: dict, logger) -> pd.DataFrame:
    field_hints = cfg["covariates"].get("field_hints", {})
    pc_count = int(cfg["covariates"].get("genetic_pc_count", 20))
    header = read_header_only(basket)
    logger.info(f"basket has {len(header)} columns")

    wanted: dict[str, str] = {}  # basket_col -> friendly_name
    # singletons
    name_map = {
        "age_at_imaging":   ("age",            field_hints.get("age_at_imaging", [])),
        "sex":              ("sex",            field_hints.get("sex", [])),
        "imaging_centre":   ("site",           field_hints.get("imaging_centre", [])),
        "head_size_scaling":("head_size",      field_hints.get("head_size_scaling", [])),
        "head_motion_dmri": ("motion_dmri",    field_hints.get("head_motion_dmri", [])),
        "array":            ("array",          field_hints.get("array", [])),
        "acquisition_date": ("acq_date",       field_hints.get("acquisition_date", [])),
    }
    for friendly_meta, (friendly, hints) in name_map.items():
        cols = _find_basket_cols(header, hints)
        if cols:
            wanted[cols[0]] = friendly
            logger.info(f"  {friendly:14s} ← {cols[0]} (from hints {hints})")
        else:
            logger.warning(f"  {friendly:14s} ← NOT FOUND (tried {hints})")

    # PCs
    pc_prefix = field_hints.get("pc_prefix", "22009-0.")
    pc_patterns = _pc_pattern(pc_prefix, pc_count)
    pc_cols_found = 0
    for i, pat in enumerate(pc_patterns, start=1):
        match = next((c for c in header if pat.match(c)), None)
        if match:
            wanted[match] = f"PC{i}"
            pc_cols_found += 1
    logger.info(f"  PCs found: {pc_cols_found}/{pc_count} (prefix '{pc_prefix}')")

    cols_to_read = ["f.eid"] + list(wanted.keys())
    logger.info(f"reading {len(cols_to_read)} covariate cols from basket")
    df = safe_read_table(basket, usecols=cols_to_read)
    df = standardize_eid_column(df, prefer="f.eid")
    df = df.rename(columns=wanted)
    # Coerce types
    for c in df.columns:
        if c == "eid":
            continue
        if c in ("sex", "site", "array"):
            df[c] = df[c].astype("category")
        else:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def _derive_extra_covars(df: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    if "age" in df.columns and cfg["covariates"].get("include_age2", True):
        df["age2"] = df["age"] ** 2
    if "age" in df.columns and "sex" in df.columns and cfg["covariates"].get("include_age_by_sex", True):
        sex_num = pd.to_numeric(df["sex"], errors="coerce")
        df["age_sex"] = df["age"] * sex_num
    return df


def _missingness_report(df: pd.DataFrame) -> pd.DataFrame:
    n = len(df)
    rows = []
    for c in df.columns:
        if c == "eid":
            continue
        miss = int(df[c].isna().sum())
        rows.append({"covariate": c, "missing_n": miss,
                     "missing_pct": (miss / n * 100.0) if n else float("nan"),
                     "n_total": n})
    return pd.DataFrame(rows).sort_values("missing_pct", ascending=False)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    cfg = load_config(args.config)
    dry = args.dry_run or bool(cfg.get("project", {}).get("dry_run", True))
    inputs = cfg["inputs"]

    logs_dir = resolve_under_repo(cfg, cfg["project"]["logs_dir"])
    inter_dir = resolve_under_repo(cfg, cfg["project"]["intermediate_dir"])
    results_dir = resolve_under_repo(cfg, cfg["project"]["output_dir"])
    inter_dir.mkdir(parents=True, exist_ok=True)
    log = setup_logger("build_covariates", logs_dir / "run_log.txt")

    out_tsv = inter_dir / "covariates.tsv"
    out_csv = results_dir / "covariate_missingness.csv"

    pre = inputs.get("genetic_covariates_file") or ""
    basket = inputs.get("basket_tab") or ""

    if pre and Path(pre).is_file():
        if dry:
            log.info(f"DRY RUN — would load SQC/covariate file {pre}")
            return 0
        log.info(f"loading pre-built covariates: {pre}")
        df = _load_sqc_or_generic(pre, cfg, log)
        # SQC has only genetic covariates; pull imaging-specific covariates
        # (site, head_size, motion_dmri, acq_date) from basket if available.
        if basket and Path(basket).is_file():
            log.info("augmenting SQC with imaging covariates from basket")
            img = _load_imaging_covariates_from_basket(basket, cfg, log)
            df = df.merge(img, on="eid", how="left")
            extra = _load_extra_covariates_from_basket(basket, cfg, log)
            df = df.merge(extra, on="eid", how="left")
    elif basket and Path(basket).is_file():
        if dry:
            log.info(f"DRY RUN — would scan basket header and load chosen covariate cols")
            # still scan header to validate availability
            header = read_header_only(basket)
            log.info(f"  basket cols = {len(header)}")
            return 0
        df = _load_basket_covariates(basket, cfg, log)
        extra = _load_extra_covariates_from_basket(basket, cfg, log)
        df = df.merge(extra, on="eid", how="left")
    else:
        log.error("no covariate source: provide inputs.genetic_covariates_file or inputs.basket_tab.")
        return 1

    df = _derive_extra_covars(df, cfg)

    # optional join: imaging confounds file
    conf = inputs.get("imaging_confounds_file") or ""
    if conf and Path(conf).is_file():
        log.info(f"merging imaging confounds: {conf}")
        cdf = safe_read_table(conf)
        cdf = standardize_eid_column(cdf)
        df = df.merge(cdf, on="eid", how="left", suffixes=("", "_conf"))

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv, sep="\t", index=False)
    log.info(f"wrote {out_tsv} ({len(df):,} rows, {len(df.columns)} cols)")

    miss = _missingness_report(df)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    miss.to_csv(out_csv, index=False)
    log.info(f"wrote {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
