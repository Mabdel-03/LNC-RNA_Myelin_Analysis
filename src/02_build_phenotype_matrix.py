#!/usr/bin/env python
"""Step 02 — build phenotype_manifest.csv and a wide per-eid phenotype TSV.

Source preference:
  1) inputs.phenotype_file (user-supplied wide CSV with eid + IDP columns).
  2) inputs.basket_tab (UKB basket, columns f.<id>.<inst>.<arr>) — primary path.
  3) inputs.mri_pheno_tsv_dirs (per-field TSVs already extracted) — used to
     ENRICH the wide matrix when a primary IDP isn't in the basket subset
     selected, OR as the sole source if the basket is unavailable.

The manifest classifies each phenotype into primary / secondary / exploratory
using a built-in UKB IDP field-id map (FA, MD, L1-L3, NODDI ICVF/ISOVF/ODI,
WMH) and the term lists in config. If inputs.data_dictionary_file is given,
descriptions are pulled from there for authoritative labels.

Outputs:
  results/phenotype_manifest.csv     (one row per IDP, columns per spec)
  results/intermediate/phenotypes_wide.tsv  (eid + selected IDPs)
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
    parse_ukb_field_id,
    read_header_only,
    resolve_under_repo,
    safe_read_table,
    setup_logger,
    standardize_eid_column,
)


# Authoritative UKB IDP field-id ranges (verified against the UKB Showcase
# Data Dictionary). TBSS = 48 fields per metric; weighted-mean = 27. Order
# within each block: FA, MD, MO, L1, L2, L3, ICVF, OD, ISOVF (note OD
# precedes ISOVF — easy to swap by accident).
# Each entry: (range_start, range_end_inclusive, metric, modality, panel_default)
KNOWN_IDP_RANGES = [
    # ---- Diffusion TBSS (skeleton) ----
    (25056, 25103, "FA",     "dMRI_TBSS",   "primary"),
    (25104, 25151, "MD",     "dMRI_TBSS",   "primary"),
    (25152, 25199, "MO",     "dMRI_TBSS",   "exploratory"),
    (25200, 25247, "L1",     "dMRI_TBSS",   "primary"),
    (25248, 25295, "L2",     "dMRI_TBSS",   "primary"),
    (25296, 25343, "L3",     "dMRI_TBSS",   "primary"),
    (25344, 25391, "ICVF",   "dMRI_TBSS",   "primary"),
    (25392, 25439, "OD",     "dMRI_TBSS",   "primary"),
    (25440, 25487, "ISOVF",  "dMRI_TBSS",   "primary"),
    # ---- Diffusion weighted-mean per tract (probtrackx) ----
    (25488, 25514, "FA",     "dMRI_weighted_mean", "primary"),
    (25515, 25541, "MD",     "dMRI_weighted_mean", "primary"),
    (25542, 25568, "MO",     "dMRI_weighted_mean", "exploratory"),
    (25569, 25595, "L1",     "dMRI_weighted_mean", "primary"),
    (25596, 25622, "L2",     "dMRI_weighted_mean", "primary"),
    (25623, 25649, "L3",     "dMRI_weighted_mean", "primary"),
    (25650, 25676, "ICVF",   "dMRI_weighted_mean", "primary"),
    (25677, 25703, "OD",     "dMRI_weighted_mean", "primary"),
    (25704, 25730, "ISOVF",  "dMRI_weighted_mean", "primary"),
]

# Single-field IDPs (id → (metric, modality, panel_default, label_hint))
KNOWN_SINGLE_IDPS = {
    25781: ("WMH_volume",       "T2_FLAIR", "secondary",  "Total volume of white matter hyperintensities"),
    25000: ("head_size_scaling","T1",       "covariate",  "Head size scaling factor"),
    25741: ("rfMRI_motion",     "rfMRI",    "exploratory","Mean rfMRI head motion"),
    25742: ("tfMRI_motion",     "tfMRI",    "exploratory","Mean tfMRI head motion"),
    25746: ("dMRI_outlier_slices","dMRI",   "exploratory","Number of dMRI outlier slices detected and corrected"),
}


# Description-based classifier — preferred when a dictionary is loaded.
# Parses canonical UKB IDP titles:
#   TBSS:  "Mean <METRIC> in <tract> on <METRIC> skeleton"
#   WM:    "Weighted-mean <METRIC> in tract <tract>"
_TBSS_RE = re.compile(
    r"^Mean\s+(\w+)\s+in\s+(.+?)\s+on\s+\w+\s+skeleton\s*(\([^)]+\))?\s*$",
    re.IGNORECASE)
_WM_RE = re.compile(r"^Weighted-mean\s+(\w+)\s+in\s+tract\s+(.+)$", re.IGNORECASE)
_METRIC_PANEL_DEFAULT = {
    "FA": "primary", "MD": "primary", "L1": "primary", "L2": "primary",
    "L3": "primary", "ICVF": "primary", "ISOVF": "primary", "OD": "primary",
    "MO": "exploratory",
}


def _classify_from_description(desc: str) -> tuple[str | None, str | None, str | None, str | None]:
    """Returns (metric, modality, panel_default, tract_name) by parsing the showcase Field text.
    All None if the description doesn't match any IDP pattern."""
    if not desc:
        return None, None, None, None
    text = str(desc).strip()
    m = _TBSS_RE.match(text)
    if m:
        metric = m.group(1).upper()
        tract = m.group(2).strip()
        side = (m.group(3) or "").strip("() ").lower()
        if side:
            tract = f"{tract} ({side})"
        panel = _METRIC_PANEL_DEFAULT.get(metric, "exploratory")
        return metric, "dMRI_TBSS", panel, tract
    m = _WM_RE.match(text)
    if m:
        metric = m.group(1).upper()
        tract = m.group(2).strip()
        panel = _METRIC_PANEL_DEFAULT.get(metric, "exploratory")
        return metric, "dMRI_weighted_mean", panel, tract
    return None, None, None, None


def _slugify(text: str, max_len: int = 40) -> str:
    """Lowercase, replace non-alphanumeric with '_', collapse repeats, trim."""
    s = re.sub(r"[^A-Za-z0-9]+", "_", str(text)).strip("_").lower()
    return s[:max_len] if max_len else s


def _basket_col_to_field(col: str) -> tuple[int | None, int | None, int | None]:
    """Parse 'f.25526.2.0' → (25526, 2, 0). Returns (None, None, None) if not matchable."""
    m = re.match(r"^f\.(\d+)\.(\d+)\.(\d+)$", col)
    if not m:
        return None, None, None
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def _classify_field(field_id: int) -> tuple[str | None, str | None, str | None]:
    """Returns (metric, modality, panel_default) from the known IDP map, or
    (None, None, None) if field_id is not in any known IDP block."""
    for lo, hi, metric, modality, panel in KNOWN_IDP_RANGES:
        if lo <= field_id <= hi:
            return metric, modality, panel
    if field_id in KNOWN_SINGLE_IDPS:
        m, mod, panel, _ = KNOWN_SINGLE_IDPS[field_id]
        return m, mod, panel
    return None, None, None


def _load_data_dictionary(path: str | None) -> dict[int, str]:
    """Optional. Expects a TSV/CSV with at least field_id + description columns.

    Handles the UKB Data_Dictionary_Showcase.csv quirk of using backslash-escaped
    inner quotes (\"...\") rather than RFC-standard doubled quotes.
    """
    if not path:
        return {}
    p = Path(path)
    if not p.is_file():
        return {}
    suf = p.suffix.lower()
    if suf == ".csv":
        df = pd.read_csv(p, engine="python", escapechar="\\", quotechar='"',
                         na_values=["", "NA", "N/A"])
    else:
        df = safe_read_table(p)
    # try common column names — order matters: prefer exact "field_id"/"fieldid" for id;
    # 'Field' in Data_Dictionary_Showcase.csv is the description, not the id.
    fid_col = None
    desc_col = None
    for c in df.columns:
        cl = c.lower()
        if fid_col is None and cl in ("field_id", "fieldid", "id"):
            fid_col = c
    for c in df.columns:
        cl = c.lower()
        if desc_col is None and cl in ("field", "description", "title", "name", "field_name"):
            desc_col = c
    if fid_col is None or desc_col is None:
        return {}
    df[fid_col] = pd.to_numeric(df[fid_col], errors="coerce")
    df = df.dropna(subset=[fid_col])
    return dict(zip(df[fid_col].astype(int), df[desc_col].astype(str)))


def _build_manifest_from_basket(basket_path: str,
                                 dict_map: dict[int, str],
                                 cfg: dict,
                                 logger) -> tuple[pd.DataFrame, list[str]]:
    """Scan basket header → manifest rows + list of columns to read."""
    logger.info(f"scanning basket header: {basket_path}")
    header = read_header_only(basket_path)
    rows = []
    cols_to_read = ["f.eid"]
    for col in header:
        fid, inst, arr = _basket_col_to_field(col)
        if fid is None:
            continue
        # imaging visit = instance 2 (and 3 for return scan); we take 2.
        if inst != 2:
            continue
        # only first array element if multiple
        if arr != 0:
            continue
        desc = dict_map.get(fid, "")
        # Prefer description-based classification (more robust than field-id ranges).
        metric, modality, panel, tract_name = _classify_from_description(desc)
        classifier = "dict"
        if metric is None:
            # fall back to the field-id range map
            metric, modality, panel = _classify_field(fid)
            tract_name = None
            classifier = "range"
        if metric is None:
            continue
        # Build a column-friendly name. Prefer tract slug from description; else field-id suffix.
        if metric in ("FA", "MD", "MO", "L1", "L2", "L3", "ICVF", "ISOVF", "OD"):
            if tract_name:
                wide_name = f"{modality}_{metric}_{_slugify(tract_name)}"
            else:
                wide_name = f"{modality}_{metric}_f{fid}"
            tract_disp = tract_name if tract_name else _slugify(f"f{fid}")
        else:
            wide_name = metric  # e.g. WMH_volume
            tract_disp = tract_name or ""
        # Pull a fallback description from the single-IDP map if dictionary missing
        if not desc and fid in KNOWN_SINGLE_IDPS:
            desc = KNOWN_SINGLE_IDPS[fid][3]
        rows.append({
            "field_id": fid,
            "instance": inst,
            "array": arr,
            "basket_col": col,
            "column_name": wide_name,
            "description": desc,
            "modality_guess": modality,
            "metric_guess": metric,
            "tract_or_region_guess": tract_disp,
            "panel": "secondary" if panel == "secondary" else
                     ("primary" if panel == "primary" else "exploratory"),
            "transform": "log1p_then_rank_inverse" if (panel == "secondary"
                          and (metric == "WMH_volume" or "WMH" in metric)) else "rank_inverse_normal",
            "include": True if panel in ("primary", "secondary") else False,
            "notes": f"classified_by={classifier}",
            "source": "basket",
        })
        cols_to_read.append(col)
    df = pd.DataFrame(rows)
    # enforce uniqueness of column_name — if two field_ids collide on the slug,
    # disambiguate by appending the field_id
    if not df.empty:
        dup_mask = df["column_name"].duplicated(keep=False)
        if dup_mask.any():
            logger.warning(f"  {dup_mask.sum()} column-name collisions; appending field_id to disambiguate")
            df.loc[dup_mask, "column_name"] = df.loc[dup_mask].apply(
                lambda r: f"{r['column_name']}_f{int(r['field_id'])}", axis=1)
    if df.empty:
        logger.warning("no IDP columns matched in basket — check basket file / field map")
    else:
        # apply config term filters (primary/secondary/exploratory) for description-based promotion
        prim_terms = [t.lower() for t in cfg["phenotypes"].get("primary_terms", [])]
        sec_terms = [t.lower() for t in cfg["phenotypes"].get("secondary_terms", [])]
        regions = [r.lower() for r in cfg["phenotypes"].get("prioritize_regions", [])]
        if not df["description"].fillna("").str.strip().eq("").all():
            d = df["description"].fillna("").str.lower()
            pri_hits = d.str.contains("|".join(re.escape(t) for t in prim_terms)) if prim_terms else False
            sec_hits = d.str.contains("|".join(re.escape(t) for t in sec_terms)) if sec_terms else False
            df.loc[pri_hits, "panel"] = "primary"
            df.loc[sec_hits & (df["panel"] != "primary"), "panel"] = "secondary"
            if regions:
                df["region_priority"] = d.str.contains("|".join(re.escape(r) for r in regions)).astype(int)
            else:
                df["region_priority"] = 0
        else:
            df["region_priority"] = 0
    logger.info(f"manifest rows: {len(df)} (primary={sum(df['panel']=='primary')}, "
                f"secondary={sum(df['panel']=='secondary')}, exploratory={sum(df['panel']=='exploratory')})")
    return df, cols_to_read


def _read_basket_columns(basket_path: str, cols_to_read: list[str],
                          dry_run: bool, logger) -> pd.DataFrame:
    if dry_run:
        logger.info(f"DRY RUN — would read {len(cols_to_read)-1} IDP cols + eid from basket")
        # return a tiny header-only frame for downstream merge sanity
        return pd.DataFrame({"eid": pd.Series([], dtype="Int64")})
    logger.info(f"reading {len(cols_to_read)} basket columns from {basket_path}")
    df = safe_read_table(basket_path, usecols=cols_to_read)
    df = standardize_eid_column(df, prefer="f.eid")
    return df


def _enrich_from_local_tsvs(wide: pd.DataFrame, manifest: pd.DataFrame,
                             tsv_dirs: list[str], logger) -> pd.DataFrame:
    """Optional enrichment: pull additional per-field TSVs not already in wide."""
    if wide.empty:
        return wide
    existing_fields = set(int(x) for x in manifest["field_id"].tolist() if pd.notna(x))
    for d in tsv_dirs or []:
        dp = Path(d)
        if not dp.is_dir():
            continue
        for tsv in dp.rglob("phenotype_*.tsv"):
            fid = parse_ukb_field_id(tsv.name)
            if not fid or int(fid) in existing_fields:
                continue
            try:
                df = pd.read_csv(tsv, sep="\t")
            except Exception as exc:
                logger.warning(f"could not parse {tsv}: {exc}")
                continue
            df = standardize_eid_column(df, prefer="IID")
            if "Phenotype" not in df.columns:
                continue
            metric, modality, panel = _classify_field(int(fid))
            if metric is None:
                continue
            wide_name = f"{modality}_{metric}_f{fid}"
            logger.info(f"  + enriching with {wide_name} from {tsv}")
            wide = wide.merge(df[["eid", "Phenotype"]].rename(
                columns={"Phenotype": wide_name}), on="eid", how="left")
    return wide


def _derive_rd(wide: pd.DataFrame, manifest: pd.DataFrame,
                logger) -> tuple[pd.DataFrame, pd.DataFrame]:
    """For each modality where L2 and L3 columns are present and no RD column,
    create RD = (L2+L3)/2 as a new column + manifest row."""
    new_rows = []
    for modality in ("dMRI_TBSS", "dMRI_weighted_mean"):
        l2_rows = manifest[(manifest["modality_guess"] == modality)
                           & (manifest["metric_guess"] == "L2")]
        l3_rows = manifest[(manifest["modality_guess"] == modality)
                           & (manifest["metric_guess"] == "L3")]
        rd_rows = manifest[(manifest["modality_guess"] == modality)
                           & (manifest["metric_guess"] == "RD")]
        if l2_rows.empty or l3_rows.empty or not rd_rows.empty:
            continue
        # Match by tract index suffix in column_name
        def _tract_key(name: str) -> str:
            return name.rsplit("_", 1)[-1]
        l2_map = dict(zip(l2_rows["column_name"].map(_tract_key), l2_rows["column_name"]))
        l3_map = dict(zip(l3_rows["column_name"].map(_tract_key), l3_rows["column_name"]))
        for k in sorted(set(l2_map) & set(l3_map)):
            rd_col = f"{modality}_RD_{k}"
            if rd_col in wide.columns:
                continue
            wide[rd_col] = (pd.to_numeric(wide[l2_map[k]], errors="coerce") +
                            pd.to_numeric(wide[l3_map[k]], errors="coerce")) / 2.0
            new_rows.append({
                "field_id": np.nan,
                "instance": 2,
                "array": 0,
                "basket_col": "",
                "column_name": rd_col,
                "description": f"Derived RD=(L2+L3)/2 for {modality} tract {k}",
                "modality_guess": modality,
                "metric_guess": "RD",
                "tract_or_region_guess": k,
                "panel": "primary",
                "transform": "rank_inverse_normal",
                "include": True,
                "notes": "derived",
                "source": "derived",
                "region_priority": 0,
            })
    if new_rows:
        logger.info(f"derived {len(new_rows)} RD phenotypes")
        manifest = pd.concat([manifest, pd.DataFrame(new_rows)], ignore_index=True)
    return wide, manifest


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
    log = setup_logger("build_phenotypes", logs_dir / "run_log.txt")

    manifest_path = results_dir / "phenotype_manifest.csv"
    wide_path = inter_dir / "phenotypes_wide.tsv"

    pheno_file = inputs.get("phenotype_file") or ""
    basket = inputs.get("basket_tab") or ""
    dict_map = _load_data_dictionary(inputs.get("data_dictionary_file"))

    # ---- branch A: user-supplied wide phenotype CSV (we don't classify here)
    if pheno_file and Path(pheno_file).is_file():
        log.info(f"using user-supplied phenotype file: {pheno_file}")
        if dry:
            cols = read_header_only(pheno_file)
            log.info(f"  → {len(cols)} columns")
            manifest = pd.DataFrame([{
                "field_id": np.nan, "column_name": c, "description": "",
                "modality_guess": "", "metric_guess": "", "tract_or_region_guess": "",
                "panel": "primary", "transform": "rank_inverse_normal",
                "include": True, "notes": "user-supplied wide csv",
                "source": "user", "region_priority": 0,
            } for c in cols if c.lower() not in ("eid", "iid", "fid")])
            manifest.to_csv(manifest_path, index=False)
            log.info(f"wrote {manifest_path} (dry-run, no wide tsv)")
            return 0
        wide = safe_read_table(pheno_file)
        wide = standardize_eid_column(wide)
    elif basket and Path(basket).is_file():
        manifest, cols_to_read = _build_manifest_from_basket(basket, dict_map, cfg, log)
        wide = _read_basket_columns(basket, cols_to_read, dry, log)
        # rename basket cols → friendly wide names
        if not dry and not wide.empty and not manifest.empty:
            ren = dict(zip(manifest["basket_col"], manifest["column_name"]))
            ren = {k: v for k, v in ren.items() if k and k in wide.columns}
            wide = wide.rename(columns=ren)
            # keep only id + classified cols
            keep = ["eid"] + [c for c in manifest["column_name"] if c in wide.columns]
            wide = wide[keep]
    else:
        log.error("no phenotype source: provide inputs.phenotype_file OR inputs.basket_tab.")
        return 1

    if not dry:
        # derive RD where possible
        wide, manifest = _derive_rd(wide, manifest, log)
        # optional enrichment from local per-field TSVs
        wide = _enrich_from_local_tsvs(wide, manifest, inputs.get("mri_pheno_tsv_dirs", []), log)
        # write wide
        wide_path.parent.mkdir(parents=True, exist_ok=True)
        wide.to_csv(wide_path, sep="\t", index=False)
        log.info(f"wrote {wide_path} ({len(wide):,} rows, {len(wide.columns)} cols)")

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_path, index=False)
    log.info(f"wrote {manifest_path} ({len(manifest)} rows)")

    if cfg["phenotypes"].get("require_manual_manifest_approval", False):
        # require user to set 'approved' col True
        if "approved" not in manifest.columns or not manifest["approved"].fillna(False).any():
            log.warning("phenotype manifest requires manual approval; add column "
                        "'approved=True' to the rows you want to include, then rerun "
                        "step 02 (or set require_manual_manifest_approval=false).")
            return 0 if dry else 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
