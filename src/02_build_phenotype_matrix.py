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
    # ---- Exploratory SWI / T2* / grey-white contrast contexts ----
    (24467, 24484, "QSM",    "SWI_deep_gm", "exploratory"),
    (25026, 25039, "T2star", "SWI_T2star",  "exploratory"),
    (26989, 27058, "GWC",    "T1_GWC",      "exploratory"),
]

# Single-field IDPs (id → (metric, modality, panel_default, label_hint))
KNOWN_SINGLE_IDPS = {
    25781: ("WMH_volume",       "T2_FLAIR", "secondary",  "Total volume of white matter hyperintensities"),
    24485: ("WMH_count",        "T2_FLAIR", "secondary",  "White matter hyperintensity count"),
    24486: ("WMH_mean_volume",  "T2_FLAIR", "secondary",  "Mean white matter hyperintensity volume"),
    25007: ("WM_volume_normalized","T1",    "exploratory","Volume of white matter normalized for head size"),
    25008: ("WM_volume_raw",    "T1",       "exploratory","Volume of white matter"),
    25000: ("head_size_scaling","T1",       "covariate",  "Head size scaling factor"),
    25741: ("rfMRI_motion",     "rfMRI",    "exploratory","Mean rfMRI head motion"),
    25742: ("tfMRI_motion",     "tfMRI",    "exploratory","Mean tfMRI head motion"),
    25746: ("dMRI_outlier_slices","dMRI",   "exploratory","Number of dMRI outlier slices detected and corrected"),
}


TBSS_OFFSETS = {
    "FA": 0, "MD": 48, "MO": 96, "L1": 144, "L2": 192, "L3": 240,
    "ICVF": 288, "OD": 336, "ISOVF": 384,
}
WEIGHTED_OFFSETS = {
    "FA": 0, "MD": 27, "MO": 54, "L1": 81, "L2": 108, "L3": 135,
    "ICVF": 162, "OD": 189, "ISOVF": 216,
}


def _as_int_set(values) -> set[int]:
    return {int(v) for v in (values or []) if pd.notna(v)}


def _field_ranges_to_set(ranges) -> set[int]:
    out: set[int] = set()
    for item in ranges or []:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            continue
        lo, hi = int(item[0]), int(item[1])
        out.update(range(lo, hi + 1))
    return out


def _field_id_base(field_id: int | float | None,
                   modality: str | None,
                   metric: str | None) -> int | float:
    if field_id is None or pd.isna(field_id):
        return np.nan
    metric_u = str(metric or "").upper()
    if modality == "dMRI_TBSS" and metric_u in TBSS_OFFSETS:
        return int(field_id) - TBSS_OFFSETS[metric_u]
    if modality == "dMRI_weighted_mean" and metric_u in WEIGHTED_OFFSETS:
        return int(field_id) - WEIGHTED_OFFSETS[metric_u]
    return int(field_id)


def _focused_metadata(field_id: int | float | None,
                      modality: str | None,
                      metric: str | None,
                      cfg: dict) -> dict:
    """Return raw-p analysis tier metadata for one manifest row."""
    fa = cfg.get("focused_analysis", {})
    metric_u = str(metric or "").upper()
    base = _field_id_base(field_id, modality, metric_u)
    primary_bases = _as_int_set(fa.get("tbss_primary_fa_base_ids"))
    replication_bases = _as_int_set(fa.get("weighted_replication_fa_base_ids"))
    primary_metrics = set(str(x).upper() for x in fa.get("primary_metrics", ["FA", "MD", "RD"]))
    directional_metrics = set(str(x).upper() for x in fa.get("directional_metrics", ["L1"]))
    mech_metrics = set(str(x).upper() for x in fa.get("mechanistic_metrics", ["ICVF", "ISOVF"]))
    control_ids = _as_int_set(fa.get("controls_field_ids", [25781]))
    exploratory_ids = _field_ranges_to_set(fa.get("exploratory_field_ranges", []))
    structural_ids = _as_int_set(fa.get("structural_field_ids", []))
    include_unfocused = bool(fa.get("include_unfocused_dmri", False))

    out = {
        "field_id_base": base,
        "analysis_tier": "unfocused",
        "confirmatory_role": "excluded",
        "include_in_lmm": False,
        "include": False,
        "family": "unfocused",
    }
    fid_int = int(field_id) if field_id is not None and pd.notna(field_id) else None
    if fid_int in control_ids or metric_u.startswith("WMH"):
        out.update({
            "analysis_tier": "controls_pathology",
            "confirmatory_role": "control",
            "include_in_lmm": True,
            "include": True,
            "family": "controls",
        })
        return out
    if fid_int in structural_ids:
        out.update({
            "analysis_tier": "exploratory_structural",
            "confirmatory_role": "exploratory",
            "include": True,
            "family": "exploratory_structural",
        })
        return out
    if fid_int in exploratory_ids:
        out.update({
            "analysis_tier": "exploratory_context",
            "confirmatory_role": "exploratory",
            "include": True,
            "family": "exploratory_context",
        })
        return out
    if modality == "dMRI_TBSS" and pd.notna(base) and int(base) in primary_bases:
        if metric_u in primary_metrics:
            out.update({
                "analysis_tier": "primary_skeleton_tensor",
                "confirmatory_role": "primary",
                "include_in_lmm": True,
                "include": True,
                "family": "primary_skeleton_tensor",
            })
        elif metric_u in directional_metrics:
            out.update({
                "analysis_tier": "directional_tensor_context",
                "confirmatory_role": "directional_context",
                "include_in_lmm": True,
                "include": True,
                "family": "directional_tensor_context",
            })
        elif metric_u in mech_metrics:
            out.update({
                "analysis_tier": "secondary_mechanistic",
                "confirmatory_role": "mechanistic_support",
                "include_in_lmm": True,
                "include": True,
                "family": "secondary_mechanistic",
            })
        else:
            out.update({
                "analysis_tier": "exploratory_focused_dmri",
                "confirmatory_role": "exploratory",
                "include": include_unfocused,
                "family": "exploratory_focused_dmri",
            })
        return out
    if modality == "dMRI_weighted_mean" and pd.notna(base) and int(base) in replication_bases:
        if metric_u in primary_metrics:
            out.update({
                "analysis_tier": "replication_weighted_tensor",
                "confirmatory_role": "replication",
                "include_in_lmm": True,
                "include": True,
                "family": "replication_weighted_tensor",
            })
        elif metric_u in mech_metrics:
            out.update({
                "analysis_tier": "replication_mechanistic",
                "confirmatory_role": "replication_support",
                "include_in_lmm": True,
                "include": True,
                "family": "replication_mechanistic",
            })
        elif metric_u in directional_metrics:
            out.update({
                "analysis_tier": "replication_directional_context",
                "confirmatory_role": "replication_context",
                "include_in_lmm": True,
                "include": True,
                "family": "replication_directional_context",
            })
        else:
            out.update({
                "analysis_tier": "exploratory_focused_dmri",
                "confirmatory_role": "exploratory",
                "include": include_unfocused,
                "family": "exploratory_focused_dmri",
            })
        return out
    if str(modality or "").startswith("dMRI_"):
        out.update({
            "analysis_tier": "unfocused_dmri",
            "confirmatory_role": "background",
            "include": include_unfocused,
            "family": "unfocused_dmri",
        })
    return out


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


def _classify_family(description: str, metric: str, panel: str,
                     primary_terms: list[str], controls_terms: list[str]) -> str:
    """Map a phenotype row to its focused analysis family.

    Single IDPs default to `exploratory` (the new PheWAS family). Myelin-sensitive
    MRI (MTR/MTsat/MWF/qT1/T2*/R2*) → `primary`. WMH and any matched control
    term → `controls`. Composites + ROI PCs are added later (step 04b) as
    `secondary`.
    """
    desc_l = (description or "").lower()
    metric_l = (metric or "").upper()
    if primary_terms:
        if any(t.lower() in desc_l for t in primary_terms):
            return "primary"
    if controls_terms:
        if any(t.lower() in desc_l for t in controls_terms):
            return "controls"
    if "WMH" in metric_l or "WHITE_MATTER_HYPERINTENS" in metric_l.upper():
        return "controls"
    return "exploratory"


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
    fam_primary_terms = list(cfg.get("families", {}).get("primary_terms", []))
    fam_controls_terms = list(cfg.get("families", {}).get("controls_terms", []))
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
            if fid in KNOWN_SINGLE_IDPS:
                wide_name = metric  # e.g. WMH_volume
            else:
                suffix = _slugify(desc, 44) if desc else f"f{fid}"
                wide_name = f"{modality}_{metric}_{suffix}"
            tract_disp = tract_name or (desc if desc else "")
        # Pull a fallback description from the single-IDP map if dictionary missing
        if not desc and fid in KNOWN_SINGLE_IDPS:
            desc = KNOWN_SINGLE_IDPS[fid][3]
        panel_norm = ("secondary" if panel == "secondary" else
                      ("primary" if panel == "primary" else "exploratory"))
        legacy_family = _classify_family(desc, metric, panel_norm,
                                         fam_primary_terms, fam_controls_terms)
        focus = _focused_metadata(fid, modality, metric, cfg)
        if not cfg.get("focused_analysis", {}).get("enabled", True):
            focus["family"] = legacy_family
            focus["include"] = True if panel in ("primary", "secondary") else False
            focus["include_in_lmm"] = focus["family"] in {"primary", "secondary", "controls"}
        row = {
            "field_id": fid,
            "field_id_base": focus["field_id_base"],
            "instance": inst,
            "array": arr,
            "basket_col": col,
            "column_name": wide_name,
            "description": desc,
            "modality_guess": modality,
            "metric_guess": metric,
            "tract_or_region_guess": tract_disp,
            "panel": panel_norm,
            "family": focus["family"],
            "analysis_tier": focus["analysis_tier"],
            "confirmatory_role": focus["confirmatory_role"],
            "include_in_lmm": focus["include_in_lmm"],
            "transform": "log1p_then_rank_inverse" if (metric.startswith("WMH") or "WMH" in metric) else "rank_inverse_normal",
            "include": focus["include"],
            "notes": f"classified_by={classifier}",
            "source": "basket",
        }
        rows.append(row)
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
        # Match by FULL post-metric tract suffix. Column names look like
        # "dMRI_TBSS_L2_acoustic_radiation_left" — split off the
        # "{modality}_{metric}_" prefix and use everything after as the key.
        # (The earlier .rsplit("_",1)[-1] collapsed every "*_left" tract
        # into one row, dropping us from ~75 RD to ~16.)
        l2_prefix = f"{modality}_L2_"
        l3_prefix = f"{modality}_L3_"
        def _strip(name: str, prefix: str) -> str:
            return name[len(prefix):] if name.startswith(prefix) else name
        l2_map = {_strip(c, l2_prefix): c for c in l2_rows["column_name"]}
        l3_map = {_strip(c, l3_prefix): c for c in l3_rows["column_name"]}
        for k in sorted(set(l2_map) & set(l3_map)):
            rd_col = f"{modality}_RD_{k}"
            if rd_col in wide.columns:
                continue
            wide[rd_col] = (pd.to_numeric(wide[l2_map[k]], errors="coerce") +
                            pd.to_numeric(wide[l3_map[k]], errors="coerce")) / 2.0
            # Reuse the L2 row's tract_or_region_guess so the derived RD row
            # keeps the human-readable tract name (not the slugified suffix).
            tract_disp = l2_rows.loc[l2_rows["column_name"] == l2_map[k],
                                       "tract_or_region_guess"].iloc[0]
            l2_row = l2_rows.loc[l2_rows["column_name"] == l2_map[k]].iloc[0]
            base = l2_row.get("field_id_base", np.nan)
            if pd.isna(base):
                base = _field_id_base(l2_row.get("field_id", np.nan), modality, "L2")
            focus = {
                "field_id_base": base,
                "analysis_tier": l2_row.get("analysis_tier", "unfocused"),
                "confirmatory_role": l2_row.get("confirmatory_role", "excluded"),
                "include_in_lmm": bool(l2_row.get("include_in_lmm", False)),
                "include": bool(l2_row.get("include", False)),
                "family": l2_row.get("family", "unfocused"),
            }
            if modality == "dMRI_TBSS":
                if str(l2_row.get("analysis_tier")) in {"primary_skeleton_tensor", "directional_tensor_context", "secondary_mechanistic", "exploratory_focused_dmri"}:
                    focus.update({
                        "analysis_tier": "primary_skeleton_tensor",
                        "confirmatory_role": "primary",
                        "include_in_lmm": True,
                        "include": True,
                        "family": "primary_skeleton_tensor",
                    })
            elif modality == "dMRI_weighted_mean":
                if str(l2_row.get("analysis_tier")) in {"replication_weighted_tensor", "replication_directional_context", "replication_mechanistic", "exploratory_focused_dmri"}:
                    focus.update({
                        "analysis_tier": "replication_weighted_tensor",
                        "confirmatory_role": "replication",
                        "include_in_lmm": True,
                        "include": True,
                        "family": "replication_weighted_tensor",
                    })
            new_rows.append({
                "field_id": np.nan,
                "field_id_base": focus["field_id_base"],
                "instance": 2,
                "array": 0,
                "basket_col": "",
                "column_name": rd_col,
                "description": f"Derived RD=(L2+L3)/2 for {modality} tract {tract_disp}",
                "modality_guess": modality,
                "metric_guess": "RD",
                "tract_or_region_guess": tract_disp,
                "panel": "primary",
                "family": focus["family"],
                "analysis_tier": focus["analysis_tier"],
                "confirmatory_role": focus["confirmatory_role"],
                "include_in_lmm": focus["include_in_lmm"],
                "transform": "rank_inverse_normal",
                "include": focus["include"],
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
