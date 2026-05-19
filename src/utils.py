"""Shared helpers for the rs2546890 → MRI IDP association pipeline.

All scripts import from here. Keep dependencies to the standard scientific
Python stack (pandas, numpy, scipy, statsmodels) so the module loads even
when matplotlib/seaborn are missing on a thin environment.
"""
from __future__ import annotations

import datetime as _dt
import logging
import os
import platform
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import yaml
from scipy import stats


# ---------------------------------------------------------------------------
# Config + paths
# ---------------------------------------------------------------------------

def load_config(path: str | os.PathLike) -> dict:
    """Read a YAML config file. Raises if path is missing."""
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"config file not found: {p}")
    with p.open("r") as fh:
        cfg = yaml.safe_load(fh) or {}
    # repo root = parent of the config file
    cfg.setdefault("_repo_root", str(p.resolve().parent))
    return cfg


def resolve_under_repo(cfg: dict, rel: str | os.PathLike) -> Path:
    """Treat relative paths in config as relative to the repo root."""
    p = Path(rel)
    if p.is_absolute():
        return p
    return Path(cfg["_repo_root"]) / p


def timestamped_outdir(base: str | os.PathLike, overwrite: bool) -> Path:
    """Return base if overwrite=True, else base/YYYYmmdd_HHMMSS/."""
    base = Path(base)
    base.mkdir(parents=True, exist_ok=True)
    if overwrite:
        return base
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out = base / stamp
    out.mkdir(parents=True, exist_ok=True)
    return out


def setup_logger(name: str, logfile: Optional[str | os.PathLike] = None,
                 level: int = logging.INFO) -> logging.Logger:
    logger = logging.getLogger(name)
    if logger.handlers:
        return logger
    logger.setLevel(level)
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                            "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if logfile is not None:
        Path(logfile).parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile, mode="a")
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


# ---------------------------------------------------------------------------
# Tool detection + version logging
# ---------------------------------------------------------------------------

def tool_available(name: str) -> Optional[str]:
    """Return the resolved path to a tool on PATH, or None."""
    return shutil.which(name)


def write_version_log(path: str | os.PathLike) -> None:
    """Snapshot interpreter + key package versions to a text file."""
    lines = [
        f"timestamp: {_dt.datetime.now().isoformat(timespec='seconds')}",
        f"python: {sys.version.splitlines()[0]}",
        f"platform: {platform.platform()}",
        f"executable: {sys.executable}",
        f"cwd: {os.getcwd()}",
    ]
    for pkg in ("pandas", "numpy", "scipy", "statsmodels", "sklearn",
                "matplotlib", "seaborn", "yaml", "pyarrow", "zstandard"):
        try:
            mod = __import__(pkg)
            ver = getattr(mod, "__version__", "unknown")
            lines.append(f"{pkg}: {ver}")
        except Exception as exc:  # noqa: BLE001
            lines.append(f"{pkg}: NOT INSTALLED ({exc.__class__.__name__})")
    for tool in ("plink2", "qctool", "bgenix", "regenie", "bcftools"):
        loc = tool_available(tool)
        lines.append(f"{tool}: {loc or 'not on PATH'}")
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Table IO
# ---------------------------------------------------------------------------

def _sniff_sep(path: Path) -> str:
    suf = path.suffix.lower()
    if suf in (".csv",):
        return ","
    if suf in (".tsv", ".tab", ".txt"):
        return "\t"
    # Sniff first line
    with path.open("r") as fh:
        head = fh.readline()
    if "\t" in head:
        return "\t"
    if "," in head:
        return ","
    return r"\s+"


def safe_read_table(path: str | os.PathLike,
                    usecols: Optional[Sequence[str]] = None,
                    chunksize: Optional[int] = None,
                    nrows: Optional[int] = None,
                    dtype=None,
                    sep: Optional[str] = None) -> pd.DataFrame | Iterable[pd.DataFrame]:
    """Read CSV/TSV/whitespace-delim with sensible defaults.

    - Sniffs separator from extension when sep is None.
    - na_values includes UKB's typical sentinels.
    - Returns a DataFrame, or an iterator if chunksize is set.
    """
    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(p)
    if sep is None:
        sep = _sniff_sep(p)
    na = ["", "NA", "N/A", "n/a", "NaN", "nan", ".", "-9", "-999"]
    kw = dict(sep=sep, na_values=na, low_memory=False)
    if usecols is not None:
        kw["usecols"] = list(usecols)
    if dtype is not None:
        kw["dtype"] = dtype
    if nrows is not None:
        kw["nrows"] = nrows
    if chunksize is not None:
        kw["chunksize"] = chunksize
    return pd.read_csv(p, **kw)


def read_header_only(path: str | os.PathLike, sep: Optional[str] = None) -> List[str]:
    """Return just the column names of a table without loading body."""
    p = Path(path)
    if sep is None:
        sep = _sniff_sep(p)
    with p.open("r") as fh:
        line = fh.readline().rstrip("\n").rstrip("\r")
    return line.split(sep) if sep != r"\s+" else line.split()


# ---------------------------------------------------------------------------
# UKB conventions
# ---------------------------------------------------------------------------

_EID_CANDIDATES = ("eid", "EID", "IID", "iid", "FID", "fid",
                   "f.eid", "ID_1", "ID")


def standardize_eid_column(df: pd.DataFrame, prefer: Optional[str] = None) -> pd.DataFrame:
    """Rename the first detected ID column to 'eid' and cast to int64 where possible."""
    if "eid" in df.columns:
        chosen = "eid"
    else:
        chosen = prefer if (prefer and prefer in df.columns) else None
        if chosen is None:
            for c in _EID_CANDIDATES:
                if c in df.columns:
                    chosen = c
                    break
        if chosen is None:
            raise KeyError(f"no recognised id column in {list(df.columns)[:10]}")
        df = df.rename(columns={chosen: "eid"})
    # Try numeric; fall back to string if it has non-numeric noise
    try:
        df["eid"] = pd.to_numeric(df["eid"], errors="raise").astype("Int64")
    except Exception:
        df["eid"] = df["eid"].astype(str)
    return df


def parse_ukb_field_id(colname: str) -> Optional[str]:
    """Extract '25526' from 'f.25526.2.0', '25526-2.0', 'phenotype_25526' etc."""
    m = re.search(r"(?:^|[^0-9])(\d{3,6})(?:-|\.|_)", colname + "-")
    if m:
        return m.group(1)
    m = re.search(r"phenotype_(\d{3,6})", colname)
    if m:
        return m.group(1)
    return None


def find_columns_by_metadata_terms(meta: pd.DataFrame,
                                   terms: Sequence[str],
                                   prioritize_regions: Optional[Sequence[str]] = None,
                                   description_col: str = "description",
                                   id_col: str = "field_id") -> pd.DataFrame:
    """Filter a metadata table (field_id, description) by case-insensitive term match.

    Returns a copy with an added `priority` column (0=base, +1 per region hit).
    """
    if description_col not in meta.columns:
        raise KeyError(f"meta must have a '{description_col}' column")
    desc = meta[description_col].fillna("").str.lower()
    pat = "|".join(re.escape(t.lower()) for t in terms) if terms else ""
    mask = desc.str.contains(pat, regex=True, na=False) if pat else pd.Series([False] * len(meta))
    out = meta.loc[mask].copy()
    out["priority"] = 0
    if prioritize_regions:
        for region in prioritize_regions:
            hit = out[description_col].fillna("").str.lower().str.contains(
                re.escape(region.lower()), regex=True)
            out.loc[hit, "priority"] += 1
    return out


# ---------------------------------------------------------------------------
# Phenotype transforms / stats helpers
# ---------------------------------------------------------------------------

def inverse_rank_normalize(x: pd.Series | np.ndarray, c: float = 3.0 / 8.0) -> pd.Series:
    """Blom-style rank-based inverse normal transform, NA-preserving."""
    x = pd.Series(x).astype(float).copy()
    n_obs = x.notna().sum()
    if n_obs == 0:
        return x
    ranks = x.rank(method="average")
    quantiles = (ranks - c) / (n_obs - 2 * c + 1)
    out = pd.Series(stats.norm.ppf(quantiles), index=x.index)
    out[x.isna()] = np.nan
    return out


def trim_outliers_z(x: pd.Series, z_thresh: float = 6.0) -> pd.Series:
    """Set |z|>z_thresh to NaN on the raw scale."""
    x = pd.Series(x).astype(float).copy()
    mu, sd = x.mean(skipna=True), x.std(skipna=True, ddof=1)
    if sd == 0 or not np.isfinite(sd):
        return x
    z = (x - mu) / sd
    x[(z.abs() > z_thresh)] = np.nan
    return x


def bh_fdr(pvals: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted q-values. NaN-safe."""
    p = np.asarray(pvals, dtype=float)
    n_total = len(p)
    q = np.full(n_total, np.nan)
    valid = ~np.isnan(p)
    pv = p[valid]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(1, n + 1))
    # enforce monotonicity from the right
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out_valid = np.empty(n)
    out_valid[order] = adj
    q[valid] = out_valid
    return q


# ---------------------------------------------------------------------------
# Design matrix + model fits
# ---------------------------------------------------------------------------

def make_design_matrix(df: pd.DataFrame,
                       cont_cols: Sequence[str],
                       cat_cols: Sequence[str],
                       add_intercept: bool = True,
                       drop_first: bool = True) -> pd.DataFrame:
    """Build a numeric design matrix from continuous + categorical cols."""
    pieces: list[pd.DataFrame] = []
    if add_intercept:
        pieces.append(pd.DataFrame({"const": 1.0}, index=df.index))
    if cont_cols:
        pieces.append(df[list(cont_cols)].astype(float))
    for c in cat_cols:
        if c not in df.columns:
            continue
        dummies = pd.get_dummies(df[c].astype("category"),
                                 prefix=c, drop_first=drop_first, dtype=float)
        pieces.append(dummies)
    return pd.concat(pieces, axis=1)


def fit_ols_model(y: pd.Series,
                  exposure: pd.Series,
                  covariates: pd.DataFrame,
                  exposure_name: str = "exposure") -> dict:
    """OLS of y on [exposure, covariates+const]. Returns a flat result dict.

    NaN rows in any of y/exposure/covariates are dropped jointly.
    """
    import statsmodels.api as sm

    frame = pd.concat([y.rename("__y__"), exposure.rename(exposure_name),
                       covariates], axis=1).dropna()
    n = len(frame)
    if n < 30:
        return {"status": "too_few_obs", "n": n}
    X = frame.drop(columns="__y__")
    Y = frame["__y__"]
    # rank check
    try:
        rank = np.linalg.matrix_rank(X.values)
    except Exception:
        rank = X.shape[1]
    if rank < X.shape[1]:
        # drop perfectly collinear cols
        _, idx = np.linalg.qr(X.values, mode="reduced")
        # simpler: rely on statsmodels' pinv; just warn
        pass
    try:
        model = sm.OLS(Y.values, X.values, hasconst=True).fit()
    except Exception as exc:  # noqa: BLE001
        return {"status": f"fit_error:{exc.__class__.__name__}", "n": n}
    cols = list(X.columns)
    if exposure_name not in cols:
        return {"status": "exposure_missing_post_design", "n": n}
    i = cols.index(exposure_name)
    return {
        "status": "ok",
        "n": int(n),
        "beta": float(model.params[i]),
        "se": float(model.bse[i]),
        "t": float(model.tvalues[i]),
        "p": float(model.pvalues[i]),
        "r2": float(model.rsquared),
        "rank": int(rank),
        "ncols": int(X.shape[1]),
    }


def fit_genotypic_model(y: pd.Series,
                        genotype: pd.Series,
                        covariates: pd.DataFrame,
                        reference: str = "GG") -> dict:
    """Fit y ~ C(genotype) + covariates. genotype is a string column ('AA','AG','GG').

    Returns AG-vs-ref and AA-vs-ref contrasts plus 2-df Wald p.
    """
    import statsmodels.api as sm

    g = pd.Series(genotype).astype("category")
    # ensure reference level
    cats = list(g.cat.categories)
    if reference not in cats:
        return {"status": f"reference_{reference}_absent", "categories": cats}
    g = g.cat.reorder_categories([reference] + [c for c in cats if c != reference],
                                 ordered=False)
    dummies = pd.get_dummies(g, prefix="geno", drop_first=True, dtype=float)
    frame = pd.concat([y.rename("__y__"), dummies, covariates], axis=1).dropna()
    n = len(frame)
    if n < 30:
        return {"status": "too_few_obs", "n": n}
    X = frame.drop(columns="__y__")
    Y = frame["__y__"]
    try:
        model = sm.OLS(Y.values, X.values, hasconst=True).fit()
    except Exception as exc:  # noqa: BLE001
        return {"status": f"fit_error:{exc.__class__.__name__}", "n": n}
    cols = list(X.columns)
    out: dict = {"status": "ok", "n": int(n), "reference": reference}
    for level in [c for c in cats if c != reference]:
        name = f"geno_{level}"
        if name in cols:
            i = cols.index(name)
            out[f"beta_{level}_vs_{reference}"] = float(model.params[i])
            out[f"se_{level}_vs_{reference}"] = float(model.bse[i])
            out[f"p_{level}_vs_{reference}"] = float(model.pvalues[i])
    # 2-df Wald: joint test of both genotype dummies
    geno_idx = [cols.index(c) for c in cols if c.startswith("geno_")]
    if len(geno_idx) >= 1:
        R = np.zeros((len(geno_idx), len(cols)))
        for r, i in enumerate(geno_idx):
            R[r, i] = 1.0
        try:
            wald = model.wald_test(R, use_f=False)
            out["wald_df"] = int(len(geno_idx))
            out["wald_chi2"] = float(np.asarray(wald.statistic).ravel()[0])
            out["wald_p"] = float(np.asarray(wald.pvalue).ravel()[0])
        except Exception as exc:  # noqa: BLE001
            out["wald_status"] = f"wald_error:{exc.__class__.__name__}"
    return out


# ---------------------------------------------------------------------------
# Hard-call assignment from imputed dosage
# ---------------------------------------------------------------------------

def dosage_to_hardcall(dosage_A: pd.Series,
                       effect_allele: str = "A",
                       other_allele: str = "G",
                       threshold: float = 0.9) -> pd.Series:
    """Assign hard genotypes from imputed dosage iff dosage is within `threshold`
    of an integer. Otherwise NaN.

    Returns strings: 'AA' (2 effect copies), 'AG', 'GG'.
    """
    d = pd.to_numeric(dosage_A, errors="coerce")
    aa_label = effect_allele * 2
    ag_label = "".join(sorted([effect_allele, other_allele]))
    gg_label = other_allele * 2
    out = pd.Series(np.nan, index=d.index, dtype=object)
    out[(d - 0).abs() <= (1 - threshold)] = gg_label
    out[(d - 1).abs() <= (1 - threshold)] = ag_label
    out[(d - 2).abs() <= (1 - threshold)] = aa_label
    return out


# ---------------------------------------------------------------------------
# Subprocess helpers
# ---------------------------------------------------------------------------

def run_cmd(cmd: Sequence[str] | str, logger: Optional[logging.Logger] = None,
            check: bool = True, env: Optional[dict] = None) -> subprocess.CompletedProcess:
    if isinstance(cmd, str):
        printable = cmd
        shell = True
    else:
        printable = " ".join(str(x) for x in cmd)
        shell = False
    if logger:
        logger.info(f"$ {printable}")
    proc = subprocess.run(cmd, shell=shell, check=False, env=env,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)
    if logger and proc.stdout:
        for line in proc.stdout.rstrip().splitlines():
            logger.info(f"  | {line}")
    if logger and proc.stderr:
        for line in proc.stderr.rstrip().splitlines():
            logger.info(f"  ! {line}")
    if check and proc.returncode != 0:
        raise RuntimeError(f"command failed (rc={proc.returncode}): {printable}")
    return proc


__all__ = [
    "bh_fdr",
    "dosage_to_hardcall",
    "find_columns_by_metadata_terms",
    "fit_genotypic_model",
    "fit_ols_model",
    "inverse_rank_normalize",
    "load_config",
    "make_design_matrix",
    "parse_ukb_field_id",
    "read_header_only",
    "resolve_under_repo",
    "run_cmd",
    "safe_read_table",
    "set_seed",
    "setup_logger",
    "standardize_eid_column",
    "timestamped_outdir",
    "tool_available",
    "trim_outliers_z",
    "write_version_log",
]
