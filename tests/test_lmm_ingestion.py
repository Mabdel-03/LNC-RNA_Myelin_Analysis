"""parse_regenie_output + parse_bolt_stats on small text fixtures.

We write the fixtures inline as tempfiles so the test is self-contained and
doesn't require any real REGENIE/BOLT binary to be installed.
"""
import gzip
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from utils import parse_bolt_stats, parse_regenie_output


REGENIE_TEXT = """\
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
5 158759890 5:158759890:T:C T C 0.234 0.95 39000 ADD 0.005 0.013 0.15 0.20 NA
5 158759900 5:158759900:A:G G A 0.516 1.00 40000 ADD 0.024 0.011 4.75 1.30 NA
5 158759910 5:158759910:G:A G A 0.110 0.91 40000 ADD -0.018 0.014 1.65 0.55 NA
"""


def test_parse_regenie_target_extraction(tmp_path):
    p = tmp_path / "regenie_step2_PHENO.regenie.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(REGENIE_TEXT)
    df = parse_regenie_output(p, target_variant_id="5:158759900:A:G")
    assert len(df) == 1
    row = df.iloc[0]
    assert row["chrom"] == 5 or str(row["chrom"]) == "5"
    assert int(row["pos"]) == 158759900
    assert abs(row["beta"] - 0.024) < 1e-9
    assert abs(row["se"] - 0.011) < 1e-9
    # log10p was given; p should be derived
    assert "p" in df.columns
    assert abs(row["p"] - 10 ** (-1.30)) < 1e-9


def test_parse_regenie_no_target_returns_all(tmp_path):
    p = tmp_path / "all.regenie"
    p.write_text(REGENIE_TEXT)
    df = parse_regenie_output(p, target_variant_id=None)
    assert len(df) == 3


BOLT_TEXT = """\
SNP CHR BP GENPOS ALLELE1 ALLELE0 A1FREQ INFO BETA SE CHISQ_BOLT_LMM P_BOLT_LMM
rs9999999 5 158759890 0.5 C T 0.234 0.95 0.005 0.013 0.15 0.7
rs2546890 5 158759900 0.5 A G 0.516 1.00 0.024 0.011 4.75 0.029
rs8888888 5 158759910 0.5 A G 0.110 0.91 -0.018 0.014 1.65 0.20
"""


def test_parse_bolt_target_by_rsid(tmp_path):
    p = tmp_path / "bolt_PHENO.stats"
    p.write_text(BOLT_TEXT)
    df = parse_bolt_stats(p, target_rsid="rs2546890")
    assert len(df) == 1
    row = df.iloc[0]
    assert abs(row["beta"] - 0.024) < 1e-9
    assert abs(row["p"] - 0.029) < 1e-9
    assert np.isfinite(row["log10p"])


def test_parse_bolt_target_by_chrpos(tmp_path):
    p = tmp_path / "bolt_PHENO.stats"
    p.write_text(BOLT_TEXT)
    df = parse_bolt_stats(p, target_chrpos=(5, 158759900))
    assert len(df) == 1
    assert df.iloc[0]["variant_id"] == "rs2546890"
