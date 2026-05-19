#!/usr/bin/env bash
set -euo pipefail

# Primary: plink2 on pgen
plink2 --pfile /net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/imp/ukb_imp vzs \
       --snp rs2546890 \
       --export A \
       --out /net/bmc-lab4/data/kellis/users/mabdel03/files/Isolation_Genetics/lnc_rna_mri/results/variant/rs2546890

# Coordinate fallback (if rsID miss)
plink2 --pfile /net/bmc-lab5/data/kellis/group/tanigawa/data/ukb21942/imp/ukb_imp vzs \
       --chr 5 --from-bp 159332892 --to-bp 159332892 \
       --export A \
       --out /net/bmc-lab4/data/kellis/users/mabdel03/files/Isolation_Genetics/lnc_rna_mri/results/variant/rs2546890_bycoord

