#!/usr/bin/env bash
#SBATCH -J lncrna_ols
#SBATCH -p kellis
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o logs/sbatch/lncrna_ols_%j.out
#SBATCH -e logs/sbatch/lncrna_ols_%j.err
#
# OLS-only pipeline for rs2546890 × UKBB MRI IDPs.
# Runs 00 → 06 with engine=ols. Includes composites + ROI PCA + per-family FDR.
# Wall-clock budget: ~30 min on the imaging subsample (39K). 8h is generous safety.
#
# Submit from the repo root:
#   sbatch scripts/sbatch/pipeline_ols.sbatch.sh
#
# Override engine on the command line (e.g. switch to regenie+ols) without editing config:
#   sbatch --export=ENGINE_OVERRIDE=regenie+ols scripts/sbatch/pipeline_ols.sbatch.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${REPO_ROOT}"
mkdir -p logs/sbatch

# Source the cluster's conda + activate analysis env (paths match the prior
# UKBB-SI-Genetics setup).
if command -v module >/dev/null 2>&1; then
    module load miniconda3/v4 || true
fi
if [ -f /home/software/conda/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/software/conda/miniconda3/etc/profile.d/conda.sh
elif command -v conda >/dev/null 2>&1; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
fi
conda activate /home/mabdel03/data/conda_envs/Python_Analysis

# Append GWAS_env bin for plink2 (active env doesn't have it; appending
# avoids shadowing Python_Analysis's `python`).
export PATH="${PATH}:/home/mabdel03/data/conda_envs/GWAS_env/bin"

# Optional engine override (default = whatever config.yaml says).
if [ -n "${ENGINE_OVERRIDE:-}" ]; then
    # Patch the engine line on the fly (don't edit config.yaml in place).
    TMP_CFG="$(mktemp --suffix=.yaml)"
    awk -v eng="${ENGINE_OVERRIDE}" '
        /^  engine:/ { print "  engine: \"" eng "\""; next } { print }
    ' config.yaml > "${TMP_CFG}"
    CFG_ARG="--config ${TMP_CFG}"
    trap 'rm -f "${TMP_CFG}"' EXIT
else
    CFG_ARG="--config config.yaml"
fi

echo "[$(date -Iseconds)] starting pipeline on $(hostname)"
echo "  engine: ${ENGINE_OVERRIDE:-from config.yaml}"
echo "  threads: ${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-4}}"
echo "  job id: ${SLURM_JOB_ID:-N/A}"

bash scripts/run_all.sh ${CFG_ARG}

echo "[$(date -Iseconds)] pipeline done"
