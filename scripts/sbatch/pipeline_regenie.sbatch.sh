#!/usr/bin/env bash
#SBATCH -J lncrna_regenie
#SBATCH -p kellis
#SBATCH -n 16
#SBATCH --mem=200G
#SBATCH -t 48:00:00
#SBATCH -o logs/sbatch/lncrna_regenie_%j.out
#SBATCH -e logs/sbatch/lncrna_regenie_%j.err
#
# REGENIE+OLS pipeline for rs2546890 × UKBB MRI IDPs.
# Runs the full 00 → 06 chain with engine=regenie+ols.
# Per-phenotype REGENIE step 1 is the wall-clock bottleneck: ~15-25 min/pheno
# on this HM3 BED + 39K imaging subjects with 16 threads. 48h budget covers
# ~100-300 phenotypes; reduce lmm.max_phenotypes or scope_families in config
# if you want faster turnaround.
#
# Submit from the repo root:
#   sbatch scripts/sbatch/pipeline_regenie.sbatch.sh
#
# Skip OLS (only run REGENIE LMM):
#   sbatch --export=ENGINE_OVERRIDE=regenie scripts/sbatch/pipeline_regenie.sbatch.sh

set -eo pipefail   # NOTE: NOT -u — conda activate scripts reference unbound vars

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"
mkdir -p logs/sbatch

if command -v module >/dev/null 2>&1; then
    module load miniconda3/v4 || true
fi
if [ -f /home/software/conda/miniconda3/bin/condainit ]; then
    source /home/software/conda/miniconda3/bin/condainit
elif [ -f /home/software/conda/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/software/conda/miniconda3/etc/profile.d/conda.sh
fi
conda activate /home/mabdel03/data/conda_envs/Python_Analysis
export PATH="${PATH}:/home/mabdel03/data/conda_envs/GWAS_env/bin"

# REGENIE thread count from SLURM allocation
N_THREADS="${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-16}}"

# IMPORTANT: keep the tempfile in the repo root — the pipeline derives
# `_repo_root` from the config-file parent dir, so a /tmp tempfile would
# route ALL output writes to /tmp/results/.
TMP_CFG="$(mktemp --tmpdir=. --suffix=.yaml runtime_config.XXXXXX)"
if [ -n "${ENGINE_OVERRIDE:-}" ]; then
    awk -v eng="${ENGINE_OVERRIDE}" -v thr="${N_THREADS}" '
        /^  engine:/    { print "  engine: \"" eng "\""; next }
        /^  threads_step1:/ { print "  threads_step1: " thr; next }
        /^  threads_step2:/ { print "  threads_step2: " thr; next }
        { print }
    ' config.yaml > "${TMP_CFG}"
else
    # Patch threads to match SLURM allocation but keep engine from config
    awk -v thr="${N_THREADS}" '
        /^  threads_step1:/ { print "  threads_step1: " thr; next }
        /^  threads_step2:/ { print "  threads_step2: " thr; next }
        { print }
    ' config.yaml > "${TMP_CFG}"
fi
CFG_ARG="--config ${TMP_CFG}"
trap 'rm -f "${TMP_CFG}"' EXIT

echo "[$(date -Iseconds)] starting REGENIE pipeline on $(hostname)"
echo "  engine: ${ENGINE_OVERRIDE:-from config.yaml (default: regenie+ols)}"
echo "  threads: ${N_THREADS}"
echo "  job id: ${SLURM_JOB_ID:-N/A}"

bash scripts/run_all.sh ${CFG_ARG}

echo "[$(date -Iseconds)] REGENIE pipeline done"
