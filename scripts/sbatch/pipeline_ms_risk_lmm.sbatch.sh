#!/usr/bin/env bash
#SBATCH -J lncrna_ms_risk_lmm
#SBATCH -p kellis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o logs/sbatch/lncrna_ms_risk_lmm_%j.out
#SBATCH -e logs/sbatch/lncrna_ms_risk_lmm_%j.err
#
# REGENIE binary-trait LMM for MS (G35) risk on rs2546890.
# Step-1 cache lives at results/lmm/ms_step1_*, isolated from the QT pipeline.
#
# Submit from the repo root:
#   sbatch scripts/sbatch/pipeline_ms_risk_lmm.sbatch.sh

set -eo pipefail   # NOTE: NOT -u — conda activate scripts reference unbound vars

REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"
mkdir -p logs/sbatch results/lmm results/lmm_inputs_ms

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

# Patch step-1 thread count from SLURM allocation into runtime tempfile.
# Keep the tempfile in the repo root so config-relative output paths still
# resolve under this checkout.
CFG=$(mktemp --tmpdir=. --suffix=.yaml runtime_config.XXXXXX)
cp config.yaml "${CFG}"
trap 'rm -f "${CFG}"' EXIT
N_CPU="${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-16}}}"
export PLINK2_THREADS="${PLINK2_THREADS:-${N_CPU}}"
export OLS_BLAS_THREADS="${OLS_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${OLS_BLAS_THREADS}"
export OPENBLAS_NUM_THREADS="${OLS_BLAS_THREADS}"
export MKL_NUM_THREADS="${OLS_BLAS_THREADS}"
export NUMEXPR_NUM_THREADS="${OLS_BLAS_THREADS}"
export VECLIB_MAXIMUM_THREADS="${OLS_BLAS_THREADS}"
# only update threads_step1; threads_step2 isn't critical for 1-variant step 2
sed -i "s/^  threads_step1: .*/  threads_step1: ${N_CPU}/" "${CFG}" || true

echo "[$(date -Iseconds)] starting MS-risk LMM pipeline (cfg=${CFG}, threads=${N_CPU})"

# 07: derive MS phenotype (or reuse if exists) + OLS-side logistic
python src/07_run_ms_risk.py --config "${CFG}" --reuse-phenotype

# 07b: REGENIE step 1 (BT, Firth) + step 2 for rs2546890
python src/07b_run_ms_risk_lmm.py --config "${CFG}"

echo "[$(date -Iseconds)] DONE"
