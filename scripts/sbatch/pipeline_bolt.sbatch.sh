#!/usr/bin/env bash
#SBATCH -J lncrna_bolt
#SBATCH -p kellis
#SBATCH -n 32
#SBATCH --mem=150G
#SBATCH -t 48:00:00
#SBATCH -o logs/sbatch/lncrna_bolt_%j.out
#SBATCH -e logs/sbatch/lncrna_bolt_%j.err
#
# BOLT-LMM (+ OLS) pipeline for rs2546890 × UKBB MRI IDPs.
# Runs 00 → 06 with engine=bolt+ols (BOLT for composites + ROI PCs only by
# default; flip lmm.bolt.include_single_idps: true in config to expand).
#
# BOLT processes one phenotype per invocation, so cost scales linearly. For
# the default scope (composites + PCs ≈ 316 phenotypes) at ~30-90 min each
# on 32 threads, 48h covers ~30-40 phenotypes. Set lmm.max_phenotypes (in
# config) to constrain scope.
#
# Submit from the repo root:
#   sbatch scripts/sbatch/pipeline_bolt.sbatch.sh

set -eo pipefail   # NOTE: NOT -u — conda activate scripts reference unbound vars

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
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
# bolt_lmm conda env has BOLT-LMM binary + plink2
export PATH="${PATH}:/home/mabdel03/data/conda_envs/bolt_lmm/bin:/home/mabdel03/data/conda_envs/GWAS_env/bin"

N_THREADS="${SLURM_CPUS_ON_NODE:-${SLURM_NTASKS:-32}}"

# Force engine=bolt+ols + patch BOLT threads to SLURM allocation
TMP_CFG="$(mktemp --suffix=.yaml)"
awk -v thr="${N_THREADS}" '
    /^  engine:/      { print "  engine: \"bolt+ols\""; next }
    /^    threads:/   { print "    threads: " thr; next }
    { print }
' config.yaml > "${TMP_CFG}"
trap 'rm -f "${TMP_CFG}"' EXIT

echo "[$(date -Iseconds)] starting BOLT-LMM pipeline on $(hostname)"
echo "  engine: bolt+ols"
echo "  threads: ${N_THREADS}"
echo "  job id: ${SLURM_JOB_ID:-N/A}"

bash scripts/run_all.sh --config "${TMP_CFG}"

echo "[$(date -Iseconds)] BOLT-LMM pipeline done"
