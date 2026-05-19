#!/usr/bin/env bash
#SBATCH -J lncrna_ols
#SBATCH -p kellis
#SBATCH -n 4
#SBATCH --mem=300G
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

set -eo pipefail   # NOTE: NOT -u — conda activate scripts reference unbound vars

# SLURM copies the submit script to a spool dir on the compute node, so
# BASH_SOURCE no longer points at the original — use SLURM_SUBMIT_DIR (the
# cwd from which sbatch was invoked, which must be the repo root) and
# fall back to BASH_SOURCE for interactive runs.
REPO_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${REPO_ROOT}"
mkdir -p logs/sbatch

# Source the cluster's conda init (Luria-specific: condainit, not conda.sh).
# Matches the UKBB-SI-Genetics SLURM pattern.
if command -v module >/dev/null 2>&1; then
    module load miniconda3/v4 || true
fi
if [ -f /home/software/conda/miniconda3/bin/condainit ]; then
    source /home/software/conda/miniconda3/bin/condainit
elif [ -f /home/software/conda/miniconda3/etc/profile.d/conda.sh ]; then
    source /home/software/conda/miniconda3/etc/profile.d/conda.sh
fi
conda activate /home/mabdel03/data/conda_envs/Python_Analysis

# Append GWAS_env bin for plink2 (active env doesn't have it; appending
# avoids shadowing Python_Analysis's `python`).
export PATH="${PATH}:/home/mabdel03/data/conda_envs/GWAS_env/bin"

# Optional engine override (default = whatever config.yaml says).
# IMPORTANT: keep the tempfile in the repo root, not /tmp — the pipeline
# derives `_repo_root` from the config-file parent dir, so a /tmp tempfile
# routes ALL output writes to /tmp/results/.
if [ -n "${ENGINE_OVERRIDE:-}" ]; then
    TMP_CFG="$(mktemp --tmpdir=. --suffix=.yaml runtime_config.XXXXXX)"
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
