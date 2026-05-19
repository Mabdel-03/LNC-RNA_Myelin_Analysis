#!/usr/bin/env bash
# Orchestrator: runs 00 → 06 sequentially, exits on first failure.
# Engine selection is read from config.yaml: models.engine ∈
#   ols | regenie | bolt | regenie+ols | bolt+ols | all
# Usage:
#   bash scripts/run_all.sh --config config.yaml [--dry-run] [--test-mode] [--with-tests]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${REPO_ROOT}/config.yaml"
EXTRA_ARGS=()
RUN_TESTS=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --dry-run) EXTRA_ARGS+=("--dry-run"); shift ;;
    --test-mode) EXTRA_ARGS+=("--test-mode"); shift ;;
    --with-tests) RUN_TESTS=1; shift ;;
    -h|--help)
      sed -n '2,7p' "${BASH_SOURCE[0]}"
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -f "${CONFIG}" ]]; then
  echo "config not found: ${CONFIG}" >&2
  exit 1
fi

cd "${REPO_ROOT}"
mkdir -p logs results figures results/intermediate results/variant results/lmm results/lmm_inputs

LOG="logs/run_log.txt"
echo "==============================================================" >> "${LOG}"
echo "[$(date -Iseconds)] run_all.sh START  config=${CONFIG}  args=${EXTRA_ARGS[*]+${EXTRA_ARGS[*]}}" >> "${LOG}"

# Activate conda env if available and configured
CONDA_ENV="$(awk -F'"' '/conda_env:/ {print $2; exit}' "${CONFIG}" 2>/dev/null || true)"
if [[ -n "${CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" || true
  # Match by basename of the env path (or exact name for system envs)
  if conda env list | awk '{print $1}' | sed 's:.*/::' | grep -qx "${CONDA_ENV}"; then
    conda activate "${CONDA_ENV}"
    echo "[$(date -Iseconds)] activated conda env: ${CONDA_ENV}" >> "${LOG}"
  else
    echo "[$(date -Iseconds)] conda env '${CONDA_ENV}' not found; using current python" >> "${LOG}"
  fi
fi

PY="${PYTHON:-python}"

# Append the GWAS_env conda bin so plink2 is visible to step 01 / 05c even
# when the active env is Python_Analysis (which lacks plink2). Appending (not
# prepending) means GWAS_env's `python` doesn't shadow the active env's Python.
GWAS_BIN="/home/mabdel03/data/conda_envs/GWAS_env/bin"
[[ -d "${GWAS_BIN}" ]] && export PATH="${PATH}:${GWAS_BIN}"

# Strip --test-mode out before forwarding (only 05a honors it via config flag)
FILTER_ARGS=()
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for a in "${EXTRA_ARGS[@]}"; do
    [[ "$a" == "--test-mode" ]] && continue
    FILTER_ARGS+=("$a")
  done
fi

run_step() {
  local step="$1"
  echo "[$(date -Iseconds)] >>> ${step}" | tee -a "${LOG}"
  "${PY}" "${step}" --config "${CONFIG}" ${FILTER_ARGS[@]+"${FILTER_ARGS[@]}"} 2>&1 | tee -a "${LOG}"
  local rc=${PIPESTATUS[0]}
  if [[ ${rc} -ne 0 ]]; then
    echo "[$(date -Iseconds)] STEP FAILED (rc=${rc}): ${step}" | tee -a "${LOG}"
    exit ${rc}
  fi
}

# Read engine from config (models.engine: "regenie+ols")
ENGINE="$(awk -F'"' '/^  engine:/ {print $2; exit}' "${CONFIG}" 2>/dev/null || true)"
ENGINE="${ENGINE:-ols}"
echo "[$(date -Iseconds)] engine=${ENGINE}" | tee -a "${LOG}"

# 00 → 04: input inspection, variant extract, phenotype matrix, covariates, QC join
run_step "src/00_inspect_inputs.py"
run_step "src/01_extract_variant.py"
run_step "src/02_build_phenotype_matrix.py"
run_step "src/03_build_covariates.py"
run_step "src/04_merge_qc_sample.py"

# 04b: derive composites + ROI PCA
run_step "src/04b_derive_composites_and_pca.py"

# 04c: prepare LMM input files (only if any LMM engine is going to run)
case "${ENGINE}" in
  regenie|bolt|regenie+ols|bolt+ols|all)
    run_step "src/04c_prepare_lmm_inputs.py"
    ;;
esac

# 05: route based on engine
case "${ENGINE}" in
  ols)
    run_step "src/05a_run_ols.py"
    ;;
  regenie)
    run_step "src/05b_run_regenie.py"
    run_step "src/05d_ingest_lmm_results.py"
    ;;
  bolt)
    run_step "src/05c_run_bolt.py"
    run_step "src/05d_ingest_lmm_results.py"
    ;;
  regenie+ols)
    run_step "src/05a_run_ols.py"
    run_step "src/05b_run_regenie.py"
    run_step "src/05d_ingest_lmm_results.py"
    ;;
  bolt+ols)
    run_step "src/05a_run_ols.py"
    run_step "src/05c_run_bolt.py"
    run_step "src/05d_ingest_lmm_results.py"
    ;;
  all)
    run_step "src/05a_run_ols.py"
    run_step "src/05b_run_regenie.py"
    run_step "src/05c_run_bolt.py"
    run_step "src/05d_ingest_lmm_results.py"
    ;;
  *)
    echo "unknown engine: ${ENGINE}" >&2; exit 1 ;;
esac

# 06: plots + reports
run_step "src/06_plots_and_report.py"

# Optional: run synthetic-data tests
if [[ ${RUN_TESTS} -eq 1 ]]; then
  echo "[$(date -Iseconds)] >>> pytest tests/" | tee -a "${LOG}"
  "${PY}" -m pytest tests/ -v 2>&1 | tee -a "${LOG}"
fi

echo "[$(date -Iseconds)] run_all.sh DONE" | tee -a "${LOG}"
