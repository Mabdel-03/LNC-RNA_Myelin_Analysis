#!/usr/bin/env bash
# Orchestrator: runs 00 → 06 sequentially, exits on first failure.
# Usage:
#   bash scripts/run_all.sh --config config.yaml [--dry-run] [--test-mode]
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG="${REPO_ROOT}/config.yaml"
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --dry-run) EXTRA_ARGS+=("--dry-run"); shift ;;
    --test-mode) EXTRA_ARGS+=("--test-mode"); shift ;;  # forwarded only to steps that honor it
    -h|--help)
      sed -n '2,4p' "${BASH_SOURCE[0]}"
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ ! -f "${CONFIG}" ]]; then
  echo "config not found: ${CONFIG}" >&2
  exit 1
fi

cd "${REPO_ROOT}"
mkdir -p logs results figures results/intermediate results/variant

LOG="logs/run_log.txt"
echo "==============================================================" >> "${LOG}"
echo "[$(date -Iseconds)] run_all.sh START  config=${CONFIG}  args=${EXTRA_ARGS[*]+${EXTRA_ARGS[*]}}" >> "${LOG}"

# Activate conda env if available and configured
CONDA_ENV="$(awk -F'"' '/conda_env:/ {print $2; exit}' "${CONFIG}" 2>/dev/null || true)"
if [[ -n "${CONDA_ENV:-}" ]] && command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh" || true
  if conda env list | awk '{print $1}' | grep -qx "${CONDA_ENV}"; then
    conda activate "${CONDA_ENV}"
    echo "[$(date -Iseconds)] activated conda env: ${CONDA_ENV}" >> "${LOG}"
  else
    echo "[$(date -Iseconds)] conda env '${CONDA_ENV}' not found; using current python" >> "${LOG}"
  fi
fi

PY="${PYTHON:-python}"

# 01-04 honor --dry-run (no --test-mode forwarding needed).
# 05 has its own test_mode flag in config (we leave it to read config).
FILTER_ARGS=()
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for a in "${EXTRA_ARGS[@]}"; do
    [[ "$a" == "--test-mode" ]] && continue
    FILTER_ARGS+=("$a")
  done
fi

for step in \
  "src/00_inspect_inputs.py" \
  "src/01_extract_variant.py" \
  "src/02_build_phenotype_matrix.py" \
  "src/03_build_covariates.py" \
  "src/04_merge_qc_sample.py" \
  "src/05_run_association.py" \
  "src/06_plots_and_report.py"
do
  echo "[$(date -Iseconds)] >>> ${step}" | tee -a "${LOG}"
  "${PY}" "${step}" --config "${CONFIG}" ${FILTER_ARGS[@]+"${FILTER_ARGS[@]}"} 2>&1 | tee -a "${LOG}"
  rc=${PIPESTATUS[0]}
  if [[ ${rc} -ne 0 ]]; then
    echo "[$(date -Iseconds)] STEP FAILED (rc=${rc}): ${step}" | tee -a "${LOG}"
    exit ${rc}
  fi
done

echo "[$(date -Iseconds)] run_all.sh DONE" | tee -a "${LOG}"
