#!/bin/bash
#SBATCH --job-name=package_build_benchmark
#SBATCH --partition=bigmem
#SBATCH --output=package_build_benchmark_%j.out
#SBATCH --error=package_build_benchmark_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=900G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=ALL

set -euo pipefail

: "${CONFIG:?Set CONFIG to an absolute benchmark YAML path}"
: "${PYTHON:?Set PYTHON to the PACKAGE environment Python executable}"

PACKAGE_ROOT=${PACKAGE_ROOT:-${SLURM_SUBMIT_DIR}}
RESULTS_DIR=${RESULTS_DIR:-${SLURM_SUBMIT_DIR}/benchmark/build_${SLURM_JOB_ID}}
ALLOW_OVERWRITE=${ALLOW_OVERWRITE:-0}

mkdir -p "${RESULTS_DIR}"
cd "${PACKAGE_ROOT}"

DB_PATH=$("${PYTHON}" -c \
  'import sys; from PACKAGE.config import load_config; print(load_config(sys.argv[1]).output_path)' \
  "${CONFIG}")

if [[ -e "${DB_PATH}" && "${ALLOW_OVERWRITE}" != "1" ]]; then
  echo "Refusing to overwrite existing benchmark output: ${DB_PATH}" >&2
  echo "Use a new output_file in the benchmark YAML." >&2
  exit 2
fi

RESOURCE_LOG=${RESULTS_DIR}/build_resources_${SLURM_JOB_ID}.txt
REPORT=${RESULTS_DIR}/build_report_${SLURM_JOB_ID}.json
GIT_COMMIT=$(git rev-parse HEAD)
STARTED_AT=$(date --iso-8601=seconds)
START_EPOCH=$(date +%s)

echo "PACKAGE root: ${PACKAGE_ROOT}"
echo "Git commit: ${GIT_COMMIT}"
echo "Python: ${PYTHON}"
echo "Config: ${CONFIG}"
echo "Database: ${DB_PATH}"
echo "Results: ${RESULTS_DIR}"
"${PYTHON}" --version
"${PYTHON}" -m PACKAGE.cli --version

/usr/bin/time -v -o "${RESOURCE_LOG}" \
  "${PYTHON}" -m PACKAGE.cli build --config "${CONFIG}"

FINISHED_AT=$(date --iso-8601=seconds)
FINISHED_EPOCH=$(date +%s)
ELAPSED_SECONDS=$((FINISHED_EPOCH - START_EPOCH))

"${PYTHON}" -m PACKAGE.benchmark.build_report \
  --config "${CONFIG}" \
  --resource-log "${RESOURCE_LOG}" \
  --out "${REPORT}" \
  --started-at "${STARTED_AT}" \
  --finished-at "${FINISHED_AT}" \
  --elapsed-seconds "${ELAPSED_SECONDS}" \
  --git-commit "${GIT_COMMIT}"

"${PYTHON}" -m PACKAGE.cli info "${DB_PATH}"
echo "Build benchmark complete: ${REPORT}"
