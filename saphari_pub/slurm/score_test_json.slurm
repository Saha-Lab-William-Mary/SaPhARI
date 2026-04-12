#!/bin/bash
#SBATCH --job-name=score_test_json
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/sciclone/home/ncnallapaneni/scr10/saphari_pub
#SBATCH --output=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%j.out
#SBATCH --error=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%j.err

set -euo pipefail

PROJ="/sciclone/home/ncnallapaneni/scr10/saphari_pub"
SCRIPT="${PROJ}/score.py"
ROOT="${PROJ}/output/test_json"

DETAILS="${ROOT}/score_details.tsv"
FAMILY_SUMMARY="${ROOT}/family_summary.tsv"
OVERALL_SUMMARY="${ROOT}/overall_summary.tsv"

mkdir -p "${PROJ}/slurm/slurm_logs"

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: score script not found at ${SCRIPT}"
    exit 1
fi

if [[ ! -d "${ROOT}" ]]; then
    echo "ERROR: output root not found at ${ROOT}"
    exit 1
fi

python3 "${SCRIPT}" "${ROOT}" \
  --write-details "${DETAILS}" \
  --write-family-summary "${FAMILY_SUMMARY}" \
  --write-overall-summary "${OVERALL_SUMMARY}"
