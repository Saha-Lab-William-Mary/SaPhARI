#!/bin/bash
#SBATCH --job-name=saphari_extract
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --chdir=/sciclone/home/ncnallapaneni/scr10/saphari_pub
#SBATCH --output=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%j.out
#SBATCH --error=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%j.err

set -euo pipefail

PROJ="/sciclone/home/ncnallapaneni/scr10/saphari_pub"
SCRIPT="${PROJ}/extract_clusters.py"
ANNOTATED_ROOT="${PROJ}/annotated"
OUT_DIR="${PROJ}/output/extraction"
JSON_NAME="Families_final.json"

mkdir -p "${OUT_DIR}"
mkdir -p "${PROJ}/slurm/slurm_logs"

if [[ ! -f "${SCRIPT}" ]]; then
    echo "ERROR: script not found at ${SCRIPT}"
    exit 1
fi

if [[ ! -d "${ANNOTATED_ROOT}" ]]; then
    echo "ERROR: annotated root not found at ${ANNOTATED_ROOT}"
    exit 1
fi

python3 "${SCRIPT}" \
  --annotated-root "${ANNOTATED_ROOT}" \
  --out-dir "${OUT_DIR}" \
  --emit-families-json "${JSON_NAME}"
