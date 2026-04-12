#!/bin/bash
#SBATCH --job-name=saphari_final
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ncnallapaneni@wm.edu
#SBATCH --array=0-5%4
#SBATCH --chdir=/sciclone/home/ncnallapaneni/scr10/saphari_pub
#SBATCH --output=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%A_%a.out
#SBATCH --error=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%A_%a.err

set -euo pipefail

PROJ="/sciclone/home/ncnallapaneni/scr10/saphari_pub"
IMG="${PROJ}/SaPhARI/containers/saphari_0.1.0.sif"
FAMILIES_JSON="${PROJ}/output/extraction/Families_final.json"

DB_HOST="/sciclone/scr20/ncnallapaneni/create_db/files/SaPhARI_11_18_25.dmnd"
DB_CONT="/scr20/ncnallapaneni/create_db/files/SaPhARI_11_18_25.dmnd"

FAMS=(
  CFPICI
  EPIP
  MFMR
  P4
  PHANIE
  PICI
)

TASK_ID="${SLURM_ARRAY_TASK_ID}"

if (( TASK_ID < 0 || TASK_ID >= ${#FAMS[@]} )); then
    echo "ERROR: SLURM_ARRAY_TASK_ID=$TASK_ID is out of range for ${#FAMS[@]} folders"
    exit 1
fi

FAM="${FAMS[$TASK_ID]}"
DATA_DIR="/work/annotated/${FAM}/annotated"
OUT="${PROJ}/output/test_json/out_${FAM}"

mkdir -p "${PROJ}/slurm/slurm_logs"
mkdir -p "${OUT}/work" "${OUT}/.nextflow/plugins" "${OUT}/tmp"

if [[ ! -f "${IMG}" ]]; then
    echo "ERROR: image not found at ${IMG}"
    exit 1
fi

if [[ ! -f "${FAMILIES_JSON}" ]]; then
    echo "ERROR: families JSON not found at ${FAMILIES_JSON}"
    exit 1
fi

if [[ ! -f "${DB_HOST}" ]]; then
    echo "ERROR: DB not found at ${DB_HOST}"
    exit 1
fi

if [[ ! -d "${PROJ}/annotated/${FAM}/annotated" ]]; then
    echo "ERROR: annotated directory missing for ${FAM} at ${PROJ}/annotated/${FAM}/annotated"
    exit 1
fi

export SINGULARITYENV_NXF_HOME=/out
export SINGULARITYENV_NXF_ASSETS=/out/.nextflow/assets
export SINGULARITYENV_NXF_TEMP=/out/tmp
export SINGULARITYENV_TMPDIR=/out/tmp

singularity exec --home "${OUT}" \
  -B "${PROJ}":/work \
  -B "${OUT}":/out \
  -B /sciclone/scr20:/scr20 \
  "${IMG}" \
  saphari \
    --data_dir "${DATA_DIR}" \
    --db "${DB_CONT}" \
    --out_dir /out \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --memory "64 GB" \
    --families_file "/work/output/extraction/$(basename "${FAMILIES_JSON}")" \
    --reanalyze
