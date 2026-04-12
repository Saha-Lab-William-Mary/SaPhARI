#!/bin/bash
#SBATCH --job-name=saphari_fasta
#SBATCH --array=1-9
#SBATCH --cpus-per-task=64
#SBATCH --mem=128G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ncnallapaneni@wm.edu
#SBATCH --chdir=/sciclone/home/ncnallapaneni/scr10/saphari_pub
#SBATCH --output=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%A_%a.out
#SBATCH --error=/sciclone/home/ncnallapaneni/scr10/saphari_pub/slurm/slurm_logs/%x.%A_%a.err

set -euo pipefail

PROJ="/sciclone/home/ncnallapaneni/scr10/saphari_pub"
IMG="${PROJ}/SaPhARI/containers/saphari_0.1.0.sif"
FAMILIES_JSON="${PROJ}/output/extraction/Families_final.json"

DB_HOST="/sciclone/scr20/ncnallapaneni/create_db/files/SaPhARI_11_18_25.dmnd"
DB_CONT="/scr20/ncnallapaneni/create_db/files/SaPhARI_11_18_25.dmnd"

DATA_HOST_ROOT="${PROJ}/assembled_genomes/"
DATA_CONT_ROOT="/work/assembled_genomes/"
OUT_ROOT="${PROJ}/output/new_genomes_fasta"

CHUNK=$(printf "%03d" "${SLURM_ARRAY_TASK_ID}")
DATA_HOST="${DATA_HOST_ROOT}/${CHUNK}"
DATA_CONT="${DATA_CONT_ROOT}/${CHUNK}"
OUT="${OUT_ROOT}/${CHUNK}"

mkdir -p "${PROJ}/slurm/slurm_logs"
mkdir -p "${OUT}/work" "${OUT}/.nextflow/plugins" "${OUT}/tmp"

echo "Running array task ${SLURM_ARRAY_TASK_ID} on chunk ${CHUNK}"
echo "Input: ${DATA_HOST}"
echo "Output: ${OUT}"

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

if [[ ! -d "${DATA_HOST}" ]]; then
    echo "ERROR: FASTA folder not found at ${DATA_HOST}"
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
    --workflow m_n_t_a \
    --data_dir "${DATA_CONT}" \
    --db "${DB_CONT}" \
    --out_dir /out \
    --cpus "${SLURM_CPUS_PER_TASK}" \
    --memory "128 GB" \
    --families_file "/work/output/extraction/$(basename "${FAMILIES_JSON}")"
