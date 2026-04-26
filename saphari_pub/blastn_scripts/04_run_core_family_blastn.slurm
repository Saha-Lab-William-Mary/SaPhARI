#!/bin/bash
#SBATCH --job-name=core_family_blastn
#SBATCH --output=/sciclone/scr10/ncnallapaneni/saphari_pub/slurm/logs/core_family_blastn_%A_%a.out
#SBATCH --error=/sciclone/scr10/ncnallapaneni/saphari_pub/slurm/logs/core_family_blastn_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --array=0-5

set -euo pipefail

REAL="/sciclone/scr10/ncnallapaneni"
SIF="$REAL/saphari_0.1.0.sif"

PROJECT="$REAL/saphari_pub"
QUERY_DIR="$PROJECT/output/blast_core_queries"
DB_DIR="$REAL/satellite_blastdb"
OUT_DIR="$PROJECT/output/blast_core_results"
TMP_DIR="$PROJECT/output/blast_core_tmp"

mkdir -p "$OUT_DIR" "$TMP_DIR" "$PROJECT/slurm/logs"

families=(CFPICI PICI EPIP MFMR P4 PHANIE)
fam="${families[$SLURM_ARRAY_TASK_ID]}"

query="$QUERY_DIR/${fam}_core_queries.fasta"
db="$DB_DIR/$fam"
out="$OUT_DIR/${fam}_core_vs_refs.tsv"
tmp="$TMP_DIR/${fam}_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$tmp"

echo "Family: $fam"
echo "Query: $query"
echo "DB: $db"
echo "Out: $out"

if [[ ! -s "$query" ]]; then
    echo "No core query file for $fam or file is empty. Exiting."
    exit 0
fi

if [[ ! -f "${db}.nsq" ]]; then
    echo "BLAST DB for $fam not found at $db"
    exit 1
fi

singularity exec --bind "$REAL":"$REAL" "$SIF" blastn \
    -query "$query" \
    -db "$db" \
    -out "$out" \
    -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
    -evalue 1e-10 \
    -num_threads "$SLURM_CPUS_PER_TASK" \
    -max_target_seqs 10 \
    -word_size 11 \
    -dust no

echo "Finished $fam"
