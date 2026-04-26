#!/bin/bash
set -euo pipefail

REAL="/sciclone/scr10/ncnallapaneni"
SIF="$REAL/saphari_0.1.0.sif"

REF_DIR="$REAL/satellite_fasta_db"
DB_DIR="$REAL/satellite_blastdb"

mkdir -p "$REF_DIR" "$DB_DIR"

families=(CFPICI PICI EPIP MFMR P4 PHANIE)

for fam in "${families[@]}"; do
    ref="$REF_DIR/${fam}_refs.fasta"
    db="$DB_DIR/$fam"

    echo "============================================================"
    echo "Family: $fam"
    echo "Reference FASTA: $ref"
    echo "BLAST DB prefix: $db"
    echo "============================================================"

    if [[ ! -s "$ref" ]]; then
        echo "ERROR: missing or empty reference FASTA: $ref" >&2
        exit 1
    fi

    singularity exec --bind "$REAL":"$REAL" "$SIF" \
        makeblastdb \
        -in "$ref" \
        -dbtype nucl \
        -out "$db"

    echo
done

echo "Done building family BLAST databases."
