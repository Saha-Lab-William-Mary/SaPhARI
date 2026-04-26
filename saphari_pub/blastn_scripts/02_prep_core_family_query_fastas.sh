#!/bin/bash
set -euo pipefail

PROJECT="/sciclone/scr10/ncnallapaneni/saphari_pub"
HIT_DIR="$PROJECT/output/core_region_fastas"
OUT_DIR="$PROJECT/output/blast_core_queries"

mkdir -p "$OUT_DIR"

families=(CFPICI PICI EPIP MFMR P4 PHANIE)

for fam in "${families[@]}"; do
    out="$OUT_DIR/${fam}_core_queries.fasta"
    rm -f "$out"

    if [[ -d "$HIT_DIR/$fam" ]]; then
        find "$HIT_DIR/$fam" -maxdepth 1 -type f -name "*.fasta" -print0 \
            | sort -z \
            | xargs -0 cat > "$out" || true
    fi

    if [[ -s "$out" ]]; then
        echo "$fam -> $(grep -c '^>' "$out") core query sequences"
    else
        echo "$fam -> no core query fasta written"
    fi
done
