#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/bacillus_amr"

INPUT_DIR="$PROJECT/08_genomic_context/windows_10kb"
OUTPUT_DIR="$PROJECT/08_genomic_context/prodigal_meta"
LOG="$PROJECT/logs/18_prodigal_meta_contexts.log"

mkdir -p "$OUTPUT_DIR"

{
    echo "BACILLUS AMR PROJECT"
    echo "STEP 18 - PRODIGAL META PREDICTION OF AMR CONTEXT WINDOWS"
    echo
    echo "Date: $(date --iso-8601=seconds)"
    echo "Prodigal: $(prodigal -v 2>&1 | head -1)"
    echo "Input windows: $(find "$INPUT_DIR" -maxdepth 1 -name '*.fna' | wc -l)"
    echo
} > "$LOG"

count=0

for fasta in "$INPUT_DIR"/*.fna; do

    count=$((count + 1))

    base=$(basename "$fasta" .fna)

    gff="$OUTPUT_DIR/${base}.gff"
    faa="$OUTPUT_DIR/${base}.faa"
    ffn="$OUTPUT_DIR/${base}.ffn"

    echo "[$count/41] $base" | tee -a "$LOG"

    prodigal \
        -i "$fasta" \
        -a "$faa" \
        -d "$ffn" \
        -f gff \
        -o "$gff" \
        -p meta \
        -g 11 \
        -q

done

echo >> "$LOG"

echo "GFF files: $(find "$OUTPUT_DIR" -maxdepth 1 -name '*.gff' | wc -l)" | tee -a "$LOG"
echo "FAA files: $(find "$OUTPUT_DIR" -maxdepth 1 -name '*.faa' | wc -l)" | tee -a "$LOG"
echo "FFN files: $(find "$OUTPUT_DIR" -maxdepth 1 -name '*.ffn' | wc -l)" | tee -a "$LOG"

TOTAL_CDS=$(grep -h -v '^#' "$OUTPUT_DIR"/*.gff | awk '$3=="CDS"{n++} END{print n+0}')

echo "Total CDS predicted: $TOTAL_CDS" | tee -a "$LOG"
echo "STEP 18 COMPLETE" | tee -a "$LOG"
