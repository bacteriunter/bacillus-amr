#!/usr/bin/env bash
set -euo pipefail

PROJECT="$HOME/bacillus_amr"

INPUT_DIR="$PROJECT/08_genomic_context/windows_10kb"
OUTPUT_DIR="$PROJECT/08_genomic_context/prokka"
LOG="$PROJECT/logs/17_prokka_contexts.log"

mkdir -p "$OUTPUT_DIR"

{
    echo "BACILLUS AMR PROJECT"
    echo "STEP 17 - PROKKA ANNOTATION OF AMR CONTEXT WINDOWS"
    echo
    echo "Date: $(date --iso-8601=seconds)"
    echo "Prokka: $(prokka --version 2>&1)"
    echo "Input windows: $(find "$INPUT_DIR" -maxdepth 1 -name '*.fna' | wc -l)"
    echo
} > "$LOG"

count=0

for fasta in "$INPUT_DIR"/*.fna; do

    count=$((count + 1))

    base=$(basename "$fasta" .fna)
    outdir="$OUTPUT_DIR/$base"

    echo "[$count/41] $base" | tee -a "$LOG"

    # Remove incomplete previous result, if present
    if [[ -d "$outdir" && ! -s "$outdir/${base}.gff" ]]; then
        rm -rf "$outdir"
    fi

    # Skip completed annotations
    if [[ -s "$outdir/${base}.gff" ]]; then
        echo "  already completed - skipping" | tee -a "$LOG"
        continue
    fi

    prokka "$fasta" \
        --outdir "$outdir" \
        --prefix "$base" \
        --kingdom Bacteria \
        --genus Bacillus \
        --cpus 2 \
        >> "$LOG" 2>&1

done

echo >> "$LOG"
echo "Completed GFF files: $(find "$OUTPUT_DIR" -name '*.gff' | wc -l)" | tee -a "$LOG"
echo "Completed FAA files: $(find "$OUTPUT_DIR" -name '*.faa' | wc -l)" | tee -a "$LOG"
echo "Completed TSV files: $(find "$OUTPUT_DIR" -name '*.tsv' | wc -l)" | tee -a "$LOG"

echo "STEP 17 COMPLETE" | tee -a "$LOG"
