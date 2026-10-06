#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/bacillus_amr"
GENOME_DIR="${PROJECT_DIR}/01_genomes/ncbi_dataset/ncbi_dataset/data"
OUT="${PROJECT_DIR}/02_qc/seqkit_genome_stats.tsv"
LOG="${PROJECT_DIR}/logs/07_genome_qc_seqkit.log"

cd "$PROJECT_DIR"

mkdir -p 02_qc logs

{
    echo "BACILLUS AMR PROJECT"
    echo "STEP 07 - STRUCTURAL GENOME QC WITH SEQKIT"
    echo "Started: $(date --iso-8601=seconds)"
    echo

    printf "Assembly_Accession\tFile\tNum_contigs\tTotal_length\tMin_length\tAverage_length\tMax_length\tN50\tGC_percent\n" \
        > "$OUT"

    count=0

    while IFS= read -r accession; do

        fasta=$(find "${GENOME_DIR}/${accession}" \
            -maxdepth 1 \
            -type f \
            -name "*_genomic.fna" \
            -print -quit)

        if [[ -z "${fasta}" ]]; then
            echo "WARNING: FASTA not found for ${accession}"
            continue
        fi

        stats=$(seqkit stats \
            -T \
            -a \
            "$fasta" | tail -n 1)

        file=$(printf '%s\n' "$stats" | cut -f1)
        num=$(printf '%s\n' "$stats" | cut -f4)
        sum_len=$(printf '%s\n' "$stats" | cut -f5)
        min_len=$(printf '%s\n' "$stats" | cut -f6)
        avg_len=$(printf '%s\n' "$stats" | cut -f7)
        max_len=$(printf '%s\n' "$stats" | cut -f8)
        n50=$(printf '%s\n' "$stats" | cut -f13)
        gc=$(printf '%s\n' "$stats" | cut -f18)

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$accession" \
            "$file" \
            "$num" \
            "$sum_len" \
            "$min_len" \
            "$avg_len" \
            "$max_len" \
            "$n50" \
            "$gc" >> "$OUT"

        count=$((count + 1))

        if (( count % 100 == 0 )); then
            echo "Processed: ${count}"
        fi

    done < 01_genomes/selected_accessions.txt

    echo
    echo "Genomes processed: ${count}"
    echo "Output rows: $(( $(wc -l < "$OUT") - 1 ))"
    echo "Finished: $(date --iso-8601=seconds)"

} 2>&1 | tee "$LOG"

