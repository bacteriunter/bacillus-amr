#!/usr/bin/env bash

set -euo pipefail

PROJECT="$HOME/bacillus_amr"
ACCESSIONS="$PROJECT/02_qc/qc_passed_accessions.txt"
GENOME_ROOT="$PROJECT/01_genomes/ncbi_dataset/ncbi_dataset/data"
OUTDIR="$PROJECT/05_amr/amrfinder_raw"
LOG="$PROJECT/logs/10_amrfinder.log"

mkdir -p "$OUTDIR"
mkdir -p "$(dirname "$LOG")"

DB_VERSION=$(amrfinder --database_version 2>&1 | grep "^Database version:" | sed 's/^Database version: //')

{
    echo "BACILLUS AMR PROJECT"
    echo "STEP 10 - AMRFinderPlus"
    echo "Started: $(date --iso-8601=seconds)"
    echo "AMRFinderPlus version: $(amrfinder --version 2>&1 | tail -n 1)"
    echo "Database: $DB_VERSION"
    echo
} > "$LOG"

TOTAL=$(wc -l < "$ACCESSIONS")
COUNT=0
FAILED=0

while IFS= read -r ACC
do
    COUNT=$((COUNT + 1))

    GENOME=$(find "$GENOME_ROOT/$ACC" \
        -type f \
        -name "*_genomic.fna" \
        -print -quit)

    OUT="$OUTDIR/${ACC}.tsv"

    if [[ -z "${GENOME:-}" ]]; then
        echo "[$COUNT/$TOTAL] $ACC - GENOME NOT FOUND" | tee -a "$LOG"
        FAILED=$((FAILED + 1))
        continue
    fi

    if [[ -s "$OUT" ]]; then
        echo "[$COUNT/$TOTAL] $ACC - already completed" | tee -a "$LOG"
        continue
    fi

    echo "[$COUNT/$TOTAL] $ACC" | tee -a "$LOG"

    if amrfinder \
        -n "$GENOME" \
        --plus \
        -o "$OUT" \
        >> "$LOG" 2>&1
    then
        echo "    OK" >> "$LOG"
    else
        echo "    FAILED" | tee -a "$LOG"
        rm -f "$OUT"
        FAILED=$((FAILED + 1))
    fi

done < "$ACCESSIONS"

{
    echo
    echo "Finished: $(date --iso-8601=seconds)"
    echo "Total accessions: $TOTAL"
    echo "Failed: $FAILED"
    echo "Output files: $(find "$OUTDIR" -type f -name '*.tsv' | wc -l)"
} >> "$LOG"

echo
echo "AMRFinderPlus finished."
echo "Total: $TOTAL"
echo "Failed: $FAILED"
echo "Results: $OUTDIR"
echo "Log: $LOG"
