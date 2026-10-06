#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/bacillus_amr"
INPUT="$PROJECT_DIR/00_metadata/bacillus_refseq_summary.jsonl"
OUTPUT="$PROJECT_DIR/00_metadata/bacillus_master_metadata.tsv"
LOG="$PROJECT_DIR/logs/02_build_master_metadata.log"

{
    echo "BACILLUS AMR PROJECT"
    echo "Step 02: Build master metadata table"
    echo "Date: $(date -Iseconds)"
    echo
    echo "Input: $INPUT"
    echo "Output: $OUTPUT"
    echo

    dataformat tsv genome \
        --inputfile "$INPUT" \
        --fields \
accession,organism-name,organism-tax-id,organism-infraspecific-strain,assminfo-level,assminfo-refseq-category,assminfo-release-date,assminfo-submitter,assminfo-bioproject,assminfo-biosample-accession,assminfo-biosample-strain,assminfo-biosample-host,assminfo-biosample-isolation-source,assminfo-biosample-geo-loc-name,assminfo-biosample-collection-date,assminfo-biosample-description-title,assminfo-biosample-description-comment \
        > "$OUTPUT"

    echo "Rows including header:"
    wc -l "$OUTPUT"

    echo
    echo "Output size:"
    ls -lh "$OUTPUT"

    echo
    echo "SHA256:"
    sha256sum "$OUTPUT"

} 2>&1 | tee "$LOG"
