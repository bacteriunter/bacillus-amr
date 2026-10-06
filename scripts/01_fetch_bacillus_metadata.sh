#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/bacillus_amr"
OUTDIR="$PROJECT_DIR/00_metadata"
LOGDIR="$PROJECT_DIR/logs"

OUTFILE="$OUTDIR/bacillus_refseq_summary.jsonl"
LOGFILE="$LOGDIR/01_fetch_bacillus_metadata.log"

mkdir -p "$OUTDIR" "$LOGDIR"

{
    echo "BACILLUS AMR PROJECT"
    echo "Step 01: Retrieve Bacillus RefSeq assembly metadata"
    echo "Date: $(date -Iseconds)"
    echo
    echo "Taxon: Bacillus (NCBI Taxonomy ID: 1386)"
    echo "Assembly source: RefSeq"
    echo

    datasets summary genome taxon 1386 \
        --assembly-source RefSeq \
        --as-json-lines \
        > "$OUTFILE"

    echo
    echo "Records retrieved:"
    wc -l "$OUTFILE"

    echo
    echo "Output file:"
    ls -lh "$OUTFILE"

    echo
    echo "SHA256:"
    sha256sum "$OUTFILE"

} 2>&1 | tee "$LOGFILE"
