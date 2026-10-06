#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$HOME/bacillus_amr"
OUT="$PROJECT_DIR/00_metadata/software_versions.txt"

{
    echo "BACILLUS AMR PROJECT"
    echo "Generated: $(date -Iseconds)"
    echo
    echo "Conda environment: ${CONDA_DEFAULT_ENV:-NA}"
    echo
    echo "=== NCBI Datasets ==="
    datasets --version
    echo
    echo "=== Dataformat ==="
    echo "Executable: $(which dataformat)"
    echo "Version output: $(dataformat version 2>&1)"
    echo
    echo "=== SeqKit ==="
    seqkit version
    echo
    echo "=== Python ==="
    python --version
} > "$OUT"

echo "Written: $OUT"
