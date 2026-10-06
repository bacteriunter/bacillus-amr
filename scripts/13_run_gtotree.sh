#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/miniconda3/envs/bacillus_phylo/bin:$PATH"

PROJECT="$HOME/bacillus_amr"

GTOTREE="$HOME/miniconda3/envs/bacillus_phylo/bin/GToTree"
HMM="$HOME/miniconda3/envs/bacillus_phylo/share/gtotree/hmm_sets/Firmicutes.hmm"

INPUT="$PROJECT/06_phylogeny/protein_faa_unique_paths.txt"
OUTDIR="$PROJECT/06_phylogeny/gtotree_firmicutes"
LOG="$PROJECT/logs/13_gtotree.log"

{
    echo "BACILLUS AMR PROJECT"
    echo "STEP 13 - GTOTREE FIRMICUTES PHYLOGENOMIC MARKERS"
    echo
    echo "Date: $(date --iso-8601=seconds)"
    echo "GToTree: $("$GTOTREE" -v 2>&1)"
    echo "HMM set: $HMM"
    echo "HMM profiles: $(grep -c '^NAME' "$HMM")"
    echo "Input proteomes: $(wc -l < "$INPUT")"
    echo
} > "$LOG"

"$GTOTREE" \
    -A "$INPUT" \
    -H "$HMM" \
    -o "$OUTDIR" \
    -N \
    2>&1 | tee -a "$LOG"
