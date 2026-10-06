#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

PROJECT = Path.home() / "bacillus_amr"

LOCI_FILE = (
    PROJECT / "08_genomic_context" / "sporadic_amr_loci.tsv"
)

GENOME_ROOT = (
    PROJECT / "01_genomes/ncbi_dataset/ncbi_dataset/data"
)

OUT_DIR = PROJECT / "08_genomic_context" / "windows_10kb"
OUT_TABLE = (
    PROJECT / "08_genomic_context" /
    "sporadic_amr_context_windows.tsv"
)
LOG_FILE = (
    PROJECT / "logs" /
    "16_extract_amr_context_windows.txt"
)

FLANK = 10_000

OUT_DIR.mkdir(parents=True, exist_ok=True)

loci = pd.read_csv(LOCI_FILE, sep="\t")

records_out = []
metadata = []


# ============================================================
# PROCESS EACH PHYSICAL AMR LOCUS
# ============================================================

for i, row in loci.iterrows():

    accession = row["accession"]
    determinant = row["determinant"]
    contig_id = row["Contig id"]

    amr_start = int(row["Start"])
    amr_stop = int(row["Stop"])

    gene_left = min(amr_start, amr_stop)
    gene_right = max(amr_start, amr_stop)

    # --------------------------------------------------------
    # Locate genomic FASTA
    # --------------------------------------------------------

    fasta_files = list(
        (GENOME_ROOT / accession).glob("*_genomic.fna")
    )

    if len(fasta_files) != 1:
        raise RuntimeError(
            f"{accession}: expected 1 genomic FASTA, "
            f"found {len(fasta_files)}"
        )

    fasta = fasta_files[0]

    # --------------------------------------------------------
    # Find target contig
    # --------------------------------------------------------

    target = None

    for rec in SeqIO.parse(fasta, "fasta"):
        if rec.id == contig_id:
            target = rec
            break

    if target is None:
        raise RuntimeError(
            f"{accession}: contig {contig_id} not found"
        )

    contig_length = len(target.seq)

    # --------------------------------------------------------
    # Determine available flanking sequence
    #
    # Coordinates from AMRFinderPlus are 1-based.
    # --------------------------------------------------------

    left_available = gene_left - 1
    right_available = contig_length - gene_right

    left_recovered = min(FLANK, left_available)
    right_recovered = min(FLANK, right_available)

    window_start = gene_left - left_recovered
    window_stop = gene_right + right_recovered

    # Python slicing:
    # convert 1-based inclusive coordinates to 0-based slicing
    window_seq = target.seq[
        window_start - 1 : window_stop
    ]

    left_complete = left_available >= FLANK
    right_complete = right_available >= FLANK

    if left_complete and right_complete:
        context_status = "complete"
    elif not left_complete and not right_complete:
        context_status = "truncated_both"
    elif not left_complete:
        context_status = "truncated_left"
    else:
        context_status = "truncated_right"

    # --------------------------------------------------------
    # Unique locus identifier
    # Important because one genome contains two
    # aph(3')-IIIa physical loci.
    # --------------------------------------------------------

    locus_id = (
        f"L{i+1:03d}_"
        f"{accession}_"
        f"{determinant}_"
        f"{contig_id}_"
        f"{gene_left}-{gene_right}"
    )

    # Make filename-safe version
    safe_determinant = (
        determinant
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("'", "")
    )

    file_id = (
        f"L{i+1:03d}_"
        f"{accession}_"
        f"{safe_determinant}"
    )

    # --------------------------------------------------------
    # FASTA record
    # --------------------------------------------------------

    out_record = SeqRecord(
        window_seq,
        id=file_id,
        description=(
            f"locus={locus_id} "
            f"contig={contig_id} "
            f"window={window_start}-{window_stop} "
            f"AMR={gene_left}-{gene_right} "
            f"strand={row['Strand']} "
            f"context={context_status}"
        )
    )

    records_out.append(out_record)

    # Individual FASTA
    SeqIO.write(
        out_record,
        OUT_DIR / f"{file_id}.fna",
        "fasta"
    )

    metadata.append({
        "locus_id": locus_id,
        "accession": accession,
        "species": row["species"],
        "determinant": determinant,
        "contig_id": contig_id,
        "contig_length_bp": contig_length,
        "amr_start": gene_left,
        "amr_stop": gene_right,
        "strand": row["Strand"],
        "target_flank_bp": FLANK,
        "left_available_bp": left_available,
        "right_available_bp": right_available,
        "left_recovered_bp": left_recovered,
        "right_recovered_bp": right_recovered,
        "window_start": window_start,
        "window_stop": window_stop,
        "window_length_bp": len(window_seq),
        "context_status": context_status,
        "window_file": f"{file_id}.fna"
    })


# ============================================================
# WRITE COMBINED FASTA
# ============================================================

combined_fasta = (
    PROJECT / "08_genomic_context" /
    "sporadic_amr_context_10kb.fna"
)

SeqIO.write(
    records_out,
    combined_fasta,
    "fasta"
)


# ============================================================
# WRITE METADATA
# ============================================================

meta = pd.DataFrame(metadata)

meta.to_csv(
    OUT_TABLE,
    sep="\t",
    index=False
)


# ============================================================
# VALIDATION
# ============================================================

if len(meta) != len(loci):
    raise RuntimeError(
        f"Expected {len(loci)} loci, extracted {len(meta)}"
    )

if meta["locus_id"].duplicated().any():
    raise RuntimeError("Duplicate locus IDs detected.")

if (meta["window_length_bp"] <= 0).any():
    raise RuntimeError("Invalid zero-length window detected.")


# ============================================================
# SUMMARY
# ============================================================

status_counts = (
    meta["context_status"]
    .value_counts()
    .to_dict()
)

log = [
    "BACILLUS AMR PROJECT",
    "STEP 16 - EXTRACT ±10 KB AMR CONTEXT WINDOWS",
    "",
    f"Physical AMR loci: {len(meta)}",
    f"Target flank: {FLANK} bp per side",
    f"Complete contexts: "
    f"{status_counts.get('complete', 0)}",
    f"Truncated left: "
    f"{status_counts.get('truncated_left', 0)}",
    f"Truncated right: "
    f"{status_counts.get('truncated_right', 0)}",
    f"Truncated both: "
    f"{status_counts.get('truncated_both', 0)}",
    "",
    f"Minimum window length: "
    f"{meta['window_length_bp'].min()} bp",
    f"Maximum window length: "
    f"{meta['window_length_bp'].max()} bp",
    "",
    f"Combined FASTA: {combined_fasta}",
    f"Metadata table: {OUT_TABLE}",
    f"Individual FASTA directory: {OUT_DIR}"
]

LOG_FILE.write_text(
    "\n".join(log) + "\n"
)

print("\n".join(log))
