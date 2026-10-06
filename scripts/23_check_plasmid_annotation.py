#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
import re

INPUT = Path("08_genomic_context/sporadic_amr_loci.tsv")
GENOME_ROOT = Path("01_genomes/ncbi_dataset/ncbi_dataset/data")

OUT = Path("08_genomic_context/amr_locus_plasmid_annotation.tsv")
LOG = Path("logs/23_check_plasmid_annotation.txt")

df = pd.read_csv(INPUT, sep="\t")

records = []

for i, row in df.iterrows():

    accession = str(row["accession"])
    contig = str(row["Contig id"])

    # Reconstruct the same locus numbering used previously: L001-L041
    locus_number = f"L{i + 1:03d}"

    genome_dir = GENOME_ROOT / accession
    fasta_files = list(genome_dir.glob("*_genomic.fna"))

    if len(fasta_files) != 1:
        raise RuntimeError(
            f"{accession}: expected 1 genomic FASTA, "
            f"found {len(fasta_files)}"
        )

    fasta = fasta_files[0]
    header = None

    with open(fasta) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue

            seqid = line[1:].split()[0]

            if seqid == contig:
                header = line[1:].strip()
                break

    if header is None:
        raise RuntimeError(
            f"{accession}: contig {contig} not found"
        )

    # Explicit NCBI annotation only.
    # Absence of "plasmid" is NOT interpreted as proof
    # of chromosomal localization.
    is_plasmid = bool(
        re.search(r"\bplasmid\b", header, flags=re.IGNORECASE)
    )

    records.append({
        "locus_number": locus_number,
        "accession": accession,
        "species": row["species"],
        "determinant": row["determinant"],
        "contig": contig,
        "start": int(row["Start"]),
        "stop": int(row["Stop"]),
        "strand": row["Strand"],
        "ncbi_header": header,
        "ncbi_explicit_plasmid": is_plasmid
    })


out = pd.DataFrame(records)

out.to_csv(
    OUT,
    sep="\t",
    index=False
)

n = len(out)
n_plasmid = int(out["ncbi_explicit_plasmid"].sum())

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 23 - NCBI PLASMID ANNOTATION OF "
        "SPORADIC AMR LOCI\n\n"
    )

    f.write(f"Physical AMR loci: {n}\n")
    f.write(
        f"Explicitly annotated as plasmid by NCBI: "
        f"{n_plasmid}\n"
    )
    f.write(
        f"Without explicit plasmid annotation: "
        f"{n - n_plasmid}\n"
    )


print("\nExplicitly plasmid-associated loci:\n")

plasmids = out[out["ncbi_explicit_plasmid"]]

if len(plasmids):

    print(
        plasmids[
            [
                "locus_number",
                "accession",
                "species",
                "determinant",
                "contig",
                "ncbi_header"
            ]
        ].to_string(index=False)
    )

else:
    print("None")

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
