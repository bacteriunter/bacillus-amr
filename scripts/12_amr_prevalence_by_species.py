#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

AMR = PROJECT / "05_amr" / "amrfinder_amr_hits.tsv"
TAX = PROJECT / "03_taxonomy" / "qc_passed_taxonomy.tsv"

OUT_MATRIX = PROJECT / "05_amr" / "amr_presence_absence.tsv"
OUT_PREV = PROJECT / "05_amr" / "amr_prevalence_by_species.tsv"
OUT_GENOME = PROJECT / "05_amr" / "amr_genome_summary.tsv"
LOG = PROJECT / "logs" / "12_amr_prevalence_summary.txt"

amr = pd.read_csv(AMR, sep="\t", dtype=str)
tax = pd.read_csv(TAX, sep="\t", dtype=str)

# Keep one occurrence per genome/determinant for presence-absence.
# Multiple copies of the same determinant in one genome count as presence = 1.
pairs = (
    amr[
        ["Assembly_Accession", "Element symbol"]
    ]
    .dropna()
    .drop_duplicates()
)

# Genome x AMR determinant matrix
matrix = (
    pairs.assign(Presence=1)
    .pivot(
        index="Assembly_Accession",
        columns="Element symbol",
        values="Presence"
    )
    .fillna(0)
    .astype(int)
)

# Add genomes without AMR hits
all_accessions = tax["Assembly Accession"].drop_duplicates()

matrix = (
    matrix
    .reindex(all_accessions, fill_value=0)
)

matrix.index.name = "Assembly_Accession"
matrix.reset_index().to_csv(
    OUT_MATRIX,
    sep="\t",
    index=False
)

# Number of distinct AMR determinants per genome
genome_summary = pd.DataFrame({
    "Assembly_Accession": matrix.index,
    "AMR_determinant_count": matrix.sum(axis=1).values
})

genome_summary = genome_summary.merge(
    tax[
        [
            "Assembly Accession",
            "Species_normalized",
            "Taxonomic_resolution"
        ]
    ],
    left_on="Assembly_Accession",
    right_on="Assembly Accession",
    how="left"
)

genome_summary.drop(
    columns=["Assembly Accession"],
    inplace=True
)

genome_summary.to_csv(
    OUT_GENOME,
    sep="\t",
    index=False
)

# Long presence/absence table with taxonomy
long = (
    matrix
    .reset_index()
    .melt(
        id_vars="Assembly_Accession",
        var_name="Element_symbol",
        value_name="Presence"
    )
)

long = long.merge(
    tax[
        [
            "Assembly Accession",
            "Species_normalized"
        ]
    ],
    left_on="Assembly_Accession",
    right_on="Assembly Accession",
    how="left"
)

long.drop(
    columns=["Assembly Accession"],
    inplace=True
)

# Number of genomes in each normalized species
species_n = (
    tax["Species_normalized"]
    .value_counts()
    .rename_axis("Species_normalized")
    .reset_index(name="Species_genome_count")
)

# Presence counts by species and determinant
prev = (
    long.groupby(
        ["Species_normalized", "Element_symbol"],
        as_index=False
    )["Presence"]
    .sum()
)

prev = prev.merge(
    species_n,
    on="Species_normalized",
    how="left"
)

prev["Prevalence_percent"] = (
    prev["Presence"]
    / prev["Species_genome_count"]
    * 100
)

# Retain only determinant/species combinations with >=1 presence
prev = prev[
    prev["Presence"] > 0
].copy()

prev = prev.sort_values(
    [
        "Species_normalized",
        "Prevalence_percent",
        "Presence"
    ],
    ascending=[True, False, False]
)

prev.to_csv(
    OUT_PREV,
    sep="\t",
    index=False
)

n_genomes = len(matrix)
n_elements = matrix.shape[1]

with_amr = (
    genome_summary["AMR_determinant_count"] > 0
).sum()

without_amr = (
    genome_summary["AMR_determinant_count"] == 0
).sum()

with open(LOG, "w") as f:
    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 12 - AMR PRESENCE/ABSENCE AND SPECIES PREVALENCE\n\n")

    f.write(f"Genomes: {n_genomes}\n")
    f.write(f"Unique AMR determinants: {n_elements}\n")
    f.write(f"Genomes with >=1 AMR determinant: {with_amr}\n")
    f.write(f"Genomes without AMR determinants: {without_amr}\n")

print("Genomes:", n_genomes)
print("Unique AMR determinants:", n_elements)
print("Genomes with >=1 AMR determinant:", with_amr)
print("Genomes without AMR determinants:", without_amr)

print("\nAMR DETERMINANTS PER GENOME:")
print(
    genome_summary["AMR_determinant_count"]
    .describe()
    .to_string()
)

print("\nMOST PREVALENT DETERMINANTS OVERALL:")
overall = (
    matrix.sum()
    .sort_values(ascending=False)
    .head(20)
)

for element, n in overall.items():
    print(
        f"{element}\t{int(n)}\t"
        f"{n/n_genomes*100:.2f}%"
    )
