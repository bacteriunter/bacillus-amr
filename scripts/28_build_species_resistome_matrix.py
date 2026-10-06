#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PREVALENCE = Path(
    "05_amr/amr_prevalence_by_species.tsv"
)

GENOME_SUMMARY = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT = Path(
    "09_analysis/species_resistome_prevalence_matrix.tsv"
)

LONG_OUT = Path(
    "09_analysis/species_resistome_prevalence_long.tsv"
)

LOG = Path(
    "logs/28_species_resistome_matrix.txt"
)

# ============================================================
# 1. Determine eligible species: species-level, n >= 10
# ============================================================

genomes = pd.read_csv(
    GENOME_SUMMARY,
    sep="\t"
)

species_only = genomes[
    genomes["Taxonomic_resolution"] == "species"
].copy()

counts = (
    species_only
    .groupby("Species_normalized")
    .size()
)

eligible = sorted(
    counts[counts >= 10].index
)

# ============================================================
# 2. Load prevalence table
# ============================================================

prev = pd.read_csv(
    PREVALENCE,
    sep="\t"
)

prev = prev[
    prev["Species_normalized"].isin(eligible)
].copy()

# ============================================================
# 3. Build complete species × determinant grid
# ============================================================

all_determinants = sorted(
    prev["Element_symbol"].unique()
)

grid = pd.MultiIndex.from_product(
    [eligible, all_determinants],
    names=[
        "Species_normalized",
        "Element_symbol"
    ]
).to_frame(index=False)

long = grid.merge(
    prev[
        [
            "Species_normalized",
            "Element_symbol",
            "Presence",
            "Species_genome_count",
            "Prevalence_percent"
        ]
    ],
    on=[
        "Species_normalized",
        "Element_symbol"
    ],
    how="left"
)

# Missing combination = determinant not detected
long["Presence"] = (
    long["Presence"]
    .fillna(0)
    .astype(int)
)

long["Prevalence_percent"] = (
    long["Prevalence_percent"]
    .fillna(0.0)
)

long["Species_genome_count"] = (
    long["Species_normalized"]
    .map(counts)
    .astype(int)
)

long.to_csv(
    LONG_OUT,
    sep="\t",
    index=False
)

# ============================================================
# 4. Wide prevalence matrix
# ============================================================

matrix = (
    long.pivot(
        index="Species_normalized",
        columns="Element_symbol",
        values="Prevalence_percent"
    )
)

matrix.to_csv(
    OUT,
    sep="\t"
)

# ============================================================
# 5. Summary
# ============================================================

detected_per_species = (
    (matrix > 0)
    .sum(axis=1)
    .sort_values(ascending=False)
)

detected_species_per_gene = (
    (matrix > 0)
    .sum(axis=0)
    .sort_values(ascending=False)
)

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 28 - SPECIES RESISTOME "
        "PREVALENCE MATRIX\n\n"
    )

    f.write(
        f"Eligible species: {len(matrix)}\n"
    )

    f.write(
        f"Determinants represented: "
        f"{matrix.shape[1]}\n"
    )

    f.write(
        f"Matrix dimensions: "
        f"{matrix.shape[0]} x {matrix.shape[1]}\n\n"
    )

    f.write(
        "DETERMINANTS DETECTED PER SPECIES\n"
    )

    for species, n in detected_per_species.items():
        f.write(f"{species}\t{n}\n")

    f.write(
        "\nNUMBER OF SPECIES CONTAINING EACH DETERMINANT\n"
    )

    for gene, n in detected_species_per_gene.items():
        f.write(f"{gene}\t{n}\n")

print(
    f"Matrix dimensions: "
    f"{matrix.shape[0]} species x "
    f"{matrix.shape[1]} determinants"
)

print("\nDeterminants detected per species:")
print(detected_per_species.to_string())

print("\nNumber of species containing each determinant:")
print(detected_species_per_gene.to_string())

print(f"\nMatrix: {OUT}")
print(f"Long:   {LONG_OUT}")
print(f"Log:    {LOG}")
