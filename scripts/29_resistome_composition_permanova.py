#!/usr/bin/env python3

import numpy as np
import pandas as pd

from pathlib import Path
from scipy.spatial.distance import pdist, squareform

from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.stats.ordination import pcoa


PA_FILE = Path(
    "05_amr/amr_presence_absence.tsv"
)

GENOME_FILE = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT_PCOA = Path(
    "09_analysis/resistome_pcoa_coordinates.tsv"
)

OUT_PERMANOVA = Path(
    "09_analysis/resistome_permanova.tsv"
)

LOG = Path(
    "logs/29_resistome_composition_permanova.txt"
)


# ============================================================
# 1. Load data
# ============================================================

pa = pd.read_csv(
    PA_FILE,
    sep="\t"
)

genomes = pd.read_csv(
    GENOME_FILE,
    sep="\t"
)


# ============================================================
# 2. Select species-level genomes and species n >= 10
# ============================================================

species_only = genomes[
    genomes["Taxonomic_resolution"] == "species"
].copy()

species_counts = (
    species_only
    .groupby("Species_normalized")
    .size()
)

eligible_species = (
    species_counts[
        species_counts >= 10
    ]
    .index
)

metadata = species_only[
    species_only["Species_normalized"].isin(
        eligible_species
    )
][
    [
        "Assembly_Accession",
        "Species_normalized"
    ]
].copy()


# ============================================================
# 3. Merge with binary AMR matrix
# ============================================================

data = metadata.merge(
    pa,
    on="Assembly_Accession",
    how="inner",
    validate="one_to_one"
)

if len(data) != 1082:
    raise RuntimeError(
        f"Expected 1082 genomes, found {len(data)}"
    )

amr_cols = [
    col for col in pa.columns
    if col != "Assembly_Accession"
]

X = data[amr_cols].to_numpy(dtype=int)

ids = data[
    "Assembly_Accession"
].tolist()

groups = data[
    "Species_normalized"
].tolist()


# ============================================================
# 4. Jaccard distance
# ============================================================

dist_condensed = pdist(
    X,
    metric="jaccard"
)

dist_square = squareform(
    dist_condensed
)

dm = DistanceMatrix(
    dist_square,
    ids=ids
)


# ============================================================
# 5. PERMANOVA
# ============================================================

permanova_result = permanova(
    dm,
    grouping=groups,
    permutations=999
)

permanova_df = pd.DataFrame(
    {
        "parameter": permanova_result.index,
        "value": permanova_result.values
    }
)

permanova_df.to_csv(
    OUT_PERMANOVA,
    sep="\t",
    index=False
)


# ============================================================
# 6. Pseudo-R2
#
# scikit-bio reports pseudo-F but not R2 directly.
# For one-factor PERMANOVA:
#
# R2 = F * (k - 1) /
#      [F * (k - 1) + (N - k)]
#
# where:
# N = number of genomes
# k = number of groups
# ============================================================

N = len(data)
k = len(eligible_species)

F = float(
    permanova_result["test statistic"]
)

pseudo_r2 = (
    F * (k - 1)
    /
    (
        F * (k - 1)
        + (N - k)
    )
)


# ============================================================
# 7. PCoA
# ============================================================

ordination = pcoa(
    dm,
    method="eigh"
)

coords = ordination.samples.copy()

coords.insert(
    0,
    "Assembly_Accession",
    coords.index
)

coords = coords.reset_index(drop=True)

coords = coords.merge(
    metadata,
    on="Assembly_Accession",
    how="left",
    validate="one_to_one"
)

coords.to_csv(
    OUT_PCOA,
    sep="\t",
    index=False
)


# ============================================================
# 8. Explained variation
# ============================================================

prop = ordination.proportion_explained

axis1 = float(prop.iloc[0]) * 100
axis2 = float(prop.iloc[1]) * 100
axis3 = float(prop.iloc[2]) * 100


# ============================================================
# 9. Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 29 - RESISTOME COMPOSITION "
        "BY SPECIES\n\n"
    )

    f.write(
        f"Genomes analyzed: {N}\n"
    )

    f.write(
        f"Species analyzed: {k}\n"
    )

    f.write(
        f"AMR determinants: {len(amr_cols)}\n"
    )

    f.write(
        "Distance: Jaccard\n"
    )

    f.write(
        "PERMANOVA permutations: 999\n\n"
    )

    f.write("PERMANOVA\n")

    for key, value in permanova_result.items():
        f.write(
            f"{key}: {value}\n"
        )

    f.write(
        f"Pseudo-R2: {pseudo_r2:.6f}\n\n"
    )

    f.write("PCoA\n")
    f.write(
        f"Axis 1 explained: {axis1:.4f}%\n"
    )
    f.write(
        f"Axis 2 explained: {axis2:.4f}%\n"
    )
    f.write(
        f"Axis 3 explained: {axis3:.4f}%\n"
    )


# ============================================================
# 10. Console output
# ============================================================

print("\nPERMANOVA:")
print(permanova_result)

print(
    f"\nPseudo-R2 = {pseudo_r2:.6f}"
)

print("\nPCoA explained variation:")
print(f"Axis 1 = {axis1:.4f}%")
print(f"Axis 2 = {axis2:.4f}%")
print(f"Axis 3 = {axis3:.4f}%")

print(f"\nPERMANOVA: {OUT_PERMANOVA}")
print(f"PCoA:      {OUT_PCOA}")
print(f"Log:       {LOG}")
