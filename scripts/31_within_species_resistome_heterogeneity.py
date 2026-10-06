#!/usr/bin/env python3

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from pathlib import Path


PA_FILE = Path(
    "05_amr/amr_presence_absence.tsv"
)

GENOME_FILE = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT = Path(
    "09_analysis/within_species_resistome_heterogeneity.tsv"
)

LOG = Path(
    "logs/31_within_species_resistome_heterogeneity.txt"
)


# ============================================================
# 1. Load data
# ============================================================

pa = pd.read_csv(PA_FILE, sep="\t")
genomes = pd.read_csv(GENOME_FILE, sep="\t")


# ============================================================
# 2. Species-level genomes with n >= 10
# ============================================================

species_only = genomes[
    genomes["Taxonomic_resolution"] == "species"
].copy()

counts = (
    species_only
    .groupby("Species_normalized")
    .size()
)

eligible = counts[counts >= 10].index

metadata = species_only[
    species_only["Species_normalized"].isin(eligible)
][
    ["Assembly_Accession", "Species_normalized"]
].copy()


# ============================================================
# 3. Merge with binary resistome matrix
# ============================================================

data = metadata.merge(
    pa,
    on="Assembly_Accession",
    how="inner",
    validate="one_to_one"
)

amr_cols = [
    c for c in pa.columns
    if c != "Assembly_Accession"
]


# ============================================================
# 4. Pairwise within-species Jaccard distances
# ============================================================

records = []

for species, group in data.groupby("Species_normalized"):

    X = group[amr_cols].to_numpy(dtype=int)

    distances = pdist(
        X,
        metric="jaccard"
    )

    # pdist may return NaN for comparisons between
    # two all-zero resistome profiles.
    valid = distances[
        np.isfinite(distances)
    ]

    total_pairs = len(distances)
    valid_pairs = len(valid)
    undefined_pairs = total_pairs - valid_pairs

    if valid_pairs > 0:

        record = {
            "species": species,
            "n_genomes": len(group),
            "total_pairwise_comparisons": total_pairs,
            "valid_jaccard_comparisons": valid_pairs,
            "undefined_zero_zero_comparisons": undefined_pairs,
            "mean_within_species_jaccard": np.mean(valid),
            "sd_within_species_jaccard": np.std(valid, ddof=1),
            "median_within_species_jaccard": np.median(valid),
            "Q1": np.quantile(valid, 0.25),
            "Q3": np.quantile(valid, 0.75),
            "min": np.min(valid),
            "max": np.max(valid)
        }

    else:

        record = {
            "species": species,
            "n_genomes": len(group),
            "total_pairwise_comparisons": total_pairs,
            "valid_jaccard_comparisons": 0,
            "undefined_zero_zero_comparisons": undefined_pairs,
            "mean_within_species_jaccard": np.nan,
            "sd_within_species_jaccard": np.nan,
            "median_within_species_jaccard": np.nan,
            "Q1": np.nan,
            "Q3": np.nan,
            "min": np.nan,
            "max": np.nan
        }

    records.append(record)


result = pd.DataFrame(records)

result = result.sort_values(
    "mean_within_species_jaccard",
    ascending=False,
    na_position="last"
)

result.to_csv(
    OUT,
    sep="\t",
    index=False
)


# ============================================================
# 5. Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 31 - WITHIN-SPECIES RESISTOME "
        "HETEROGENEITY\n\n"
    )

    f.write(
        "Metric: pairwise Jaccard distance "
        "among genomes within each species\n"
    )

    f.write(
        "Species included: n >= 10 genomes\n\n"
    )

    f.write(
        result.round(6).to_csv(
            sep="\t",
            index=False
        )
    )


print(
    result[
        [
            "species",
            "n_genomes",
            "mean_within_species_jaccard",
            "median_within_species_jaccard",
            "Q1",
            "Q3",
            "min",
            "max",
            "undefined_zero_zero_comparisons"
        ]
    ].round(4).to_string(index=False)
)

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")

