#!/usr/bin/env python3

import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp
from pathlib import Path

INPUT = Path("05_amr/amr_genome_summary.tsv")

OUT_SUMMARY = Path(
    "09_analysis/amr_burden_by_species_summary.tsv"
)

OUT_DUNN = Path(
    "09_analysis/amr_burden_dunn_bh.tsv"
)

LOG = Path(
    "logs/27_amr_burden_by_species.txt"
)

# ============================================================
# 1. Load data
# ============================================================

df = pd.read_csv(INPUT, sep="\t")

# Species-level assignments only
species_df = df[
    df["Taxonomic_resolution"] == "species"
].copy()

# ============================================================
# 2. Species with n >= 10
# ============================================================

species_counts = (
    species_df
    .groupby("Species_normalized")
    .size()
)

eligible_species = (
    species_counts[species_counts >= 10]
    .index
)

x = species_df[
    species_df["Species_normalized"].isin(
        eligible_species
    )
].copy()

# ============================================================
# 3. Descriptive statistics
# ============================================================

summary = (
    x.groupby("Species_normalized")[
        "AMR_determinant_count"
    ]
    .agg(
        n="count",
        mean="mean",
        sd="std",
        median="median",
        min="min",
        max="max"
    )
)

quartiles = (
    x.groupby("Species_normalized")[
        "AMR_determinant_count"
    ]
    .quantile([0.25, 0.75])
    .unstack()
)

quartiles.columns = ["Q1", "Q3"]

summary = summary.join(quartiles)

summary = summary[
    [
        "n",
        "mean",
        "sd",
        "median",
        "Q1",
        "Q3",
        "min",
        "max"
    ]
]

summary = summary.sort_values(
    "mean",
    ascending=False
)

summary.to_csv(
    OUT_SUMMARY,
    sep="\t"
)

# ============================================================
# 4. Kruskal-Wallis
# ============================================================

groups = [
    group["AMR_determinant_count"].values
    for _, group in
    x.groupby("Species_normalized")
]

H, p_kw = kruskal(*groups)

# ============================================================
# 5. Dunn post hoc with BH correction
# ============================================================

dunn = sp.posthoc_dunn(
    x,
    val_col="AMR_determinant_count",
    group_col="Species_normalized",
    p_adjust="fdr_bh"
)

dunn.to_csv(
    OUT_DUNN,
    sep="\t"
)

# ============================================================
# 6. Count significant pairwise comparisons
# ============================================================

species = list(dunn.columns)

sig_pairs = []

for i in range(len(species)):
    for j in range(i + 1, len(species)):

        p = dunn.iloc[i, j]

        if p < 0.05:
            sig_pairs.append(
                (
                    species[i],
                    species[j],
                    p
                )
            )

# ============================================================
# 7. Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 27 - AMR DETERMINANT BURDEN "
        "BY SPECIES\n\n"
    )

    f.write(
        f"Eligible species (n >= 10): "
        f"{len(eligible_species)}\n"
    )

    f.write(
        f"Genomes analyzed: {len(x)}\n\n"
    )

    f.write("KRUSKAL-WALLIS\n")
    f.write(f"H = {H:.6f}\n")
    f.write(f"p = {p_kw:.6e}\n\n")

    f.write(
        "DUNN POST HOC WITH "
        "BENJAMINI-HOCHBERG CORRECTION\n"
    )

    f.write(
        f"Significant pairwise comparisons "
        f"(q < 0.05): {len(sig_pairs)}\n"
    )

    for sp1, sp2, p in sig_pairs:
        f.write(
            f"{sp1}\t{sp2}\t{p:.6e}\n"
        )

# ============================================================
# 8. Console output
# ============================================================

print("\nAMR burden by species:\n")
print(summary.round(3).to_string())

print("\nKruskal-Wallis:")
print(f"H = {H:.6f}")
print(f"p = {p_kw:.6e}")

print(
    "\nSignificant Dunn comparisons "
    f"(BH-adjusted q < 0.05): {len(sig_pairs)}"
)

print(f"\nSummary: {OUT_SUMMARY}")
print(f"Dunn:    {OUT_DUNN}")
print(f"Log:     {LOG}")
