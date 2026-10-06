#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# Files
# ============================================================

BURDEN_FILE = Path(
    "09_analysis/amr_burden_by_species_summary.tsv"
)

HETEROGENEITY_FILE = Path(
    "09_analysis/within_species_resistome_heterogeneity.tsv"
)

OUT = Path(
    "11_tables/Table2_species_AMR_burden_heterogeneity.tsv"
)

LOG = Path(
    "logs/37_table2_species_AMR_burden_heterogeneity.txt"
)


# ============================================================
# Load
# ============================================================

burden = pd.read_csv(
    BURDEN_FILE,
    sep="\t"
)

hetero = pd.read_csv(
    HETEROGENEITY_FILE,
    sep="\t"
)


# ============================================================
# Inspect/validate required columns
# ============================================================

required_hetero = {
    "species",
    "n_genomes",
    "mean_within_species_jaccard",
    "sd_within_species_jaccard",
    "median_within_species_jaccard",
    "Q1",
    "Q3",
    "min",
    "max"
}

missing = required_hetero - set(
    hetero.columns
)

if missing:
    raise RuntimeError(
        f"Missing heterogeneity columns: {missing}"
    )

if len(hetero) != 16:
    raise RuntimeError(
        f"Expected 16 species in heterogeneity table, "
        f"found {len(hetero)}"
    )


# ============================================================
# Detect burden column names
# ============================================================

print("Burden columns:")
print(burden.columns.tolist())


def find_column(candidates):

    for c in candidates:
        if c in burden.columns:
            return c

    raise RuntimeError(
        "Could not identify burden column among: "
        + ", ".join(candidates)
    )


species_col = find_column([
    "Species_normalized",
    "species"
])

n_col = find_column([
    "n_genomes",
    "n",
    "N"
])

mean_col = find_column([
    "mean",
    "mean_AMR_determinants",
    "mean_amr_determinants"
])

sd_col = find_column([
    "sd",
    "std",
    "SD",
    "sd_AMR_determinants"
])

median_col = find_column([
    "median",
    "median_AMR_determinants"
])

q1_col = find_column([
    "Q1",
    "q1"
])

q3_col = find_column([
    "Q3",
    "q3"
])

min_col = find_column([
    "min",
    "minimum"
])

max_col = find_column([
    "max",
    "maximum"
])


# ============================================================
# Standardize burden table
# ============================================================

b = burden[[
    species_col,
    n_col,
    mean_col,
    sd_col,
    median_col,
    q1_col,
    q3_col,
    min_col,
    max_col
]].copy()

b.columns = [
    "species",
    "n_genomes",
    "burden_mean",
    "burden_sd",
    "burden_median",
    "burden_Q1",
    "burden_Q3",
    "burden_min",
    "burden_max"
]


# ============================================================
# Merge
# ============================================================

merged = b.merge(
    hetero,
    on=["species", "n_genomes"],
    how="inner",
    validate="one_to_one"
)

if len(merged) != 16:
    raise RuntimeError(
        f"Expected 16 species after merge, "
        f"found {len(merged)}"
    )


# ============================================================
# Order by AMR burden
# ============================================================

merged = merged.sort_values(
    "burden_mean",
    ascending=False
).reset_index(drop=True)


# ============================================================
# Build manuscript table
# ============================================================

table = pd.DataFrame()

table["Species"] = merged["species"]

table["n"] = merged[
    "n_genomes"
].astype(int)

table[
    "AMR determinants per genome, mean ± SD"
] = (
    merged["burden_mean"].map(
        lambda x: f"{x:.2f}"
    )
    + " ± "
    + merged["burden_sd"].map(
        lambda x: f"{x:.2f}"
    )
)

table[
    "AMR determinants per genome, median (IQR)"
] = (
    merged["burden_median"].map(
        lambda x: f"{x:.2f}"
    )
    + " ("
    + merged["burden_Q1"].map(
        lambda x: f"{x:.2f}"
    )
    + "–"
    + merged["burden_Q3"].map(
        lambda x: f"{x:.2f}"
    )
    + ")"
)

table[
    "AMR determinant range"
] = (
    merged["burden_min"].map(
        lambda x: f"{x:.0f}"
    )
    + "–"
    + merged["burden_max"].map(
        lambda x: f"{x:.0f}"
    )
)

table[
    "Mean within-species Jaccard distance"
] = merged[
    "mean_within_species_jaccard"
].map(
    lambda x: f"{x:.3f}"
)

table[
    "Median within-species Jaccard distance (IQR)"
] = (
    merged[
        "median_within_species_jaccard"
    ].map(
        lambda x: f"{x:.3f}"
    )
    + " ("
    + merged["Q1"].map(
        lambda x: f"{x:.3f}"
    )
    + "–"
    + merged["Q3"].map(
        lambda x: f"{x:.3f}"
    )
    + ")"
)


# ============================================================
# Save
# ============================================================

table.to_csv(
    OUT,
    sep="\t",
    index=False
)


# ============================================================
# Log
# ============================================================

with open(LOG, "w") as f:

    f.write(
        "BACILLUS AMR PROJECT\n"
    )

    f.write(
        "TABLE 2 - SPECIES-LEVEL AMR BURDEN AND "
        "WITHIN-SPECIES RESISTOME HETEROGENEITY\n\n"
    )

    f.write(
        table.to_string(index=False)
    )

    f.write("\n\n")

    f.write(
        "Species included: species-level assignments "
        "represented by >=10 genomes.\n"
    )

    f.write(
        "Within-species heterogeneity: pairwise Jaccard "
        "distance among binary AMR determinant profiles.\n"
    )


print(table.to_string(index=False))

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
