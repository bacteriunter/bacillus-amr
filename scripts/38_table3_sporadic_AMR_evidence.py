#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# Files
# ============================================================

INPUT = Path(
    "09_analysis/sporadic_amr_integrated_evidence.tsv"
)

OUT = Path(
    "11_tables/Table3_sporadic_AMR_evidence.tsv"
)

LOG = Path(
    "logs/38_table3_sporadic_AMR_evidence.txt"
)


# ============================================================
# Load
# ============================================================

df = pd.read_csv(
    INPUT,
    sep="\t"
)


# ============================================================
# Validation
# ============================================================

required = {
    "species",
    "determinant",
    "species_n_global",
    "carrier_n_global",
    "prevalence_pct_global",
    "min_distance_to_negative",
    "physical_loci",
    "complete_contexts",
    "truncated_contexts",
    "loci_with_mobility_evidence",
    "total_IS_transposases",
    "total_recombinase_resolvase",
    "minimum_mobility_distance_bp",
    "explicitly_plasmid_loci",
    "sporadic_distribution",
    "local_mobility_detected",
    "same_species_negative_comparator"
}

missing = required - set(df.columns)

if missing:
    raise RuntimeError(
        f"Missing required columns: {missing}"
    )

if len(df) != 21:
    raise RuntimeError(
        f"Expected 21 species-determinant combinations, "
        f"found {len(df)}"
    )

if not df["sporadic_distribution"].all():
    raise RuntimeError(
        "Not all rows are classified as sporadic."
    )

if not df[
    "same_species_negative_comparator"
].all():
    raise RuntimeError(
        "At least one case lacks a same-species "
        "negative phylogenetic comparator."
    )


# ============================================================
# Order rows
#
# Mobility-positive combinations first.
# Within each group: descending prevalence.
# ============================================================

df = df.sort_values(
    by=[
        "local_mobility_detected",
        "prevalence_pct_global",
        "species",
        "determinant"
    ],
    ascending=[
        False,
        False,
        True,
        True
    ]
).reset_index(drop=True)


# ============================================================
# Formatting helpers
# ============================================================

def yes_no(value):
    return "Yes" if bool(value) else "No"


def format_distance(value):

    if pd.isna(value):
        return "NA"

    return f"{value:.6f}"


def format_mobility_distance(value):

    if pd.isna(value):
        return "—"

    return f"{int(round(value)):,}"


# ============================================================
# Build manuscript table
# ============================================================

table = pd.DataFrame()

table["Species"] = df["species"]

table["AMR determinant"] = df[
    "determinant"
]

table[
    "Global prevalence, carriers/n (%)"
] = (
    df["carrier_n_global"].astype(int).astype(str)
    + "/"
    + df["species_n_global"].astype(int).astype(str)
    + " ("
    + df["prevalence_pct_global"].map(
        lambda x: f"{x:.2f}"
    )
    + ")"
)

table[
    "Minimum phylogenetic distance to negative comparator"
] = df[
    "min_distance_to_negative"
].map(
    format_distance
)

table[
    "Physical loci"
] = df[
    "physical_loci"
].astype(int)

table[
    "Complete/truncated contexts"
] = (
    df["complete_contexts"]
    .astype(int)
    .astype(str)
    + "/"
    + df["truncated_contexts"]
    .astype(int)
    .astype(str)
)

table[
    "Loci with local mobility evidence"
] = (
    df["loci_with_mobility_evidence"]
    .astype(int)
    .astype(str)
    + "/"
    + df["physical_loci"]
    .astype(int)
    .astype(str)
)

table[
    "IS-associated proteins"
] = df[
    "total_IS_transposases"
].astype(int)

table[
    "Recombinase/resolvase"
] = df[
    "total_recombinase_resolvase"
].map(
    lambda x: yes_no(x > 0)
)

table[
    "Minimum mobility distance (bp)"
] = df[
    "minimum_mobility_distance_bp"
].map(
    format_mobility_distance
)

table[
    "Explicit NCBI plasmid annotation"
] = df[
    "explicitly_plasmid_loci"
].map(
    lambda x: yes_no(x > 0)
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
# Summary
# ============================================================

n_mobility_cases = int(
    df["local_mobility_detected"].sum()
)

n_no_mobility_cases = (
    len(df) - n_mobility_cases
)

n_is_cases = int(
    (df["total_IS_transposases"] > 0).sum()
)

n_recomb_cases = int(
    (
        df["total_recombinase_resolvase"] > 0
    ).sum()
)

n_plasmid_cases = int(
    (
        df["explicitly_plasmid_loci"] > 0
    ).sum()
)


# ============================================================
# Log
# ============================================================

with open(LOG, "w") as f:

    f.write(
        "BACILLUS AMR PROJECT\n"
    )

    f.write(
        "TABLE 3 - SPORADIC AMR DETERMINANTS AND "
        "EVIDENCE CONSISTENT WITH POTENTIAL "
        "HORIZONTAL ACQUISITION\n\n"
    )

    f.write(
        f"Species-determinant combinations: "
        f"{len(df)}\n"
    )

    f.write(
        f"Combinations with local mobility evidence: "
        f"{n_mobility_cases}\n"
    )

    f.write(
        f"Combinations without detected local mobility: "
        f"{n_no_mobility_cases}\n"
    )

    f.write(
        f"Combinations with IS-associated proteins: "
        f"{n_is_cases}\n"
    )

    f.write(
        f"Combinations with recombinase/resolvase: "
        f"{n_recomb_cases}\n"
    )

    f.write(
        f"Combinations with explicit NCBI plasmid annotation: "
        f"{n_plasmid_cases}\n\n"
    )

    f.write(
        "All combinations have within-species prevalence "
        ">0–5% among species represented by >=10 genomes.\n"
    )

    f.write(
        "Phylogenetic distance is the minimum patristic "
        "distance from a carrier to a determinant-negative "
        "genome of the same species in the targeted "
        "87-genome phylogeny.\n"
    )

    f.write(
        "Complete/truncated contexts refer to recovery of "
        "the target +/-10 kb genomic-context window.\n"
    )

    f.write(
        "Absence of detected local mobility evidence is not "
        "interpreted as evidence of absence, particularly "
        "for truncated genomic contexts.\n"
    )

    f.write(
        "No explicit NCBI plasmid annotation is interpreted "
        "only as absence of explicit plasmid annotation, "
        "not as proof of chromosomal localization.\n\n"
    )

    f.write(
        table.to_string(index=False)
    )

    f.write("\n")


# ============================================================
# Console
# ============================================================

print(table.to_string(index=False))

print(
    f"\nSpecies-determinant combinations: {len(df)}"
)

print(
    "With local mobility evidence:",
    n_mobility_cases
)

print(
    "Without detected local mobility:",
    n_no_mobility_cases
)

print(
    "With IS-associated proteins:",
    n_is_cases
)

print(
    "With recombinase/resolvase:",
    n_recomb_cases
)

print(
    "With explicit NCBI plasmid annotation:",
    n_plasmid_cases
)

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
