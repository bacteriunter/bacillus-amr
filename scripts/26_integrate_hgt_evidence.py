#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

CONTEXT = Path(
    "08_genomic_context/sporadic_amr_locus_evidence.tsv"
)

PHYLO = Path(
    "09_analysis/sporadic_amr_phylogenetic_discontinuity_summary.tsv"
)

OUT = Path(
    "09_analysis/sporadic_amr_integrated_evidence.tsv"
)

LOG = Path(
    "logs/26_integrate_hgt_evidence.txt"
)

# ============================================================
# 1. Load data
# ============================================================

context = pd.read_csv(CONTEXT, sep="\t")
phylo = pd.read_csv(PHYLO, sep="\t")

# ============================================================
# 2. Summarize genomic-context evidence
#    by species-determinant combination
# ============================================================

context_summary = (
    context
    .groupby(
        ["species", "determinant"],
        as_index=False
    )
    .agg(
        physical_loci=("locus_number", "count"),
        complete_contexts=(
            "context_status",
            lambda x: (x == "complete").sum()
        ),
        truncated_contexts=(
            "context_status",
            lambda x: (x != "complete").sum()
        ),
        loci_with_mobility_evidence=(
            "mobility_evidence_detected",
            "sum"
        ),
        total_IS_transposases=(
            "n_IS_transposases",
            "sum"
        ),
        total_recombinase_resolvase=(
            "n_recombinase_resolvase",
            "sum"
        ),
        minimum_mobility_distance_bp=(
            "nearest_mobility_distance_bp",
            "min"
        ),
        explicitly_plasmid_loci=(
            "ncbi_explicit_plasmid",
            "sum"
        )
    )
)

# Convert count-like columns to integers
for col in [
    "physical_loci",
    "complete_contexts",
    "truncated_contexts",
    "loci_with_mobility_evidence",
    "total_IS_transposases",
    "total_recombinase_resolvase",
    "explicitly_plasmid_loci"
]:
    context_summary[col] = (
        context_summary[col]
        .fillna(0)
        .astype(int)
    )

# ============================================================
# 3. Merge with phylogenetic evidence
# ============================================================

final = phylo.merge(
    context_summary,
    on=["species", "determinant"],
    how="left",
    validate="one_to_one"
)

if len(final) != 21:
    raise RuntimeError(
        f"Expected 21 species-determinant cases, "
        f"found {len(final)}"
    )

# ============================================================
# 4. Objective evidence flags
#
# These are descriptive flags, NOT HGT classifications.
# ============================================================

final["sporadic_distribution"] = True

final["local_mobility_detected"] = (
    final["loci_with_mobility_evidence"] > 0
)

final["IS_detected"] = (
    final["total_IS_transposases"] > 0
)

final["recombinase_resolvase_detected"] = (
    final["total_recombinase_resolvase"] > 0
)

final["explicit_plasmid_annotation"] = (
    final["explicitly_plasmid_loci"] > 0
)

final["same_species_negative_comparator"] = (
    final["negative_n_tree"] > 0
)

# ============================================================
# 5. Number of independent evidence dimensions
#
# Distributional evidence is common to all 21 cases because
# these cases were selected using the >0-5% prevalence rule.
#
# Mobility and plasmid evidence are counted separately.
# Phylogenetic distance itself remains continuous and is NOT
# thresholded here.
# ============================================================

final["n_context_mobility_dimensions"] = (
    final["IS_detected"].astype(int)
    + final["recombinase_resolvase_detected"].astype(int)
    + final["explicit_plasmid_annotation"].astype(int)
)

# ============================================================
# 6. Sort
# ============================================================

final = final.sort_values(
    [
        "local_mobility_detected",
        "minimum_mobility_distance_bp",
        "prevalence_pct_global"
    ],
    ascending=[False, True, True],
    na_position="last"
)

# ============================================================
# 7. Save
# ============================================================

final.to_csv(
    OUT,
    sep="\t",
    index=False
)

# ============================================================
# 8. Log
# ============================================================

n_cases = len(final)

n_mobility = int(
    final["local_mobility_detected"].sum()
)

n_is = int(
    final["IS_detected"].sum()
)

n_recomb = int(
    final["recombinase_resolvase_detected"].sum()
)

n_plasmid = int(
    final["explicit_plasmid_annotation"].sum()
)

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 26 - INTEGRATED EVIDENCE FOR "
        "SPORADIC AMR DETERMINANTS\n\n"
    )

    f.write(
        f"Species-determinant cases: {n_cases}\n"
    )

    f.write(
        f"Cases with local mobility evidence: "
        f"{n_mobility}\n"
    )

    f.write(
        f"Cases with IS-associated proteins: "
        f"{n_is}\n"
    )

    f.write(
        f"Cases with recombinase/resolvase evidence: "
        f"{n_recomb}\n"
    )

    f.write(
        f"Cases with explicit NCBI plasmid annotation: "
        f"{n_plasmid}\n"
    )

print(
    final[
        [
            "species",
            "determinant",
            "prevalence_pct_global",
            "carrier_n_global",
            "min_distance_to_negative",
            "median_distance_to_negative",
            "physical_loci",
            "complete_contexts",
            "truncated_contexts",
            "loci_with_mobility_evidence",
            "total_IS_transposases",
            "total_recombinase_resolvase",
            "minimum_mobility_distance_bp",
            "explicitly_plasmid_loci"
        ]
    ].to_string(index=False)
)

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
