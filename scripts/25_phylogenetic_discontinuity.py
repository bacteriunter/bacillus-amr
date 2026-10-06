#!/usr/bin/env python3

import pandas as pd
from Bio import Phylo
from pathlib import Path

TREE_FILE = Path(
    "06_phylogeny/iqtree_subset87/final_tree.treefile"
)
SUBSET_FILE = Path(
    "06_phylogeny/phylogenetic_subset.tsv"
)
CARRIERS_FILE = Path(
    "05_amr/sporadic_amr_carriers.tsv"
)
PA_FILE = Path(
    "05_amr/amr_presence_absence.tsv"
)

OUT = Path(
    "09_analysis/sporadic_amr_phylogenetic_discontinuity.tsv"
)
SUMMARY = Path(
    "09_analysis/sporadic_amr_phylogenetic_discontinuity_summary.tsv"
)
LOG = Path(
    "logs/25_phylogenetic_discontinuity.txt"
)

# ============================================================
# 1. Load data
# ============================================================

tree = Phylo.read(TREE_FILE, "newick")

subset = pd.read_csv(
    SUBSET_FILE,
    sep="\t"
)

carriers = pd.read_csv(
    CARRIERS_FILE,
    sep="\t"
)

pa = pd.read_csv(
    PA_FILE,
    sep="\t"
)

pa = pa.set_index("Assembly_Accession")

# ============================================================
# 2. Tree tips
# ============================================================

tips = {
    tip.name
    for tip in tree.get_terminals()
}

if len(tips) != 87:
    raise RuntimeError(
        f"Expected 87 tree tips, found {len(tips)}"
    )

# ============================================================
# 3. Identify accession/species columns in subset
# ============================================================

subset_acc_col = "Assembly_Accession"

if subset_acc_col not in subset.columns:
    raise RuntimeError(
        "Assembly_Accession not found in phylogenetic subset"
    )

if "Species_normalized" not in subset.columns:
    raise RuntimeError(
        "Species_normalized not found in phylogenetic subset"
    )

subset = subset[
    subset[subset_acc_col].isin(tips)
].copy()

species_map = dict(
    zip(
        subset[subset_acc_col],
        subset["Species_normalized"]
    )
)

# ============================================================
# 4. Validate all tips have species assignment
# ============================================================

missing_species = tips - set(species_map)

if missing_species:
    raise RuntimeError(
        "Tree tips without species assignment: "
        + ", ".join(sorted(missing_species))
    )

# ============================================================
# 5. Analyze each sporadic carrier
# ============================================================

records = []

for _, row in carriers.iterrows():

    accession = row["accession"]
    species = row["species"]
    determinant = row["determinant"]

    if accession not in tips:
        raise RuntimeError(
            f"Carrier absent from tree: {accession}"
        )

    if determinant not in pa.columns:
        raise RuntimeError(
            f"Determinant absent from presence/absence "
            f"matrix: {determinant}"
        )

    # All sampled genomes of same species in reduced tree
    same_species = sorted([
        acc
        for acc in tips
        if species_map[acc] == species
    ])

    # Determine actual AMR state from complete matrix
    positive = [
        acc for acc in same_species
        if acc in pa.index
        and int(pa.loc[acc, determinant]) == 1
    ]

    negative = [
        acc for acc in same_species
        if acc in pa.index
        and int(pa.loc[acc, determinant]) == 0
    ]

    # Nearest negative
    nearest_negative = None
    distance_negative = None

    if negative:
        distances = [
            (other, tree.distance(accession, other))
            for other in negative
        ]

        nearest_negative, distance_negative = min(
            distances,
            key=lambda x: x[1]
        )

    # Nearest OTHER positive
    other_positive = [
        acc for acc in positive
        if acc != accession
    ]

    nearest_positive = None
    distance_positive = None

    if other_positive:
        distances = [
            (other, tree.distance(accession, other))
            for other in other_positive
        ]

        nearest_positive, distance_positive = min(
            distances,
            key=lambda x: x[1]
        )

    records.append({
        "accession": accession,
        "species": species,
        "determinant": determinant,
        "species_n_global": int(row["species_n"]),
        "carrier_n_global": int(row["carrier_n"]),
        "prevalence_pct_global": row["prevalence_pct"],
        "species_n_tree": len(same_species),
        "positive_n_tree": len(positive),
        "negative_n_tree": len(negative),
        "nearest_negative": nearest_negative,
        "distance_to_nearest_negative": distance_negative,
        "nearest_positive": nearest_positive,
        "distance_to_nearest_positive": distance_positive
    })

result = pd.DataFrame(records)

result.to_csv(
    OUT,
    sep="\t",
    index=False
)

# ============================================================
# 6. Species-determinant summary
# ============================================================

summary = (
    result.groupby(
        ["species", "determinant"],
        as_index=False
    )
    .agg(
        species_n_global=("species_n_global", "first"),
        carrier_n_global=("carrier_n_global", "first"),
        prevalence_pct_global=("prevalence_pct_global", "first"),
        species_n_tree=("species_n_tree", "first"),
        positive_n_tree=("positive_n_tree", "first"),
        negative_n_tree=("negative_n_tree", "first"),
        carriers_analyzed=("accession", "nunique"),
        min_distance_to_negative=(
            "distance_to_nearest_negative", "min"
        ),
        median_distance_to_negative=(
            "distance_to_nearest_negative", "median"
        ),
        max_distance_to_negative=(
            "distance_to_nearest_negative", "max"
        )
    )
)

summary.to_csv(
    SUMMARY,
    sep="\t",
    index=False
)

# ============================================================
# 7. Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 25 - PHYLOGENETIC DISCONTINUITY OF "
        "SPORADIC AMR DETERMINANTS\n\n"
    )

    f.write(f"Tree tips: {len(tips)}\n")
    f.write(
        f"Carrier-determinant observations: "
        f"{len(result)}\n"
    )
    f.write(
        f"Species-determinant cases: "
        f"{len(summary)}\n"
    )

    f.write(
        "Carriers with same-species negative comparator: "
        f"{result['distance_to_nearest_negative'].notna().sum()}\n"
    )

    f.write(
        "Carriers with another same-species positive comparator: "
        f"{result['distance_to_nearest_positive'].notna().sum()}\n"
    )


print("\nSpecies-determinant summary:\n")

print(
    summary[
        [
            "species",
            "determinant",
            "prevalence_pct_global",
            "species_n_tree",
            "positive_n_tree",
            "negative_n_tree",
            "min_distance_to_negative",
            "median_distance_to_negative",
            "max_distance_to_negative"
        ]
    ].to_string(index=False)
)

print(f"\nOutput:  {OUT}")
print(f"Summary: {SUMMARY}")
print(f"Log:     {LOG}")
