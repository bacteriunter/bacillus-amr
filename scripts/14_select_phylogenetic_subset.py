#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

AMR_FILE = PROJECT / "05_amr" / "amr_presence_absence.tsv"
TAX_FILE = PROJECT / "03_taxonomy" / "qc_passed_taxonomy.tsv"
CARRIERS_FILE = PROJECT / "05_amr" / "sporadic_amr_carrier_summary.tsv"

OUT_FILE = PROJECT / "06_phylogeny" / "phylogenetic_subset.tsv"
OUT_ACCESSIONS = PROJECT / "06_phylogeny" / "phylogenetic_subset_accessions.txt"
LOG_FILE = PROJECT / "logs" / "14_phylogenetic_subset.txt"

MIN_SPECIES_N = 10
N_CONTEXT = 3


# ============================================================
# LOAD DATA
# ============================================================

amr = pd.read_csv(AMR_FILE, sep="\t")
tax = pd.read_csv(TAX_FILE, sep="\t")
carriers = pd.read_csv(CARRIERS_FILE, sep="\t")

ACC = "Assembly_Accession"
TAX_ACC = "Assembly Accession"
SPECIES = "Species_normalized"

amr_cols = [c for c in amr.columns if c != ACC]


# ============================================================
# MERGE AMR + TAXONOMY
# ============================================================

df = amr.merge(
    tax[[TAX_ACC, SPECIES]],
    left_on=ACC,
    right_on=TAX_ACC,
    how="left"
)

df = df.drop(columns=[TAX_ACC])

if df[SPECIES].isna().any():
    raise ValueError("Some genomes lack normalized taxonomy.")


# ============================================================
# ELIGIBLE SPECIES
# ============================================================

species_counts = (
    df.groupby(SPECIES)[ACC]
    .nunique()
)

eligible_species = species_counts[
    (species_counts >= MIN_SPECIES_N) &
    (species_counts.index != "Bacillus sp.")
].index.tolist()

if len(eligible_species) != 16:
    raise ValueError(
        f"Expected 16 eligible species, found {len(eligible_species)}."
    )


# ============================================================
# AMR PROFILE
# ============================================================

df["AMR_profile"] = (
    df[amr_cols]
    .astype(int)
    .astype(str)
    .agg("".join, axis=1)
)


# ============================================================
# SPORADIC CARRIERS
# ============================================================

carrier_set = set(carriers["accession"])

carrier_rows = df[
    df[ACC].isin(carrier_set)
].copy()

if carrier_rows[ACC].nunique() != 39:
    raise ValueError(
        f"Expected 39 sporadic carriers, found "
        f"{carrier_rows[ACC].nunique()}."
    )

carrier_rows["selection_type"] = "sporadic_AMR_carrier"

carrier_info = carriers[
    ["accession", "n_sporadic_determinants", "sporadic_determinants"]
].copy()

carrier_rows = carrier_rows.merge(
    carrier_info,
    left_on=ACC,
    right_on="accession",
    how="left"
).drop(columns=["accession"])


# ============================================================
# SELECT 3 NON-CARRIERS PER SPECIES
#
# Rule:
# 1. Exclude all sporadic carriers.
# 2. Sort deterministically by accession.
# 3. Select one genome from each distinct AMR profile first.
# 4. If fewer than 3 profiles are available, fill remaining
#    positions with the lowest-accession remaining genomes.
# ============================================================

context_rows = []

for species in sorted(eligible_species):

    sub = df[
        (df[SPECIES] == species) &
        (~df[ACC].isin(carrier_set))
    ].copy()

    sub = sub.sort_values(ACC)

    if len(sub) < N_CONTEXT:
        raise ValueError(
            f"{species} has only {len(sub)} available non-carriers."
        )

    # One representative per distinct AMR profile
    profile_reps = (
        sub.drop_duplicates(
            subset="AMR_profile",
            keep="first"
        )
        .sort_values(ACC)
    )

    selected = profile_reps.head(N_CONTEXT).copy()

    # Complete to N_CONTEXT if species has <3 distinct profiles
    if len(selected) < N_CONTEXT:

        already = set(selected[ACC])

        remaining = sub[
            ~sub[ACC].isin(already)
        ].sort_values(ACC)

        needed = N_CONTEXT - len(selected)

        selected = pd.concat(
            [selected, remaining.head(needed)],
            ignore_index=True
        )

    if len(selected) != N_CONTEXT:
        raise ValueError(
            f"Could not select {N_CONTEXT} context genomes for {species}."
        )

    selected["selection_type"] = "intraspecific_context"
    selected["n_sporadic_determinants"] = 0
    selected["sporadic_determinants"] = ""

    context_rows.append(selected)


context_rows = pd.concat(
    context_rows,
    ignore_index=True
)


# ============================================================
# COMBINE
# ============================================================

keep_cols = [
    ACC,
    SPECIES,
    "selection_type",
    "AMR_profile",
    "n_sporadic_determinants",
    "sporadic_determinants"
]

subset = pd.concat(
    [
        carrier_rows[keep_cols],
        context_rows[keep_cols]
    ],
    ignore_index=True
)

subset = subset.sort_values(
    [SPECIES, "selection_type", ACC]
).reset_index(drop=True)


# ============================================================
# VALIDATION
# ============================================================

n_total = subset[ACC].nunique()
n_carriers = (
    subset["selection_type"] == "sporadic_AMR_carrier"
).sum()
n_context = (
    subset["selection_type"] == "intraspecific_context"
).sum()

duplicates = subset[ACC].duplicated().sum()

if duplicates != 0:
    raise ValueError(
        f"Duplicate accessions detected: {duplicates}"
    )

if n_carriers != 39:
    raise ValueError(
        f"Expected 39 carriers, found {n_carriers}."
    )

if n_context != 48:
    raise ValueError(
        f"Expected 48 context genomes, found {n_context}."
    )

if n_total != 87:
    raise ValueError(
        f"Expected 87 total genomes, found {n_total}."
    )


# ============================================================
# SAVE
# ============================================================

subset.to_csv(
    OUT_FILE,
    sep="\t",
    index=False
)

subset[ACC].to_csv(
    OUT_ACCESSIONS,
    index=False,
    header=False
)


# ============================================================
# SUMMARY
# ============================================================

summary = (
    subset.groupby(SPECIES)
    .agg(
        total_selected=(ACC, "nunique"),
        sporadic_carriers=(
            "selection_type",
            lambda x: (x == "sporadic_AMR_carrier").sum()
        ),
        context_genomes=(
            "selection_type",
            lambda x: (x == "intraspecific_context").sum()
        ),
        AMR_profiles_represented=("AMR_profile", "nunique")
    )
    .reset_index()
    .sort_values(
        ["total_selected", SPECIES],
        ascending=[False, True]
    )
)

log = [
    "BACILLUS AMR PROJECT",
    "STEP 14 - PHYLOGENETIC SUBSET SELECTION",
    "",
    f"Eligible species: {len(eligible_species)}",
    f"Sporadic AMR carriers: {n_carriers}",
    f"Context genomes per species: {N_CONTEXT}",
    f"Context genomes total: {n_context}",
    f"Final phylogenetic subset: {n_total}",
    "",
    "Selection rule:",
    "All sporadic AMR carriers retained.",
    "Three non-carrier context genomes selected per eligible species.",
    "Distinct AMR profiles prioritized; accession order used",
    "for deterministic selection and completion.",
    "",
    summary.to_string(index=False)
]

LOG_FILE.write_text(
    "\n".join(log) + "\n"
)

print("\n".join(log))
