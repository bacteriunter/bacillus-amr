#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

AMR_FILE = PROJECT / "05_amr" / "amr_presence_absence.tsv"
TAX_FILE = PROJECT / "03_taxonomy" / "qc_passed_taxonomy.tsv"

OUT_CASES = PROJECT / "05_amr" / "sporadic_amr_cases.tsv"
OUT_CARRIERS = PROJECT / "05_amr" / "sporadic_amr_carriers.tsv"
OUT_SUMMARY = PROJECT / "05_amr" / "sporadic_amr_carrier_summary.tsv"
LOG_FILE = PROJECT / "logs" / "13_sporadic_amr_carriers.txt"

MIN_SPECIES_N = 10
MAX_PREVALENCE = 5.0


# ============================================================
# LOAD DATA
# ============================================================

amr = pd.read_csv(AMR_FILE, sep="\t")
tax = pd.read_csv(TAX_FILE, sep="\t")

print("AMR matrix:", amr.shape)
print("Taxonomy:", tax.shape)
print()

print("AMR first columns:")
print(amr.columns[:10].tolist())
print()

# Detect accession column
accession_candidates = [
    "Assembly_Accession",
    "Assembly Accession",
    "accession",
    "Accession"
]

amr_acc = next((x for x in accession_candidates if x in amr.columns), None)
tax_acc = next((x for x in accession_candidates if x in tax.columns), None)

if amr_acc is None:
    raise ValueError("Accession column not found in AMR matrix.")

if tax_acc is None:
    raise ValueError("Accession column not found in taxonomy table.")

species_col = "Species_normalized"

if species_col not in tax.columns:
    raise ValueError(f"{species_col} not found in taxonomy table.")


# ============================================================
# MERGE AMR + TAXONOMY
# ============================================================

df = amr.merge(
    tax[[tax_acc, species_col]],
    left_on=amr_acc,
    right_on=tax_acc,
    how="left"
)

if tax_acc != amr_acc:
    df = df.drop(columns=[tax_acc])

if df[species_col].isna().any():
    missing = df.loc[df[species_col].isna(), amr_acc].tolist()
    raise ValueError(
        f"{len(missing)} genomes have no normalized taxonomy."
    )


# ============================================================
# SPECIES ELIGIBLE FOR WITHIN-SPECIES INFERENCE
# ============================================================

species_counts = (
    df.groupby(species_col)[amr_acc]
    .nunique()
    .sort_values(ascending=False)
)

eligible_species = species_counts[
    (species_counts >= MIN_SPECIES_N) &
    (species_counts.index != "Bacillus sp.")
].index.tolist()

print(f"Species with n >= {MIN_SPECIES_N}: {len(eligible_species)}")
print()


# ============================================================
# DETECT AMR COLUMNS
# ============================================================

metadata_cols = {
    amr_acc,
    species_col,
    "total_amr_determinants",
    "amr_count",
    "n_amr"
}

amr_cols = []

for col in amr.columns:
    if col == amr_acc or col in metadata_cols:
        continue

    values = pd.to_numeric(amr[col], errors="coerce")

    non_na = values.dropna()

    if len(non_na) > 0 and set(non_na.unique()).issubset({0, 1}):
        amr_cols.append(col)

print(f"Binary AMR determinant columns detected: {len(amr_cols)}")
print()


# ============================================================
# IDENTIFY SPORADIC SPECIES-DETERMINANT CASES
# prevalence >0 and <=5%
# ============================================================

cases = []
carriers = []

for species in eligible_species:

    sub = df[df[species_col] == species].copy()
    n_species = sub[amr_acc].nunique()

    for gene in amr_cols:

        values = pd.to_numeric(sub[gene], errors="coerce").fillna(0)

        carrier_mask = values == 1
        n_carriers = int(carrier_mask.sum())

        if n_carriers == 0:
            continue

        prevalence = 100.0 * n_carriers / n_species

        if prevalence <= MAX_PREVALENCE:

            cases.append({
                "species": species,
                "determinant": gene,
                "species_n": n_species,
                "carrier_n": n_carriers,
                "prevalence_pct": prevalence
            })

            carrier_accessions = sub.loc[
                carrier_mask, amr_acc
            ].tolist()

            for accession in carrier_accessions:
                carriers.append({
                    "accession": accession,
                    "species": species,
                    "determinant": gene,
                    "species_n": n_species,
                    "carrier_n": n_carriers,
                    "prevalence_pct": prevalence
                })


cases_df = pd.DataFrame(cases).sort_values(
    ["species", "prevalence_pct", "determinant"]
)

carriers_df = pd.DataFrame(carriers).sort_values(
    ["species", "accession", "determinant"]
)


# ============================================================
# UNIQUE-CARRIER SUMMARY
# ============================================================

summary = (
    carriers_df.groupby(["accession", "species"])
    .agg(
        n_sporadic_determinants=("determinant", "nunique"),
        sporadic_determinants=(
            "determinant",
            lambda x: ";".join(sorted(set(x)))
        )
    )
    .reset_index()
    .sort_values(
        ["species", "accession"]
    )
)


# ============================================================
# SAVE
# ============================================================

cases_df.to_csv(OUT_CASES, sep="\t", index=False)
carriers_df.to_csv(OUT_CARRIERS, sep="\t", index=False)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)


# ============================================================
# LOG
# ============================================================

unique_determinants = cases_df["determinant"].nunique()
unique_carriers = summary["accession"].nunique()
multi_carriers = (summary["n_sporadic_determinants"] > 1).sum()

log = [
    "BACILLUS AMR PROJECT",
    "STEP 13 - SPORADIC AMR CARRIERS",
    "",
    f"Minimum species size: {MIN_SPECIES_N}",
    f"Sporadic prevalence threshold: >0 to <= {MAX_PREVALENCE}%",
    "",
    f"Eligible species: {len(eligible_species)}",
    f"AMR determinants evaluated: {len(amr_cols)}",
    f"Sporadic species-determinant cases: {len(cases_df)}",
    f"Unique sporadic determinants: {unique_determinants}",
    f"Unique carrier genomes: {unique_carriers}",
    f"Carriers with >1 sporadic determinant: {multi_carriers}",
]

LOG_FILE.write_text("\n".join(log) + "\n")

print("\n".join(log))
print()
print("=== UNIQUE CARRIERS ===")
print(summary.to_string(index=False))
