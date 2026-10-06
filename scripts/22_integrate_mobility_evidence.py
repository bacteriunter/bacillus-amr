#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

IS_FILE = Path("08_genomic_context/amr_IS_proximity.tsv")

OUT = Path("08_genomic_context/mobility_evidence.tsv")
SUMMARY = Path("08_genomic_context/mobility_evidence_summary.tsv")
LOG = Path("logs/22_integrate_mobility_evidence.txt")

# ============================================================
# 1. IS-associated proteins detected in Steps 20-21
# ============================================================

df = pd.read_csv(IS_FILE, sep="\t")

is_hits = df[df["IS_detected"] == True].copy()

records = []

for _, r in is_hits.iterrows():

    records.append({
        "locus_id": r["locus_id"],
        "assembly_accession": r["assembly_accession"],
        "determinant": r["determinant"],
        "context_status": r["context_status"],
        "mobility_class": "insertion_sequence_transposase",
        "protein_id": r["IS_protein_id"],
        "reference": r["IS_reference"],
        "annotation": r["IS_family_annotation"],
        "protein_start_local": int(r["IS_start_local"]),
        "protein_stop_local": int(r["IS_stop_local"]),
        "protein_strand": r["IS_strand"],
        "distance_to_AMR_bp": int(r["distance_bp"]),
        "relative_position": r["relative_position"],
        "evidence_source": "Prokka bacterial IS database / BLASTP"
    })


# ============================================================
# 2. L024 transposon-associated recombinase/resolvase
#
# Supported by convergent Swiss-Prot hits:
#   bin3 - Putative transposon Tn552 DNA-invertase
#   tnpR - Transposon Tn3 resolvase
#   tnpR - Transposon gamma-delta resolvase
#   hin  - DNA-invertase
#
# Conservative functional assignment:
# transposon-associated recombinase/resolvase
#
# Coordinates from Prodigal:
# 9908-10615, strand -
#
# AMRFinder fosB coordinates in local window:
# 7064-7477
#
# Edge-to-edge distance:
# 9908 - 7477 - 1 = 2430 bp
# ============================================================

l024 = df[df["locus_id"].str.startswith("L024_")]

if len(l024) == 0:
    raise RuntimeError("L024 not found in AMR/IS proximity table")

r = l024.iloc[0]

records.append({
    "locus_id": r["locus_id"],
    "assembly_accession": r["assembly_accession"],
    "determinant": r["determinant"],
    "context_status": r["context_status"],
    "mobility_class": "transposon_associated_recombinase_resolvase",
    "protein_id": "L024_GCF_002568725.1_fosB_13",
    "reference": "Swiss-Prot convergent hits: bin3/tnpR/hin",
    "annotation": "transposon-associated recombinase/resolvase",
    "protein_start_local": 9908,
    "protein_stop_local": 10615,
    "protein_strand": "-",
    "distance_to_AMR_bp": 2430,
    "relative_position": "downstream",
    "evidence_source": "Prokka bacterial Swiss-Prot database / BLASTP"
})


mobility = pd.DataFrame(records)

mobility = mobility.sort_values(
    ["locus_id", "distance_to_AMR_bp", "mobility_class"]
)

mobility.to_csv(OUT, sep="\t", index=False)


# ============================================================
# 3. One-row-per-AMR-locus summary
# ============================================================

all_loci = (
    df[
        [
            "locus_id",
            "assembly_accession",
            "determinant",
            "context_status"
        ]
    ]
    .drop_duplicates()
    .copy()
)

counts = (
    mobility.groupby("locus_id")
    .agg(
        n_mobility_proteins=("protein_id", "count"),
        n_IS_transposases=(
            "mobility_class",
            lambda x: (x == "insertion_sequence_transposase").sum()
        ),
        n_recombinase_resolvase=(
            "mobility_class",
            lambda x: (
                x == "transposon_associated_recombinase_resolvase"
            ).sum()
        ),
        nearest_mobility_distance_bp=("distance_to_AMR_bp", "min")
    )
    .reset_index()
)

summary = all_loci.merge(
    counts,
    on="locus_id",
    how="left"
)

for col in [
    "n_mobility_proteins",
    "n_IS_transposases",
    "n_recombinase_resolvase"
]:
    summary[col] = summary[col].fillna(0).astype(int)

summary["mobility_evidence_detected"] = (
    summary["n_mobility_proteins"] > 0
)

summary.to_csv(SUMMARY, sep="\t", index=False)


# ============================================================
# 4. Log
# ============================================================

n_loci = len(summary)
n_positive = summary["mobility_evidence_detected"].sum()

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 22 - INTEGRATED LOCAL MOBILITY EVIDENCE\n\n")

    f.write(f"Physical AMR loci: {n_loci}\n")
    f.write(
        f"IS-associated proteins: "
        f"{(mobility['mobility_class'] == 'insertion_sequence_transposase').sum()}\n"
    )
    f.write(
        f"Transposon-associated recombinase/resolvase proteins: "
        f"{(mobility['mobility_class'] == 'transposon_associated_recombinase_resolvase').sum()}\n"
    )
    f.write(
        f"Total mobility-associated proteins: {len(mobility)}\n"
    )
    f.write(
        f"AMR loci with >=1 local mobility-associated protein: {n_positive}\n"
    )
    f.write(
        f"AMR loci without detected local mobility-associated protein: "
        f"{n_loci - n_positive}\n"
    )

print(summary.to_string(index=False))

print(f"\nOutput:  {OUT}")
print(f"Summary: {SUMMARY}")
print(f"Log:     {LOG}")
