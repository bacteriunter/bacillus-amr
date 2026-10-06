#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

META = PROJECT / "00_metadata" / "plant_associated_selected.tsv"
PASS = PROJECT / "02_qc" / "qc_passed_accessions.txt"

OUT = PROJECT / "03_taxonomy" / "qc_passed_taxonomy.tsv"
SUMMARY = PROJECT / "03_taxonomy" / "species_counts.tsv"
LOG = PROJECT / "logs" / "09_taxonomy_normalization.txt"

(PROJECT / "03_taxonomy").mkdir(exist_ok=True)

meta = pd.read_csv(META, sep="\t", dtype=str)

passed = set(
    pd.read_csv(PASS, header=None, dtype=str)[0]
)

df = meta[
    meta["Assembly Accession"].isin(passed)
].copy()


def normalize_species(name):
    if pd.isna(name):
        return "Unresolved"

    name = name.strip()
    parts = name.split()

    # Bacillus sp. strain/accession
    if len(parts) >= 2 and parts[0] == "Bacillus" and parts[1] == "sp.":
        return "Bacillus sp."

    # Named Bacillus species
    if len(parts) >= 2 and parts[0] == "Bacillus":
        return f"{parts[0]} {parts[1]}"

    return "Unresolved"


df["Species_normalized"] = (
    df["Organism Name"]
    .apply(normalize_species)
)

df["Taxonomic_resolution"] = "species"

df.loc[
    df["Species_normalized"] == "Bacillus sp.",
    "Taxonomic_resolution"
] = "genus"

df.loc[
    df["Species_normalized"] == "Unresolved",
    "Taxonomic_resolution"
] = "unresolved"


# Preserve complete original metadata plus normalized taxonomy
df.to_csv(
    OUT,
    sep="\t",
    index=False
)


counts = (
    df.groupby(
        ["Species_normalized", "Taxonomic_resolution"]
    )
    .size()
    .reset_index(name="Genome_count")
    .sort_values(
        "Genome_count",
        ascending=False
    )
)

counts.to_csv(
    SUMMARY,
    sep="\t",
    index=False
)


with open(LOG, "w") as f:
    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 09 - TAXONOMIC NAME NORMALIZATION\n\n")

    f.write(f"QC-passed genomes: {len(df)}\n")
    f.write(
        f"Species-level genomes: "
        f"{(df['Taxonomic_resolution'] == 'species').sum()}\n"
    )
    f.write(
        f"Genus-level Bacillus sp.: "
        f"{(df['Taxonomic_resolution'] == 'genus').sum()}\n"
    )
    f.write(
        f"Unresolved: "
        f"{(df['Taxonomic_resolution'] == 'unresolved').sum()}\n"
    )
    f.write(
        f"Normalized named species: "
        f"{df.loc[df['Taxonomic_resolution'] == 'species', 'Species_normalized'].nunique()}\n"
    )


print("QC-passed genomes:", len(df))
print(
    "Species-level:",
    (df["Taxonomic_resolution"] == "species").sum()
)
print(
    "Bacillus sp.:",
    (df["Taxonomic_resolution"] == "genus").sum()
)
print(
    "Unresolved:",
    (df["Taxonomic_resolution"] == "unresolved").sum()
)
print(
    "Normalized named species:",
    df.loc[
        df["Taxonomic_resolution"] == "species",
        "Species_normalized"
    ].nunique()
)

print("\nTOP 20:")
print(counts.head(20).to_string(index=False))
