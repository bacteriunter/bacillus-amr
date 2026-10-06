#!/usr/bin/env python3

import pandas as pd
import re
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

INPUT = PROJECT / "00_metadata" / "bacillus_master_metadata.tsv"
OUTPUT = PROJECT / "00_metadata" / "plant_associated_candidates.tsv"

df = pd.read_csv(INPUT, sep="\t", dtype=str).fillna("")

# ---------------------------------------------------------
# Fields used as metadata evidence
# ---------------------------------------------------------

fields = {
    "host": "Assembly BioSample Host",
    "source": "Assembly BioSample Isolation source",
    "title": "Assembly BioSample Description Title",
    "comment": "Assembly BioSample Description Comment",
}

# ---------------------------------------------------------
# Broad candidate vocabulary
# This is intentionally sensitive rather than definitive.
# Final inclusion will be audited subsequently.
# ---------------------------------------------------------

plant_terms = re.compile(
    r"\b(?:"
    r"plant|plants|"
    r"root|roots|rootzone|root-soil|"
    r"rhizo\w*|"
    r"leaf|leaves|"
    r"seed|seeds|seedling|seedlings|"
    r"stem|stems|"
    r"endoph\w*|"
    r"phyllo\w*|"
    r"flower|flowers|"
    r"anther|"
    r"nodule|nodules|"
    r"maize|corn|"
    r"wheat|"
    r"rice|"
    r"soybean|"
    r"tomato|"
    r"potato|"
    r"arabidopsis|"
    r"tobacco|"
    r"cotton|"
    r"cucumber|"
    r"pepper|"
    r"banana|"
    r"citrus|"
    r"alfalfa"
    r")\b",
    re.IGNORECASE
)

# ---------------------------------------------------------
# Detect evidence independently in each field
# ---------------------------------------------------------

for short, col in fields.items():
    df[f"plant_match_{short}"] = (
        df[col]
        .str.contains(plant_terms, na=False)
    )

match_cols = [
    "plant_match_host",
    "plant_match_source",
    "plant_match_title",
    "plant_match_comment",
]

df["plant_match_count"] = df[match_cols].sum(axis=1)

df["plant_candidate"] = df["plant_match_count"] > 0

# Record which metadata fields generated the candidate
def evidence_fields(row):
    evidence = []

    for short in fields:
        if row[f"plant_match_{short}"]:
            evidence.append(short)

    return ";".join(evidence)

df["plant_evidence_fields"] = df.apply(evidence_fields, axis=1)

# ---------------------------------------------------------
# Export candidates only
# ---------------------------------------------------------

candidates = df[df["plant_candidate"]].copy()

candidates.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

# ---------------------------------------------------------
# Summary
# ---------------------------------------------------------

print("BACILLUS AMR PROJECT")
print("Step 03: Identify plant-associated metadata candidates")
print()
print(f"Total assemblies:          {len(df)}")
print(f"Plant candidates:          {len(candidates)}")
print()
print("Evidence by metadata field:")
for col in match_cols:
    print(f"  {col}: {df[col].sum()}")

print()
print("Number of matching evidence fields:")
print(candidates["plant_match_count"].value_counts().sort_index())

print()
print(f"Output: {OUTPUT}")
