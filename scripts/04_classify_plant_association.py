#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Classify preliminary Bacillus plant-associated candidates according to
metadata evidence.

This step does NOT remove assemblies.

Plant-associated is defined strictly as explicit evidence of:
- plant or plant tissue/organ
- plant surface
- endosphere/endophyte
- rhizoplane
- rhizosphere

Agricultural soil, crop fields, compost, plant residues, foods, fermented
products, and processed plant materials are not sufficient by themselves.
"""

from pathlib import Path
import pandas as pd
import re


PROJECT_DIR = Path.home() / "bacillus_amr"

INPUT = (
    PROJECT_DIR
    / "00_metadata"
    / "plant_associated_candidates.tsv"
)

OUTPUT = (
    PROJECT_DIR
    / "00_metadata"
    / "plant_association_classification.tsv"
)


# ============================================================
# DIRECT PLANT-ASSOCIATION EVIDENCE
# ============================================================

DIRECT_PATTERN = re.compile(
    r"\b(?:"
    r"rhizosphere|rhizospheric|rhizophere|rhisosphere|"
    r"rhizoplane|rhizobacterium|"
    r"phyllosphere|phylloplane|"
    r"endophyte|endophytes|endophytic|"
    r"endosphere|"
    r"plant tissue|plant tissues|"
    r"plant root|plant roots|"
    r"plant leaf|plant leaves|"
    r"plant surface|"
    r"root|roots|rootzone|"
    r"leaf|leaves|foliar|"
    r"stem|stems|"
    r"seed|seeds|seedling|seedlings|"
    r"flower|flowers|anther|"
    r"nodule|nodules"
    r")\b",
    flags=re.IGNORECASE,
)


EVIDENCE_COLUMNS = {
    "host": "Assembly BioSample Host",
    "source": "Assembly BioSample Isolation source",
    "title": "Assembly BioSample Description Title",
    "comment": "Assembly BioSample Description Comment",
}


df = pd.read_csv(INPUT, sep="\t", dtype=str).fillna("")


# ============================================================
# IDENTIFY DIRECT EVIDENCE IN EACH METADATA FIELD
# ============================================================

direct_columns = []

for label, column in EVIDENCE_COLUMNS.items():

    outcol = f"direct_match_{label}"

    df[outcol] = (
        df[column]
        .astype(str)
        .str.contains(DIRECT_PATTERN, na=False)
    )

    direct_columns.append(outcol)


df["direct_match_count"] = (
    df[direct_columns]
    .astype(int)
    .sum(axis=1)
)


df["direct_evidence_fields"] = df.apply(
    lambda row: ";".join(
        label
        for label in EVIDENCE_COLUMNS
        if row[f"direct_match_{label}"]
    ),
    axis=1,
)


df["direct_plant_evidence"] = (
    df["direct_match_count"] > 0
)


# ============================================================
# SAVE AUDIT TABLE
# ============================================================

df.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ============================================================
# REPORT
# ============================================================

print("BACILLUS AMR PROJECT")
print("Step 04A: Detect strict direct plant-association evidence")
print()

print(f"Candidate assemblies:       {len(df)}")
print(
    f"Direct evidence detected:   "
    f"{df['direct_plant_evidence'].sum()}"
)
print(
    f"No direct evidence yet:     "
    f"{(~df['direct_plant_evidence']).sum()}"
)

print("\nDirect evidence by metadata field:")

for col in direct_columns:
    print(f"  {col}: {df[col].sum()}")

print("\nNumber of direct evidence fields:")
print(
    df["direct_match_count"]
    .value_counts()
    .sort_index()
)

print()
print(f"Output: {OUTPUT}")
