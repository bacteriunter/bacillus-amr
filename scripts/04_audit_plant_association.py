#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
04_audit_plant_association.py

Audit trail for plant-association classification.

This script summarizes the evidence used to evaluate candidate
plant-associated Bacillus assemblies.

It does NOT modify, include, or exclude assemblies.

Inputs
------
00_metadata/plant_association_classification.tsv
00_metadata/manual_exclusions_step04.tsv

Output
------
logs/04_plant_association_audit.txt
"""

from pathlib import Path
import pandas as pd


# ============================================================
# PATHS
# ============================================================

INPUT = Path("00_metadata/plant_association_classification.tsv")
EXCLUSIONS = Path("00_metadata/manual_exclusions_step04.tsv")
OUTPUT = Path("logs/04_plant_association_audit.txt")


# ============================================================
# LOAD DATA
# ============================================================

df = pd.read_csv(INPUT, sep="\t", dtype=str).fillna("")

if EXCLUSIONS.exists():
    exclusions = pd.read_csv(EXCLUSIONS, sep="\t", dtype=str).fillna("")
else:
    exclusions = pd.DataFrame()


# ============================================================
# BASIC GROUPS
# ============================================================

direct = df["direct_plant_evidence"].str.lower().eq("true")
no_direct = ~direct

host_only = (
    no_direct &
    df["plant_evidence_fields"].eq("host")
)

core_emissing = (
    host_only &
    df["Assembly BioSample Isolation source"]
    .str.strip()
    .str.lower()
    .isin(["core", "emissingophyte"])
)


# ============================================================
# REPORT
# ============================================================

lines = []

lines.append("BACILLUS AMR PROJECT")
lines.append("PLANT-ASSOCIATION METADATA AUDIT")
lines.append("=" * 80)

lines.append("\n1. CANDIDATE CLASSIFICATION")
lines.append("-" * 80)
lines.append(f"Candidate assemblies: {len(df)}")
lines.append(f"Direct plant evidence detected: {direct.sum()}")
lines.append(f"No direct plant evidence detected: {no_direct.sum()}")

lines.append("\nDirect evidence match count:")
lines.append(
    df.loc[direct, "direct_match_count"]
    .value_counts()
    .sort_index()
    .to_string()
)

lines.append("\n2. EVIDENCE FIELDS AMONG RECORDS WITHOUT DIRECT EVIDENCE")
lines.append("-" * 80)
lines.append(
    df.loc[no_direct, "plant_evidence_fields"]
    .value_counts()
    .to_string()
)

lines.append("\n3. HOST-ONLY CANDIDATES")
lines.append("-" * 80)
lines.append(f"Host-only candidates: {host_only.sum()}")

lines.append("\nTop 50 hosts:")
lines.append(
    df.loc[host_only, "Assembly BioSample Host"]
    .str.strip()
    .replace("", "<EMPTY>")
    .value_counts()
    .head(50)
    .to_string()
)

lines.append("\nTop 50 isolation sources:")
lines.append(
    df.loc[host_only, "Assembly BioSample Isolation source"]
    .str.strip()
    .replace("", "<EMPTY>")
    .value_counts()
    .head(50)
    .to_string()
)

lines.append("\n4. CORE / EMISSINGOPHYTE BLOCK")
lines.append("-" * 80)
lines.append(f"Records: {core_emissing.sum()}")

lines.append("\nIsolation source distribution:")
lines.append(
    df.loc[core_emissing, "Assembly BioSample Isolation source"]
    .value_counts()
    .to_string()
)

lines.append("\nHost distribution:")
lines.append(
    df.loc[core_emissing, "Assembly BioSample Host"]
    .value_counts()
    .to_string()
)

lines.append("\nBioProject distribution:")
lines.append(
    df.loc[core_emissing, "Assembly BioProject Accession"]
    .value_counts()
    .to_string()
)

lines.append("\nTop 20 organisms:")
lines.append(
    df.loc[core_emissing, "Organism Name"]
    .value_counts()
    .head(20)
    .to_string()
)

lines.append("\nRepresentative records:")
representatives = (
    df.loc[core_emissing]
    .groupby(
        [
            "Assembly BioSample Isolation source",
            "Assembly BioSample Host"
        ],
        group_keys=False
    )
    .head(3)
)

rep_cols = [
    "Assembly Accession",
    "Organism Name",
    "Assembly BioProject Accession",
    "Assembly BioSample Accession",
    "Assembly BioSample Host",
    "Assembly BioSample Isolation source",
    "Assembly BioSample Description Title",
    "Assembly BioSample Description Comment"
]

lines.append(
    representatives[rep_cols].to_string(index=False)
)


# ============================================================
# MANUAL EXCLUSIONS
# ============================================================

lines.append("\n5. MANUALLY REVIEWED EXCLUSIONS")
lines.append("-" * 80)

if len(exclusions):
    lines.append(f"Confirmed exclusions: {len(exclusions)}")
    lines.append(exclusions.to_string(index=False))
else:
    lines.append("No manual exclusion file found.")


# ============================================================
# DOCUMENTED RESCUES
# ============================================================

lines.append("\n6. DOCUMENTED RESCUES FROM CONTEXTUAL REVIEW")
lines.append("-" * 80)

rescues = [
    (
        "GCF_002214765.1",
        "Retained: Isolation source is Sesame leaf and BioSample comment "
        "explicitly states that the sample was isolated from a sesame leaf."
    ),
    (
        "GCF_004294525.1",
        "Retained: Isolation source is Maize rhizosphere; food-related wording "
        "occurs only in another metadata field."
    ),
    (
        "GCF_041950475.1",
        "Retained: Isolation source is tobacco leaves; fermentation-related "
        "wording occurs only in another metadata field."
    ),
]

for accession, reason in rescues:
    lines.append(f"{accession}\t{reason}")


# ============================================================
# METHODOLOGICAL RULE
# ============================================================

lines.append("\n7. CURRENT OPERATIONAL DEFINITION")
lines.append("-" * 80)

lines.append(
    "Include only isolates whose metadata explicitly supports direct "
    "association with a plant, plant tissue/organ, plant surface, "
    "endosphere/endophyte, rhizoplane, or rhizosphere."
)

lines.append(
    "Agricultural soil, crop fields, compost, plant residues, leaf litter, "
    "foods, fermented products, and processed plant materials are not "
    "sufficient by themselves."
)

lines.append(
    "Plant organs such as leaf, root, stem, seed, flower, or nodule are "
    "acceptable when they represent the biological material from which "
    "the isolate was obtained."
)

lines.append(
    "Exclusion terms are interpreted in metadata context rather than "
    "applied globally across all fields."
)

lines.append(
    "No assembly is automatically excluded solely because another metadata "
    "field contains food- or fermentation-related terminology when the "
    "isolation source itself provides explicit direct plant association."
)


# ============================================================
# WRITE REPORT
# ============================================================

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

with open(OUTPUT, "w", encoding="utf-8") as handle:
    handle.write("\n".join(lines))
    handle.write("\n")

print(f"Audit written to: {OUTPUT}")
print(f"Candidate assemblies: {len(df)}")
print(f"Direct evidence: {direct.sum()}")
print(f"No direct evidence: {no_direct.sum()}")
print(f"Host-only candidates: {host_only.sum()}")
print(f"Core/emissingophyte records: {core_emissing.sum()}")
print(f"Manual exclusions recorded: {len(exclusions)}")
