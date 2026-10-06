#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
05_finalize_plant_selection.py

Final metadata-based classification of candidate plant-associated
Bacillus assemblies.

Each candidate is assigned to:
    INCLUDE
    EXCLUDE
    UNRESOLVED

Inputs
------
00_metadata/plant_association_classification.tsv
00_metadata/manual_exclusions_step04.tsv
00_metadata/manual_inclusions_step04.tsv

Outputs
-------
00_metadata/plant_association_final_classification.tsv
00_metadata/plant_associated_selected.tsv
00_metadata/plant_association_unresolved.tsv
logs/05_plant_selection_summary.txt

Important
---------
This script performs metadata selection only.
No genome-quality filtering or taxonomic filtering is performed here.
"""

from pathlib import Path
import pandas as pd
import re


# ============================================================
# PATHS
# ============================================================

CLASSIFICATION = Path(
    "00_metadata/plant_association_classification.tsv"
)

MANUAL_EXCLUSIONS = Path(
    "00_metadata/manual_exclusions_step04.tsv"
)

PROJECT_DECISIONS = Path(
    "00_metadata/manual_inclusions_step04.tsv"
)

OUT_ALL = Path(
    "00_metadata/plant_association_final_classification.tsv"
)

OUT_SELECTED = Path(
    "00_metadata/plant_associated_selected.tsv"
)

OUT_UNRESOLVED = Path(
    "00_metadata/plant_association_unresolved.tsv"
)

OUT_LOG = Path(
    "logs/05_plant_selection_summary.txt"
)


# ============================================================
# LOAD
# ============================================================

df = pd.read_csv(
    CLASSIFICATION,
    sep="\t",
    dtype=str
).fillna("")

manual_exclusions = pd.read_csv(
    MANUAL_EXCLUSIONS,
    sep="\t",
    dtype=str
).fillna("")

project_decisions = pd.read_csv(
    PROJECT_DECISIONS,
    sep="\t",
    dtype=str
).fillna("")


# ============================================================
# VALIDATION
# ============================================================

if df["Assembly Accession"].duplicated().any():
    raise ValueError(
        "Duplicate Assembly Accession values detected "
        "in candidate table."
    )

manual_exclusion_accessions = set(
    manual_exclusions["Assembly Accession"]
)

missing_exclusions = (
    manual_exclusion_accessions -
    set(df["Assembly Accession"])
)

if missing_exclusions:
    raise ValueError(
        "Manual exclusion accessions absent from candidate table: "
        + ", ".join(sorted(missing_exclusions))
    )


# ============================================================
# EXPLICIT DIRECT-ASSOCIATION PATTERN
# ============================================================

# Direct biological association accepted under the fixed study
# definition. Plant organs are accepted when they represent the
# biological isolation material.

DIRECT_PATTERN = re.compile(
    r"(?:"
    r"\brhizosphere\b|"
    r"\brhizospheric\b|"
    r"\brhizophere\b|"
    r"\brhisosphere\b|"
    r"\brhizoplane\b|"
    r"\brhizobacterium\b|"
    r"\bphyllosphere\b|"
    r"\bphylloplane\b|"
    r"\bendophyte\b|"
    r"\bendophytes\b|"
    r"\bendophytic\b|"
    r"\bendosphere\b|"
    r"\bplant tissue\b|"
    r"\bplant tissues\b|"
    r"\bplant root\b|"
    r"\bplant roots\b|"
    r"\bplant leaf\b|"
    r"\bplant leaves\b|"
    r"\bplant surface\b|"
    r"\broot surface\b|"
    r"\brootzone\b|"
    r"\broot\b|"
    r"\broots\b|"
    r"\bleaf\b|"
    r"\bleaves\b|"
    r"\bfoliar\b|"
    r"\bstem\b|"
    r"\bstems\b|"
    r"\bseed\b|"
    r"\bseeds\b|"
    r"\bseedling\b|"
    r"\bseedlings\b|"
    r"\bflower\b|"
    r"\bflowers\b|"
    r"\banther\b|"
    r"\bnodule\b|"
    r"\bnodules\b|"
    r"\bpetiole\b|"
    r"\bxylem sap\b"
    r")",
    flags=re.IGNORECASE
)


# ============================================================
# CLEARLY INCOMPATIBLE ISOLATION MATRICES
# ============================================================

# Applied ONLY to Isolation source.
# Manual reviewed decisions take precedence.

INCOMPATIBLE_SOURCE_PATTERN = re.compile(
    r"(?:"
    r"\bfood\b|"
    r"\bfermented\b|"
    r"\bfermentation\b|"
    r"\bfermentated\b|"
    r"\bpaste\b|"
    r"\bpuree\b|"
    r"\bpowder\b|"
    r"\bflour\b|"
    r"\bsilage\b|"
    r"\bstew\b|"
    r"\bdessert\b|"
    r"\bseed oil\b|"
    r"\bsalad leaves\b|"
    r"\bleaf litter\b|"
    r"\bcompost\b|"
    r"\bcomposting\b|"
    r"\bwastewater\b|"
    r"\bwaste water\b|"
    r"\bsewage\b|"
    r"\btreatment plant\b|"
    r"\bmanganese nodule\b|"
    r"\bgrainbin dust\b|"
    r"\bgastrointestinal tract\b"
    r")",
    flags=re.IGNORECASE
)


# ============================================================
# PROJECT-LEVEL DOCUMENTED INCLUSIONS
# ============================================================

documented_projects = set(
    project_decisions.loc[
        project_decisions["Decision"]
        .str.upper()
        .eq("INCLUDE"),
        "BioProject"
    ]
)

unresolved_projects = set(
    project_decisions.loc[
        project_decisions["Decision"]
        .str.upper()
        .eq("UNRESOLVED"),
        "BioProject"
    ]
)


# ============================================================
# CLASSIFY
# ============================================================

decisions = []
reasons = []
evidence_types = []

for _, row in df.iterrows():

    accession = row["Assembly Accession"]
    bioproject = row["Assembly BioProject Accession"]

    source = row[
        "Assembly BioSample Isolation source"
    ].strip()

    host = row[
        "Assembly BioSample Host"
    ].strip()

    title = row[
        "Assembly BioSample Description Title"
    ].strip()

    comment = row[
        "Assembly BioSample Description Comment"
    ].strip()

    # --------------------------------------------------------
    # 1. Manually reviewed exclusions have highest priority
    # --------------------------------------------------------

    if accession in manual_exclusion_accessions:

        reason = manual_exclusions.loc[
            manual_exclusions["Assembly Accession"]
            .eq(accession),
            "Reason"
        ].iloc[0]

        decisions.append("EXCLUDE")
        reasons.append(
            f"Manually reviewed exclusion: {reason}"
        )
        evidence_types.append("manual_exclusion")
        continue

    # --------------------------------------------------------
    # 2. Explicit direct association in Isolation source
    # --------------------------------------------------------

    if (
        source and
        DIRECT_PATTERN.search(source) and
        not INCOMPATIBLE_SOURCE_PATTERN.search(source)
    ):

        decisions.append("INCLUDE")
        reasons.append(
            "Isolation source explicitly indicates direct "
            "plant association or plant biological material."
        )
        evidence_types.append("isolation_source")
        continue

    # --------------------------------------------------------
    # 3. Documented project-level inclusion
    # --------------------------------------------------------

    if bioproject in documented_projects:

        project_row = project_decisions.loc[
            (project_decisions["BioProject"] == bioproject) &
            (
                project_decisions["Decision"]
                .str.upper()
                .eq("INCLUDE")
            )
        ].iloc[0]

        metadata_value = (
            project_row["Metadata_value"]
            .strip()
            .lower()
        )

        if (
            source.lower() == metadata_value and
            host
        ):
            decisions.append("INCLUDE")
            reasons.append(
                "Documented project-level plant association: "
                + project_row["Reason"]
            )
            evidence_types.append(
                "documented_project"
            )
            continue

    # --------------------------------------------------------
    # 4. Explicitly unresolved project metadata
    # --------------------------------------------------------

    if bioproject in unresolved_projects:

        project_row = project_decisions.loc[
            (project_decisions["BioProject"] == bioproject) &
            (
                project_decisions["Decision"]
                .str.upper()
                .eq("UNRESOLVED")
            )
        ].iloc[0]

        metadata_value = (
            project_row["Metadata_value"]
            .strip()
            .lower()
        )

        if source.lower() == metadata_value:

            decisions.append("UNRESOLVED")
            reasons.append(
                "Project-level metadata remains ambiguous: "
                + project_row["Reason"]
            )
            evidence_types.append(
                "unresolved_project"
            )
            continue

    # --------------------------------------------------------
    # 5. Clearly incompatible isolation source
    # --------------------------------------------------------

    if (
        source and
        INCOMPATIBLE_SOURCE_PATTERN.search(source)
    ):

        decisions.append("EXCLUDE")
        reasons.append(
            "Isolation source represents a non-direct or "
            "processed matrix incompatible with the study "
            "inclusion definition."
        )
        evidence_types.append(
            "incompatible_source"
        )
        continue

    # --------------------------------------------------------
    # 6. Direct evidence elsewhere, but source unavailable
    # --------------------------------------------------------

    other_text = " | ".join(
        x for x in [host, title, comment] if x
    )

    source_missing = (
        not source or
        source.lower() in {
            "missing",
            "not collected",
            "not applicable",
            "na",
            "n/a"
        }
    )

    if (
        source_missing and
        DIRECT_PATTERN.search(other_text)
    ):

        decisions.append("INCLUDE")
        reasons.append(
            "Isolation source is unavailable, but other "
            "BioSample metadata explicitly indicates direct "
            "plant association or plant biological material."
        )
        evidence_types.append(
            "other_explicit_metadata"
        )
        continue

    # --------------------------------------------------------
    # 7. Remaining candidates are unresolved
    # --------------------------------------------------------

    decisions.append("UNRESOLVED")
    reasons.append(
        "Candidate metadata does not provide sufficiently "
        "explicit evidence to demonstrate direct plant "
        "association under the fixed inclusion criteria."
    )
    evidence_types.append(
        "insufficient_metadata"
    )


# ============================================================
# ADD FINAL COLUMNS
# ============================================================

df["plant_final_decision"] = decisions
df["plant_final_reason"] = reasons
df["plant_final_evidence"] = evidence_types


# ============================================================
# SANITY CHECK
# ============================================================

valid = {"INCLUDE", "EXCLUDE", "UNRESOLVED"}

observed = set(df["plant_final_decision"])

if not observed.issubset(valid):
    raise ValueError(
        f"Unexpected final decisions: {observed}"
    )

if len(df) != 2550:
    raise ValueError(
        f"Expected 2550 candidate assemblies; found {len(df)}."
    )


# ============================================================
# EXPORT
# ============================================================

selected = df[
    df["plant_final_decision"] == "INCLUDE"
].copy()

unresolved = df[
    df["plant_final_decision"] == "UNRESOLVED"
].copy()

OUT_ALL.parent.mkdir(parents=True, exist_ok=True)
OUT_LOG.parent.mkdir(parents=True, exist_ok=True)

df.to_csv(
    OUT_ALL,
    sep="\t",
    index=False
)

selected.to_csv(
    OUT_SELECTED,
    sep="\t",
    index=False
)

unresolved.to_csv(
    OUT_UNRESOLVED,
    sep="\t",
    index=False
)


# ============================================================
# SUMMARY
# ============================================================

summary = []

summary.append("BACILLUS AMR PROJECT")
summary.append("FINAL PLANT-ASSOCIATION METADATA SELECTION")
summary.append("=" * 70)

summary.append(
    f"\nCandidate assemblies evaluated: {len(df)}"
)

summary.append("\nFinal decisions:")
summary.append(
    df["plant_final_decision"]
    .value_counts()
    .to_string()
)

summary.append("\nEvidence categories:")
summary.append(
    df["plant_final_evidence"]
    .value_counts()
    .to_string()
)

summary.append(
    "\nSelected assemblies by BioProject (top 20):"
)

summary.append(
    selected["Assembly BioProject Accession"]
    .value_counts()
    .head(20)
    .to_string()
)

summary.append(
    "\nSelected assemblies by organism (top 30):"
)

summary.append(
    selected["Organism Name"]
    .value_counts()
    .head(30)
    .to_string()
)

summary.append(
    "\nIMPORTANT: UNRESOLVED records are retained in the "
    "audit table but are not included in the analytical "
    "dataset."
)

with open(
    OUT_LOG,
    "w",
    encoding="utf-8"
) as handle:

    handle.write("\n".join(summary))
    handle.write("\n")


# ============================================================
# TERMINAL OUTPUT
# ============================================================

print(f"Candidates evaluated: {len(df)}")

print("\nFinal decisions:")
print(
    df["plant_final_decision"]
    .value_counts()
    .to_string()
)

print("\nEvidence categories:")
print(
    df["plant_final_evidence"]
    .value_counts()
    .to_string()
)

print(
    f"\nSelected dataset: {len(selected)} assemblies"
)

print(
    f"Unresolved: {len(unresolved)} assemblies"
)

print("\nOutputs:")
print(f"  {OUT_ALL}")
print(f"  {OUT_SELECTED}")
print(f"  {OUT_UNRESOLVED}")
print(f"  {OUT_LOG}")
