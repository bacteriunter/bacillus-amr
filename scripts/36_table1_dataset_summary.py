#!/usr/bin/env python3

from pathlib import Path
import pandas as pd


# ============================================================
# Files
# ============================================================

CANDIDATES_FILE = Path(
    "00_metadata/plant_association_classification.tsv"
)

QC_FILE = Path(
    "02_qc/seqkit_genome_stats.tsv"
)

TAXONOMY_FILE = Path(
    "03_taxonomy/qc_passed_taxonomy.tsv"
)

AMR_HITS_FILE = Path(
    "05_amr/amrfinder_amr_hits.tsv"
)

AMR_SUMMARY_FILE = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT = Path(
    "11_tables/Table1_dataset_QC_AMR_summary.tsv"
)

LOG = Path(
    "logs/36_table1_dataset_summary.txt"
)


# ============================================================
# Load
# ============================================================

candidates = pd.read_csv(
    CANDIDATES_FILE,
    sep="\t"
)

qc = pd.read_csv(
    QC_FILE,
    sep="\t"
)

taxonomy = pd.read_csv(
    TAXONOMY_FILE,
    sep="\t"
)

amr_hits = pd.read_csv(
    AMR_HITS_FILE,
    sep="\t"
)

amr = pd.read_csv(
    AMR_SUMMARY_FILE,
    sep="\t"
)


# ============================================================
# Validation
# ============================================================

if len(candidates) != 2550:
    raise RuntimeError(
        f"Expected 2550 plant-associated candidates, "
        f"found {len(candidates)}"
    )

if len(qc) != 1374:
    raise RuntimeError(
        f"Expected 1374 genomes entering structural QC, "
        f"found {len(qc)}"
    )

if len(taxonomy) != 1312:
    raise RuntimeError(
        f"Expected 1312 QC-passed genomes, "
        f"found {len(taxonomy)}"
    )

if len(amr) != 1312:
    raise RuntimeError(
        f"Expected 1312 genomes in AMR summary, "
        f"found {len(amr)}"
    )


# ============================================================
# QC-passed structural statistics
# ============================================================

passed_accessions = set(
    taxonomy["Assembly Accession"]
)

qc_pass = qc[
    qc["Assembly_Accession"].isin(
        passed_accessions
    )
].copy()

if len(qc_pass) != 1312:
    raise RuntimeError(
        f"Expected 1312 QC-passed genomes after merge, "
        f"found {len(qc_pass)}"
    )


# ============================================================
# Taxonomic summary
# ============================================================

species_level = (
    taxonomy["Taxonomic_resolution"] == "species"
).sum()

genus_level = (
    taxonomy["Taxonomic_resolution"] != "species"
).sum()

named_species = (
    taxonomy.loc[
        taxonomy["Taxonomic_resolution"] == "species",
        "Species_normalized"
    ]
    .nunique()
)


# ============================================================
# AMR summary
# ============================================================

genomes_with_amr = (
    amr["AMR_determinant_count"] > 0
).sum()

genomes_without_amr = (
    amr["AMR_determinant_count"] == 0
).sum()

unique_determinants = (
    amr_hits["Element symbol"]
    .nunique()
)


# ============================================================
# Helper
# ============================================================

rows = []


def add(section, metric, value):

    rows.append({
        "Section": section,
        "Metric": metric,
        "Value": value
    })


# ============================================================
# Dataset selection
# ============================================================

add(
    "Dataset selection",
    "Plant-associated candidate assemblies",
    f"{len(candidates):,}"
)

add(
    "Dataset selection",
    "Assemblies selected after plant-association assessment",
    f"{len(qc):,}"
)

add(
    "Dataset selection",
    "Assemblies retained after structural QC",
    f"{len(qc_pass):,}"
)

add(
    "Dataset selection",
    "Assemblies excluded by structural QC",
    f"{len(qc) - len(qc_pass):,}"
)


# ============================================================
# Genome quality and taxonomy
# ============================================================

add(
    "Genome quality and taxonomy",
    "Contigs per genome, median (range)",
    (
        f"{qc_pass['Num_contigs'].median():.0f} "
        f"({qc_pass['Num_contigs'].min():.0f}–"
        f"{qc_pass['Num_contigs'].max():.0f})"
    )
)

add(
    "Genome quality and taxonomy",
    "Genome length, median Mb (range)",
    (
        f"{qc_pass['Total_length'].median()/1e6:.2f} "
        f"({qc_pass['Total_length'].min()/1e6:.2f}–"
        f"{qc_pass['Total_length'].max()/1e6:.2f})"
    )
)

add(
    "Genome quality and taxonomy",
    "N50, median kb (range)",
    (
        f"{qc_pass['N50'].median()/1e3:.1f} "
        f"({qc_pass['N50'].min()/1e3:.1f}–"
        f"{qc_pass['N50'].max()/1e3:.1f})"
    )
)

add(
    "Genome quality and taxonomy",
    "GC content, median % (range)",
    (
        f"{qc_pass['GC_percent'].median():.2f} "
        f"({qc_pass['GC_percent'].min():.2f}–"
        f"{qc_pass['GC_percent'].max():.2f})"
    )
)

add(
    "Genome quality and taxonomy",
    "Species-level assignments",
    f"{species_level:,}"
)

add(
    "Genome quality and taxonomy",
    "Genus-level Bacillus sp. assignments",
    f"{genus_level:,}"
)

add(
    "Genome quality and taxonomy",
    "Named species represented",
    f"{named_species:,}"
)


# ============================================================
# Global resistome
# ============================================================

add(
    "Global resistome",
    "Genomes with ≥1 AMR determinant",
    (
        f"{genomes_with_amr:,} "
        f"({100*genomes_with_amr/len(amr):.1f}%)"
    )
)

add(
    "Global resistome",
    "Genomes without detected AMR determinants",
    (
        f"{genomes_without_amr:,} "
        f"({100*genomes_without_amr/len(amr):.1f}%)"
    )
)

add(
    "Global resistome",
    "AMR hits",
    f"{len(amr_hits):,}"
)

add(
    "Global resistome",
    "Unique AMR determinants",
    f"{unique_determinants:,}"
)

add(
    "Global resistome",
    "AMR determinants per genome, mean ± SD",
    (
        f"{amr['AMR_determinant_count'].mean():.2f} ± "
        f"{amr['AMR_determinant_count'].std():.2f}"
    )
)

add(
    "Global resistome",
    "AMR determinants per genome, median (range)",
    (
        f"{amr['AMR_determinant_count'].median():.0f} "
        f"({amr['AMR_determinant_count'].min():.0f}–"
        f"{amr['AMR_determinant_count'].max():.0f})"
    )
)


# ============================================================
# Save
# ============================================================

table = pd.DataFrame(rows)

table.to_csv(
    OUT,
    sep="\t",
    index=False
)


# ============================================================
# Log
# ============================================================

with open(LOG, "w") as f:

    f.write(
        "BACILLUS AMR PROJECT\n"
        "TABLE 1 - DATASET, QC AND GLOBAL AMR SUMMARY\n\n"
    )

    f.write(
        table.to_string(index=False)
    )

    f.write("\n")


print(table.to_string(index=False))

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
