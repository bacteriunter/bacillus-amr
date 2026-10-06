#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

HITS_FILE = PROJECT / "05_amr" / "amrfinder_amr_hits.tsv"
CARRIERS_FILE = PROJECT / "05_amr" / "sporadic_amr_carriers.tsv"

OUT_FILE = PROJECT / "08_genomic_context" / "sporadic_amr_loci.tsv"
LOG_FILE = PROJECT / "logs" / "15_sporadic_amr_loci.txt"

OUT_FILE.parent.mkdir(parents=True, exist_ok=True)

hits = pd.read_csv(HITS_FILE, sep="\t")
carriers = pd.read_csv(CARRIERS_FILE, sep="\t")

# ------------------------------------------------------------
# Strict match:
# accession + determinant
# vs.
# Assembly_Accession + Element symbol
# ------------------------------------------------------------

loci = carriers.merge(
    hits,
    left_on=["accession", "determinant"],
    right_on=["Assembly_Accession", "Element symbol"],
    how="left",
    validate="one_to_many"
)

# Check unmatched candidate cases
unmatched = loci["Contig id"].isna()

if unmatched.any():
    print("\nERROR: candidate cases without matching AMRFinder locus:")
    print(
        loci.loc[
            unmatched,
            ["accession", "species", "determinant"]
        ].to_string(index=False)
    )
    raise SystemExit(1)

# Number of physical loci per genome-determinant case
locus_counts = (
    loci.groupby(["accession", "determinant"])
    .size()
    .reset_index(name="n_loci")
)

loci = loci.merge(
    locus_counts,
    on=["accession", "determinant"],
    how="left"
)

# Select useful fields
keep = [
    "accession",
    "species",
    "determinant",
    "species_n",
    "carrier_n",
    "prevalence_pct",
    "n_loci",
    "Contig id",
    "Start",
    "Stop",
    "Strand",
    "Element name",
    "Class",
    "Subclass",
    "Method",
    "% Coverage of reference",
    "% Identity to reference",
    "Closest reference accession",
    "Closest reference name"
]

loci = loci[keep].sort_values(
    ["species", "determinant", "accession", "Contig id", "Start"]
)

loci.to_csv(
    OUT_FILE,
    sep="\t",
    index=False
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

n_genomes = loci["accession"].nunique()

n_cases = (
    loci[["accession", "determinant"]]
    .drop_duplicates()
    .shape[0]
)

n_loci = len(loci)

multi = locus_counts[
    locus_counts["n_loci"] > 1
]

log = [
    "BACILLUS AMR PROJECT",
    "STEP 15 - MAP SPORADIC AMR LOCI",
    "",
    f"Unique carrier genomes: {n_genomes}",
    f"Genome-determinant cases: {n_cases}",
    f"Physical AMR loci: {n_loci}",
    f"Cases with >1 physical locus: {len(multi)}",
    ""
]

if len(multi):
    log.append("MULTI-LOCUS CASES")
    log.append(multi.to_string(index=False))
else:
    log.append("No multi-locus cases detected.")

LOG_FILE.write_text("\n".join(log) + "\n")

print("\n".join(log))
