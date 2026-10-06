#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

INDIR = PROJECT / "05_amr" / "amrfinder_raw"
OUT_ALL = PROJECT / "05_amr" / "amrfinder_all_hits.tsv"
OUT_AMR = PROJECT / "05_amr" / "amrfinder_amr_hits.tsv"
OUT_SUMMARY = PROJECT / "05_amr" / "amrfinder_type_summary.tsv"
LOG = PROJECT / "logs" / "11_amrfinder_summary.txt"

files = sorted(INDIR.glob("*.tsv"))

tables = []

for file in files:
    df = pd.read_csv(file, sep="\t", dtype=str)

    # Empty AMRFinder output: header only
    if df.empty:
        continue

    df.insert(0, "Assembly_Accession", file.stem)
    tables.append(df)

if tables:
    all_hits = pd.concat(tables, ignore_index=True)
else:
    raise RuntimeError("No AMRFinder hits were found.")

# Preserve complete AMRFinderPlus output
all_hits.to_csv(
    OUT_ALL,
    sep="\t",
    index=False
)

# AMR only
amr = all_hits[
    all_hits["Type"].eq("AMR")
].copy()

amr.to_csv(
    OUT_AMR,
    sep="\t",
    index=False
)

# Summary by AMRFinder Type
type_summary = (
    all_hits["Type"]
    .value_counts(dropna=False)
    .rename_axis("Type")
    .reset_index(name="Hit_count")
)

type_summary.to_csv(
    OUT_SUMMARY,
    sep="\t",
    index=False
)

total_files = len(files)
genomes_any_hit = all_hits["Assembly_Accession"].nunique()
genomes_amr = amr["Assembly_Accession"].nunique()

with open(LOG, "w") as f:
    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 11 - AMRFinderPlus SUMMARY\n\n")

    f.write(f"AMRFinder result files: {total_files}\n")
    f.write(f"Total hits (--plus): {len(all_hits)}\n")
    f.write(f"Genomes with any hit: {genomes_any_hit}\n")
    f.write(f"AMR hits: {len(amr)}\n")
    f.write(f"Genomes with >=1 AMR hit: {genomes_amr}\n")
    f.write(
        f"Genomes without AMR hits: "
        f"{total_files - genomes_amr}\n\n"
    )

    f.write("Hits by Type:\n")
    f.write(type_summary.to_string(index=False))
    f.write("\n")


print("AMRFinder result files:", total_files)
print("Total hits (--plus):", len(all_hits))
print("Genomes with any hit:", genomes_any_hit)
print("AMR hits:", len(amr))
print("Genomes with >=1 AMR hit:", genomes_amr)
print("Genomes without AMR hits:", total_files - genomes_amr)

print("\nHITS BY TYPE:")
print(type_summary.to_string(index=False))

print("\nTOP 20 AMR ELEMENTS:")
print(
    amr["Element symbol"]
    .value_counts()
    .head(20)
    .to_string()
)
