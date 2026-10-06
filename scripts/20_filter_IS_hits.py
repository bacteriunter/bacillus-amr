#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

INFILE = Path("08_genomic_context/context_IS_best_hits.tsv")
OUT_ALL = Path("08_genomic_context/context_IS_best_hits_classified.tsv")
OUT_ACCEPTED = Path("08_genomic_context/context_IS_hits_accepted.tsv")
LOG = Path("logs/20_filter_IS_hits.txt")

cols = [
    "protein_id",
    "IS_reference",
    "pident",
    "alignment_length",
    "query_length",
    "subject_length",
    "query_coverage_pct",
    "evalue",
    "bitscore",
    "reference_annotation",
]

df = pd.read_csv(INFILE, sep="\t", names=cols)

# Conservative operational criterion for IS-associated proteins
df["IS_hit_accepted"] = (
    (df["pident"] >= 40.0)
    & (df["query_coverage_pct"] >= 70.0)
    & (df["evalue"] <= 1e-10)
)

df["classification"] = df["IS_hit_accepted"].map(
    {True: "accepted", False: "rejected"}
)

df.to_csv(OUT_ALL, sep="\t", index=False)

accepted = df[df["IS_hit_accepted"]].copy()
accepted.to_csv(OUT_ACCEPTED, sep="\t", index=False)

with open(LOG, "w") as f:
    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 20 - FILTER IS PROTEIN HITS\n\n")
    f.write("Operational acceptance criteria:\n")
    f.write("Identity >= 40%\n")
    f.write("Query coverage >= 70%\n")
    f.write("E-value <= 1e-10\n\n")
    f.write(f"Proteins with BLAST hit: {len(df)}\n")
    f.write(f"Accepted IS-associated proteins: {len(accepted)}\n")
    f.write(f"Rejected weak/partial hits: {(~df['IS_hit_accepted']).sum()}\n")

print(df[
    [
        "protein_id",
        "IS_reference",
        "pident",
        "query_coverage_pct",
        "evalue",
        "classification",
    ]
].to_string(index=False))

print(f"\nAccepted: {len(accepted)}")
print(f"Rejected: {(~df['IS_hit_accepted']).sum()}")
print(f"\nOutput: {OUT_ACCEPTED}")
print(f"Log:    {LOG}")

