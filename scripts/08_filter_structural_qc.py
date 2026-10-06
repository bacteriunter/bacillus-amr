#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

PROJECT = Path.home() / "bacillus_amr"

INPUT = PROJECT / "02_qc" / "seqkit_genome_stats.tsv"
OUTPUT = PROJECT / "02_qc" / "structural_qc_classification.tsv"
PASSED = PROJECT / "02_qc" / "qc_passed_accessions.txt"
FAILED = PROJECT / "02_qc" / "qc_failed_accessions.txt"
LOG = PROJECT / "logs" / "08_structural_qc_summary.txt"

MAX_CONTIGS = 500
MIN_N50 = 20000

df = pd.read_csv(INPUT, sep="\t")

df["fail_contigs"] = df["Num_contigs"] > MAX_CONTIGS
df["fail_n50"] = df["N50"] < MIN_N50

df["QC_status"] = "PASS"
df.loc[df["fail_contigs"] | df["fail_n50"], "QC_status"] = "FAIL"

df["QC_reason"] = ""

df.loc[
    df["fail_contigs"] & ~df["fail_n50"],
    "QC_reason"
] = ">500_contigs"

df.loc[
    ~df["fail_contigs"] & df["fail_n50"],
    "QC_reason"
] = "N50<20000"

df.loc[
    df["fail_contigs"] & df["fail_n50"],
    "QC_reason"
] = ">500_contigs;N50<20000"

df.to_csv(OUTPUT, sep="\t", index=False)

passed = df.loc[df["QC_status"] == "PASS", "Assembly_Accession"]
failed = df.loc[df["QC_status"] == "FAIL", "Assembly_Accession"]

passed.to_csv(PASSED, index=False, header=False)
failed.to_csv(FAILED, index=False, header=False)

n_total = len(df)
n_pass = len(passed)
n_fail = len(failed)

only_contigs = (
    df["fail_contigs"] & ~df["fail_n50"]
).sum()

only_n50 = (
    ~df["fail_contigs"] & df["fail_n50"]
).sum()

both = (
    df["fail_contigs"] & df["fail_n50"]
).sum()

with open(LOG, "w") as f:
    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 08 - STRUCTURAL GENOME QC\n\n")

    f.write("Criteria:\n")
    f.write("FAIL if Num_contigs > 500 OR N50 < 20000 bp\n\n")

    f.write(f"Total genomes: {n_total}\n")
    f.write(f"PASS: {n_pass}\n")
    f.write(f"FAIL: {n_fail}\n")
    f.write(f"PASS percentage: {n_pass/n_total*100:.2f}%\n")
    f.write(f"FAIL percentage: {n_fail/n_total*100:.2f}%\n\n")

    f.write(f"Only >500 contigs: {only_contigs}\n")
    f.write(f"Only N50 <20000: {only_n50}\n")
    f.write(f"Both criteria: {both}\n")

print(f"Total genomes: {n_total}")
print(f"QC PASS: {n_pass}")
print(f"QC FAIL: {n_fail}")
print(f"Only >500 contigs: {only_contigs}")
print(f"Only N50 <20000: {only_n50}")
print(f"Both criteria: {both}")
print(f"\nPASS list: {PASSED}")
print(f"FAIL list: {FAILED}")

