#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

WINDOWS = Path(
    "08_genomic_context/sporadic_amr_context_windows.tsv"
)
MOBILITY = Path(
    "08_genomic_context/mobility_evidence_summary.tsv"
)
PLASMID = Path(
    "08_genomic_context/amr_locus_plasmid_annotation.tsv"
)

OUT = Path(
    "08_genomic_context/sporadic_amr_locus_evidence.tsv"
)
LOG = Path(
    "logs/24_integrate_locus_evidence.txt"
)

# ============================================================
# 1. Load tables
# ============================================================

windows = pd.read_csv(WINDOWS, sep="\t")
mobility = pd.read_csv(MOBILITY, sep="\t")
plasmid = pd.read_csv(PLASMID, sep="\t")

# ============================================================
# 2. Create common locus number L001-L041
# ============================================================

windows["locus_number"] = (
    windows["locus_id"]
    .str.extract(r"^(L\d{3})", expand=False)
)

mobility["locus_number"] = (
    mobility["locus_id"]
    .str.extract(r"^(L\d{3})", expand=False)
)

# plasmid table already contains locus_number

# ============================================================
# 3. Validate uniqueness
# ============================================================

for name, df in [
    ("windows", windows),
    ("mobility", mobility),
    ("plasmid", plasmid)
]:
    if len(df) != 41:
        raise RuntimeError(
            f"{name}: expected 41 rows, found {len(df)}"
        )

    if df["locus_number"].duplicated().any():
        raise RuntimeError(
            f"{name}: duplicated locus_number detected"
        )

# ============================================================
# 4. Select relevant columns
# ============================================================

base = windows[
    [
        "locus_number",
        "locus_id",
        "accession",
        "species",
        "determinant",
        "contig_id",
        "contig_length_bp",
        "amr_start",
        "amr_stop",
        "strand",
        "left_recovered_bp",
        "right_recovered_bp",
        "window_length_bp",
        "context_status"
    ]
].copy()

mob = mobility[
    [
        "locus_number",
        "n_mobility_proteins",
        "n_IS_transposases",
        "n_recombinase_resolvase",
        "nearest_mobility_distance_bp",
        "mobility_evidence_detected"
    ]
].copy()

pla = plasmid[
    [
        "locus_number",
        "ncbi_explicit_plasmid",
        "ncbi_header"
    ]
].copy()

# ============================================================
# 5. Merge
# ============================================================

final = (
    base
    .merge(
        mob,
        on="locus_number",
        how="left",
        validate="one_to_one"
    )
    .merge(
        pla,
        on="locus_number",
        how="left",
        validate="one_to_one"
    )
)

if len(final) != 41:
    raise RuntimeError(
        f"Final table expected 41 loci, found {len(final)}"
    )

# ============================================================
# 6. Add conservative evidence descriptor
# ============================================================

def classify(row):

    if row["mobility_evidence_detected"]:
        if (
            row["n_IS_transposases"] > 0
            and row["n_recombinase_resolvase"] > 0
        ):
            return "IS + recombinase/resolvase"

        if row["n_IS_transposases"] > 0:
            return "IS-associated protein"

        if row["n_recombinase_resolvase"] > 0:
            return "recombinase/resolvase"

    return "no local mobility protein detected"


final["local_mobility_evidence"] = final.apply(
    classify,
    axis=1
)

# ============================================================
# 7. Save
# ============================================================

final.to_csv(
    OUT,
    sep="\t",
    index=False
)

# ============================================================
# 8. Summary
# ============================================================

n_total = len(final)

n_mobility = int(
    final["mobility_evidence_detected"].sum()
)

n_plasmid = int(
    final["ncbi_explicit_plasmid"].sum()
)

n_complete = int(
    (final["context_status"] == "complete").sum()
)

n_truncated = n_total - n_complete

combined = int(
    (
        (final["n_IS_transposases"] > 0)
        &
        (final["n_recombinase_resolvase"] > 0)
    ).sum()
)

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 24 - INTEGRATED SPORADIC AMR "
        "LOCUS EVIDENCE\n\n"
    )

    f.write(f"Physical AMR loci: {n_total}\n")
    f.write(f"Complete ±10 kb contexts: {n_complete}\n")
    f.write(f"Truncated contexts: {n_truncated}\n")

    f.write(
        f"Loci with local mobility evidence: "
        f"{n_mobility}\n"
    )

    f.write(
        f"Loci with both IS and "
        f"recombinase/resolvase evidence: "
        f"{combined}\n"
    )

    f.write(
        f"Loci explicitly annotated as plasmid "
        f"by NCBI: {n_plasmid}\n"
    )


print(
    final[
        [
            "locus_number",
            "species",
            "determinant",
            "context_status",
            "n_IS_transposases",
            "n_recombinase_resolvase",
            "nearest_mobility_distance_bp",
            "local_mobility_evidence",
            "ncbi_explicit_plasmid"
        ]
    ].to_string(index=False)
)

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
