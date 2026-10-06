#!/usr/bin/env python3

import re
from pathlib import Path
import pandas as pd

WINDOWS = Path("08_genomic_context/sporadic_amr_context_windows.tsv")
IS_HITS = Path("08_genomic_context/context_IS_hits_accepted.tsv")
GFF_DIR = Path("08_genomic_context/prodigal_meta")

OUT = Path("08_genomic_context/amr_IS_proximity.tsv")
SUMMARY = Path("08_genomic_context/amr_IS_proximity_summary.tsv")
LOG = Path("logs/21_map_IS_proximity_to_AMR.txt")


def parse_gff(gff_file):
    records = []

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) != 9 or parts[2] != "CDS":
                continue

            seqid = parts[0]
            start = int(parts[3])
            stop = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            match = re.search(r"ID=1_(\d+)", attributes)

            if not match:
                continue

            cds_number = int(match.group(1))

            records.append({
                "seqid": seqid,
                "cds_number": cds_number,
                "cds_start": start,
                "cds_stop": stop,
                "cds_strand": strand
            })

    return pd.DataFrame(records)


windows = pd.read_csv(WINDOWS, sep="\t")
hits = pd.read_csv(IS_HITS, sep="\t")

# Map accepted IS protein IDs to locus/window ID.
#
# Protein IDs have structure:
# L002_GCF_002566845.1_abaF_15
#
# The final integer is the Prodigal CDS number.
hits["cds_number"] = (
    hits["protein_id"]
    .str.extract(r"_(\d+)$")[0]
    .astype(int)
)

hits["locus_prefix"] = hits["protein_id"].str.extract(r"^(L\d+)")[0]

rows = []

for _, w in windows.iterrows():

    locus_id = w["locus_id"]
    locus_prefix = locus_id.split("_")[0]

    # Convert original AMR coordinates to coordinates within
    # the extracted local window.
    amr_start_local = int(w["amr_start"]) - int(w["window_start"]) + 1
    amr_stop_local = int(w["amr_stop"]) - int(w["window_start"]) + 1

    gff_files = list(GFF_DIR.glob(f"{locus_prefix}_*.gff"))

    if len(gff_files) != 1:
        raise RuntimeError(
            f"{locus_id}: expected exactly one GFF, found {len(gff_files)}"
        )

    cds = parse_gff(gff_files[0])

    locus_hits = hits[hits["locus_prefix"] == locus_prefix]

    # Keep AMR loci with no accepted IS hit.
    if locus_hits.empty:

        rows.append({
            "locus_id": locus_id,
            "assembly_accession": w["accession"],
            "determinant": w["determinant"],
            "context_status": w["context_status"],
            "amr_start_local": amr_start_local,
            "amr_stop_local": amr_stop_local,
            "IS_detected": False,
            "IS_protein_id": pd.NA,
            "IS_reference": pd.NA,
            "IS_family_annotation": pd.NA,
            "IS_start_local": pd.NA,
            "IS_stop_local": pd.NA,
            "IS_strand": pd.NA,
            "distance_bp": pd.NA,
            "relative_position": pd.NA,
            "pident": pd.NA,
            "query_coverage_pct": pd.NA,
            "evalue": pd.NA
        })

        continue

    for _, h in locus_hits.iterrows():

        match = cds[cds["cds_number"] == h["cds_number"]]

        if len(match) != 1:
            raise RuntimeError(
                f"{h['protein_id']}: expected one CDS in GFF, found {len(match)}"
            )

        c = match.iloc[0]

        is_start = int(c["cds_start"])
        is_stop = int(c["cds_stop"])

        # Minimum edge-to-edge distance.
        # 0 means the intervals overlap.
        if is_stop < amr_start_local:
            distance = amr_start_local - is_stop - 1
            relative_position = "upstream"

        elif is_start > amr_stop_local:
            distance = is_start - amr_stop_local - 1
            relative_position = "downstream"

        else:
            distance = 0
            relative_position = "overlap"

        rows.append({
            "locus_id": locus_id,
            "assembly_accession": w["accession"],
            "determinant": w["determinant"],
            "context_status": w["context_status"],
            "amr_start_local": amr_start_local,
            "amr_stop_local": amr_stop_local,
            "IS_detected": True,
            "IS_protein_id": h["protein_id"],
            "IS_reference": h["IS_reference"],
            "IS_family_annotation": h["reference_annotation"],
            "IS_start_local": is_start,
            "IS_stop_local": is_stop,
            "IS_strand": c["cds_strand"],
            "distance_bp": distance,
            "relative_position": relative_position,
            "pident": h["pident"],
            "query_coverage_pct": h["query_coverage_pct"],
            "evalue": h["evalue"]
        })


result = pd.DataFrame(rows)

result.to_csv(OUT, sep="\t", index=False)


# One-row-per-AMR-locus summary
summary_rows = []

for locus_id, group in result.groupby("locus_id", sort=False):

    first = group.iloc[0]
    detected = group[group["IS_detected"] == True]

    if detected.empty:
        n_is = 0
        min_distance = pd.NA
        nearest_is = pd.NA
    else:
        n_is = len(detected)
        nearest = detected.loc[detected["distance_bp"].astype(float).idxmin()]
        min_distance = int(nearest["distance_bp"])
        nearest_is = nearest["IS_reference"]

    summary_rows.append({
        "locus_id": locus_id,
        "assembly_accession": first["assembly_accession"],
        "determinant": first["determinant"],
        "context_status": first["context_status"],
        "n_IS_proteins": n_is,
        "nearest_IS": nearest_is,
        "min_distance_bp": min_distance
    })

summary = pd.DataFrame(summary_rows)

summary.to_csv(SUMMARY, sep="\t", index=False)


n_loci = summary["locus_id"].nunique()
n_with_is = (summary["n_IS_proteins"] > 0).sum()
n_without_is = (summary["n_IS_proteins"] == 0).sum()

distances = pd.to_numeric(
    summary["min_distance_bp"],
    errors="coerce"
).dropna()


with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write("STEP 21 - IS PROXIMITY TO SPORADIC AMR LOCI\n\n")

    f.write(f"Physical AMR loci: {n_loci}\n")
    f.write(f"Accepted IS-associated proteins: {result['IS_detected'].sum()}\n")
    f.write(f"AMR loci with >=1 accepted IS protein: {n_with_is}\n")
    f.write(f"AMR loci without accepted IS protein: {n_without_is}\n")

    if len(distances):
        f.write("\nNearest IS distance to AMR locus (bp):\n")
        f.write(distances.describe().to_string())
        f.write("\n")


print(summary.to_string(index=False))

print(f"\nPhysical AMR loci: {n_loci}")
print(f"Loci with >=1 accepted IS protein: {n_with_is}")
print(f"Loci without accepted IS protein: {n_without_is}")

print(f"\nOutput: {OUT}")
print(f"Summary: {SUMMARY}")
print(f"Log: {LOG}")
