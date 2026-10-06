#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

PROJECT = Path.home() / "bacillus_amr"

META = PROJECT / "08_genomic_context/sporadic_amr_context_windows.tsv"
GFF_DIR = PROJECT / "08_genomic_context/prodigal_meta"

OUT = PROJECT / "08_genomic_context/amr_cds_overlap_validation.tsv"
LOG = PROJECT / "logs/19_amr_cds_overlap_validation.txt"

meta = pd.read_csv(META, sep="\t")

results = []

for _, row in meta.iterrows():

    window_file = row["window_file"]
    base = Path(window_file).stem
    gff = GFF_DIR / f"{base}.gff"

    # Coordenadas AMR dentro de la ventana
    amr_start = int(row["amr_start"]) - int(row["window_start"]) + 1
    amr_stop = int(row["amr_stop"]) - int(row["window_start"]) + 1

    cds_list = []

    with open(gff) as fh:
        for line in fh:

            if line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")

            if len(parts) < 9 or parts[2] != "CDS":
                continue

            cds_start = int(parts[3])
            cds_stop = int(parts[4])
            strand = parts[6]
            attributes = parts[8]

            # Intersección en bp
            overlap = max(
                0,
                min(amr_stop, cds_stop) -
                max(amr_start, cds_start) + 1
            )

            if overlap > 0:

                amr_length = amr_stop - amr_start + 1
                cds_length = cds_stop - cds_start + 1

                cds_list.append({
                    "cds_start": cds_start,
                    "cds_stop": cds_stop,
                    "cds_strand": strand,
                    "attributes": attributes,
                    "overlap_bp": overlap,
                    "amr_coverage_pct": 100 * overlap / amr_length,
                    "cds_coverage_pct": 100 * overlap / cds_length
                })

    if cds_list:

        # CDS con mayor solapamiento
        best = max(cds_list, key=lambda x: x["overlap_bp"])

        results.append({
            "locus_id": row["locus_id"],
            "accession": row["accession"],
            "species": row["species"],
            "determinant": row["determinant"],
            "window_file": window_file,
            "amr_window_start": amr_start,
            "amr_window_stop": amr_stop,
            "amr_strand": row["strand"],
            "cds_start": best["cds_start"],
            "cds_stop": best["cds_stop"],
            "cds_strand": best["cds_strand"],
            "overlap_bp": best["overlap_bp"],
            "amr_coverage_pct": best["amr_coverage_pct"],
            "cds_coverage_pct": best["cds_coverage_pct"],
            "strand_match": row["strand"] == best["cds_strand"],
            "n_overlapping_cds": len(cds_list),
            "cds_attributes": best["attributes"]
        })

    else:

        results.append({
            "locus_id": row["locus_id"],
            "accession": row["accession"],
            "species": row["species"],
            "determinant": row["determinant"],
            "window_file": window_file,
            "amr_window_start": amr_start,
            "amr_window_stop": amr_stop,
            "amr_strand": row["strand"],
            "cds_start": pd.NA,
            "cds_stop": pd.NA,
            "cds_strand": pd.NA,
            "overlap_bp": 0,
            "amr_coverage_pct": 0,
            "cds_coverage_pct": 0,
            "strand_match": False,
            "n_overlapping_cds": 0,
            "cds_attributes": ""
        })

out = pd.DataFrame(results)

out.to_csv(OUT, sep="\t", index=False)

with open(LOG, "w") as fh:

    fh.write("BACILLUS AMR PROJECT\n")
    fh.write("STEP 19 - AMR/CDS OVERLAP VALIDATION\n\n")

    fh.write(f"Physical AMR loci: {len(out)}\n")
    fh.write(
        f"Loci overlapping >=1 predicted CDS: "
        f"{(out['n_overlapping_cds'] > 0).sum()}\n"
    )
    fh.write(
        f"Loci without predicted CDS overlap: "
        f"{(out['n_overlapping_cds'] == 0).sum()}\n"
    )
    fh.write(
        f"Strand matches: "
        f"{out['strand_match'].sum()}\n"
    )

    fh.write("\nAMR coverage by predicted CDS (%):\n")
    fh.write(out["amr_coverage_pct"].describe().to_string())
    fh.write("\n")

print(out[
    [
        "determinant",
        "amr_coverage_pct",
        "cds_coverage_pct",
        "strand_match",
        "n_overlapping_cds"
    ]
].to_string(index=False))

print()
print(f"Output: {OUT}")
print(f"Log:    {LOG}")
