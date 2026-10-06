#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from scipy.spatial.distance import pdist, squareform
from skbio import DistanceMatrix
from skbio.stats.distance import permdisp


PA_FILE = Path(
    "05_amr/amr_presence_absence.tsv"
)

GENOME_FILE = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT = Path(
    "09_analysis/resistome_permdisp.tsv"
)

LOG = Path(
    "logs/30_resistome_permdisp.txt"
)


# ============================================================
# 1. Load data
# ============================================================

pa = pd.read_csv(PA_FILE, sep="\t")
genomes = pd.read_csv(GENOME_FILE, sep="\t")


# ============================================================
# 2. Select species-level genomes with n >= 10
# ============================================================

species_only = genomes[
    genomes["Taxonomic_resolution"] == "species"
].copy()

counts = (
    species_only
    .groupby("Species_normalized")
    .size()
)

eligible = counts[counts >= 10].index

metadata = species_only[
    species_only["Species_normalized"].isin(eligible)
][
    ["Assembly_Accession", "Species_normalized"]
].copy()


# ============================================================
# 3. Merge with AMR presence/absence
# ============================================================

data = metadata.merge(
    pa,
    on="Assembly_Accession",
    how="inner",
    validate="one_to_one"
)

if len(data) != 1082:
    raise RuntimeError(
        f"Expected 1082 genomes, found {len(data)}"
    )

amr_cols = [
    c for c in pa.columns
    if c != "Assembly_Accession"
]

X = data[amr_cols].to_numpy(dtype=int)

ids = data["Assembly_Accession"].tolist()
groups = data["Species_normalized"].tolist()


# ============================================================
# 4. Jaccard distance
# ============================================================

dist = squareform(
    pdist(X, metric="jaccard")
)

dm = DistanceMatrix(
    dist,
    ids=ids
)


# ============================================================
# 5. PERMDISP
# ============================================================

result = permdisp(
    dm,
    grouping=groups,
    permutations=999,
    test="centroid"
)

result_df = pd.DataFrame({
    "parameter": result.index,
    "value": result.values
})

result_df.to_csv(
    OUT,
    sep="\t",
    index=False
)


# ============================================================
# 6. Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "STEP 30 - MULTIVARIATE DISPERSION "
        "OF RESISTOME COMPOSITION\n\n"
    )

    f.write("Distance: Jaccard\n")
    f.write("Grouping factor: Species_normalized\n")
    f.write("Center: centroid\n")
    f.write("Permutations: 999\n\n")

    for key, value in result.items():
        f.write(f"{key}: {value}\n")


print("\nPERMDISP:")
print(result)

print(f"\nOutput: {OUT}")
print(f"Log:    {LOG}")
