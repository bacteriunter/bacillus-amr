#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# ============================================================
# Files
# ============================================================

INPUT = Path(
    "09_analysis/species_resistome_prevalence_matrix.tsv"
)

OUT_PDF = Path(
    "10_figures/Figure2_resistome_prevalence_heatmap.pdf"
)

OUT_PNG = Path(
    "10_figures/Figure2_resistome_prevalence_heatmap.png"
)

LOG = Path(
    "logs/33_figure2_resistome_heatmap.txt"
)


# ============================================================
# Load prevalence matrix
# ============================================================

df = pd.read_csv(
    INPUT,
    sep="\t"
)

species_col = "Species_normalized"

matrix = df.set_index(species_col)

if matrix.shape != (16, 33):
    raise RuntimeError(
        f"Expected 16 x 33 matrix, found {matrix.shape}"
    )


# ============================================================
# Hierarchical clustering
#
# Descriptive ordering only:
# Euclidean distance + average linkage
# ============================================================

species_dist = pdist(
    matrix.values,
    metric="euclidean"
)

species_linkage = linkage(
    species_dist,
    method="average"
)

species_order_idx = leaves_list(
    species_linkage
)

determinant_dist = pdist(
    matrix.T.values,
    metric="euclidean"
)

determinant_linkage = linkage(
    determinant_dist,
    method="average"
)

determinant_order_idx = leaves_list(
    determinant_linkage
)


species_order = (
    matrix.index[
        species_order_idx
    ].tolist()
)

determinant_order = (
    matrix.columns[
        determinant_order_idx
    ].tolist()
)

ordered = matrix.loc[
    species_order,
    determinant_order
]


# ============================================================
# Species labels in italics
# ============================================================

species_labels = []

for species in species_order:

    epithet = species.replace(
        "Bacillus ",
        ""
    )

    species_labels.append(
        rf"$\it{{B.\ {epithet}}}$"
    )


# ============================================================
# Plot
# ============================================================

sns.set_theme(
    context="paper",
    style="white"
)

fig, ax = plt.subplots(
    figsize=(14.5, 7.5)
)

hm = sns.heatmap(
    ordered,
    cmap="viridis",
    vmin=0,
    vmax=100,
    linewidths=0.25,
    linecolor="white",
    cbar_kws={
        "label": "Prevalence (%)",
        "shrink": 0.75
    },
    ax=ax
)


# ============================================================
# Axis formatting
# ============================================================

ax.set_xlabel(
    "AMR determinant",
    fontsize=11
)

ax.set_ylabel("")

ax.set_xticklabels(
    determinant_order,
    rotation=60,
    ha="right",
    fontsize=8
)

ax.set_yticklabels(
    species_labels,
    rotation=0,
    fontsize=9
)

ax.tick_params(
    axis="both",
    length=0
)

fig.tight_layout()


# ============================================================
# Save
# ============================================================

fig.savefig(
    OUT_PDF,
    bbox_inches="tight"
)

fig.savefig(
    OUT_PNG,
    dpi=600,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# Log
# ============================================================

with open(LOG, "w") as f:

    f.write("BACILLUS AMR PROJECT\n")
    f.write(
        "FIGURE 2 - SPECIES RESISTOME "
        "PREVALENCE HEATMAP\n\n"
    )

    f.write(
        f"Species: {matrix.shape[0]}\n"
    )

    f.write(
        f"AMR determinants: {matrix.shape[1]}\n"
    )

    f.write(
        "Values: within-species prevalence (%)\n"
    )

    f.write(
        "Species ordering: hierarchical clustering, "
        "Euclidean distance, average linkage\n"
    )

    f.write(
        "Determinant ordering: hierarchical clustering, "
        "Euclidean distance, average linkage\n"
    )

    f.write(
        "Clustering purpose: descriptive visualization only\n\n"
    )

    f.write("SPECIES ORDER\n")

    for species in species_order:
        f.write(
            f"{species}\n"
        )

    f.write(
        "\nDETERMINANT ORDER\n"
    )

    for determinant in determinant_order:
        f.write(
            f"{determinant}\n"
        )


print(
    f"Matrix plotted: "
    f"{matrix.shape[0]} species x "
    f"{matrix.shape[1]} determinants"
)

print("\nSpecies order:")
for x in species_order:
    print(x)

print("\nDeterminant order:")
print(", ".join(determinant_order))

print(f"\nPDF: {OUT_PDF}")
print(f"PNG: {OUT_PNG}")
print(f"Log: {LOG}")
