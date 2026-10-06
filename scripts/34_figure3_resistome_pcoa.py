#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# Files
# ============================================================

INPUT = Path(
    "09_analysis/resistome_pcoa_coordinates.tsv"
)

OUT_PDF = Path(
    "10_figures/Figure3_resistome_PCoA.pdf"
)

OUT_PNG = Path(
    "10_figures/Figure3_resistome_PCoA.png"
)

LOG = Path(
    "logs/34_figure3_resistome_pcoa.txt"
)


# ============================================================
# Load only required columns
# ============================================================

df = pd.read_csv(
    INPUT,
    sep="\t",
    usecols=[
        "Assembly_Accession",
        "PC1",
        "PC2",
        "Species_normalized"
    ]
)

if len(df) != 1082:
    raise RuntimeError(
        f"Expected 1082 genomes, found {len(df)}"
    )


# ============================================================
# Species order
# ============================================================

species_order = (
    df["Species_normalized"]
    .value_counts()
    .index
    .tolist()
)


# ============================================================
# Labels
# ============================================================

label_map = {}

for species in species_order:

    epithet = species.replace(
        "Bacillus ",
        ""
    )

    n = (
        df["Species_normalized"] == species
    ).sum()

    label_map[species] = (
        rf"$\it{{B.\ {epithet}}}$"
        + f" (n={n})"
    )


# ============================================================
# Plot
# ============================================================

sns.set_theme(
    context="paper",
    style="ticks"
)

palette = sns.color_palette(
    "tab20",
    n_colors=len(species_order)
)

fig, ax = plt.subplots(
    figsize=(10.5, 7.5)
)


for species, color in zip(
    species_order,
    palette
):

    sub = df[
        df["Species_normalized"] == species
    ]

    # Individual genomes
    ax.scatter(
        sub["PC1"],
        sub["PC2"],
        s=22,
        alpha=0.60,
        edgecolors="none",
        color=color,
        label=label_map[species]
    )

    # Species centroid
    centroid_x = sub["PC1"].mean()
    centroid_y = sub["PC2"].mean()

    ax.scatter(
        centroid_x,
        centroid_y,
        s=115,
        marker="o",
        color=color,
        edgecolors="black",
        linewidths=1.2,
        zorder=10
    )


# ============================================================
# Axis formatting
# ============================================================

ax.set_xlabel(
    "PCoA axis 1 (26.02%)",
    fontsize=11
)

ax.set_ylabel(
    "PCoA axis 2 (19.11%)",
    fontsize=11
)

ax.axhline(
    0,
    linewidth=0.5,
    color="0.8",
    zorder=0
)

ax.axvline(
    0,
    linewidth=0.5,
    color="0.8",
    zorder=0
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.tick_params(
    labelsize=9
)


# ============================================================
# Legend
# ============================================================

legend = ax.legend(
    title="Species",
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    frameon=False,
    fontsize=8,
    title_fontsize=9,
    markerscale=1.2
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
        "FIGURE 3 - PCoA OF RESISTOME COMPOSITION\n\n"
    )

    f.write(
        f"Genomes plotted: {len(df)}\n"
    )

    f.write(
        f"Species plotted: {len(species_order)}\n"
    )

    f.write(
        "Distance used in original ordination: Jaccard\n"
    )

    f.write(
        "PCoA axis 1: 26.0206%\n"
    )

    f.write(
        "PCoA axis 2: 19.1088%\n"
    )

    f.write(
        "PCoA axes 1+2: 45.1294%\n\n"
    )

    f.write(
        "Species centroids displayed as enlarged markers "
        "with black outlines.\n"
    )

    f.write(
        "Centroids are shown for visualization only.\n\n"
    )

    f.write(
        "PERMANOVA: pseudo-F = 200.311429; "
        "R2 = 0.738127; p = 0.001; "
        "999 permutations\n"
    )

    f.write(
        "PERMDISP: F = 28.473100; "
        "p = 0.001; "
        "999 permutations\n\n"
    )

    f.write("SPECIES\n")

    for species in species_order:

        n = (
            df["Species_normalized"] == species
        ).sum()

        f.write(
            f"{species}\t{n}\n"
        )


print(f"Genomes plotted: {len(df)}")
print(f"Species plotted: {len(species_order)}")

print(f"\nPDF: {OUT_PDF}")
print(f"PNG: {OUT_PNG}")
print(f"Log: {LOG}")
