#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# ============================================================
# Files
# ============================================================

INPUT = Path(
    "05_amr/amr_genome_summary.tsv"
)

OUT_PDF = Path(
    "10_figures/Figure1_AMR_burden_by_species.pdf"
)

OUT_PNG = Path(
    "10_figures/Figure1_AMR_burden_by_species.png"
)

LOG = Path(
    "logs/32_figure1_amr_burden.txt"
)


# ============================================================
# Load data
# ============================================================

df = pd.read_csv(
    INPUT,
    sep="\t"
)

df = df[
    df["Taxonomic_resolution"] == "species"
].copy()


# ============================================================
# Species with n >= 10
# ============================================================

counts = (
    df.groupby("Species_normalized")
    .size()
)

eligible = counts[
    counts >= 10
].index

df = df[
    df["Species_normalized"].isin(eligible)
].copy()


# ============================================================
# Order species by mean AMR burden
# ============================================================

means = (
    df.groupby("Species_normalized")
    ["AMR_determinant_count"]
    .mean()
    .sort_values(
        ascending=False
    )
)

order = means.index.tolist()


# ============================================================
# Labels: scientific names in italics
# ============================================================

label_map = {}

for species in order:

    n = counts.loc[species]

    epithet = species.replace(
        "Bacillus ",
        ""
    )

    label_map[species] = (
        rf"$\it{{B.\ {epithet}}}$"
        + f"\n(n={n})"
    )

# ============================================================
# Plot
# ============================================================

sns.set_theme(
    context="paper",
    style="ticks"
)

fig, ax = plt.subplots(
    figsize=(10.5, 6.5)
)


sns.boxplot(
    data=df,
    x="Species_normalized",
    y="AMR_determinant_count",
    order=order,
    width=0.58,
    showfliers=False,
    ax=ax
)


sns.stripplot(
    data=df,
    x="Species_normalized",
    y="AMR_determinant_count",
    order=order,
    jitter=0.22,
    size=1.8,
    alpha=0.35,
    ax=ax
)


# ============================================================
# Axis formatting
# ============================================================

ax.set_xlabel("")

ax.set_ylabel(
    "AMR determinants per genome",
    fontsize=11
)

ax.set_xticks(
    range(len(order))
)

ax.set_xticklabels(
    [
        label_map[s]
        for s in order
    ],
    rotation=55,
    ha="right",
    fontsize=8
)

ax.tick_params(
    axis="y",
    labelsize=9
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.grid(False)

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
        "FIGURE 1 - AMR BURDEN BY SPECIES\n\n"
    )

    f.write(
        f"Genomes plotted: {len(df)}\n"
    )

    f.write(
        f"Species plotted: {len(order)}\n"
    )

    f.write(
        "Inclusion criterion: species-level "
        "assignment with n >= 10 genomes\n\n"
    )

    f.write(
        "Species order (descending mean burden):\n"
    )

    for species in order:

        sub = df[
            df["Species_normalized"] == species
        ]

        f.write(
            f"{species}\t"
            f"n={len(sub)}\t"
            f"mean="
            f"{sub['AMR_determinant_count'].mean():.4f}\t"
            f"median="
            f"{sub['AMR_determinant_count'].median():.4f}\n"
        )

    f.write(
        "\nStatistics associated with figure:\n"
    )

    f.write(
        "Kruskal-Wallis H = 876.402007\n"
    )

    f.write(
        "p = 3.953134e-177\n"
    )

    f.write(
        "Dunn post hoc with Benjamini-Hochberg "
        "correction: 81/120 comparisons q < 0.05\n"
    )


print(f"Genomes plotted: {len(df)}")
print(f"Species plotted: {len(order)}")

print("\nSpecies order:")
for species in order:
    print(
        species,
        f"mean={means.loc[species]:.4f}",
        f"n={counts.loc[species]}"
    )

print(f"\nPDF: {OUT_PDF}")
print(f"PNG: {OUT_PNG}")
print(f"Log: {LOG}")
