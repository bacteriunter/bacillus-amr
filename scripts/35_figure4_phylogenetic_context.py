#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo


# ============================================================
# Files
# ============================================================

TREE_FILE = Path(
    "06_phylogeny/iqtree_subset87/final_tree.treefile"
)

SUBSET_FILE = Path(
    "06_phylogeny/phylogenetic_subset.tsv"
)

EVIDENCE_FILE = Path(
    "09_analysis/sporadic_amr_integrated_evidence.tsv"
)

LOCUS_EVIDENCE_FILE = Path(
    "08_genomic_context/sporadic_amr_locus_evidence.tsv"
)

OUT_PDF = Path(
    "10_figures/Figure4_phylogenetic_context.pdf"
)

OUT_PNG = Path(
    "10_figures/Figure4_phylogenetic_context.png"
)

LOG = Path(
    "logs/35_figure4_phylogenetic_context.txt"
)


# ============================================================
# Load data
# ============================================================

tree = Phylo.read(
    TREE_FILE,
    "newick"
)

subset = pd.read_csv(
    SUBSET_FILE,
    sep="\t"
)

evidence = pd.read_csv(
    EVIDENCE_FILE,
    sep="\t"
)

locus_evidence = pd.read_csv(
    LOCUS_EVIDENCE_FILE,
    sep="\t"
)


# ============================================================
# Validation
# ============================================================

tips = [
    tip.name
    for tip in tree.get_terminals()
]

if len(tips) != 87:
    raise RuntimeError(
        f"Expected 87 tree tips, found {len(tips)}"
    )

subset_accessions = set(
    subset["Assembly_Accession"]
)

missing = set(tips) - subset_accessions

if missing:
    raise RuntimeError(
        f"Tree tips missing from subset metadata: {missing}"
    )

n_expected_carriers = (
    subset["n_sporadic_determinants"] > 0
).sum()

if n_expected_carriers != 39:
    raise RuntimeError(
        f"Expected 39 sporadic carriers, "
        f"found {n_expected_carriers}"
    )


# ============================================================
# Metadata dictionary
# ============================================================

meta = (
    subset
    .set_index("Assembly_Accession")
    .to_dict("index")
)


# ============================================================
# Genome-level local mobility evidence
# ============================================================

mobility_positive = locus_evidence.loc[
    locus_evidence["mobility_evidence_detected"] == True
].copy()

mobility_accessions = set(
    mobility_positive["accession"]
)

if len(mobility_positive) != 10:
    raise RuntimeError(
        f"Expected 10 mobility-positive loci, "
        f"found {len(mobility_positive)}"
    )

if len(mobility_accessions) != 10:
    raise RuntimeError(
        f"Expected 10 mobility-positive genomes, "
        f"found {len(mobility_accessions)}"
    )


# ============================================================
# Determine tree coordinates
# ============================================================

depths = tree.depths()

if not max(depths.values()):
    depths = tree.depths(
        unit_branch_lengths=True
    )

terminals = tree.get_terminals()

y_positions = {
    tip: i
    for i, tip in enumerate(terminals)
}


def assign_internal_y(clade):

    if clade in y_positions:
        return y_positions[clade]

    child_y = [
        assign_internal_y(child)
        for child in clade.clades
    ]

    y_positions[clade] = (
        sum(child_y) / len(child_y)
    )

    return y_positions[clade]


assign_internal_y(tree.root)


# ============================================================
# Figure
# ============================================================

fig, ax = plt.subplots(
    figsize=(17, 17)
)

# ============================================================
# Draw tree
# ============================================================

for clade in tree.find_clades(
    order="preorder"
):

    x = depths[clade]
    y = y_positions[clade]

    if clade.clades:

        child_ys = [
            y_positions[c]
            for c in clade.clades
        ]

        # Vertical connector
        ax.plot(
            [x, x],
            [min(child_ys), max(child_ys)],
            color="0.35",
            linewidth=0.7
        )

        # Horizontal branches
        for child in clade.clades:

            child_x = depths[child]
            child_y = y_positions[child]

            ax.plot(
                [x, child_x],
                [child_y, child_y],
                color="0.35",
                linewidth=0.7
            )


# ============================================================
# Layout coordinates
# ============================================================

tree_max = max(
    depths[t]
    for t in terminals
)

species_x = tree_max + 0.12
status_x = tree_max + 0.27
determinant_x = tree_max + 0.34
mobility_x = tree_max + 0.50


# ============================================================
# Tip annotations
# ============================================================

n_carriers = 0
n_context = 0
n_mobility_marked = 0


for tip in terminals:

    accession = tip.name
    y = y_positions[tip]

    row = meta[accession]

    species = row[
        "Species_normalized"
    ]

    epithet = species.replace(
        "Bacillus ",
        ""
    )

    # --------------------------------------------------------
    # Accession
    # --------------------------------------------------------

    ax.text(
        depths[tip] + 0.003,
        y,
        accession,
        va="center",
        fontsize=5.6
    )

    # --------------------------------------------------------
    # Species
    # --------------------------------------------------------

    ax.text(
        species_x,
        y,
        rf"$\it{{B.\ {epithet}}}$",
        va="center",
        fontsize=6.0
    )

    # --------------------------------------------------------
    # Carrier
    # --------------------------------------------------------

    if row["n_sporadic_determinants"] > 0:

        n_carriers += 1

        ax.scatter(
            status_x,
            y,
            marker="o",
            s=20,
            facecolor="black",
            edgecolor="black",
            zorder=5
        )

        determinants = row[
            "sporadic_determinants"
        ]

        if pd.isna(determinants):
            determinants = ""
        else:
            determinants = str(
                determinants
            )

        ax.text(
            determinant_x,
            y,
            determinants,
            va="center",
            fontsize=5.8
        )

        # ----------------------------------------------------
        # Genome-specific local mobility evidence
        # ----------------------------------------------------

        if accession in mobility_accessions:

            n_mobility_marked += 1

            ax.scatter(
                mobility_x,
                y,
                marker="*",
                s=42,
                facecolor="black",
                edgecolor="black",
                zorder=6
            )

    # --------------------------------------------------------
    # Intraspecific context genome
    # --------------------------------------------------------

    else:

        n_context += 1

        ax.scatter(
            status_x,
            y,
            marker="o",
            s=18,
            facecolor="white",
            edgecolor="0.35",
            linewidth=0.7,
            zorder=5
        )


# ============================================================
# Final validation
# ============================================================

if n_carriers != 39:
    raise RuntimeError(
        f"Expected 39 carrier tips, "
        f"found {n_carriers}"
    )

if n_context != 48:
    raise RuntimeError(
        f"Expected 48 context genomes, "
        f"found {n_context}"
    )

if n_mobility_marked != 10:
    raise RuntimeError(
        f"Expected 10 mobility-positive carrier genomes, "
        f"found {n_mobility_marked}"
    )


# ============================================================
# Column headers
# ============================================================

header_y = len(terminals) + 1.0

ax.text(
    species_x,
    header_y,
    "Species",
    fontsize=8,
    fontweight="bold",
    ha="left"
)

ax.text(
    status_x,
    header_y,
    "Carrier",
    fontsize=8,
    fontweight="bold",
    ha="center"
)

ax.text(
    determinant_x,
    header_y,
    "Sporadic determinant",
    fontsize=8,
    fontweight="bold",
    ha="left"
)

ax.text(
    mobility_x,
    header_y,
    "Local mobility",
    fontsize=8,
    fontweight="bold",
    ha="center"
)


# ============================================================
# Legend
# ============================================================

legend_elements = [

    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=5,
        label="Sporadic AMR carrier"
    ),

    Line2D(
        [0],
        [0],
        marker="o",
        linestyle="None",
        markerfacecolor="white",
        markeredgecolor="0.35",
        markersize=5,
        label="Intraspecific context genome"
    ),

    Line2D(
        [0],
        [0],
        marker="*",
        linestyle="None",
        markerfacecolor="black",
        markeredgecolor="black",
        markersize=7,
        label="Carrier genome with local mobility evidence"
    )
]

ax.legend(
    handles=legend_elements,
    loc="lower left",
    bbox_to_anchor=(0.0, -0.035),
    frameon=False,
    fontsize=7
)


# ============================================================
# Axes
# ============================================================

ax.set_ylim(
    -1,
    len(terminals) + 2
)

ax.set_xlim(
    -0.005,
    mobility_x + 0.08
)

ax.set_xlabel(
    "Branch length (substitutions per site)",
    fontsize=9
)

ax.set_yticks([])

ax.spines["left"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

ax.tick_params(
    axis="x",
    labelsize=8
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

    f.write(
        "BACILLUS AMR PROJECT\n"
    )

    f.write(
        "FIGURE 4 - TARGETED PHYLOGENETIC CONTEXT "
        "OF SPORADIC AMR DETERMINANTS\n\n"
    )

    f.write(
        f"Tree tips: {len(terminals)}\n"
    )

    f.write(
        f"Sporadic carrier genomes: "
        f"{n_carriers}\n"
    )

    f.write(
        f"Intraspecific context genomes: "
        f"{n_context}\n"
    )

    f.write(
        f"Carrier genomes with local mobility evidence: "
        f"{n_mobility_marked}\n\n"
    )

    f.write(
        "Tree source: final ML tree from the targeted "
        "87-genome phylogenetic subset.\n"
    )

    f.write(
        "Alignment: 119 Firmicutes single-copy marker genes, "
        "23,930 amino-acid positions.\n"
    )

    f.write(
        "Model: Q.PLANT+F+I+G4.\n"
    )

    f.write(
        "Branch support analysis: 1,000 ultrafast bootstrap "
        "replicates and 1,000 SH-aLRT replicates.\n\n"
    )

    f.write(
        "Filled circles indicate sporadic AMR carriers; "
        "open circles indicate intraspecific context genomes.\n"
    )

    f.write(
        "Stars indicate carrier genomes with at least one "
        "AMR locus showing local mobility evidence in the "
        "genomic-context analysis.\n"
    )


# ============================================================
# Console summary
# ============================================================

print(
    f"Tree tips: {len(terminals)}"
)

print(
    f"Sporadic carriers: {n_carriers}"
)

print(
    f"Context genomes: {n_context}"
)

print(
    "Carrier genomes with local mobility evidence:",
    n_mobility_marked
)

print(
    f"\nPDF: {OUT_PDF}"
)

print(
    f"PNG: {OUT_PNG}"
)

print(
    f"Log: {LOG}"
)
