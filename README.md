# Genomic architecture of antimicrobial resistance in plant-associated *Bacillus*

Reproducible computational workflow and supporting data for the study:

**Genomic architecture of antimicrobial resistance in plant-associated Bacillus: conserved resistomes and signatures of potential horizontal acquisition**

## Overview

This study investigates the genomic architecture of antimicrobial resistance (AMR) determinants across plant-associated *Bacillus* lineages.

The central question is:

> How are antimicrobial resistance determinants distributed across plant-associated *Bacillus* lineages, and which genomic features distinguish conserved lineage-associated determinants from those showing signatures consistent with potential horizontal acquisition?

Plant association is used exclusively as an **inclusion criterion**. Repository-derived environmental categories are not treated as explanatory ecological variables.

The analytical framework integrates genome selection, structural quality control, AMR detection, species-level resistome analyses, targeted phylogenomic reconstruction, genomic-context analysis, and local mobility evidence.

## Study design

The workflow starts from RefSeq assemblies assigned to *Bacillus* and identifies genomes with direct evidence of plant association.

Accepted evidence includes direct association with plant tissues, surfaces, endosphere, rhizoplane, rhizosphere, roots, leaves, seeds, or equivalent explicitly plant-associated material.

Generic soil, agricultural soil, crop fields, compost, plant residues, foods, fermented products, and processed plant materials are not considered sufficient evidence by themselves.

Plant-association decisions and curated inclusions/exclusions are retained as auditable metadata files.

## Dataset

The principal dataset contains:

- 14,111 RefSeq *Bacillus* assemblies screened
- 2,550 plant-associated candidates
- 1,374 assemblies selected after metadata classification
- 1,312 genomes retained after structural quality control
- 1,159 genomes with species-level assignments
- 153 genomes retained as *Bacillus* sp.
- 44 named species
- 1,181 genomes containing at least one AMR determinant
- 3,913 AMR hits
- 36 unique AMR determinants

Structural genome QC excluded assemblies with:

- more than 500 contigs, or
- N50 < 20,000 bp

No genome-size or GC-content filter was applied.

## AMR detection

AMR determinants were identified with:

- AMRFinderPlus 4.2.7
- AMRFinderPlus database 2026-08-07.1

AMRFinderPlus was run with the `--plus` option. Complete output was retained during the analysis, but the primary resistome analyses use records classified as `AMR`. `VIRULENCE` and `STRESS` records are not combined with AMR determinants.

## Species-level analyses

Quantitative within-species analyses were restricted to named species represented by at least 10 genomes, resulting in 16 species and 1,082 genomes.

AMR burden was compared among species using the Kruskal-Wallis test followed by Dunn pairwise comparisons with Benjamini-Hochberg correction.

Resistome composition was evaluated from binary AMR presence/absence profiles using Jaccard distances. The analysis included PERMANOVA, principal coordinates analysis (PCoA), PERMDISP, and within-species pairwise Jaccard heterogeneity.

The significant PERMDISP result was retained as a biological feature of lineage-specific resistome heterogeneity and considered when interpreting PERMANOVA.

## Conserved and sporadic AMR determinants

AMR prevalence was treated as a continuous species-level property. Highly prevalent determinants were interpreted as conserved lineage-associated determinants rather than direct evidence of intrinsic resistance.

For targeted analysis of discontinuous AMR distributions, sporadic species-determinant combinations were operationally defined as determinants present in more than 0% and no more than 5% of genomes within species represented by at least 10 genomes.

This identified:

- 21 species-determinant combinations
- 14 unique sporadic determinants
- 39 carrier genomes
- 40 genome-determinant observations
- 41 physical AMR loci

## Phylogenomic context

Phylogenomic reconstruction used GToTree 1.8.16 with the Firmicutes single-copy gene HMM set containing 119 profiles.

The full phylogenomic alignment retained 1,311 of the 1,312 QC-passed genomes. One genome was excluded from phylogenetic reconstruction because only two unique single-copy genes were recovered; this genome remained in the global AMR analyses.

For targeted phylogenetic context, a deterministic subset was constructed containing all 39 sporadic AMR carrier genomes plus three non-carrier intraspecific context genomes for each of the 16 eligible species. The resulting subset contained 87 genomes and 23,930 aligned amino-acid positions.

ModelFinder selected `Q.PLANT+F+I+G4` according to AIC, AICc, and BIC.

The final maximum-likelihood phylogeny was inferred with IQ-TREE 3.1.2 using 1,000 ultrafast bootstrap replicates and 1,000 SH-aLRT replicates.

## Genomic context and mobility evidence

A ±10 kb genomic region was extracted around each sporadic AMR locus when assembly boundaries permitted.

Protein-coding sequences were predicted uniformly with Prodigal 2.6.3 using metagenomic mode (`-p meta`) and translation table 11. AMR coordinates were validated against predicted CDS coordinates before mobility analysis.

Insertion-sequence-associated proteins were identified by BLASTP against the bacterial insertion-sequence database distributed with the Prokka database installation. Accepted matches required:

- amino-acid identity ≥ 40%
- query coverage ≥ 70%
- E-value ≤ 1e-10

Additional mobility-associated proteins were evaluated against the bacterial Swiss-Prot database distributed with the same database installation. Explicit plasmid annotation in the corresponding NCBI sequence records was also evaluated.

Evidence compatible with potential horizontal acquisition was interpreted from convergent features, particularly discontinuous within-lineage occurrence and local mobility-associated genomic context. These signatures were not treated as direct proof of horizontal gene transfer.

## Main biological interpretation

The study evaluates a hierarchical model of AMR architecture in plant-associated *Bacillus*: a dominant resistome component structured by lineage and conserved within species, together with a smaller discontinuously distributed component in which some determinants additionally show local genomic signatures compatible with potential horizontal acquisition.

## Reproducible workflow

The complete workflow is defined in `Snakefile` and was validated with Snakemake 9.16.3.

The final workflow connects metadata acquisition, genome selection, structural QC, AMR detection, resistome analyses, targeted phylogenomics, genomic-context analysis, mobility evidence, statistical analyses, and generation of the final manuscript figures and tables.

Historical or superseded analyses are intentionally excluded from the final workflow, including the incomplete CheckM2 analysis, preliminary Prokka annotation of short context windows, and exploratory full-dataset IQ-TREE ModelFinder runs.

## Conda environments

Environment exports are provided in `envs/`:

- `bacillus_amr.yml`
- `amrfinder_env.yml`
- `bacillus_phylo.yml`
- `bacillus_ecogenomics.yml`
- `smk.yml`

These files represent the functional Conda environments exported during repository preparation. Some Python packages were subsequently updated; therefore, versions recorded during the statistical analyses are documented separately below.

### Recorded statistical-analysis environment

- Python 3.10.18
- pandas 2.3.2
- NumPy 2.2.6
- SciPy 1.15.2
- scikit-bio 0.7.2
- Biopython 1.85
- scikit-posthocs 0.14.0

### Principal command-line software

- NCBI Datasets CLI 18.38.0
- SeqKit 2.14.0
- AMRFinderPlus 4.2.7
- AMRFinderPlus database 2026-08-07.1
- GToTree 1.8.16
- HMMER 3.4
- MAFFT 7.526
- IQ-TREE 3.1.2
- Prodigal 2.6.3
- Snakemake 9.16.3

## Running the workflow

The workflow manager can be activated with:

    conda activate smk

The workflow can first be inspected without executing analyses:

    snakemake --snakefile Snakefile --dry-run --printshellcmds

The complete workflow can be executed with:

    snakemake --snakefile Snakefile --cores 4

The `Snakefile` invokes the analysis environments required by individual steps through `conda run`.

Some upstream NCBI resources are dynamic. The accession lists and processed metadata retained in this repository define the dataset analyzed in the study.

## Repository structure

    00_metadata/          Metadata, selection, and audit tables
    01_genomes/           Selected genome accession list
    02_qc/                Structural genome QC
    03_taxonomy/          Taxonomic normalization
    04_annotation/        Reserved annotation directory
    05_amr/               AMR summaries and prevalence data
    06_phylogeny/         Phylogenomic analyses
    07_mge/               Reserved MGE directory
    08_genomic_context/   Genomic-context and mobility analyses
    09_analysis/          Statistical and integrated analyses
    10_figures/           Final manuscript figures
    11_tables/            Final manuscript tables
    scripts/              Ordered analysis scripts
    envs/                 Conda environment exports
    logs/                 Analysis logs
    Snakefile             Snakemake workflow

Large downloadable or reconstructable resources are intentionally excluded from version control through `.gitignore`, including downloaded NCBI genome packages, downloaded protein packages, raw per-genome AMRFinderPlus outputs, and selected large intermediate phylogenomic files.

## Main outputs

Four principal manuscript figures are provided:

1. AMR burden across species
2. Species-level AMR prevalence heatmap
3. Resistome composition PCoA
4. Phylogenetic context of sporadic AMR determinants

Three principal manuscript tables are provided:

1. Dataset, genome QC, and AMR summary
2. Species-level AMR burden and within-species resistome heterogeneity
3. Integrated evidence for sporadic AMR determinants

## Reproducibility and traceability

The repository preserves the ordered analytical scripts, metadata-selection criteria, audit files, intermediate analytical tables, final figures, final tables, software information, Conda environment exports, and Snakemake workflow used to document the computational analysis.

Downloaded public sequence data and other large reconstructable resources are not duplicated in version control.

## Citation

Citation information and the permanent DOI will be provided through the archived Zenodo release.
