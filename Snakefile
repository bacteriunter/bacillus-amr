# ============================================================
# Bacillus AMR project
# Reproducible Snakemake workflow
#
# Snakemake: 9.16.3
#
# Analysis environments used:
#   bacillus_amr
#   amrfinder_env
#   bacillus_phylo
#   bacillus_ecogenomics
#
# Historical/experimental analyses intentionally excluded:
#   - CheckM2 attempt
#   - Prokka short-window annotation (Step 17)
#   - exploratory full-dataset IQ-TREE ModelFinder runs
# ============================================================

import os
from pathlib import Path


# ------------------------------------------------------------
# Conda environments
# ------------------------------------------------------------

ENV_AMR = "bacillus_amr"
ENV_AMRFINDER = "amrfinder_env"
ENV_PHYLO = "bacillus_phylo"
ENV_ECO = "bacillus_ecogenomics"


# ------------------------------------------------------------
# Final manuscript outputs
# ------------------------------------------------------------

FINAL_OUTPUTS = [
    "10_figures/Figure1_AMR_burden_by_species.pdf",
    "10_figures/Figure1_AMR_burden_by_species.png",

    "10_figures/Figure2_resistome_prevalence_heatmap.pdf",
    "10_figures/Figure2_resistome_prevalence_heatmap.png",

    "10_figures/Figure3_resistome_PCoA.pdf",
    "10_figures/Figure3_resistome_PCoA.png",

    "10_figures/Figure4_phylogenetic_context.pdf",
    "10_figures/Figure4_phylogenetic_context.png",

    "11_tables/Table1_dataset_QC_AMR_summary.tsv",
    "11_tables/Table2_species_AMR_burden_heterogeneity.tsv",
    "11_tables/Table3_sporadic_AMR_evidence.tsv"
]


rule all:
    input:
        FINAL_OUTPUTS


# ============================================================
# 00. Software versions
# ============================================================

rule record_versions:
    output:
        "00_metadata/software_versions.txt"
    shell:
        """
        conda run -n {ENV_AMR} \
            bash scripts/00_record_versions.sh
        """


# ============================================================
# 01. Retrieve Bacillus RefSeq metadata
# ============================================================

rule fetch_bacillus_metadata:
    output:
        "00_metadata/bacillus_refseq_summary.jsonl"
    shell:
        """
        conda run -n {ENV_AMR} \
            bash scripts/01_fetch_bacillus_metadata.sh
        """


# ============================================================
# 02. Master metadata table
# ============================================================

rule build_master_metadata:
    input:
        "00_metadata/bacillus_refseq_summary.jsonl"
    output:
        "00_metadata/bacillus_master_metadata.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            bash scripts/02_build_master_metadata.sh
        """


# ============================================================
# 03. Plant-associated candidates
# ============================================================

rule identify_plant_candidates:
    input:
        "00_metadata/bacillus_master_metadata.tsv"
    output:
        "00_metadata/plant_associated_candidates.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/03_identify_plant_candidates.py
        """


# ============================================================
# 04. Plant-association classification
# ============================================================

rule classify_plant_association:
    input:
        "00_metadata/plant_associated_candidates.tsv"
    output:
        "00_metadata/plant_association_classification.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/04_classify_plant_association.py
        """


rule audit_plant_association:
    input:
        "00_metadata/plant_association_classification.tsv"
    output:
        touch("logs/04_audit_plant_association.done")
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/04_audit_plant_association.py

        touch {output}
        """


# ============================================================
# 05. Final plant-associated dataset
# ============================================================

rule finalize_plant_selection:
    input:
        classification="00_metadata/plant_association_classification.tsv",
        audit="logs/04_audit_plant_association.done"
    output:
        selected="00_metadata/plant_associated_selected.tsv",
        unresolved="00_metadata/plant_association_unresolved.tsv",
        final="00_metadata/plant_association_final_classification.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/05_finalize_plant_selection.py
        """


# ============================================================
# 06. Download selected genomes
# ============================================================

rule download_selected_genomes:
    input:
        "00_metadata/plant_associated_selected.tsv"
    output:
        accessions="01_genomes/selected_accessions.txt",
        dataset=directory("01_genomes/ncbi_dataset")
    shell:
        """
        conda run -n {ENV_AMR} \
            bash scripts/06_download_selected_genomes.sh
        """


# ============================================================
# 07. Structural genome statistics
# ============================================================

rule genome_qc_seqkit:
    input:
        dataset="01_genomes/ncbi_dataset"
    output:
        "02_qc/seqkit_genome_stats.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            bash scripts/07_genome_qc_seqkit.sh
        """


# ============================================================
# 08. Structural QC filter
# ============================================================

rule filter_structural_qc:
    input:
        "02_qc/seqkit_genome_stats.tsv"
    output:
        classification="02_qc/structural_qc_classification.tsv",
        passed="02_qc/qc_passed_accessions.txt",
        failed="02_qc/qc_failed_accessions.txt"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/08_filter_structural_qc.py
        """


# ============================================================
# 09. Taxonomic normalization
# ============================================================

rule normalize_taxonomy:
    input:
        qc="02_qc/qc_passed_accessions.txt",
        metadata="00_metadata/plant_association_final_classification.tsv"
    output:
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv",
        counts="03_taxonomy/species_counts.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/09_normalize_taxonomy.py
        """


# ============================================================
# 10. AMRFinderPlus
# ============================================================

rule run_amrfinder:
    input:
        qc="02_qc/qc_passed_accessions.txt",
        genomes="01_genomes/ncbi_dataset"
    output:
        touch("05_amr/amrfinder_complete.done")
    threads: 4
    shell:
        """
        conda run -n {ENV_AMRFINDER} \
            bash scripts/10_run_amrfinder.sh

        touch {output}
        """


# ============================================================
# 11. Consolidate AMRFinderPlus results
# ============================================================

rule summarize_amrfinder:
    input:
        "05_amr/amrfinder_complete.done"
    output:
        all_hits="05_amr/amrfinder_all_hits.tsv",
        amr_hits="05_amr/amrfinder_amr_hits.tsv",
        type_summary="05_amr/amrfinder_type_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/11_summarize_amrfinder.py
        """


# ============================================================
# 12. AMR prevalence
# ============================================================

rule amr_prevalence:
    input:
        amr="05_amr/amrfinder_amr_hits.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv"
    output:
        presence="05_amr/amr_presence_absence.tsv",
        prevalence="05_amr/amr_prevalence_by_species.tsv",
        genome_summary="05_amr/amr_genome_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/12_amr_prevalence_by_species.py
        """


# ============================================================
# 13. Sporadic AMR carriers
# ============================================================

rule identify_sporadic_carriers:
    input:
        prevalence="05_amr/amr_prevalence_by_species.tsv",
        presence="05_amr/amr_presence_absence.tsv"
    output:
        cases="05_amr/sporadic_amr_cases.tsv",
        carriers="05_amr/sporadic_amr_carriers.tsv",
        summary="05_amr/sporadic_amr_carrier_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/13_identify_sporadic_amr_carriers.py
        """


# ============================================================
# 13b. GToTree phylogenomics
# ============================================================

rule run_gtotree:
    input:
        qc="02_qc/qc_passed_accessions.txt",
        genomes="01_genomes/ncbi_dataset"
    output:
        alignment="06_phylogeny/gtotree_firmicutes/Aligned_SCGs.faa"
    threads: 4
    shell:
        """
        conda run -n {ENV_PHYLO} \
            bash scripts/13_run_gtotree.sh
        """


# ============================================================
# 14. Targeted phylogenetic subset
# ============================================================

rule select_phylogenetic_subset:
    input:
        carriers="05_amr/sporadic_amr_carriers.tsv",
        presence="05_amr/amr_presence_absence.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv"
    output:
        table="06_phylogeny/phylogenetic_subset.tsv",
        accessions="06_phylogeny/phylogenetic_subset_accessions.txt"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/14_select_phylogenetic_subset.py
        """


# ============================================================
# Extract 87-genome alignment
# ============================================================

rule extract_phylogenetic_subset_alignment:
    input:
        alignment="06_phylogeny/gtotree_firmicutes/Aligned_SCGs.faa",
        accessions="06_phylogeny/phylogenetic_subset_accessions.txt"
    output:
        "06_phylogeny/Aligned_SCGs_subset87.faa"
    shell:
        """
        conda run -n {ENV_ECO} \
            seqkit grep \
            -f {input.accessions} \
            {input.alignment} \
            > {output}
        """


# ============================================================
# ModelFinder
# ============================================================

rule phylogenetic_model_selection:
    input:
        "06_phylogeny/Aligned_SCGs_subset87.faa"
    output:
        iqtree="06_phylogeny/iqtree_subset87/modelfinder.iqtree"
    threads: 4
    shell:
        """
        mkdir -p 06_phylogeny/iqtree_subset87

        conda run -n {ENV_ECO} \
            iqtree \
            -s {input} \
            -m TESTONLY \
            -T {threads} \
            --prefix 06_phylogeny/iqtree_subset87/modelfinder
        """


# ============================================================
# Final targeted ML phylogeny
# ============================================================

rule final_phylogenetic_tree:
    input:
        alignment="06_phylogeny/Aligned_SCGs_subset87.faa",
        model="06_phylogeny/iqtree_subset87/modelfinder.iqtree"
    output:
        tree="06_phylogeny/iqtree_subset87/final_tree.treefile",
        iqtree="06_phylogeny/iqtree_subset87/final_tree.iqtree"
    threads: 4
    shell:
        """
        conda run -n {ENV_ECO} \
            iqtree \
            -s {input.alignment} \
            -m Q.PLANT+F+I+G4 \
            -B 1000 \
            -alrt 1000 \
            -T {threads} \
            --prefix 06_phylogeny/iqtree_subset87/final_tree
        """


# ============================================================
# 15. Map sporadic AMR loci
# ============================================================

rule map_sporadic_amr_loci:
    input:
        hits="05_amr/amrfinder_amr_hits.tsv",
        carriers="05_amr/sporadic_amr_carriers.tsv"
    output:
        "08_genomic_context/sporadic_amr_loci.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/15_map_sporadic_amr_loci.py
        """


# ============================================================
# 16. Extract +/-10 kb genomic contexts
# ============================================================

rule extract_context_windows:
    input:
        loci="08_genomic_context/sporadic_amr_loci.tsv",
        genomes="01_genomes/ncbi_dataset"
    output:
        fasta="08_genomic_context/sporadic_amr_context_10kb.fna",
        table="08_genomic_context/sporadic_amr_context_windows.tsv",
        windows=directory("08_genomic_context/windows_10kb")
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/16_extract_amr_context_windows.py
        """


# ============================================================
# 18. Uniform Prodigal meta annotation
# Step 17 Prokka run intentionally excluded
# ============================================================

rule prodigal_meta_contexts:
    input:
        windows="08_genomic_context/windows_10kb"
    output:
        proteins="08_genomic_context/context_all_proteins.faa",
        done=touch("08_genomic_context/prodigal_meta_complete.done")
    shell:
        """
        conda run -n {ENV_ECO} \
            bash scripts/18_prodigal_meta_contexts.sh

        cat 08_genomic_context/prodigal_meta/*.faa \
            > {output.proteins}

        touch {output.done}
        """


# ============================================================
# 19. Validate AMR/CDS overlap
# ============================================================

rule validate_amr_cds_overlap:
    input:
        loci="08_genomic_context/sporadic_amr_loci.tsv",
        prodigal="08_genomic_context/prodigal_meta_complete.done"
    output:
        "08_genomic_context/amr_cds_overlap_validation.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/19_validate_amr_cds_overlap.py
        """


# ============================================================
# IS database BLASTP
# ============================================================

rule blast_IS:
    input:
        "08_genomic_context/context_all_proteins.faa"
    output:
        "08_genomic_context/context_IS_blastp_raw.tsv"
    threads: 4
    shell:
        """
        conda run -n {ENV_ECO} bash -c '
        blastp \
          -query {input} \
          -db "$CONDA_PREFIX/db/kingdom/Bacteria/IS" \
          -evalue 1e-10 \
          -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle" \
          -max_target_seqs 5 \
          -num_threads {threads} \
          > {output}
        '
        """


# ============================================================
# Swiss-Prot mobility search
# ============================================================

rule blast_sprot:
    input:
        "08_genomic_context/context_all_proteins.faa"
    output:
        "08_genomic_context/context_sprot_blastp_raw.tsv"
    threads: 4
    shell:
        """
        conda run -n {ENV_ECO} bash -c '
        blastp \
          -query {input} \
          -db "$CONDA_PREFIX/db/kingdom/Bacteria/sprot" \
          -evalue 1e-10 \
          -outfmt "6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle" \
          -max_target_seqs 5 \
          -num_threads {threads} \
          > {output}
        '
        """


# ============================================================
# 20. Filter IS hits
# ============================================================

rule filter_IS_hits:
    input:
        "08_genomic_context/context_IS_blastp_raw.tsv"
    output:
        classified="08_genomic_context/context_IS_hits_classified.tsv",
        accepted="08_genomic_context/context_IS_hits_accepted.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/20_filter_IS_hits.py
        """


# ============================================================
# 21. IS proximity
# ============================================================

rule map_IS_proximity:
    input:
        accepted="08_genomic_context/context_IS_hits_accepted.tsv",
        overlap="08_genomic_context/amr_cds_overlap_validation.tsv"
    output:
        proximity="08_genomic_context/amr_IS_proximity.tsv",
        summary="08_genomic_context/amr_IS_proximity_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/21_map_IS_proximity_to_AMR.py
        """


# ============================================================
# 22. Integrate mobility evidence
# ============================================================

rule integrate_mobility:
    input:
        proximity="08_genomic_context/amr_IS_proximity.tsv",
        sprot="08_genomic_context/context_sprot_blastp_raw.tsv"
    output:
        evidence="08_genomic_context/mobility_evidence.tsv",
        summary="08_genomic_context/mobility_evidence_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/22_integrate_mobility_evidence.py
        """


# ============================================================
# 23. Explicit NCBI plasmid annotation
# ============================================================

rule check_plasmid_annotation:
    input:
        loci="08_genomic_context/sporadic_amr_loci.tsv",
        genomes="01_genomes/ncbi_dataset"
    output:
        "08_genomic_context/amr_locus_plasmid_annotation.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/23_check_plasmid_annotation.py
        """


# ============================================================
# 24. Integrated locus evidence
# ============================================================

rule integrate_locus_evidence:
    input:
        contexts="08_genomic_context/sporadic_amr_context_windows.tsv",
        mobility="08_genomic_context/mobility_evidence.tsv",
        plasmid="08_genomic_context/amr_locus_plasmid_annotation.tsv"
    output:
        "08_genomic_context/sporadic_amr_locus_evidence.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/24_integrate_locus_evidence.py
        """


# ============================================================
# 25. Phylogenetic discontinuity
# ============================================================

rule phylogenetic_discontinuity:
    input:
        tree="06_phylogeny/iqtree_subset87/final_tree.treefile",
        subset="06_phylogeny/phylogenetic_subset.tsv",
        presence="05_amr/amr_presence_absence.tsv",
        carriers="05_amr/sporadic_amr_carriers.tsv"
    output:
        detail="09_analysis/sporadic_amr_phylogenetic_discontinuity.tsv",
        summary="09_analysis/sporadic_amr_phylogenetic_discontinuity_summary.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/25_phylogenetic_discontinuity.py
        """


# ============================================================
# 26. Integrated evidence
# ============================================================

rule integrate_hgt_evidence:
    input:
        phylo="09_analysis/sporadic_amr_phylogenetic_discontinuity_summary.tsv",
        context="08_genomic_context/sporadic_amr_locus_evidence.tsv"
    output:
        "09_analysis/sporadic_amr_integrated_evidence.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/26_integrate_hgt_evidence.py
        """


# ============================================================
# 27. Species-level AMR burden
# ============================================================

rule amr_burden:
    input:
        "05_amr/amr_genome_summary.tsv"
    output:
        summary="09_analysis/amr_burden_by_species_summary.tsv",
        dunn="09_analysis/amr_burden_dunn_bh.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/27_amr_burden_by_species.py
        """


# ============================================================
# 28. Species resistome prevalence matrix
# ============================================================

rule species_resistome_matrix:
    input:
        prevalence="05_amr/amr_prevalence_by_species.tsv"
    output:
        matrix="09_analysis/species_resistome_prevalence_matrix.tsv",
        long="09_analysis/species_resistome_prevalence_long.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/28_build_species_resistome_matrix.py
        """


# ============================================================
# 29. PERMANOVA + PCoA
# ============================================================

rule resistome_permanova_pcoa:
    input:
        presence="05_amr/amr_presence_absence.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv"
    output:
        permanova="09_analysis/resistome_permanova.tsv",
        pcoa="09_analysis/resistome_pcoa_coordinates.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/29_resistome_composition_permanova.py
        """


# ============================================================
# 30. PERMDISP
# ============================================================

rule resistome_permdisp:
    input:
        presence="05_amr/amr_presence_absence.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv"
    output:
        "09_analysis/resistome_permdisp.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/30_resistome_permdisp.py
        """


# ============================================================
# 31. Within-species heterogeneity
# ============================================================

rule within_species_heterogeneity:
    input:
        presence="05_amr/amr_presence_absence.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv"
    output:
        "09_analysis/within_species_resistome_heterogeneity.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/31_within_species_resistome_heterogeneity.py
        """


# ============================================================
# Figures
# ============================================================

rule figure1:
    input:
        burden="09_analysis/amr_burden_by_species_summary.tsv",
        genomes="05_amr/amr_genome_summary.tsv"
    output:
        pdf="10_figures/Figure1_AMR_burden_by_species.pdf",
        png="10_figures/Figure1_AMR_burden_by_species.png"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/32_figure1_amr_burden.py
        """


rule figure2:
    input:
        "09_analysis/species_resistome_prevalence_matrix.tsv"
    output:
        pdf="10_figures/Figure2_resistome_prevalence_heatmap.pdf",
        png="10_figures/Figure2_resistome_prevalence_heatmap.png"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/33_figure2_resistome_heatmap.py
        """


rule figure3:
    input:
        pcoa="09_analysis/resistome_pcoa_coordinates.tsv",
        permanova="09_analysis/resistome_permanova.tsv",
        permdisp="09_analysis/resistome_permdisp.tsv"
    output:
        pdf="10_figures/Figure3_resistome_PCoA.pdf",
        png="10_figures/Figure3_resistome_PCoA.png"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/34_figure3_resistome_pcoa.py
        """


rule figure4:
    input:
        tree="06_phylogeny/iqtree_subset87/final_tree.treefile",
        subset="06_phylogeny/phylogenetic_subset.tsv",
        evidence="09_analysis/sporadic_amr_integrated_evidence.tsv",
        loci="08_genomic_context/sporadic_amr_locus_evidence.tsv"
    output:
        pdf="10_figures/Figure4_phylogenetic_context.pdf",
        png="10_figures/Figure4_phylogenetic_context.png"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/35_figure4_phylogenetic_context.py
        """


# ============================================================
# Tables
# ============================================================

rule table1:
    input:
        candidates="00_metadata/plant_association_classification.tsv",
        qc="02_qc/seqkit_genome_stats.tsv",
        taxonomy="03_taxonomy/qc_passed_taxonomy.tsv",
        hits="05_amr/amrfinder_amr_hits.tsv",
        summary="05_amr/amr_genome_summary.tsv"
    output:
        "11_tables/Table1_dataset_QC_AMR_summary.tsv"
    shell:
        """
        conda run -n {ENV_AMR} \
            python scripts/36_table1_dataset_summary.py
        """


rule table2:
    input:
        burden="09_analysis/amr_burden_by_species_summary.tsv",
        heterogeneity="09_analysis/within_species_resistome_heterogeneity.tsv"
    output:
        "11_tables/Table2_species_AMR_burden_heterogeneity.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/37_table2_species_burden_heterogeneity.py
        """


rule table3:
    input:
        "09_analysis/sporadic_amr_integrated_evidence.tsv"
    output:
        "11_tables/Table3_sporadic_AMR_evidence.tsv"
    shell:
        """
        conda run -n {ENV_ECO} \
            python scripts/38_table3_sporadic_AMR_evidence.py
        """
