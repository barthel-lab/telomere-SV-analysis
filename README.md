Notebooks

This directory contains the Jupyter notebooks and scripts used to perform the telomere-SV analysis and generate data for the project figures. All processed genomics data is hosted on Synapse:
Project Link:https://www.synapse.org/Synapse:syn72016012/datasets/

Directory Structure

Fig 1: TCGA Cohort Overview & Telomere Metrics

    TCGA_GLASS_paired_samples_only_updated.ipynb: Comparative analysis of paired TCGA and GLASS datasets.

    TCGA_GLASS_unpaired_samples.ipynb: Analysis of unpaired cohort samples.

    TCGA_telomere_associationplots_clean.ipynb: Telomere association plotting for TCGA

    long-reads_cov_n50.ipynb: Sequencing coverage and N50 quality metrics from our long reads cohort.
    
Fig 2: Structural Variant & Clinical Integration

    clinical_data_heatmap.ipynb: Integration of clinical metadata with genomic findings.

    severus/: Sub-analysis of SVs using the Severus caller (BND subtypes, barplots, and heatmaps).

    sniffles/: Comparative SV analysis using the Sniffles caller.

    CNVs/: Visualization of Copy Number Variants via barplots and heatmaps.

Fig 3: ecDNA & Genomic Repeats

    decoil.ipynb: Analysis using the Decoil tool for ecDNA reconstruction.

    HGs_extraction.ipynb: Extraction of Gigh-Gains from Copy number data.
    
Fig 4: Telomere Associations & Genomic Instability

    telomers_analysis_v2.ipynb: Primary telomere length analysis.

    generate_figures.py: Final figure generation (A–E) and statistical modeling.
    
Setup environment for generate_figures.py

    conda create -n telomere_env python=3.12 -y

    conda activate telomere_env

    pip install pandas numpy matplotlib seaborn scipy statsmodels

Fig 5: Distance & Spatial Analysis

    TCGA_GLASS_distance_analysis.ipynb: Distance analysis across TCGA and GLASS cohorts.

     severus/: Spatial SV and somatic distance analysis using Severus.

     sniffles/: Spatial SV analysis and line plots using Sniffles.
     
Additional Analyses

    aneuploidy_plot.ipynb: Chromosomal arm-level alterations.

    survival_analysis.ipynb: Clinical survival modeling and genomic metrics.
    
Usage

Most notebooks depend on processed data located in the ../data directory and output results to ../output or ../plots.# telomere-SV-analysis
