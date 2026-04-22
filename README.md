Notebooks

This directory contains the Jupyter notebooks and scripts used to perform the telomere-SV analysis and generate data for the project figures.
Directory Structure
Fig 1: TCGA Cohort Overview & Telomere Metrics

    CNVs_tcga_glass_analysis_v2.ipynb: Comparative analysis of TCGA and GLASS datasets.

    cov_n50.ipynb: Sequencing coverage and N50 quality metrics.

Fig 2: Structural Variant & Clinical Integration

    clinical_data_heatmap.ipynb: Integration of clinical metadata with genomic findings.

    severus/: Sub-analysis of SVs using the Severus caller (BND subtypes, barplots, and heatmaps).

    sniffles/: Comparative SV analysis using the Sniffles caller.

    CNVs/: Visualization of Copy Number Variants via barplots and heatmaps.

Fig 3: ecDNA & Genomic Repeats

    decoil.ipynb: Analysis using the Decoil tool for ecDNA reconstruction.

    HARs_extraction.ipynb & HARs_repeats_density.ipynb: Analysis of HARs.

Fig 4: Telomere Associations & Genomic Instability

    telomers_analysis_v2.ipynb: Primary telomere length analysis.

    generate_figures.py: Final figure generation (A–E) and statistical modeling.
    
Fig 5: Distance & Spatial Analysis

    severus/ & sniffles/: Spatial SV analysis

Additional Analyses

    aneuploidy/: Chromosomal arm-level alterations and survival analysis.

    tcga_analysis/: Validation and comparative workflows using TCGA glioma datasets (IDH subtypes, ecDNA, and patient aneuploidy).

Usage

Most notebooks depend on processed data located in the ../data directory and output results to ../output or ../plots.# telomere-SV-analysis
