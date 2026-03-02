# Multi-Tool SV Calling Snakemake Pipeline

## Overview
This Snakemake pipeline enables structural variant (SV) calling from **long-read sequencing data** across multiple samples.  
It integrates **three SV callers** – [Sniffles2](https://github.com/fritzsedlazeck/Sniffles), [Severus](https://github.com/chhylp123/severus), and [SAVANA](https://github.com/fcrisci/savana) – to detect **germline** and **somatic** SVs.

The pipeline supports:
- **Tumor/Normal paired analysis**  
- **Tumor-only mode** (no control available)  

By default, all tools run in parallel (if resources allow), each producing results in its own output directory.

---

## Key Features
- **Multiple SV Callers**
  - *Sniffles2* → germline + optional mosaic mode  
  - *Severus* → somatic detection (tumor/normal or tumor-only)  
  - *SAVANA* → somatic detection + copy number analysis (optional)  

- **Parallel Execution**  
  Orchestrated by Snakemake, with configurable **threads** and **memory** for HPC use.

- **Flexible Tool Selection**  
  Run all tools or only a subset (e.g., `snakemake --cores 8 sniffles_only`).

- **Reference Genome Support**  
  Works with **T2T-CHM13**, **GRCh38**, or **GRCh37/hg19**.  
  Maps automatically to required resources (FASTA, repeats, PoN, known variants).

- **Simple Configuration**  
  Controlled via a single `config.yaml` file (inputs, reference, tool options, output paths).  

---

## Tools & Modes

### Sniffles2
- Germline mode (default).  
- Optional `--mosaic` mode for rare/low-frequency variants.  
- Uses tandem repeat annotations for improved accuracy.  
- **Output:** VCF  

### Severus
- Designed for **somatic SVs** in long reads.  
- Modes:
  - Tumor/Normal → outputs **somatic-only** + **unfiltered** VCFs  
  - Tumor-only → outputs all SVs  
- Uses VNTR BED and Panel of Normals for filtering.  
- **Output:** VCF  

### SAVANA
- Detects **somatic SVs** and **copy number alterations**.  
- Modes:
  - Tumor/Normal → somatic calls  
  - Tumor-only → all SVs  
- Optional CNA calling with SNP VCF (or 1000 Genomes resource).  
- **Output:** VCF + BEDPE  

---

## Workflow & Outputs

- **Default target runs:**
  - Sniffles2 (standard + mosaic)  
  - Severus (somatic + unfiltered)  
  - SAVANA (somatic + tumor-only)  

- **Selective runs:**  
  ```bash
  snakemake --cores 16 severus_only savana_only




results/
├── sniffles/
│   ├── standard/
│   └── mosaic/
├── severus/
│   ├── somatic/
│   └── unfiltered/
└── savana/
    ├── somatic/
    └── tumor_only/


## Configuration

Edit the `config.yaml` file to specify:

- **Input BAMs** and sample roles (tumor/normal)
- **Reference genome**: `t2t`, `hg38`, or `hg19`
- **Output directory** path
- **Tool options** (e.g., enable/disable Sniffles2 mosaic mode)
- **Resource allocation** (threads, memory per job)

---

## Requirements

- [Snakemake](https://snakemake.readthedocs.io/) ≥ 6.0  
- Conda environments or pre-installed binaries for:  
  - [Sniffles2](https://github.com/fritzsedlazeck/Sniffles)  
  - [Severus](https://github.com/chhylp123/severus)  
  - [SAVANA](https://github.com/fcrisci/savana)  
- HPC or local machine (default: **8 threads**, **16 GB RAM per job**)  

---

## Citation

If you use this pipeline, please cite:

- **Sniffles2**: Sedlazeck et al., *Nature Biotechnology* (2018 & updates)  
- **Severus**: Cheng et al., GitHub Project  
- **SAVANA**: Crisci et al., GitHub Project  


