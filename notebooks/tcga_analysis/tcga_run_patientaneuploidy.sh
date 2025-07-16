#!/bin/bash
#SBATCH --job-name=patient_aneuploidy          # Job name
#SBATCH --output=patient_aneuploidy.out          # Standard output file
#SBATCH --error=patient_aneuploidy.err           # Standard error file
#SBATCH --ntasks=1                             # Number of tasks (processes)
#SBATCH --cpus-per-task=1                      # Number of CPU cores per task
#SBATCH --mem=500G                             # Memory allocation (500 GB)
#SBATCH --time=24:00:00                        # Maximum runtime, adjust as needed

# Activate your Conda environment where the R libraries are installed
source /tgen_labs/barthel/software/miniforge3/etc/profile.d/conda.sh
conda activate aneuploidy

# Change directory to where your R script and data are located, if needed
#cd /scratch/mjehangir/Glioma_project/public_glioma_data/

# Run the R script using Rscript command with the required arguments.
Rscript Maria_aneuploidy.r 

