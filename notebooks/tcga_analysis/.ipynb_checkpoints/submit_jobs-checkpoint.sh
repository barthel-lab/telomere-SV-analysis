#!/bin/bash

# Directory containing the .seg.txt files
input_dir="/scratch/mjehangir/Glioma_project/public_glioma_data/split_by_prefix"

# List all .seg.txt files in the directory
input_files=($(ls $input_dir/*.seg))

# Loop through each file and submit an individual job for each
for input_file in "${input_files[@]}"; do
  file_name=$(basename $input_file)

  # Submit a job for each input file
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=aneuploidy_${file_name}
#SBATCH --output=aneuploidy_${file_name}.out
#SBATCH --error=aneuploidy_${file_name}.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --time=5:00:00

# Activate the Conda environment
source /tgen_labs/barthel/software/miniforge3/etc/profile.d/conda.sh
conda activate aneuploidy

# Run the R script for the current input file
Rscript Maria_aneuploidy_for_TCGA_chuknks_V3.r ${file_name}
EOF

  echo "Job submitted for $file_name"
done

