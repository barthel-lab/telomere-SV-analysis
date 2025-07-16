import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Chromosome sizes for T2T human genome (in base pairs)
chrom_sizes = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}

# Load the data from your file
file_path = 'sniffles_merged_all_samples_SVs_v6.txt'  # Adjust this path
df = pd.read_csv(file_path, sep="\t")

# Filter out unwanted chromosomes (e.g., 'chrM', 'chr', etc.)
df = df[~df['chr'].isin(['chrM', 'chr'])]

# Get the unique SV types
sv_types = df['type'].unique()

# Prepare the plot for each SV type
for sv_type in sv_types:
    # Filter the dataframe for the current SV type
    df_sv_type = df[df['type'] == sv_type]
    
    # Create an empty DataFrame to store the count of each SV type per chromosome for all samples
    sv_count_by_chr = pd.DataFrame(columns=df['filename'].unique(), index=df_sv_type['chr'].unique())
    
    # Fill the dataframe with the counts of the current SV type per chromosome for each sample
    for chrom in sv_count_by_chr.index:
        for sample in sv_count_by_chr.columns:
            count = len(df_sv_type[(df_sv_type['chr'] == chrom) & (df_sv_type['filename'] == sample)])
            sv_count_by_chr.at[chrom, sample] = count
    
    # Normalize by chromosome size
    for chrom in sv_count_by_chr.index:
        size = chrom_sizes.get(chrom, 0)  # Get the size for the chromosome
        if size > 0:
            sv_count_by_chr.loc[chrom] = sv_count_by_chr.loc[chrom] / size
    
    # Convert to float and apply log transformation to SV counts (log2 is commonly used)
    sv_count_by_chr = sv_count_by_chr.astype(float)  # Convert to float
    sv_count_by_chr_log = np.log2(sv_count_by_chr + 1)  # Adding 1 to avoid log(0)

    # Plotting the heatmap for the current SV type with log-transformed counts
    plt.figure(figsize=(14, 8))
    sns.heatmap(sv_count_by_chr_log, annot=False, cmap='Blues', cbar_kws={'label': 'Log2(SV count + 1)'}, linewidths=0.5)

    # Set the title and labels
    plt.title(f'{sv_type} SV Distribution Across Chromosomes for All Samples (Log2 scale, Normalized by Chromosome Size)')
    plt.xlabel('Samples')
    plt.ylabel('Chromosomes')

    # Save the heatmap as an image
    output_file = f'{sv_type}_SV_heatmap_all_samples_log2_normalized_no_values.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')  # Save as high-resolution PNG
    plt.close()  # Close the plot after saving


