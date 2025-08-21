import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Define the chromosome lengths for the human T2T genome (in base pairs)
chromosome_lengths = {
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

# Optionally include or exclude 'chrM' and 'chr'
exclude_chromosomes = ['chrM', 'chr']

# Load your data (adjust the file path if necessary)
file_path = 'sniffles_merged_all_samples_SVs_v6.txt'
df = pd.read_csv(file_path, sep="\t")

# Filter out unwanted chromosomes
df = df[~df['chr'].isin(exclude_chromosomes)]

# Ensure that the 'start' and 'end' columns are numeric
df['start'] = pd.to_numeric(df['start'], errors='coerce')
df['end'] = pd.to_numeric(df['end'], errors='coerce')

# Group by chromosome and sample filename, counting the number of SVs for each
sv_counts = df.groupby(['chr', 'filename']).size().reset_index(name='SV_count')

# Pivot the data to create a matrix with chromosomes as rows and samples as columns
sv_matrix = sv_counts.pivot_table(index='chr', columns='filename', values='SV_count', fill_value=0)

# Ensure the SV counts are numeric
sv_matrix = sv_matrix.apply(pd.to_numeric, errors='coerce')

# Sort chromosomes based on numerical order
# Define a function to handle the sorting of chromosomes
def sort_chromosomes(chromosome):
    if chromosome.startswith('chr'):
        # For chromosomes like chr1, chr2, ..., chr22, chrX, chrY
        if chromosome[3:].isdigit():
            return int(chromosome[3:])  # Return numeric value for numeric chromosomes
        elif chromosome == 'chrX':
            return 23  # For chrX
        elif chromosome == 'chrY':
            return 24  # For chrY
    return chromosome

# Apply the sorting function and reorder the rows in the matrix
sorted_chromosomes = sorted(sv_matrix.index, key=sort_chromosomes)
sv_matrix = sv_matrix.loc[sorted_chromosomes]

# Normalize the SV counts by chromosome length (SV density)
for chr in sv_matrix.index:
    chromosome_length = chromosome_lengths.get(chr, 0)
    if chromosome_length > 0:
        # Normalize to per million base pairs for each chromosome
        sv_matrix.loc[chr] = (sv_matrix.loc[chr] / chromosome_length) * 1e6  # Normalizing to per million bases

# Plot the heatmap with normalized values
plt.figure(figsize=(12, 8))
sns.heatmap(sv_matrix, cmap='YlGnBu', annot=True, fmt='.2f', linewidths=0.5)

# Add labels and title
plt.title("Normalized Heatmap of Structural Variants (SVs) by Chromosome and Sample")
plt.ylabel('Chromosome')
plt.xlabel('Sample')

# Save the heatmap as a PNG file
output_file = 'Normalized_SVs_heatmap_T2T_sorted.png'
plt.tight_layout()
plt.savefig(output_file)

# Show the plot
plt.show()

print(f"Heatmap saved as {output_file}")

