import pandas as pd
import matplotlib.pyplot as plt

# Load the data
data = pd.read_csv("sniffles_merged_all_samples_SVs_v6.txt", sep="\t")

# Ensure 'start' and 'end' columns are numeric (convert if needed)
data['start'] = pd.to_numeric(data['start'], errors='coerce')
data['end'] = pd.to_numeric(data['end'], errors='coerce')

# Create a new column for the SV length
data['length'] = data['end'] - data['start'] + 1

# Sort the data by sample in ascending order
sample_order = sorted(data['filename'].unique())

# Create a 2x1 grid of subplots (2 rows, 1 column)
fig, axs = plt.subplots(2, 1, figsize=(12, 12))

# --- Plot 1: Total Number of SVs in Each Sample ---
# Count the number of SVs per sample
sv_counts = data['filename'].value_counts().reindex(sample_order)

# Plot total number of SVs per sample with the sample order
axs[0].bar(sv_counts.index, sv_counts.values, color='steelblue')
axs[0].set_title("Total Number of SVs in Each Sample")
axs[0].set_xlabel("Sample")
axs[0].set_ylabel("Number of SVs")
axs[0].tick_params(axis='x', rotation=45)

# --- Plot 2: Different Types of SVs in Each Sample ---
# Group by sample and SV type to get counts of each SV type
sv_types = data.groupby(['filename', 'type']).size().unstack(fill_value=0)

# Reorder the rows (samples) in the sv_types DataFrame based on the sample_order
sv_types = sv_types.loc[sample_order]  # Reorder rows by sample_order

# Define colors for each SV type
colors = {'INS': 'blue', 'DEL': 'red', 'DUP': 'green', 'INV': 'orange'}

# Plot different types of SVs in each sample with custom colors
sv_types.plot(
    kind='bar', 
    stacked=True, 
    ax=axs[1], 
    figsize=(12, 6), 
    color=[colors.get(x, 'gray') for x in sv_types.columns]
)
axs[1].set_title("SV Types in Each Sample")
axs[1].set_xlabel("Sample")
axs[1].set_ylabel("Number of SVs")
axs[1].tick_params(axis='x', rotation=45)

# Place the legend on the right side of the second plot
axs[1].legend(title='SV Type', loc='center left', bbox_to_anchor=(1, 0.5))

# Adjust layout to prevent overlapping
plt.tight_layout()

# Save the figure to a file (e.g., PNG or PDF)
plt.savefig("sv_analysis_plots_updated.png", dpi=300)  # Save as PNG with high resolution
# Alternatively, you can save as PDF:
# plt.savefig("sv_analysis_plots_updated.pdf")

# Show the plot
plt.show()

