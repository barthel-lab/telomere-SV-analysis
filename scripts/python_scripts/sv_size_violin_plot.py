import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
data = pd.read_csv("sniffles_merged_all_samples_SVs_v6.txt", sep="\t")

# Ensure 'start' and 'end' columns are numeric
data['start'] = pd.to_numeric(data['start'], errors='coerce')
data['end'] = pd.to_numeric(data['end'], errors='coerce')

# Create a new column for the SV length
data['length'] = data['end'] - data['start'] + 1

# Filter out non-positive values from the length column
data = data[data['length'] > 0]

# Sort the data by sample in ascending order
sample_order = sorted(data['filename'].unique())

# --- Create a figure with subplots for each SV type ---
sv_types = data['type'].unique()

# Adjust figure height to make it smaller while keeping readability
fig, axs = plt.subplots(len(sv_types), 1, figsize=(12, 3 * len(sv_types)))

# If there's only one subplot, axs will not be an array, so we need to handle it
if len(sv_types) == 1:
    axs = [axs]

# Define a color palette for the SV types
palette = sns.color_palette("husl", len(sv_types))

# Loop over each SV type to create a violin plot for each
for i, sv_type in enumerate(sv_types):
    # Filter data for the current SV type
    sv_data = data[data['type'] == sv_type]
    
    # Violin plot for SV length by sample for each SV type
    sns.violinplot(
        x='filename', y='length', data=sv_data, ax=axs[i], inner="quart",
        width=0.8, palette=[palette[i]] * len(sample_order)
    )

    # Set plot labels and title
    axs[i].set_xticklabels(sample_order, rotation=90)
    axs[i].set_title(f"{sv_type} Lengths by Sample", fontweight='bold')  # Bold header
    axs[i].set_xlabel("Sample")
    axs[i].set_ylabel("SV Length")

    # Use log scale for the y-axis to better visualize the range of SV lengths
    axs[i].set_yscale('log')

# Adjust layout to prevent overlapping
plt.tight_layout()

# Save the figure to a file (e.g., PNG or PDF)
plt.savefig("sv_type_violin_plots.png", dpi=300)  # Save as PNG with high resolution
# Alternatively, you can save as PDF:
# plt.savefig("sv_type_violin_plots.pdf")

# Show the plot
plt.show()

