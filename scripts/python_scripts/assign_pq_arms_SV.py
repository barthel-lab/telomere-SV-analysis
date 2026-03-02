import pandas as pd

# Define chromosome arm boundaries from the provided data
chromosome_arms = {
    "chr1": 124048267, "chr2": 93503283, "chr3": 94076514, "chr4": 52452474,
    "chr5": 48317879, "chr6": 59672548, "chr7": 62064435, "chr8": 45270456,
    "chr9": 46267185, "chr10": 40649191, "chr11": 52743313, "chr12": 35911664,
    "chr13": 16522942, "chr14": 11400261, "chr15": 17186630, "chr16": 36838903,
    "chr17": 25689679, "chr18": 18449624, "chr19": 27792923, "chr20": 28012753,
    "chr21": 11134529, "chr22": 14249622
}

# Load your variant file
df = pd.read_csv("sniffles_merged_all_samples_SVs_v6.txt", sep="\t")

# Function to assign p/q arms
def assign_arm(row):
    centromere = chromosome_arms.get(row["chr"], None)
    if centromere:
        return "p" if row["start"] < centromere else "q"
    return "Unknown"

# Apply function
df["arm"] = df.apply(assign_arm, axis=1)

# Save the output
df.to_csv("sniffles_merged_all_samples_SVs_v7.txt", sep="\t", index=False)

print("Processing complete. Check 'variants_with_arms.txt' for results.")

