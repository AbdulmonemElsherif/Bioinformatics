import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import wilcoxon
from ast import literal_eval

# Load the CSV file with dissimilar regions
dissimilar_regions_df = pd.read_csv("dissimilar_regions.csv")

# Calculate the length of each dissimilar region
dissimilar_regions_df['Region Length'] = dissimilar_regions_df['Positions'].apply(lambda x: len(eval(x)))

# Descriptive statistics of region lengths
region_length_stats = dissimilar_regions_df['Region Length'].describe()
print(region_length_stats)

# Plot the distribution of dissimilar region lengths
plt.figure(figsize=(10, 6))
sns.histplot(dissimilar_regions_df['Region Length'], bins=20, kde=True)
plt.title('Distribution of Dissimilar Region Lengths')
plt.xlabel('Region Length')
plt.ylabel('Frequency')
plt.show()

# Plot the distribution of dissimilar region lengths with log scale
plt.figure(figsize=(10, 6))
sns.histplot(dissimilar_regions_df['Region Length'], bins=20, kde=True)
plt.yscale('log')  
plt.title('Distribution of Dissimilar Region Lengths (Log Scale)')
plt.xlabel('Region Length')
plt.ylabel('Frequency (log scale)')
plt.show()

# Base composition analysis
base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0, 'CG': 0}

for bases in dissimilar_regions_df['Delta Consensus']:
    bases_list = literal_eval(bases)
    for i, base in enumerate(bases_list):
        if base in base_counts:
            base_counts[base] += 1
            if i < len(bases_list) - 1 and base == 'C' and bases_list[i+1] == 'G':
                base_counts['CG'] += 1

# Convert to DataFrame for plotting
base_counts_df = pd.DataFrame.from_dict(base_counts, orient='index', columns=['Count'])

# Plot base composition
plt.figure(figsize=(10, 6))
sns.barplot(x=base_counts_df.index, y=base_counts_df['Count'])
plt.title('Base Composition in Dissimilar Regions')
plt.xlabel('Base')
plt.ylabel('Count')
plt.show()

# Perform Wilcoxon signed-rank test
# Example data: mismatch densities for two regions
# In practice, replace these with your actual data
genome_length = max(max(literal_eval(pos)) for pos in dissimilar_regions_df['Positions']) + 1
mismatch_density = [0] * genome_length

# Populate mismatch densities
for positions in dissimilar_regions_df['Positions']:
    for pos in literal_eval(positions):
        mismatch_density[pos] += 1

# Split mismatch densities into two halves for comparison
first_half_density = mismatch_density[:genome_length // 2]
second_half_density = mismatch_density[genome_length // 2:]

# Perform Wilcoxon signed-rank test
stat, p_value = wilcoxon(first_half_density, second_half_density)

print(f"Wilcoxon Signed-Rank Test: statistic = {stat}, p-value = {p_value}")

if p_value < 0.05:
    print("There is a significant difference between the two halves of the genome.")
else:
    print("There is no significant difference between the two halves of the genome.")
