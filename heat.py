import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ast import literal_eval

dissimilar_regions_df = pd.read_csv("dissimilar_regions.csv")


genome_length = max(max(literal_eval(pos)) for pos in dissimilar_regions_df['Positions']) + 1
mismatch_density = [0] * genome_length


for positions in dissimilar_regions_df['Positions']:
    for pos in literal_eval(positions):
        mismatch_density[pos] += 1


heatmap_data = pd.DataFrame(mismatch_density, columns=["Mismatch Density"])


window_size = 100
resampled_data = heatmap_data["Mismatch Density"].groupby(heatmap_data.index // window_size).sum()
resampled_data_df = pd.DataFrame(resampled_data, columns=["Mismatch Density"])


plt.figure(figsize=(20, 5))
sns.heatmap(resampled_data_df.T, cmap="YlGnBu", cbar=True, xticklabels=window_size)
plt.title("Mismatch Density Across the Genome (Resampled)")
plt.xlabel(f"Genome Position (aggregated over {window_size} positions)")
plt.ylabel("Mismatch Density")
plt.show()
