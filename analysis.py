import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2_contingency


dissimilar_regions_df = pd.read_csv("dissimilar_regions.csv")

dissimilar_regions_df['Region Length'] = dissimilar_regions_df['Positions'].apply(lambda x: len(eval(x)))


region_length_stats = dissimilar_regions_df['Region Length'].describe()


plt.figure(figsize=(10, 6))
sns.histplot(dissimilar_regions_df['Region Length'], bins=20, kde=True)
plt.title('Distribution of Dissimilar Region Lengths')
plt.xlabel('Region Length')
plt.ylabel('Frequency')
plt.show()





plt.figure(figsize=(10, 6))
sns.histplot(dissimilar_regions_df['Region Length'], bins=20, kde=True)
plt.yscale('log')  
plt.title('Distribution of Dissimilar Region Lengths')
plt.xlabel('Region Length')
plt.ylabel('Frequency (log scale)')
plt.show()



base_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0,'CG':0}

for bases in dissimilar_regions_df['Delta Consensus']:
    bases_list = eval(bases)
    for i, base in enumerate(bases_list):
        if base in base_counts:
           
            base_counts[base] += 1
            if i < len(bases_list) - 1 and base=='C':
                if bases_list[i+1] == 'G':
                    base_counts['CG'] += 1

        

base_counts_df = pd.DataFrame.from_dict(base_counts, orient='index', columns=['Count'])



plt.figure(figsize=(10, 6))
sns.barplot(x=base_counts_df.index, y=base_counts_df['Count'])
plt.title('Base Composition in Dissimilar Regions')
plt.xlabel('Base')
plt.ylabel('Count')
plt.show()

expected_base_composition = {'A': 0.30, 'T': 0.30, 'C': 0.20, 'G': 0.20, 'CG': 0.20}
total_bases = sum(base_counts.values())

expected_counts = {base: total_bases * proportion for base, proportion in expected_base_composition.items()}
observed_counts = list(base_counts.values())
expected_counts = list(expected_counts.values())

chi2, p, _, _ = chi2_contingency([observed_counts, expected_counts])

print(f"Chi-Square Test: chi2 = {chi2}, p-value = {p}")
