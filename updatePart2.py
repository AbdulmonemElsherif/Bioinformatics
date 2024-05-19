from Bio import AlignIO
import pandas as pd
from scipy.stats import shapiro
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from scipy.stats import mannwhitneyu
import csv

omicron_alignment = AlignIO.read("Omicron(msa).aln-clustal_num", "clustal")
delta_alignment = AlignIO.read("DELTA(msa).aln-clustal_num", "clustal")


omicron_sequences = [str(record.seq) for record in omicron_alignment]
delta_sequences = [str(record.seq) for record in delta_alignment]


def calculate_sequence_stats(sequences, prefix):
    seq_stats = []
    for i, seq in enumerate(sequences):
        total_length = len(seq)
        c_percent = (seq.count("C") / total_length) * 100
        g_percent = (seq.count("G") / total_length) * 100
        t_percent = (seq.count("T") / total_length) * 100
        a_percent = (seq.count("A") / total_length) * 100
        cg_content = ((seq.count("C") + seq.count("G")) / total_length) * 100
        seq_stats.append({"Sequence Number": f"{prefix}{i+1}", "C Percentage": c_percent,
                          "G Percentage": g_percent, "T Percentage": t_percent,
                          "A Percentage": a_percent, "CG Content": cg_content})
    return seq_stats

omicron_stats = calculate_sequence_stats(omicron_sequences, "Omicron_seq")
delta_stats = calculate_sequence_stats(delta_sequences, "Delta_seq")

omicron_df = pd.DataFrame(omicron_stats)
delta_df = pd.DataFrame(delta_stats)


combined_df = pd.concat([omicron_df, delta_df], ignore_index=True)


additional_stats = combined_df.describe()

normality_results = {}
for column in combined_df.columns[1:]:  
    data = combined_df[column]
    # Shapiro-Wilk test
    statistic, p_value = shapiro(data)
    normality_results[column] = {"Shapiro-Wilk Statistic": statistic, "p-value": p_value}
   
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
   
    sns.histplot(data, kde=True, ax=axes[0])
    axes[0].set_title(f"Distribution of {column}")
    axes[0].set_xlabel(column)
    axes[0].set_ylabel("Frequency")
  
    stats.probplot(data, dist="norm", plot=axes[1])
    axes[1].set_title(f"Q-Q Plot of {column}")
    axes[1].set_xlabel("Theoretical Quantiles")
    axes[1].set_ylabel(f"Ordered Values of {column}")
    plt.tight_layout()
    plt.show()


combined_df.to_csv("sequence_stats.csv", index=False)
additional_stats.to_csv("additional_stats.csv")
normality_results_df = pd.DataFrame(normality_results).T
normality_results_df.to_csv("normality_test_results.csv")
# Perform non-parametric tests (Mann-Whitney U test) for comparing Omicron and Delta
# Null hypothesis: there's no significant difference between Omicron and Delta sequences.
# Alternative hypothesis: there's a significant difference between Omicron and Delta sequences.

# Dictionary to store test results
test_results_nonparametric = []

for column in ["C Percentage", "G Percentage", "T Percentage", "A Percentage", "CG Content"]:
    # Perform Mann-Whitney U test
    statistic, p_value = mannwhitneyu(omicron_df[column], delta_df[column], alternative='two-sided')
    
    # Store test results
    result = {"Variable": column, "Statistic": statistic, "p-value": p_value}
    
    # Determine hypothesis based on p-value
    if p_value < 0.05:
        result["Hypothesis"] = "Reject the null hypothesis: there's a significant difference."
    else:
        result["Hypothesis"] = "Fail to reject the null hypothesis: there's no significant difference."
    
    # Append test results to the list
    test_results_nonparametric.append(result)

# Save test results to CSV file
with open("test.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Variable", "Statistic", "p-value", "Hypothesis"])
    writer.writeheader()
    writer.writerows(test_results_nonparametric)