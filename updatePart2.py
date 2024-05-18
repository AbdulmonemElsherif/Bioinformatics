from Bio import AlignIO
import pandas as pd
from scipy.stats import ttest_rel

# Read the MSA from the file
omicron_alignment = AlignIO.read("Omicron(msa).aln-clustal_num", "clustal")
delta_sequences = AlignIO.read("DELTA(msa).aln-clustal_num", "clustal")
# Extract sequences from the alignment
omicron_sequences = [str(record.seq) for record in omicron_alignment]

# Define a function to calculate the percentage of each chemical constituent and CG content for each sequence
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

# Create DataFrames for Omicron and Delta
omicron_df = pd.DataFrame(omicron_stats)
delta_df = pd.DataFrame(delta_stats)

# Concatenate the DataFrames
combined_df = pd.concat([omicron_df, delta_df], ignore_index=True)

# Calculate additional statistics
additional_stats = combined_df.describe()

# Hypothesis Testing (Paired t-test) for comparing Omicron and Delta
# Assume null hypothesis: there's no significant difference between Omicron and Delta sequences.
# Alternative hypothesis: there's a significant difference between Omicron and Delta sequences.
# We'll perform a two-tailed t-test.

hypothesis_comments = {
    "C Percentage": "Null hypothesis: there's no significant difference between Omicron and Delta in terms of the percentage of C. Alternative hypothesis: there's a significant difference between Omicron and Delta in terms of the percentage of C.",
    "G Percentage": "Null hypothesis: there's no significant difference between Omicron and Delta in terms of the percentage of G. Alternative hypothesis: there's a significant difference between Omicron and Delta in terms of the percentage of G.",
    "T Percentage": "Null hypothesis: there's no significant difference between Omicron and Delta in terms of the percentage of T. Alternative hypothesis: there's a significant difference between Omicron and Delta in terms of the percentage of T.",
    "A Percentage": "Null hypothesis: there's no significant difference between Omicron and Delta in terms of the percentage of A. Alternative hypothesis: there's a significant difference between Omicron and Delta in terms of the percentage of A.",
    "CG Content": "Null hypothesis: there's no significant difference between Omicron and Delta in terms of the CG content. Alternative hypothesis: there's a significant difference between Omicron and Delta in terms of the CG content."
}

threshold = 0.05  

test_results = {}
for column in ["C Percentage", "G Percentage", "T Percentage", "A Percentage", "CG Content"]:
    test_result = ttest_rel(omicron_df[column], delta_df[column])
    if test_result.pvalue < threshold:
        hypothesis = "Alternative hypothesis"
    else:
        hypothesis = "Null hypothesis"
    test_results[column] = {"Statistic": test_result.statistic, "p-value": test_result.pvalue, "Hypothesis": hypothesis}

# Save the combined DataFrame, additional statistics, and test results to a CSV file
combined_df.to_csv("sequence_stats.csv", index=False)
additional_stats.to_csv("additional_stats.csv")

# Save the test results with comments to the CSV file
with open("test_results.csv", "a") as f:
    for column in test_results:
        f.write(f"{column}, {test_results[column]['Statistic']}, {test_results[column]['p-value']}, "
                f"{hypothesis_comments[column]}, {test_results[column]['Hypothesis']}\n")
