from Bio import AlignIO
import pandas as pd

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

# Save the combined DataFrame and additional statistics to a CSV file
combined_df.to_csv("sequence_stats.csv", index=False)
additional_stats.to_csv("additional_stats.csv")
