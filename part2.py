from Bio import AlignIO
from collections import Counter
import pandas as pd

# Read the MSA from the file
msa = AlignIO.read("clustalo-I20240517-224204-0361-10233966-p1m(20seq).aln-clustal_num", "clustal")

# Convert the MSA to a list of strings
msa_sequences = [str(record.seq) for record in msa]

# Split the sequences into Omicron and Delta
omicron_sequences = msa_sequences[:10]
delta_sequences = msa_sequences[10:]

# Define a function to calculate the average nucleotide percentages and CG content for a list of sequences
def calculate_averages(sequences):
    # Calculate the frequency of each nucleotide and the CG content for each sequence
    nucleotide_freqs = [Counter(seq) for seq in sequences]
    nucleotide_percents = [{nucleotide: count / len(seq) * 100 for nucleotide, count in freq.items()} for seq, freq in zip(sequences, nucleotide_freqs)]
    cg_contents = [(freq.get('C', 0) + freq.get('G', 0)) / len(seq) * 100 for seq, freq in zip(sequences, nucleotide_freqs)]

    # Calculate the average percentage of each nucleotide and the average CG content
    average_nucleotide_percent = {nucleotide: sum(percent[nucleotide] for percent in nucleotide_percents if nucleotide in percent) / len(nucleotide_percents) for nucleotide in ('A', 'C', 'G', 'T')}
    average_cg_content = sum(cg_contents) / len(cg_contents)

    return average_nucleotide_percent, average_cg_content

# Calculate the averages for Omicron and Delta
omicron_averages = calculate_averages(omicron_sequences)
delta_averages = calculate_averages(delta_sequences)

# Create a DataFrame to store the averages
df = pd.DataFrame([omicron_averages, delta_averages], index=['Omicron', 'Delta'], columns=['Average Nucleotide Percent', 'Average CG Content'])

# Save the DataFrame to a CSV file
df.to_csv("averages.csv")