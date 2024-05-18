from Bio import AlignIO
from collections import Counter

# Read the MSA from the file
msa = AlignIO.read("clustalo-I20240517-224204-0361-10233966-p1m(20seq).aln-clustal_num", "clustal")

# Convert the MSA to a list of strings
msa_sequences = [str(record.seq) for record in msa]

# Calculate the frequency of each nucleotide and the CG content for each sequence in the MSA
nucleotide_freqs = [Counter(seq) for seq in msa_sequences]
nucleotide_percents = [{nucleotide: count / len(seq) * 100 for nucleotide, count in freq.items()} for seq, freq in zip(msa_sequences, nucleotide_freqs)]
cg_contents = [(freq.get('C', 0) + freq.get('G', 0)) / len(seq) * 100 for seq, freq in zip(msa_sequences, nucleotide_freqs)]

# Calculate the average percentage of each nucleotide and the average CG content
average_nucleotide_percent = {nucleotide: sum(percent[nucleotide] for percent in nucleotide_percents if nucleotide in percent) / len(nucleotide_percents) for nucleotide in ('A', 'C', 'G', 'T')}
average_cg_content = sum(cg_contents) / len(cg_contents)

print(average_nucleotide_percent)
print(f"Average CG content: {average_cg_content}%")