from collections import Counter

# Read the consensus sequence from the file
with open("clustalo-I20240517-224204-0361-10233966-p1m(20seq).aln-clustal_num_consensus.txt", "r") as file:
    consensus = file.read().replace('\n', '')

# Calculate the frequency of each nucleotide in the consensus sequence
nucleotide_freq = Counter(consensus)

# Calculate the percentage of each nucleotide
nucleotide_percent = {nucleotide: count / len(consensus) * 100 for nucleotide, count in nucleotide_freq.items()}

# Calculate the CG content
cg_content = (nucleotide_freq.get('C', 0) + nucleotide_freq.get('G', 0)) / len(consensus) * 100

print(nucleotide_freq)
print(nucleotide_percent)
print(f"CG content: {cg_content}%")