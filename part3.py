from Bio import AlignIO

# Read the consensus sequence from its file
with open("clustalo-I20240517-192453-0191-73736047-p1m.aln-clustal_num_consensus.txt", "r") as file:
    consensus = file.read().replace('\n', '')

# Read the alignment from the file
msa = AlignIO.read("clustalo-I20240517-201759-0434-51888817-p1m.aln-clustal_num", "clustal")

# Convert the MSA to a list of strings
msa_sequences = [str(record.seq) for record in msa]

# Find the maximum length among the sequences in the MSA and the consensus sequence
max_length = max(max(len(seq) for seq in msa_sequences), len(consensus))

# Pad the sequences in the MSA and the consensus sequence with '-' to make them of the same length
msa_sequences = [seq.ljust(max_length, '-') for seq in msa_sequences]
consensus = consensus.ljust(max_length, '-')

# Find the positions where the sequences in the MSA differ from the consensus sequence
dissimilar_positions = {i: [seq[i] for seq in msa_sequences if seq[i] != consensus[i]] for i in range(len(consensus)) if any(seq[i] != consensus[i] for seq in msa_sequences)}

# Extract the dissimilar regions/columns as sequences
dissimilar_regions = [''.join(seq[i] for seq in msa_sequences if seq[i] != consensus[i]) for i in range(len(consensus)) if any(seq[i] != consensus[i] for seq in msa_sequences)]

# Save the dissimilar regions in a text file
with open("dissimilar_regions.txt", "w") as file:
    for region in dissimilar_regions:
        file.write(region + '\n')