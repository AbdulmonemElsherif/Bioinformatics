from Bio import AlignIO

# Read the consensus sequence from its file
with open("clustalo-I20240517-192453-0191-73736047-p1m.aln-clustal_num_consensus.txt", "r") as file:
    consensus = file.read().replace('\n', '')

# Read the alignment from the file
msa = AlignIO.read("clustalo-I20240517-201759-0434-51888817-p1m.aln-clustal_num", "clustal")


# Convert the MSA to a list of strings
msa_sequences = [str(record.seq) for record in msa]

# Ensure that the sequences in the MSA and the consensus sequence are of the same length
if all(len(seq) == len(consensus) for seq in msa_sequences):
    # Find the positions where the sequences in the MSA differ from the consensus sequence
    dissimilar_positions = {i: [seq[i] for seq in msa_sequences if seq[i] != consensus[i]] for i in range(len(consensus)) if any(seq[i] != consensus[i] for seq in msa_sequences)}
    print(dissimilar_positions)
else:
    print("The sequences in the MSA and the consensus sequence are not of the same length.")