from Bio import AlignIO

# Load the MSA file
alignment = AlignIO.read("referenceVSomicron.aln-clustal_num", "clustal")

# Identify the reference sequence
reference = alignment[0]

dissimilar_regions = []

# Compare each sequence with the reference
for record in alignment[1:]:
    for i, (ref_residue, target_residue) in enumerate(zip(reference.seq, record.seq)):
        # If residues are not the same, store the position
        if ref_residue != target_residue:
            dissimilar_regions.append(i)

# Print the dissimilar regions
print(dissimilar_regions)