from Bio import AlignIO
from Bio.Align import AlignInfo

# Path to the alignment file
alignment_file = "D:\\BioinformaticsProj\\DELTA(msa).aln-clustal_num"

# Read the MSA from the file
alignment = AlignIO.read(alignment_file, "clustal")

# Construct the summary info object
summary_align = AlignInfo.SummaryInfo(alignment)

# Generate the consensus sequence
consensus = summary_align.dumb_consensus()

# Save the consensus sequence to a file with the same extension
consensus_file = alignment_file.replace(".aln-clustal_num", "_consensus.aln-clustal_num")
with open(consensus_file, "w") as file:
    file.write(f">Consensus\n{str(consensus)}\n")

print(f"Consensus sequence saved to {consensus_file}")
