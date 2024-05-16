from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio import SeqIO
from collections import Counter
import os

# Define input and output directories
input_dir = "data/Delta"
output_dir = "data/SingleRef"

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Iterate over input FASTA files
for filename in os.listdir(input_dir):
    if filename.endswith(".fasta"):
        # Read sequences from input FASTA file
        sequences = []
        filepath = os.path.join(input_dir, filename)
        for record in SeqIO.parse(filepath, "fasta"):
            sequences.append(str(record.seq))

        # Get the maximum sequence length
        max_length = max(len(seq) for seq in sequences)

        # Pad sequences with gaps (-) to make them all the same length
        padded_sequences = [seq.ljust(max_length, '-') for seq in sequences]

        # Construct the consensus sequence
        consensus_sequence = ""
        for i in range(max_length):
            column = [sequence[i] for sequence in padded_sequences]
            most_common_base = Counter(column).most_common(1)[0][0]
            consensus_sequence += most_common_base

        # Save the consensus sequence to a new FASTA file
        output_filename = "Consensus_" + filename
        output_filepath = os.path.join(output_dir, output_filename)
        with open(output_filepath, "w") as output_file:
            output_file.write(">Consensus\n")
            output_file.write(consensus_sequence)

        print("Consensus sequence for", filename, "saved to:", output_filepath)
