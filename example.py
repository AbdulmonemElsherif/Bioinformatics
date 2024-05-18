from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline

sequences = ["TGACTTCA", "TCACGTCA", "TGACGCCA", "TGACGTCA", "TGTCGCCA", "AGACGTCA"]

# Create SeqRecord objects for each sequence
seq_records = [SeqRecord(Seq(seq), id=f"Seq_{i}") for i, seq in enumerate(sequences)]

# Write the sequences to a FASTA file
SeqIO.write(seq_records, "sequences.fasta", "fasta")

# Run ClustalW to align the sequences
clustalw_cline = ClustalwCommandline("clustalw2", infile="sequences.fasta")
stdout, stderr = clustalw_cline()

# Read the alignment from the ClustalW output file
alignment = AlignIO.read("sequences.aln", "clustal")

# Calculate the consensus sequence
consensus = alignment.dumb_consensus()

# Print the consensus sequence
print(consensus)