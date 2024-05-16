from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

sequences = []

fasta_files = [
                "./data/Delta/EPI_ISL_14594262.fasta", 
                "./data/Delta/EPI_ISL_14594268.fasta", 
                "./data/Delta/EPI_ISL_14594269.fasta", 
                "./data/Delta/EPI_ISL_14594270.fasta", 
                "./data/Delta/EPI_ISL_14594272.fasta", 
                "./data/Delta/EPI_ISL_17047775.fasta", 
                "./data/Delta/EPI_ISL_17047776.fasta", 
                "./data/Delta/EPI_ISL_17047777.fasta", 
                "./data/Delta/EPI_ISL_17047778.fasta", 
                "./data/Delta/EPI_ISL_17047779.fasta"
                ]
for file in fasta_files:
    with open(file, "r") as f:
        sequence = Seq(f.read().strip())  
        record = SeqRecord(sequence)
        record.id = file.split(".")[0] 
        sequences.append(record)
        print(f"Sequence {record.id} length: {len(sequence)}")

alignment = MultipleSeqAlignment(sequences)
summary_align = AlignInfo.SummaryInfo(alignment)
consensus = summary_align.dumb_consensus()

print("Consensus Sequence:")
print(consensus)
