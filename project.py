from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import motifs
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# List of alignment files
files = ["clustalo-I20240517-224204-0361-10233966-p1m(20seq).aln-clustal_num"]

for file in files:
    # Read the alignment file
    alignment = AlignIO.read(file, "clustal")

    # Convert the alignment to a list of SeqRecord objects, ensuring unique names
    sequences = [SeqRecord(Seq(str(record.seq)), id=f"{record.id}_{i}") for i, record in enumerate(alignment)]

    # Create a MultipleSeqAlignment object from the sequences
    msa = MultipleSeqAlignment(sequences)

    # Create a Motif object from the alignment
    motif = motifs.create([record.seq for record in msa])

    # Calculate the consensus sequence
    consensus = motif.consensus

    # Write the consensus sequence to a file
    with open(f"{file}_consensus.txt", "w") as output_file:
        output_file.write(str(consensus))

    # Print the consensus sequence
    print(consensus)

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)

    # Construct the phylogenetic tree using UPGMA algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)

    # Draw the phylogenetic tree
    Phylo.draw(tree)