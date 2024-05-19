from Bio import AlignIO
from Bio import SeqIO
import pandas as pd


omicron_alignment_file = "D:\\BioinformaticsProj\\Omicron(msa).aln-clustal_num"
delta_consensus_file = "D:\\BioinformaticsProj\\DELTA(msa)_consensus.aln-clustal_num"


omicron_alignment = AlignIO.read(omicron_alignment_file, "clustal")

with open(delta_consensus_file, "r") as file:
    delta_consensus = SeqIO.read(file, "fasta").seq


alignment_length = len(omicron_alignment[0].seq)
omicron_df = pd.DataFrame([[seq[i] for seq in omicron_alignment] for i in range(alignment_length)])
delta_consensus_df = pd.DataFrame([list(delta_consensus)])

# Identify the positions with mismatches
dissimilar_positions = []
for i in range(alignment_length):
    omicron_column = omicron_df.iloc[i]
    consensus_base = delta_consensus_df.iloc[0, i]
    if any(base != consensus_base for base in omicron_column):
        dissimilar_positions.append(i)

# Group continuous dissimilar positions into regions
def group_continuous_positions(positions):
    grouped_positions = []
    start = positions[0]
    end = positions[0]

    for pos in positions[1:]:
        if pos == end + 1:
            end = pos
        else:
            grouped_positions.append((start, end))
            start = pos
            end = pos
    grouped_positions.append((start, end))
    return grouped_positions

grouped_positions = group_continuous_positions(dissimilar_positions)

# Extract the dissimilar regions and prepare for saving
dissimilar_data = {
    "Region": [],
    "Positions": [],
    "Delta Consensus": [],
    "Omicron Sequences": [],
    "Mismatch Details": []
}

for start, end in grouped_positions:
    region_positions = list(range(start, end + 1))
    delta_bases = [delta_consensus_df.iloc[0, pos] for pos in region_positions]
    omicron_bases = [[omicron_df.iloc[pos, i] for i in range(len(omicron_alignment))] for pos in region_positions]
    
    mismatches = []
    for i, pos in enumerate(region_positions):
        base_mismatches = [f"Seq {j+1}: {base}" for j, base in enumerate(omicron_bases[i]) if base != delta_bases[i]]
        if base_mismatches:
            mismatch_details = f"Position {pos}: Delta={delta_bases[i]}, Omicron=({', '.join(base_mismatches)})"
            mismatches.append(mismatch_details)
    
    dissimilar_data["Region"].append(f"{start}-{end}")
    dissimilar_data["Positions"].append(region_positions)
    dissimilar_data["Delta Consensus"].append(delta_bases)
    dissimilar_data["Omicron Sequences"].append(omicron_bases)
    dissimilar_data["Mismatch Details"].append("; ".join(mismatches))


dissimilar_regions_df = pd.DataFrame(dissimilar_data)


dissimilar_regions_df.to_csv("dissimilar_regions.csv", index=False)

print("Dissimilar regions saved to dissimilar_regions.csv")
