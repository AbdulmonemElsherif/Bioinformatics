from Bio import AlignIO
from collections import Counter
import pandas as pd

alignment_file = "referenceVSomicron.aln-clustal_num"

alignment = AlignIO.read(alignment_file, "clustal")


reference_seq = alignment[0].seq


def find_dominant_nucleotide(column):
    count = Counter(column)
    dominant_nucleotide = count.most_common(1)[0][0]
    return dominant_nucleotide

alignment_length = len(reference_seq)
comparison_results = []
dissimilar_positions = []

for i in range(alignment_length):
    column = [record.seq[i] for record in alignment]
    dominant_nucleotide = find_dominant_nucleotide(column[1:])  
    reference_nucleotide = reference_seq[i]
    if reference_nucleotide != dominant_nucleotide:
        dissimilar_positions.append(i)
    comparison_results.append((i, reference_nucleotide, dominant_nucleotide, reference_nucleotide != dominant_nucleotide))


def group_continuous_positions(positions):
    if not positions:
        return []
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


dissimilar_data = {
    "Region": [],
    "Positions": [],
    "Delta Consensus": [],
    "Omicron Sequences": [],
    "Mismatch Details": []
}

for start, end in grouped_positions:
    region_positions = list(range(start, end + 1))
    delta_bases = [reference_seq[pos] for pos in region_positions]
    omicron_bases = [[record.seq[pos] for record in alignment[1:]] for pos in region_positions]
    dominant_omicron_bases = [find_dominant_nucleotide(bases) for bases in omicron_bases]

    mismatches = []
    for i, pos in enumerate(region_positions):
        base_mismatches = [f"Seq {j+2}: {base}" for j, base in enumerate(omicron_bases[i]) if base != delta_bases[i]]
        if base_mismatches:
            mismatch_details = f"Position {pos}: Delta={delta_bases[i]}, Omicron=({', '.join(base_mismatches)})"
            mismatches.append(mismatch_details)

    dissimilar_data["Region"].append(f"{start}-{end}")
    dissimilar_data["Positions"].append(region_positions)
    dissimilar_data["Delta Consensus"].append(delta_bases)
    dissimilar_data["Omicron Sequences"].append(dominant_omicron_bases)
    dissimilar_data["Mismatch Details"].append("; ".join(mismatches))

dissimilar_regions_df = pd.DataFrame(dissimilar_data)
dissimilar_regions_df.to_csv("dissimilar_regions.csv", index=False)

print("Dissimilar regions saved to dissimilar_regions.csv")
