from Bio import AlignIO, SeqIO

def extract_dissimilar_regions(alignment_file, reference_sequence_file):
    with open(alignment_file, 'r') as f:
        data = f.read()  # Include the first line
    with open('temp_alignment_file.fasta', 'w') as f:
        f.write(data)

    # Check if the file is empty
    with open('temp_alignment_file.fasta', 'r') as file:
        if not file.read(1):
            print("File is empty")
            return []
        else:
            print("File is not empty")

    # Check if the file is in the correct FASTA format
    with open('temp_alignment_file.fasta', 'r') as file:
        first_line = file.readline().strip()
        if not first_line.startswith('>'):
            print("File is not in FASTA format")
            return []
        else:
            print("File is in FASTA format")

    alignment = AlignIO.read('temp_alignment_file.fasta', "fasta")
    reference_sequence_record = SeqIO.read(reference_sequence_file, "fasta")
    reference_sequence = str(reference_sequence_record.seq)
    dissimilar_regions = []

    for i in range(min(len(reference_sequence), min(len(record.seq) for record in alignment))):
        for record in alignment:
            if record.id != reference_sequence_record.id and record.seq[i] != reference_sequence[i]:
                dissimilar_regions.append((record.id, i, record.seq[i]))
                break

    return dissimilar_regions

alignment_file = "data\\Omicron\\EPI_ISL_17585694.fasta"
reference_sequence_file = "reference.fasta"
dissimilar_regions = extract_dissimilar_regions(alignment_file, reference_sequence_file)

for record_id, position, base in dissimilar_regions:
    print(f"Dissimilar region found in sequence {record_id} at position {position}: {base}")