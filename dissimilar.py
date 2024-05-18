from Bio import AlignIO, SeqIO
import tkinter as tk
from tkinter import scrolledtext

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

def display_dissimilar_regions(dissimilar_regions, reference_sequence, selected_sequence):
    print("Displaying dissimilar regions...")

    root = tk.Tk()

    ref_text_box = scrolledtext.ScrolledText(root, width=200, height=20)
    sel_text_box = scrolledtext.ScrolledText(root, width=200, height=20)

    ref_text_box.pack()
    sel_text_box.pack()

    dissimilar_positions = {position for _, position, _ in dissimilar_regions}

    for i in range(min(1000, len(reference_sequence), len(selected_sequence))):
        ref_color = 'black'
        sel_color = 'red' if i in dissimilar_positions else 'black'
        ref_text_box.insert("insert", reference_sequence[i], ref_color)
        sel_text_box.insert("insert", selected_sequence[i], sel_color)

    ref_text_box.tag_config('red', foreground='red')
    sel_text_box.tag_config('red', foreground='red')
    ref_text_box.tag_config('black', foreground='black')
    sel_text_box.tag_config('black', foreground='black')

    root.mainloop()
    
alignment_file = "data\\Omicron\\EPI_ISL_17585694.fasta"
reference_sequence_file = "reference.fasta"
dissimilar_regions = extract_dissimilar_regions(alignment_file, reference_sequence_file)
reference_sequence_record = SeqIO.read(reference_sequence_file, "fasta")
reference_sequence = str(reference_sequence_record.seq)

# Assuming the selected sequence is the first sequence in the alignment file
selected_sequence_record = SeqIO.read(alignment_file, "fasta")
selected_sequence = str(selected_sequence_record.seq)

print(f"Reference sequence: {reference_sequence[:100]}...")  # Print the first 100 characters
print(f"Selected sequence: {selected_sequence[:100]}...")  # Print the first 100 characters
print(f"Dissimilar regions: {dissimilar_regions[:10]}...")  # Print the first 10 dissimilar regions

display_dissimilar_regions(dissimilar_regions, reference_sequence, selected_sequence)