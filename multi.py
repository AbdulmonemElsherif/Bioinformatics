import subprocess

# Define the path to the Clustal Omega executable
clustal = r"D:\clustal-omega-1.2.2-win64\clustal-omega-1.2.2-win64\clustalo.exe"

# Define the command and arguments to run Clustal Omega, including --force to overwrite existing files
clustalomega = [
    clustal,
    "-i", r"D:\BioinformaticsProj\clustalo-I20240517-224204-0361-10233966-p1m(20seq).aln-clustal_num",
    "-o", r"D:\BioinformaticsProj\test.fasta",
    "--auto",
    "-v",
    "--force"
]

# Run the command using subprocess.run
subprocess.run(clustalomega, check=True)

# Print a completion message
print("Done")
