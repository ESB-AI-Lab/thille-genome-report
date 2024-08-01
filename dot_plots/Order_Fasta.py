from Bio import SeqIO
import pandas as pd

# Load lengths into a DataFrame for sorting
lengths_file = "UoQ.sorted.lengths.txt"
lengths_df = pd.read_csv(lengths_file, sep=" ", header=None, names=["id", "length"])

# Sort by length in descending order
lengths_df = lengths_df.sort_values(by="length", ascending=False).drop_duplicates(subset=['id'])

# Read the FASTA file
fasta_file = "GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# Create a new sorted list of sequences
sorted_sequences = [sequences[id] for id in lengths_df['id'] if id in sequences]

# Write the sorted sequences to a new FASTA file
with open("GCA_029852735.1_ASM2985273v1_genomic.UoQ.sorted.fasta", "w") as output_handle:
    SeqIO.write(sorted_sequences, output_handle, "fasta")

print("FASTA file has been reordered and saved as 'GCA_029852735.1_ASM2985273v1_genomic.UoQ.sorted.fasta'")
