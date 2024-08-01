import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Function to load BLAST outfmt6 files and return SNP positions
def load_snp_positions(blast_outfmt6_file):
    snp_positions = {'snp_ref': [], 'snp_pos': []}
    with open(blast_outfmt6_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            snp_ref = parts[1]  # Reference ID
            snp_pos = int(parts[8])  # SNP position in reference
            snp_positions['snp_ref'].append(snp_ref)
            snp_positions['snp_pos'].append(snp_pos)
    return pd.DataFrame(snp_positions)

# Extracting command line arguments
REF = sys.argv[1]
QRY = sys.argv[2]
CONTIGS_FILE = sys.argv[3]
PAF_FILE = sys.argv[4]
Thille_SNP_FILE = sys.argv[5]

# Load SNP data
Thille_snps = load_snp_positions(Thille_SNP_FILE)

# Step 1: Read the PAF file and parse alignment positions and scores
positions = {'ref_id': [], 'REF_Start': [], 'REF_End': [], 'QRY_id': [], 'QRY_Start': [], 'QRY_End': [], 'score': [], 'num_matches': [], 'block_length': []}
with open(PAF_FILE, 'r') as paf:
    for line in paf:
        fields = line.strip().split()
        ref_id = fields[5]  # Reference sequence name
        ref_start = int(fields[7])  # 0-based start on the reference
        ref_end = int(fields[8])  # 0-based end on the reference
        qry_id = fields[0]  # Query sequence name
        qry_start = int(fields[2])  # 0-based start on the query
        qry_end = int(fields[3])  # 0-based end on the query
        score = int(fields[11]) if len(fields) > 11 else None  # Mapping quality (if present)
        num_matches = int(fields[9])
        block_length = int(fields[10])
        positions['ref_id'].append(ref_id)
        positions['REF_Start'].append(ref_start)
        positions['REF_End'].append(ref_end)
        positions['QRY_id'].append(qry_id)
        positions['QRY_Start'].append(qry_start)
        positions['QRY_End'].append(qry_end)
        positions['score'].append(score)
        positions['num_matches'].append(num_matches)
        positions['block_length'].append(block_length)

# Step 2: Read contigs data (adjusted for the new file format) and Create DataFrame from positions dictionary
contigs = pd.read_csv(CONTIGS_FILE, sep='\t', header=None, names=['Chromosome', 'length'])
df = pd.DataFrame(positions)

# Step 4: Merge contigs lengths with the alignment positions DataFrame
df = pd.merge(df, contigs, left_on='ref_id', right_on='Chromosome', how='left').drop('Chromosome', axis=1)

# Step 5: Sort DataFrame by contig length
df.sort_values(by='length', ascending=False, inplace=True)

# Print DataFrame head
pd.set_option('display.max_columns', None)  # Shows all columns
pd.set_option('display.width', 1000)  # Sets the display width for very wide rows
pd.set_option('display.max_colwidth', None)  # Shows the full content of each column
print(df.head())

# Create the plot
plt.figure(figsize=(10, 8))

# Determine the ranges for plot
max_ref_position = max(df['REF_End'])
max_qry_position = max(df['QRY_End'])

# Plotting with a fixed size that is smaller
fixed_dot_size = 20  # You can adjust this size to make the dots smaller
sns.scatterplot(data=df, x='REF_Start', y='QRY_Start', color='blue', s=fixed_dot_size, edgecolor='none')

# Add SNP markers as distinct points on the plot with the same smaller size
sns.scatterplot(x=[max_ref_position] * len(Thille_snps), y=Thille_snps['snp_pos'], color='orange', label='Thille SNPs', s=fixed_dot_size, edgecolor='none')

# Update titles as needed
plt.title('Genome Alignment Dot Plot with SNP Markers')
plt.xlabel('Thille Position')
plt.ylabel('Gwen Position')
plt.grid(True)

# Set axis limits to the max positions to reduce white space
plt.ylim(0, max_qry_position)

# Improve the legend to handle both dot plot points and SNP markers
plt.legend()

# Create the output folder if it doesn't exist and save the plot
output_folder = 'dot_plots'
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, 'Gwen_Thille_genome_alignment_dot_plot.png')
plt.savefig(output_file, dpi=300)  # Save the plot
plt.show()  # Show the plot
