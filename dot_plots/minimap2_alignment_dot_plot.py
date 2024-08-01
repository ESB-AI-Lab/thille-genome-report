import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Step 1: Load positions and alignment scores from PAF file format
paf_data = pd.read_csv(sys.argv[4], sep='\t', header=None, usecols=[0, 2, 3, 5, 7, 8, 9, 10, 11], 
                        names=['query_id', 'query_start', 'query_end', 'reference_id', 'ref_start', 'ref_end', 'num_matches', 'block_length', 'score'])

# Load contig lengths
contig_sizes_df = pd.read_csv(sys.argv[3], sep='\t', header=None, names=['Chromosome', 'length'])
contig_sizes_series = contig_sizes_df.set_index('Chromosome')['length']

# Create a DataFrame for the PAF data
paf_df = pd.DataFrame(paf_data)

# Extract query IDs from PAF data
paf_query_ids = set(paf_df['query_id'])

# Extract contig IDs from contig sizes file
contig_ids = set(contig_sizes_df['Chromosome'])

# Filter PAF data to include only alignments with query IDs present in contig sizes file
valid_query_ids = contig_ids.intersection(paf_query_ids)
paf_df_filtered = paf_df[paf_df['query_id'].isin(valid_query_ids)]
print(paf_df_filtered)
# Recalculate query size after filtering
paf_df_filtered['query_size'] = paf_df_filtered['query_id'].map(contig_sizes_series)

# Sort the filtered PAF data by the recalculated query size
paf_df_sorted = paf_df_filtered.sort_values(by='query_size', ascending=False)

# Now we can plot the sorted data
fig, ax = plt.subplots(figsize=(10, 8))

# Scatter plot for dot plot
# X-axis: reference position, Y-axis: query position (midpoint for simplicity)
# Size and color can be adjusted for better visibility
ax.scatter(
    x=(paf_df_sorted['ref_start'] + paf_df_sorted['ref_end']) / 2,
    y=(paf_df_sorted['query_start'] + paf_df_sorted['query_end']) / 2,
    s=10,  # size of points
    c='blue',  # color of points
)

# Set labels and title
plt.title('Genome Alignment Dot Plot with SNP Markers (Filtered)')
plt.xlabel('Cumulative Genome Position (bp)')
plt.ylabel('Query Position (bp)')
plt.grid(True)
plt.legend()

# Save the plot
output_folder = 'dot_plots'
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, 'genome_alignment_dot_plot_filtered.png')
plt.savefig(output_file, dpi=300)
plt.show()
