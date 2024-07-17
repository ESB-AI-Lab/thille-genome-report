import matplotlib.pyplot as plt
import sys
import os
import numpy as np

def read_data(file_path):
    ids = []
    lengths = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                id_str, length_str = line.strip().split()
                ids.append(id_str)
                lengths.append(int(length_str))
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    return ids, lengths

def plot_cumulative_sum(files, output_file):
    plt.figure(figsize=(20, 12))
    plt.xlabel('Number of Fragments (log)')
    plt.ylabel('Assembly Size (M)')
    plt.title('Avocado Assembly Cumulative Sum Graph')
    plt.tight_layout()

    legend_labels = []

    for file_path in files:
        ids, lengths = read_data(file_path)
        cumulative_sum = [sum(lengths[:i+1]) for i in range(len(lengths))]
        num_fragments = np.arange(1, len(ids) + 1)
        base_name = os.path.basename(file_path)

        if base_name == 'UoQ_sorted_contigLengths2.txt':
            legend_label = 'UoQ Contig Assembly'
        elif base_name == 'UoQ_sorted_scaffold_lengths.txt':
            legend_label = 'UoQ Scaffolded Assembly'
        elif base_name == 'sorted_pame_lengths.txt':
            legend_label = 'Pame Hass Scaffolded Assembly'
        elif base_name == 'gwen.hapsolo.contigs.primary.lengths':
            legend_label = 'Gwen Contig Assembly'
        elif base_name == 'thille.hifiv3.asm.p_ctg.masked.lengths':
            legend_label = 'Thille Contig Assembly'
        else:
            legend_label = base_name
        
        plt.plot(num_fragments, [x / 1000000 for x in cumulative_sum], label=legend_label)
        legend_labels.append(legend_label)

    plt.xscale('log')
    plt.grid(True)
    plt.legend()

    plt.axvline(x=12, linestyle='--', color='black')
        
    plt.savefig(output_file, dpi=600, bbox_inches='tight')

if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python script.py <file_path1> <file_path2> ... <file_pathN>")
        sys.exit(1)

    files = sys.argv[1:]

    current_dir = os.getcwd()
    output_file = os.path.join(current_dir, "cumulative_sum_graph3.pdf")

    plot_cumulative_sum(files, output_file)
