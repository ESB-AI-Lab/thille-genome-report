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

    for file_path in files:
        ids, lengths = read_data(file_path)
        cumulative_sum = [sum(lengths[:i+1]) for i in range(len(lengths))]
        num_fragments = np.arange(1, len(ids) + 1)
        plt.plot(num_fragments, [x / 1000000 for x in cumulative_sum], label=os.path.basename(file_path))  # Use plt.plot() instead of plt.step()

    plt.xscale('log')
    plt.grid(True)
    plt.legend()
    
    plt.axvline(x=12, linestyle='--', color='black')

    plt.savefig(output_file, dpi=600, bbox_inches='tight')
    
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python script.py <file_path1> <file_path2> ... <file_pathN>")
        sys.exit(1)

    files = sys.argv[1:]

    current_dir = os.getcwd()
    output_file = os.path.join(current_dir, "cumulative_sum_graph.pdf")

    plot_cumulative_sum(files, output_file)
