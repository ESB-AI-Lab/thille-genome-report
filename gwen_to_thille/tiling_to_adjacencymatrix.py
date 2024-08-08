#!/usr/bin/env/python3

import sys
import os
import pandas as pd
import numpy as np

input_file = sys.argv[1] # show-tiling output file

# create lists
refs = []
queries = []
connections = []

# parse the input file
with open(input_file, 'r') as file:
    for line in file:
        if line.startswith('>'):
            # get the reference sequence name, remove '>' and split by space, add to list of ref sequences
            feature = line.split()[0][1:]
            refs.append(feature)
        else:
            # get query information, add contig ID to list of queries, add connection pair
            query_info = line.strip().split()
            contig_id = query_info[-1]
            queries.append(contig_id)
            connections.append((feature, contig_id))

    # remove duplicates
    refs = list(set(refs))
    queries = list(set(queries))

    # create adjacency matrix
    adj_matrix = pd.DataFrame(0, index=refs, columns=queries)

    # update connections in adjacency matrix
    for ref, query in connections:
        adj_matrix.at[ref, query] = 1

    output_file = 'adjacency_matrix.txt'

    cwd = os.getcwd()

    # write the adjacency matrix to a text file
    adj_matrix.to_csv(os.path.join(cwd, output_file))

    print(f"Adjacency matrix saved to {os.path.join(cwd, output_file)}")
