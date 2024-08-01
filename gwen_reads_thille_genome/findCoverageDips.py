#!/usr/bin/env/python3

import sys
import pandas as pd
import numpy as np
import os

# read name of file from command line argument
file = sys.argv[1]
chromosome = sys.argv[2]

# make pandas data structure, tab separated, given file has no header, and label first 5 columns
coverage_df = pd.read_csv(file, sep='\t', header=None, names=['Feature', 'Start', 'End', 'Position', 'Coverage'])

# select positions corresponding to given chromosome
coverage_df = coverage_df[coverage_df['Feature'] == chromosome]

# trim away the Feature column
coverage_df = coverage_df.iloc[:, 1:]

# Calculate the rolling median of 'Coverage' with a window size of 100, 500, 1000
coverage_df['RollingMedian'] = coverage_df['Coverage'].shift().rolling(window=1000).median()

# Create a new column 'Drop' where coverage drops are marked as True if it's less than 0.5 or 0.2 of the rolling median
coverage_df['Drop'] = (coverage_df['Coverage'] < 0.5 * coverage_df['RollingMedian']) | (coverage_df['Coverage'] < 0.2 * coverage_df['RollingMedian'])

# Assigns a value to 'DropValue' in the 'coverage_df' DataFrame based on the condition of 'Coverage' relative to 'RollingMedian'.
coverage_df['DropValue'] = np.where(coverage_df['Coverage'] < 0.2 * coverage_df['RollingMedian'], 0, np.where(coverage_df['Coverage'] < 0.5 * coverage_df['RollingMedian'], 1, 0))

# Specify the directory for the output files
output_dir = 'coverage_dips_1k'

# Write the positions of coverage drops to a new text file
drop_positions = coverage_df[coverage_df['Drop']]['Position'].values
drop_values = coverage_df[coverage_df['Drop']]['DropValue'].values

# Create ranges of consecutive positions
ranges = []
if len(drop_positions) > 0:
    start = end = drop_positions[0]
    drop_value = drop_values[0]
    for pos, value in zip(drop_positions[1:], drop_values[1:]):
        if pos == end + 1:
            end = pos
            drop_value = value
        else:
            ranges.append((start, end, drop_value))
            start = end = pos
            drop_value = value
    ranges.append((start, end, drop_value))

with open(os.path.join(output_dir, chromosome + 'drop_locations.w1k.txt'), 'w') as f:
    
    f.write(f'Chromosome: {chromosome}\n')
    f.write('0: Coverage < 0.2 * Rolling Median\n')
    f.write('1: 0.2 * Rolling Median <= Coverage < 0.5 * Rolling Median\n')
    f.write('Start\tEnd\tDropValue\n\n')
    
    if len(drop_positions) == 0:
        f.write("No drop positions found.\n")
    else:
        for start, end, value in ranges:
            f.write(f'{start}\t{end}\t{value}\n')
