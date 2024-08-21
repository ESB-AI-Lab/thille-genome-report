#!/usr/bin/env/python3

import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# read name of file from command line argument
title = sys.argv[1]
chromosome = sys.argv[2]
smoothing = int(sys.argv[3])
genes = sys.argv[4]
opacity = 0.5
file = title + "/coverage.txt"

print(chromosome)

# make pandas datas structure, tab separated, given file has no header, and label first 3 columns
#coverage_df = pd.read_csv(file, sep='\t', header=None, names=['Feature', 'Position', 'Coverage'])

# for gavins coverage data
coverage_df = pd.read_csv(file, sep='\t', header=None, names=['Feature', 'x1', 'x2', 'Position', 'Coverage'], usecols=[0, 3, 4])

# select positions corresponding to given chromosome
coverage_df = coverage_df[coverage_df['Feature'] == chromosome]

# trim away the Feature column
coverage_df = coverage_df.iloc[:, 1:]

# get quartiles of the 'Coverage' column
quartiles = coverage_df.iloc[:, 1].describe()

# set cutoff as third quartile
cutoff = quartiles[6]

# replace any coverage value greater than the cutoff with -10
coverage_df['Coverage'].mask(coverage_df['Coverage'].gt(cutoff), other=-10, inplace=True)

# Calculate the median of 'Coverage' excluding -10 values
median_coverage = coverage_df[coverage_df['Coverage'] != -10]['Coverage'].median()

# SMOOTHING: makes a new column 'Window' giving every 1000 rows the same value.
coverage_df = coverage_df.assign(Window = lambda x: (coverage_df['Position'] - 1) // smoothing)

# SMOOTHING: median of 'Coverage' is calculated across all rows with the same 'Window' value 
coverage_df = coverage_df.groupby(['Window']).median()

# SMOOTHING: previous operation also changed 'Position' values. Overwrite it to go like 0, 1000, 2000... or however big the window is
coverage_df = coverage_df.assign(Position = lambda x: (coverage_df.index * smoothing))

# Create a new column 'Color' where negative values are 'red' and non-negative values are 'blue'
coverage_df['Color'] = np.where(coverage_df['Coverage']<0, 'red', 'blue')

# Create a custom palette
custom_palette = {'red': sns.color_palette('colorblind')[3], 'blue': sns.color_palette('colorblind')[0]}

# Scatterplot
sns.set_style('whitegrid')
plt.figure(figsize=(20,4))
ax = sns.scatterplot(x='Position', y='Coverage', hue='Color', data=coverage_df, palette=custom_palette, s=15, linewidth=0, alpha=opacity)
ax.set_title(chromosome)

# format xticks and keep original range
xlim = plt.xlim()
ticksInfo = plt.xticks()
newLabels = []
for position in ticksInfo[0]:
    newLabels.append(str(int(position / 1000000)) + 'M')
plt.xticks(ticksInfo[0], newLabels)
plt.xlim(xlim[0], xlim[1])

# increase vertical space by 10% for legend
ylim = plt.ylim()
plt.ylim(ylim[0], ylim[1] * 1.1)

# Add a line to the graph representing the median coverage
median_line_color = sns.color_palette('colorblind')[2]
plt.axhline(y=median_coverage, color=median_line_color, linestyle='solid')

# Vertical lines for genes
genes_df = pd.read_csv(genes, sep='\t', header=None, names=['Chromosome', 'Start'], usecols=[0, 1])
genes_df = genes_df.loc[genes_df['Chromosome'] == chromosome]
for i in genes_df['Start']:
    plt.axvline(x = i, color=sns.color_palette('colorblind')[4])

# Add a legend to the graph
ax.legend(['Lower Three Quartiles', 'Upper Quartile'])

legend_elements = [Line2D([0], [0], marker='o', color='w', label='Below 3rd Quartile',
                          markerfacecolor=sns.color_palette('colorblind')[0], markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='Above 3rd Quartile',
                          markerfacecolor=sns.color_palette('colorblind')[3], markersize=10),
                   Line2D([0], [0], color=sns.color_palette('colorblind')[4], lw=2, label='Gene'),
                   Line2D([0], [0], color=sns.color_palette('colorblind')[2], lw=2, label='Median Coverage: ' + str(median_coverage))]
                

ax.legend(handles=legend_elements, loc='upper right', ncol=4)

plt.savefig(title + "/" + chromosome + '.png', dpi=600)

# write out median coverages
with open(title + '/median_coverages.txt', 'a') as f:
    f.write(str(chromosome) + '\t' + str(median_coverage) + '\n')

# Use the contig, date, and time to create a unique filename
graph_data = title + '/' + chromosome + '_graph_data.tsv'

coverage_df.to_csv(path_or_buf=graph_data, sep='\t')
