#!/usr/bin/env/python3

import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# read name of file from command line argument. 
file = sys.argv[1]
contig = sys.argv[2]
smoothing = int(sys.argv[3])
opacity = 0.5

# make pandas datas structure, tab separated, given file has no header, and label the 5 columns
coverage_df = pd.read_csv(file, sep='\t', header=None, names=['Feature', 'Start', 'End', 'Position', 'Coverage'])

# select positions corresponding to given contig
coverage_df = coverage_df[coverage_df['Feature'] == contig]

# trim away the Feature column
coverage_df = coverage_df.iloc[:, 1:]

# get quartiles of the 'Coverage' column
quartiles = coverage_df.iloc[:, 3].describe()

# set cutoff as third quartile
cutoff = quartiles[6]

# replace any coverage value greater than the cutoff with -10
coverage_df['Coverage'].mask(coverage_df['Coverage'].gt(cutoff), other=-10, inplace=True)

# Calculate the median of 'Coverage' excluding -10 values
median_coverage = coverage_df[coverage_df['Coverage'] != -10]['Coverage'].median() 

# SMOOTHING: makes a new column 'Window' giving every 30000 rows the same value.
coverage_df = coverage_df.assign(Window = lambda x: (coverage_df['Position'] - 1) // smoothing)

# SMOOTHING: median of 'Coverage' is calculated across all rows with the same 'Window' value 
coverage_df = coverage_df.groupby(['Window']).median()

# SMOOTHING: previous operation also changed 'Position' values. Overwrite it to go like 0, 1000, 2000... or however big the window is
coverage_df = coverage_df.assign(Position = lambda x: (coverage_df.index * smoothing))

# Create a new column 'Color' where negative values are 'red' and non-negative values are 'blue'
coverage_df['Color'] = np.where(coverage_df['Coverage']<0, 'red', 'blue')

# Create a custom palette
custom_palette = {'red': sns.color_palette('colorblind')[3], 'blue': sns.color_palette('colorblind')[0]}

# Minimum opacity is .002
sns.set_style('whitegrid')
plt.figure(figsize=(20,4))
ax = sns.scatterplot(x='Position', y='Coverage', hue='Color', data=coverage_df, palette=custom_palette, s=10, linewidth=0, alpha=opacity)
ax.set_title(contig)

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

# Add a legend to the graph
ax.legend(['Lower Three Quartiles', 'Upper Quartile'])

legend_elements = [Line2D([0], [0], marker='o', color='w', label='Below 3rd Quartile',
                          markerfacecolor=sns.color_palette('colorblind')[0], markersize=10),
                   Line2D([0], [0], marker='o', color='w', label='Above 3rd Quartile',
                          markerfacecolor=sns.color_palette('colorblind')[3], markersize=10),
                   Line2D([0], [0], color=sns.color_palette('colorblind')[2], lw=2, label='Median Coverage: ' + str(median_coverage))]

ax.legend(handles=legend_elements, loc='upper right', ncol=3)

# and save to graphs folder, cut off the .txt extension, add contig, add smoothing, png extension, 
str_smoothing = str(smoothing // 1000)
plt.savefig('graph_thille_onto_thille_30k_median_ex-10/' + file[0:-4]+ '_' + contig + '_smoothing' + str_smoothing + 'k.png', dpi=600)

# Use the contig, date, and time to create a unique filename
filename = f"30k_graph_data/{contig}_graph_data.txt"

# Save the 'Position' and 'Coverage' data to the file
with open(filename, 'w') as f:
    for index, row in coverage_df.iterrows():
        f.write(f"{row['Position']}\t{row['Coverage']}\n")
