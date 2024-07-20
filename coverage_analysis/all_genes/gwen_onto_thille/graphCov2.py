#!/usr/bin/env/python3

import sys

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

file = sys.argv[1]
chromosome = sys.argv[2]
smoothing = sys.argv[3]
opacity = 0.5
median=float(sys.argv[4])

# read x and y values from file
xy_data = pd.read_csv(file, sep='\t', header=None)
coverage_df = pd.DataFrame({
    'Position': xy_data.iloc[:, 0],
    'Coverage': xy_data.iloc[:, 1]
})


# read gene locations from file
genes = 'gwen_genes_thille.bed'
genes_df = pd.read_csv(genes, sep='\t', header=None, names=['Chromosome', 'Start'], usecols=[0, 1])
genes_df = genes_df.loc[genes_df['Chromosome'] == chromosome]
gene_x_values = genes_df['Start']
gene_y_values = [-20] * len(genes_df)

# Create a new column 'Color' where negative values are 'red' and non-negative values are 'blue'
#colors = np.where(coverage_df['Coverage'] == -10, 'red', 'blue')
coverage_df['Color'] = np.where(coverage_df['Coverage'] < 0, 'red', 'blue')

# Create a custom palette
custom_palette = {'red': sns.color_palette('colorblind')[3], 
                  'blue': sns.color_palette('colorblind')[0], 
                  'purple':sns.color_palette('colorblind')[4]}

# Scatterplot
sns.set_style('whitegrid')
plt.figure(figsize=(20, 4))
ax = sns.scatterplot(x='Position', y='Coverage', hue='Color', data=coverage_df, palette=custom_palette, s=15, linewidth=0, alpha=opacity)
ax.scatter(gene_x_values, gene_y_values, color=custom_palette['purple'], s=15, alpha=opacity)
ax.set_title(chromosome)

# format xticks and keep original range
xlim = plt.xlim()
ticksInfo = plt.xticks()
newLabels = []
for position in ticksInfo[0]:
    newLabels.append(str(int(position / 1000000)) + 'M')
plt.xticks(ticksInfo[0], newLabels)
plt.xlim(xlim[0], xlim[1])

# Adjust the y-axis limits to make space for the legend
y_min = -20
y_max = max(coverage_df['Coverage']) * 1.4  # Increase the maximum y-value
y_ticks = np.arange(y_min, y_max, 20)
ax.set_yticks(y_ticks)
ax.set_yticklabels([str(int(y)) for y in y_ticks])
ax.set_ylim(y_min, y_max)


# Add a line to the graph representing the median coverage
median_line_color = sns.color_palette('colorblind')[2]
plt.axhline(y=median, color=median_line_color, linestyle='solid')

# Set the y-axis limits
ax.set_ylim(-30, max(coverage_df['Coverage']) * 1.2)

# Add a legend to the graph
ax.legend(['Lower Three Quartiles', 'Upper Quartile'])
legend_elements = [
    Line2D([0], [0], marker='o', color='w', label='Below 3rd Quartile',
           markerfacecolor=sns.color_palette('colorblind')[0], markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Above 3rd Quartile',
           markerfacecolor=sns.color_palette('colorblind')[3], markersize=10),
    Line2D([0], [0], marker='o', color='w', label='Gene',
           markerfacecolor=sns.color_palette('colorblind')[4], markersize=10),
    Line2D([0], [0], color=sns.color_palette('colorblind')[2], lw=2, label='Median Coverage: ' + str(median))
]
ax.legend(handles=legend_elements, loc='upper right', ncol=4)

# save to graphs folder
name_folder = 'graphs/'
plt.savefig(name_folder + chromosome + '_gwen_genes.pdf', dpi=600)
