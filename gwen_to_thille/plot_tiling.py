#!/usr/bin/env/python3

import sys
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

input_file = sys.argv[1] # show-tiling output file (or bed file?)

# Scatter plot
# Line plot from start to stop
# ex:
# Gwen is x-axis
# Start and stop of query is on x-axis

# Y-axis is thille
# Start and stop of reference is on y-axis
# Diagonal

# >ptg000001l 82066437 bases
# ref           start in ref  stop in ref        query 
# ptg000001l	18579554	  18900465	tig00021315	+

# need the start/stop of query

tiling_df = pd.read_csv(input_file, sep='\t', header=None, names=['Reference', 'Start', 'Stop', 'Query', 'Alignment'])

# 'query' is x-axis
# 'ref' is y-axis

sns.set_style('whitegrid')
plt.figure(figsize=(20,4))
ax = sns.scatterplot(x='', y='', data=tiling_df, s=15, linewidth=0, alpha=opacity)
ax.set_title('')
