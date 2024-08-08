import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Extracting command line arguments
REF = sys.argv[1]
QRY = sys.argv[2]
LENGTHS_FILE = sys.argv[3]
PAF_FILE = sys.argv[4]

# python minimap2_alignment_dot_plot.py thille.hifiv3.asm.p_ctg.softmasked.fasta gwen.hapsolo.contigs.primary.softmasked.fasta thille.hifiv3.asm.p_ctg.masked.lengths gwen_to_thille.paf

df = pd.read_csv('gwen_to_thille.paf', sep='\t', header=None, usecols=[0,1, 2, 3, 5, 6, 7, 8, 9, 10, 11], names=['QRY_id', 'QRY_len','QRY_Start', 'QRY_End', 'ref_id', 'REF_len','REF_Start', 'REF_End', 'num_matches', 'block_length', 'score'])
df.head()

# mydf[['REF_Start', 'REF_End', 'QRY_Start', 'QRY_End']]
# sns.scatterplot(data = mydf, x = 'REF_Start', y = 'QRY_Start')
# mydf[mydf['ref_id'] == 'ptg000010l']
