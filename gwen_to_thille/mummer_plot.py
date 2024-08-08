import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

#contig to contig dot plot
#thille sorted lengths file as REF
#alignment, QRY ID, start, stop, ref ID, start, stop (may not have ref start/stop, need to find other way if not)
#use QRY ID start/stop and Ref ID start/stop
#sort by size
#largest contig closest to x=0, smallest at largest x
#sort Ref contigs on x-axis from largest length to smallest
#x-values: Ref ID start/stop
#y-values: QRY ID start/stop
#connect lines from Ref/QRY ID start to Ref ID stop
#flipped if max/min are switched

# mummer_dotplot.py thille.hifiv3.asm.p_ctg.softmasked.fasta gwen.hapsolo.contigs.primary.softmasked.fasta mummer_dotplot.py nm_gwen.hapsolo.contigs.primary.softmasked.delta

# Extracting command line arguments
# REF = sys.argv[1]
# QRY = sys.argv[2]
mummer_delta = sys.argv[1]

# generated using nucmer -l 100 -c 500 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
# -l, -c, -d, -banded, -D options

# ref qry 
# >ptg000149l tig00004821 224091 268589
# The four coordinates are the start and end in the reference and the start and end in the query respectively.
# start and end in ref, start and end in query
# next alignment is preceded by a 0
# 1 6399 67295 73693 11 11 0
# -5918
# 242
# 0
# 94406 96289 4362 6249 10 10 0
# 1189
# 1

# Step 1: Read the delta file and parse alignment positions, store lengths of reference/query       
with open(mummer_delta, 'r') as delta:
    lines = delta.readlines()
    
    alignments = []
    alignment = {}
    parsing_alignment = False

    for line in lines:
        if line.startswith('>'): # append previous alignment's info at header lines
            # header_info = line.strip().split()
            # ref = header_info[0].split('>')[1]
            # query = header_info[1]
            # print(ref)
            # print(query)
            # print(parsing_alignment)
            # if alignment:
            #     print('appending new alignment')
            #     alignments.append(alignment)
            # print(alignment)
            alignment = {}
            # parsing_alignment = False
            continue

        if parsing_alignment:
            # print("parsing")
            # 1 6399 67295 73693 11 11 0
            # first 2: ref start/stop, second 2: query start/stop
            # three digits following locations: number of errors (non-identities + indels), similarity errors (non-positive match scores), stop codons (does not apply to DNA alignments, will be "0")
            line_data = line.strip().split()
            if len(line_data) == 7:
                # print("info line")
                # alignment['%ID'] = float(line_data[3])
                alignment['Reference_Length'] = int(line_data[1]) - int(line_data[0])
                alignment['Reference_Start'] = int(line_data[0])
                alignment['Reference_End'] = int(line_data[1])
                alignment['Query_Length'] = int(line_data[3]) - int(line_data[2])
                alignment['Query_Start'] = int(line_data[2])
                alignment['Query_End'] = int(line_data[3])
                # print(alignment)
                alignments.append(alignment)
                alignment = {}
                # print(alignments)
                # print(alignment)
                
        elif line.startswith('NUCMER'):
            parsing_alignment = True
            
    # if alignment:
    #     print('appending at end')
    #     alignments.append(alignment)

# print(alignments)  

# Step 2: Create DataFrame from alignment positions
df = pd.DataFrame(alignments)
new_df = pd.DataFrame()

# try merging reference position cols and query position cols???
ref_start = list(df.get('Reference_Start'))
ref_end = list(df.get('Reference_End'))
all_ref = ref_start + ref_end
new_df = new_df.assign(Reference = all_ref)
qry_start = list(df.get('Query_Start'))
qry_end = list(df.get('Query_End'))
all_qry = qry_start + qry_end
new_df = new_df.assign(Query = all_qry)
# df['Query'] = df['Query_Start'] + df['Query_End']
# df['Reference'] = pd.concat([df['Reference_Start'], df['Reference_End']])
# print(str(df['Reference'].shape[0]))
# print(str(df['Query'].shape[0]))
# print(df.shape[0])

# Step 3: Create the plot

# Calculate the length of the contigs and sort the dataframe by it
# df['Reference_Length'] = abs(df['Reference_End'] - df['Reference_Start'])
# df.sort_values(by='Reference_Length', ascending=False, inplace=True)

plt.figure(figsize=(10, 6))
# need to plot the pairs of ref start/end AND query start/end
# wait i'm confusing myself so ref/query start are the same just on different axes??? and then ref/query end are also the same but on different axes???
# ax = sns.scatterplot(data=df, x="Reference_Start", y="Query_Start")
ax = sns.scatterplot(data=new_df, x="Reference", y="Query")
ax.set_title('Dot Plot of Gwen to Thille Alignments')
ax.set_xlabel('Reference Position (Thille)')
ax.set_ylabel('Query Position (Gwen)')

print('plot created')
plt.show()

# Step 5: Create the folder if it doesn't exist
output_folder = 'dot_plots'
os.makedirs(output_folder, exist_ok=True)

# Step 6: Save the plot as PNG in the specified folder
output_file = os.path.join(output_folder, f"gwen_onto_thille_dot_plot.png")
plt.savefig(output_file)

# Step 7: Show the plot
# plt.show()
