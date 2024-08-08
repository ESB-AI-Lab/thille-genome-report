#!/bin/bash
#SBATCH -J minimap2_dot_plot     # jobname
#SBATCH -o minimap2_dot_plot.o%A.%a   # output file name (%A is the jobid and %a is the arraytaskid)
#SBATCH -e minimap2_dot_plot.e%A.%a   # error file name (%A is the jobid and %a is the arraytaskid)
#SBATCH -a 1-1                 
#SBATCH -A mcb180013p
####SBATCH -N 1
#SBATCH -p RM-shared       # queue (partition) -- changed from RM to EM
#SBATCH --ntasks-per-node=48
###SBATCH --mem-per-cpu=20000
#SBATCH -t 48:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=kac012@ucsd.edu
#SBATCH --mail-type=ALL          # Email notifications
#SBATCH --export=ALL

# Variables for plotting
REF="thille.hifiv3.asm.p_ctg.softmasked.fasta"  # Thille genome as the reference
QRY="gwen.hapsolo.contigs.primary.softmasked.fasta"  # Gwen genome as the query
LENGTHS_FILE="thille.hifiv3.asm.p_ctg.masked.lengths" # Reference lengths file

OUTPUT_PREFIX="gwen_to_thille"                         # Prefix for output files
OUTPUT_PAF="${OUTPUT_PREFIX}.paf"                     # SAM output file name

# Run Minimap2
# PATH=/jet/home/kchene/bin/minimap2-2.24_x64-linux:$PATH
#minimap2 -cx asm5 -t 4 "$REF" "$QRY" > "$OUTPUT_PAF"

# Execute the Python script for dot plot generation
/jet/packages/anaconda3-2022.10/bin/python3 minimap2_alignment_dot_plot.py "$REF" "$QRY" "$LENGTHS_FILE" "$OUTPUT_PAF"
