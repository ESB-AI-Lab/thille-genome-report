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
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=ALL          # Email notifications
#SBATCH --export=ALL

# Variables for Minimap2
REF="GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"  # UoQ genome as the reference
QRY="thille.hifiv3.asm.p_ctg.fasta"                   # Thille genome as the query
CONTIGS_FILE="UoQ.flitered.sorted.lengths.txt" # Reference lengths file
UoQ_SNP_FILE="vanessaMSL_UoQ_blastn_5.outfmt6"  # UoQ SNP file
Thille_SNP_FILE="vanessaMSL_UoQ_blastn_5.outfmt6"  # Thille SNP file

OUTPUT_PREFIX="thille_to_UoQ"                         # Prefix for output files
OUTPUT_PAF="${OUTPUT_PREFIX}.paf"                     # PAF output file name

# Run Minimap2
#PATH=/jet/home/bsanchez/bin/minimap2-2.24_x64-linux:$PATH
#minimap2 -cx asm5 -t 4 "$REF" "$QRY" > "$OUTPUT_PAF"

# Execute the Python script for dot plot generation
/jet/packages/anaconda3-2022.10/bin/python3 minimap2_alignment_dot_plot.py "$REF" "$QRY" "$CONTIGS_FILE" "$OUTPUT_PAF" "$UoQ_SNP_FILE" "$Thille_SNP_FILE"
