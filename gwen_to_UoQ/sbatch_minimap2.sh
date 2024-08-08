#!/bin/bash
#SBATCH -J minimap2     # jobname
#SBATCH -o minimap2.o%A.%a   # output file name (%A is the jobid and %a is the arraytaskid)
#SBATCH -e minimap2.e%A.%a   # error file name (%A is the jobid and %a is the arraytaskid)
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
REF="GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"  # UoQ genome as the reference
QRY="gwen.hapsolo.contigs.primary.softmasked.fasta"  # Gwen genome as the query

OUTPUT_PREFIX="gwen_to_UoQ"                         # Prefix for output files
OUTPUT_PAF="${OUTPUT_PREFIX}.paf"                     # PAF output file name

# confirm files
echo "Ref: " $REF
echo "Query: " $QRY

# Run Minimap2
PATH=/jet/home/kchene/bin/minimap2-2.24_x64-linux:$PATH
minimap2 -cx asm5 -t 4 "$REF" "$QRY" > "$OUTPUT_PAF"
