#!/bin/bash
#SBATCH -J intersect     # jobname
#SBATCH -o intersect.o%A.%a   # output file name (%A is the jobid and %a is the arraytaskid)
#SBATCH -e intersect.e%A.%a   # error file name (%A is the jobid and %a is the arraytaskid)
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

# Variables
GENES="thille_low_cov_genes.bed"  # genes intersected against low coverage of GWEN reads mapped onto UoQ
REF="thille_to_UoQ"  # gwen to UoQ 


module load bedtools/2.30.0

# -wo: write the original A and B entries plus the number of base pairs of overlap between the two features. Only A features with overlap are reported. Restricted by -f and -r.
bedtools intersect -wo -a $GENES -b ${REF}.bed > ${REF}_intersect.txt
