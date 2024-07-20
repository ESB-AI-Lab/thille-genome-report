#!/bin/bash
#SBATCH -J 6_sbatch_coverage_analysis		# jobname
###SBATCH -o step_6.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
###SBATCH -e step_6.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-30		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

CHROMOSOME=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' jobfile_chromosomes_medians.txt)
MEDIAN=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' jobfile_chromosomes_medians.txt)

python3 findCoverageDips.py $CHROMOSOME $MEDIAN

FILE="coverage_analysis/coverage_analysis_gwen_reads_onto_UoQ_${CHROMOSOME}.bed"
sort -k3nr -k2 -o $FILE $FILE
