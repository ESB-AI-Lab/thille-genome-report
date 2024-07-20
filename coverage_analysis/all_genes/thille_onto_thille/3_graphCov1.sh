#!/bin/bash
#SBATCH -J 3_sbatch_graphCoverage		# jobname
###SBATCH -o step_4.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
###SBATCH -e step_4.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 30-30		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL




CONTIG=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' thille.hifiv3.asm.p_ctg.masked.lengths)
SMOOTHING=30000

python3 graphCoverage.py thille_coverage_output.cov.txt $CONTIG $SMOOTHING
