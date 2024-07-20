#!/bin/bash
#SBATCH -J 5_sbatch_graph_cov		# jobname
###SBATCH -o graphCov.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e graphCov.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 13-30		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL


FILE=gwen_thille_cov_2,3_delete.txt

CONTIG=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' contigs.txt)
SMOOTHING=30000

mkdir -p graphs # make the folder if it doesn't exist
mkdir -p low_coverage # make the folder if it doesn't exist
touch low_coverage/median_coverages.txt

python3 graphCov.py $FILE $CONTIG $SMOOTHING
