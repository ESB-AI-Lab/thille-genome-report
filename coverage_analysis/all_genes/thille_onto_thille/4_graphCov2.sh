#!/bin/bash

#SBATCH -J 4_sbatch_graph_cov # jobname
###SBATCH -o graphCov.o%A.%a # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e graphCov.e%A.%a # error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-30 # start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared # queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00 # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin # email me when the job starts
#SBATCH --mail-type=end # email me when the job finishes
#SBATCH --export=ALL

JOB_FILE="jobfile_graph_data.txt"
FILE=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' $JOB_FILE)
MEDIAN=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' $JOB_FILE)
CONTIG=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $0}' contigs.txt)
SMOOTHING=30000

GRAPH_DATA_DIR="../thille_onto_thille/30k_graph_data/"

python3 graphCov2.py "$GRAPH_DATA_DIR$FILE" $CONTIG $SMOOTHING $MEDIAN
