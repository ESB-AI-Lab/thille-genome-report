#!/bin/bash

#SBATCH -J 7_sbatch_graph_cov # jobname
###SBATCH -o graphCov2.o%A.%a # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e graphCov2.e%A.%a # error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-12 # start and stop of the array start-end
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

JOB_FILE=contigs.txt
FILE=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' $JOB_FILE)
MEDIAN=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' $JOB_FILE)
CONTIG=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $0}' contigs.txt)
SMOOTHING=30000

GRAPH_DATA_DIR="../../brianna/GwenReads_ThilleGenome/30k_graph_data/"

mkdir -p graphs # make the folder if it doesn't exist
mkdir -p low_coverage # make the folder if it doesn't exist

touch low_coverage/median_coverages.txt

GRAPH_FILE="${GRAPH_DATA_DIR}${FILE}_graph_data.txt"

python3 graphCov2.py "$GRAPH_FILE" $CONTIG $SMOOTHING $MEDIAN
