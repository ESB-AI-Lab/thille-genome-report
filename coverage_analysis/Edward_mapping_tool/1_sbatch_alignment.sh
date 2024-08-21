#!/bin/bash
#SBATCH -J 1_sbatch_alignment		# jobname
#SBATCH -o step_1.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e step_1.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

title=$1
reads_path=$2
reference_path=$3
minimap2_path=$4
minimap2_options=$5

PATH=${minimap2_path}:$PATH
threads=$SLURM_CPUS_ON_NODE
temp="${title}/temp"

echo "minimap2 -t ${threads} -ax ${minimap2_options} ${reference_path} ${reads_path} > ${temp}.sam"

minimap2 -t ${threads} -ax ${minimap2_options} ${reference_path} ${reads_path} > ${temp}.sam

if [[ $? -eq 0 ]]
then
    touch success.txt
else
    touch failure.txt
fi
