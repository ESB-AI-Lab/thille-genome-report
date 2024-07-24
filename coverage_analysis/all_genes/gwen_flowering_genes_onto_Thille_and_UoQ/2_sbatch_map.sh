#!/bin/bash
#SBATCH -J map		# jobname
#SBATCH -o map.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e map.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-2		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

PATH=$(cat ../path.txt)$PATH

gene_name=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' jobfile_map.txt)
other_fasta=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' jobfile_map.txt)
other_name=$(awk -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $3}' jobfile_map.txt)
temp_name="temp_${gene_name}_onto_${other_name}"

minimap2 -ax splice $other_fasta temp_$gene_name.fasta > $temp_name.sam