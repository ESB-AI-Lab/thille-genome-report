#!/bin/bash
#SBATCH -J minimap2		# jobname
#SBATCH -o minimap2.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e minimap2.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 8:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1
#change_primary_group mcb190095p

#JOBFILE="jobfile.txt"

#SEED=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

PATH=$(cat ../path.txt)$PATH
uoq_fasta="../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"

minimap2 -ax splice $uoq_fasta temp.fasta > temp.sam

wait
