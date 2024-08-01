#!/bin/bash
#SBATCH -J show-coords_UoQTHILLE			# jobname
#SBATCH -o show-coords_UoQTHILLE.o%A.%a		# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e show-coords_UoQTHILLE.e%A.%a		# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1				# start and stop of the array start-end
#SBATCH -A mcb180013p
####SBATCH -N 1
#SBATCH -p RM-shared	#RM-512  #p1,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --ntasks-per-node=64
###SBATCH --mem-per-cpu=20000
#SBATCH -t 24:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes
#SBATCH --mail-type=end         # email me when the job fails
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
#export MV2_SHOW_CPU_BINDING=1
#CPUS=$(($SLURM_CPUS_ON_NODE/2))
CPUS=$SLURM_CPUS_ON_NODE
#echo  $SLURM_CPUS_PER_TASK

echo ${CPUS}

PATH=/jet/home/bsanchez/bin/MUMmer3.23:$PATH

# Thille ASM to Hass ASM

# REFERENCE_QUERY:
REF_QRY=UoQ_thille.hifiv3.asm.p_ctg.fasta.out

# identify structural variants
show-coords -rcl ${REF_QRY}.delta > ${REF_QRY}.delta.coords

# -r	Sort output lines by reference
# -c	Include percent coverage columns in the output
# -l	Include sequence length columns in the output
