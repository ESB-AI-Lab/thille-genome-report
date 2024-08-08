#!/bin/bash
#SBATCH -J gznbam			# jobname
#SBATCH -o gznbam.o%A.%a		# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e gznbam.e%A.%a		# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 12-14 #1-10				# start and stop of the array start-end
#SBATCH -A mcb180013p
####SBATCH -N 1
#SBATCH -p RM-shared	#RM-512  #p1,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --ntasks-per-node=1
###SBATCH --mem-per-cpu=20000
#SBATCH -t 48:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=solarese@uci.edu
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

module load samtools/1.13.0

JOBFILE=bamfiles.txt
BAMFILE=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

samtools index -b $(basename ${BAMFILE} .bam)_sorted.bam $(basename ${BAMFILE} .bam)_sorted.bai
