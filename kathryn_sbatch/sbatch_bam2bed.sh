#!/bin/bash
#SBATCH -J bam2bed              # jobname
#SBATCH -o bam2bed.o%A.%a       # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bam2bed.e%A.%a       # error file name A is the jobid and a is the arraytaskid
#SBATCH -a 13 #1-11                    # start and stop of the array start-end
#SBATCH -A mcb180013p
####SBATCH -N 1
#SBATCH -p RM-shared    #RM-512  #p1,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --ntasks-per-node=1
###SBATCH --mem-per-cpu=20000
#SBATCH -t 48:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=solarese@uci.edu
##SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes
#SBATCH --mail-type=fail         # email me when the job fails
#SBATCH --export=ALL

## SLURM ENV Variables
#  SLURM_CPUS_ON_NODE
#  SLURM_CPUS_PER_TASK
#  SLURM_ARRAY_TASK_ID

###bedtools 2.27.1
module load bedtools/2.30.0

JOBFILE=bamfiles.txt
#JOBFILE=newbams4reveiwersucounts.txt
BAM=$(head -n ${SLURM_ARRAY_TASK_ID} ${JOBFILE} | tail -n 1)

mkdir -p unsortedbeds
bedtools bamtobed -i $(basename ${BAM} .bam)_sorted.bam > unsortedbeds/$(basename ${BAM} .bam).bed
