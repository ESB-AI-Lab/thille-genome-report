#!/bin/bash
#SBATCH -J sortbed              # jobname
#SBATCH -o sortbed.o%A.%a       # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e sortbed.e%A.%a       # error file name A is the jobid and a is the arraytaskid
#SBATCH -a 13-15%1                    # start and stop of the array start-end
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
BEDFILE=$(basename ${BAM} .bam).bed

#sort -k1,1 -k2,2n -T /gpool/tmp/ unsortedbeds/${BEDFILE} > sortedbeds/$(basename ${BEDFILE} .bed).sorted.bed
sort -k1,1 -k2,2n -T /tmp/ ${BEDFILE} > $(dirname ${BEDFILE})/$(basename ${BEDFILE} .bed).sorted.bed
