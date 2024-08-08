#!/bin/bash
#SBATCH -J hasscovcalc			# jobname
#SBATCH -o hasscovcalc.o%A.%a		# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e hasscovcalc.e%A.%a		# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1 #1-10				# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH -N 1
#SBATCH -p RM	#RM-512  #p1,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --ntasks-per-node=4
###SBATCH --mem-per-cpu=20000
#SBATCH -t 24:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=solarese@ucsd.edu
##SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes
#SBATCH --mail-type=fail         # email me when the job fails
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
#export MV2_SHOW_CPU_BINDING=1
CPUS=$(($SLURM_CPUS_ON_NODE/2))
#CPUS=$SLURM_CPUS_ON_NODE
#echo  $SLURM_CPUS_PER_TASK

echo ${CPUS}

###bedtools 2.30.0
module load bedtools/2.30.0

ALIGNJOBFILES=alignmentbeds.txt
#CONTIGFILE=hass.hifi.pooled.asm.p_ctg.contigs.bed
#CONTIGFILE=hass.hifi.pooled.asm.a_ctg.contigs.bed


ALIGNFILE=$(head -n ${SLURM_ARRAY_TASK_ID} ${ALIGNJOBFILES} | tail -n 1)

COVOUTFILE=$(basename ${CONTIGFILE} .bed)_$(basename ${ALIGNFILE} .sorted.bed).coverage

echo "Begin on ${MFEFILE} and ${ALIGNFILE}"
bedtools coverage -a ${CONTIGFILE} -b ${ALIGNFILE} -d -sorted > covfiles/${COVOUTFILE}
echo "Finished"
cd covfiles
[[ -f "${COVOUTFILE}.gz" ]] && rm ${COVOUTFILE}.gz
gzip ${COVOUTFILE}
