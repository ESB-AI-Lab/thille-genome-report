#!/bin/bash
#SBATCH -J mumh2t_THILLE			# jobname
#SBATCH -o mumh2t_THILLE.o%A.%a		# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e mumh2t_THILLE.e%A.%a		# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1				# start and stop of the array start-end
#SBATCH -A mcb180013p
####SBATCH -N 1
#SBATCH -p EM        # queue (partition) -- changed from RM to EM
#SBATCH --ntasks-per-node=96
###SBATCH --mem-per-cpu=20000
#SBATCH -t 72:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
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

#Thille ASM to UoQ ASM

# REFERENCE: This sets the reference for the mummer plots and the prefix for the reference
REF="GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta" # masked UoQ genome

PREFIX="UoQ"

#Here we use ls return in order to iterate through each assembly. each assembly will be run independently
#SEED=$(ls iso1*.fasta | head -n $SGE_TASK_ID | tail -n 1)

# QUERY: sets the qry as the seed from above and modifies
QRY=thille.hifiv3.asm.p_ctg.fasta.out # soft masked thille genome
#FILEEXT="${QRY#*.}"
PREFIX=${PREFIX}_$(basename ${QRY} .fasta)

#Print to stdout the $QRY AND PREFIX
echo $QRY
echo $PREFIX

nucmer -l 100 -c 500 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
nucmer --maxmatch -prefix ${PREFIX}MM ${REF} ${QRY}

# mummerplot --fat --layout --filter -p ${PREFIX}MM ${PREFIX}MM.delta \
#   -R ${REF} -Q ${QRY} --postscript
   
#Change above if png or ps files are desired
#--postscript
#--png
