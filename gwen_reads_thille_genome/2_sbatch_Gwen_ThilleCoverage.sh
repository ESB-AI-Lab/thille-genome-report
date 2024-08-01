#!/bin/bash
#SBATCH -J 2_Gwen_ThilleCoverage		# jobname
#SBATCH -o 2_Gwen_ThilleCoverage.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e 2_Gwen_ThilleCoverage.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH -p EM        # queue (partition) -- changed from RM to EM
### SBATCH --mem-per-cpu=1000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

REF=thille.hifiv3.asm.p_ctg.fasta
QRY=gwen_pacb_subreads.fastq.gz
BASENAME=$(basename $QRY .fasta)_$(basename $REF .fasta)

# Create Bed File
cat thille.hifiv3.asm.p_ctg.masked.lengths | awk '{print $1"\t1\t"$2-1}' > $BASENAME.bed  

# Bedtools extract coverages
module load bedtools/2.30.0

bedtools coverage -b $BASENAME.sorted.bam -a $BASENAME.bed  -d > gwen_thilleCoverage_output.cov.txt
