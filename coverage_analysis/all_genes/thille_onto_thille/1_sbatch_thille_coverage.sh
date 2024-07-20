#!/bin/bash
#SBATCH -J ThilleCoverage		# jobname
#SBATCH -o ThilleCoverage.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e ThilleCoverage.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

CPUS=$SLURM_CPUS_ON_NODE
echo ${CPUS}

REF=thille.hifiv3.asm.p_ctg.fasta 
QRY=thille.hifiv3.fastq.gz
BASENAME=$(basename $QRY .fasta)_$(basename $REF .fasta)

# minimap2
PATH=/jet/home/simmonsg/shared/avocado/gavin/bin/minimap2/minimap2-2.24_x64-linux:$PATH
minimap2 -t ${CPUS} -ax map-hifi $REF $QRY > $BASENAME.sam

# samtools
module load samtools/1.13.0
samtools view -S -b $BASENAME.sam > $BASENAME.bam
    samtools sort $BASENAME.bam > $BASENAME.sorted.bam
    samtools index $BASENAME.sorted.bam
    samtools faidx $REF

wait
