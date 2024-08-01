#!/bin/bash
#SBATCH -J mm2_samtools_bedtools_FAD.s		# jobname
#SBATCH -o mm2_samtools_bedtools_FAD.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e mm2_samtools_bedtools_FAD.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 8:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

REF=hass.hifi.pooled.asm.p_ctg.masked.fasta
QRY=fatty_acid_genes.fasta
BASENAME=$(basename $QRY .fasta)_$(basename $REF .fasta)

# minimap2
PATH=/jet/home/bsanchez/bin/minimap2-2.24_x64-linux:$PATH
minimap2 -ax splice -k14 $REF $QRY > $BASENAME.sam

# samtools
module load samtools/1.13.0
samtools view -S -b $BASENAME.sam > $BASENAME.bam
samtools sort $BASENAME.bam > $BASENAME.sorted.bam
samtools index $BASENAME.sorted.bam
samtools faidx $REF

# bedtools
module load bedtools
bedtools bamtobed -i $BASENAME.sorted.bam > $BASENAME.sorted.bed

#_Files Produced_
echo _Files Produced_
# fasta file: 
echo Fasta File: $REF
# bam file:
echo BAM File: $BASENAME.sorted.bam
# bai file:  
echo BAM Index File: $BASENAME.sorted.bam.bai
# fai file: 
echo Fasta Index File: $REF.fai
# bed file: 
echo BED File: $BASENAME.sorted.bed

wait
