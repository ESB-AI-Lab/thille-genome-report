#!/bin/bash
#SBATCH -J tools		# jobname
#SBATCH -o tools.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e tools.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=kac012@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

module load samtools/1.13.0
module load bedtools/2.30.0

file_name="top1pct_onto_thille"
                  
samtools view -b $file_name.sam -o $file_name.bam
samtools sort $file_name.bam > $file_name.sorted.bam
bedtools bamtobed -i $file_name.sorted.bam > $file_name.bed

wait
