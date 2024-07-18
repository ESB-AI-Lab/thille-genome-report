#!/bin/bash
#SBATCH -J 4_sbatch_get_bed		# jobname
#SBATCH -o bedtools.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e bedtools.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

module load bedtools/2.30.0

file_name="gwen_genes_thille"

bedtools bamtobed -i $file_name.sorted.bam > $file_name.bed

wait
