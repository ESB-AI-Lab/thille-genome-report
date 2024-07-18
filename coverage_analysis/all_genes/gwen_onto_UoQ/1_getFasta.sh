#!/bin/bash
#SBATCH -J getFasta		# jobname
#SBATCH -o getFasta.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e getFasta.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

module load bedtools/2.30.0

gwen_fasta="gwen.scaffolds.20kNs.fasta"
gwen_gtf="gwen.scaffolds.20kNs.gene.annotation.genes.gtf"
temp_name="gwen_genes_UoQ"

awk -F "\t" '$3 == "gene"' $gwen_gtf > $temp_name.gtf
bedtools getfasta -fi $gwen_fasta -bed $temp_name.gtf -fo $temp_name.fasta
