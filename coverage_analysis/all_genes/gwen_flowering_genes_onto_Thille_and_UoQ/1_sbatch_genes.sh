#!/bin/bash
#SBATCH -J genes		# jobname
#SBATCH -o genes.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e genes.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

module load bedtools/2.30.0

cat selected_genes.txt | xargs -I % grep -F % ../symlinks/gwen.scaffolds.20kNs.gene.annotation.gff3 | grep -i mrna > selected_genes.gtf

genome_fasta="../symlinks/gwen.scaffolds.20kNs.fasta"
gene_gtf="selected_genes.gtf"
gene_name="selected_genes"
temp_name="temp_$gene_name"

bedtools getfasta -fi $genome_fasta -bed $gene_gtf -fo $temp_name.fasta