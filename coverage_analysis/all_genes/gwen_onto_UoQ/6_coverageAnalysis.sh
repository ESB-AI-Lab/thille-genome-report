#!/bin/bash
#SBATCH -J 6_sbatch_coverage_analysis		# jobname
#SBATCH -o covAnalysis.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e covAnalysis.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-12		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

title='coverage_analysis'
gene_bed='gwen_genes_UoQ.bed'

module load bedtools/2.30.0

chromosome=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' chromosomes.txt)
median=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' chromosomes.txt)

awk -v cm=${chromosome} '{ if($1 == cm) print $3 }' coverage.txt > temp_coverage_${chromosome}.txt
python3 cov_analysis.py $title $chromosome $median

awk -v cm=${chromosome} '{ if($1 == cm) print }' ${gene_bed} > temp_genes_${chromosome}.bed
bedtools intersect -a ${chromosome}.bed -b temp_genes_${chromosome}.bed > temp_intersection_${chromosome}.bed

touch complete_${SLURM_ARRAY_TASK_ID}.txt
