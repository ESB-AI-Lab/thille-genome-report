#!/bin/bash
#SBATCH -J 6_sbatch_coverage_analysis		# jobname
#SBATCH -o step_6.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e step_6.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1-30		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

title=$1
bedtools_module=$2
gene_bed=$3

module load ${bedtools_module}

chromosome=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $1}' ${title}/median_coverages.txt)
median=$(awk -F "\t" -v id=${SLURM_ARRAY_TASK_ID} 'FNR == id {print $2}' ${title}/median_coverages.txt)

awk -v cm=${chromosome} '{ if($1 == cm) print $3 }' ${title}/coverage.txt > ${title}/temp_coverage_${chromosome}.txt
python3 coverage_analysis.py $title $chromosome $median

awk -v cm=${chromosome} '{ if($1 == cm) print }' ${gene_bed} > ${title}/temp_genes_${chromosome}.bed
bedtools intersect -a ${title}/${chromosome}.bed -b ${title}/temp_genes_${chromosome}.bed > ${title}/temp_intersection_${chromosome}.bed

touch complete_${SLURM_ARRAY_TASK_ID}.txt
