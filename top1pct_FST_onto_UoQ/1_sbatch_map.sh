#!/bin/bash
#SBATCH -J map		# jobname
#SBATCH -o map.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e map.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=kac012@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

PATH=/jet/home/kchene/bin/minimap2-2.24_x64-linux:$PATH

gwen_fasta=gwen_braker_rnaseq_hapsolo_annotation.genes.avo_typeA_v_typeB.top1pct.20k.windowed.weir.stranded.fasta
UoQ_fasta=GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta

minimap2 -ax splice $UoQ_fasta $gwen_fasta  > top1pct_onto_UoQ.sam
