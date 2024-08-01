#!/bin/bash
#SBATCH -J tiling_to_bed_UoQTHILLE		# jobname
#SBATCH -o tiling_to_bed_UoQTHILLE.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e tiling_to_bed_UoQTHILLE.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 8:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
export MV2_SHOW_CPU_BINDING=1

# Parsing show-tiling into bed formatted file
awk '/^>/{ref=substr($1,2); next} {start=$1; if(start < 0) start = -$1; print ref"\t"start-1"\t"$2"\t"$8"\t"$7}' Uo_thille.hifiv3.asm.p_ctg.fasta.out.delta.txt > thille_hass_tilting_output.bed
