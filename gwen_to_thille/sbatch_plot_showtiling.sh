#!/bin/bash
#SBATCH -J plot_tiling_GWEN_THILLE		# jobname
#SBATCH -o plot_tiling_GWEN_THILLE.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e plot_tiling_GWEN_THILLE.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 8:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=kac012@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

module load anaconda3/2022.10

tiling_file=nm_gwen.hapsolo.contigs.primary.softmasked.delta.txt

python3 plot_tiling.py $tiling_file
