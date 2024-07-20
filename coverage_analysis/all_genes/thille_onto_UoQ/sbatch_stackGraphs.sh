#!/bin/bash
#SBATCH -J sbatch_stackGraphs		# jobname
###SBATCH -o sbatch_stack_graphs.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
###SBATCH -e sbatch_stack_graphs.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

python3 stackGraphs.py
