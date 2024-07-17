#!/bin/bash
#SBATCH -J 2_sbatch_cum_sum		# jobname
###SBATCH -o step_2.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e step_2.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL



python3 cumulativesum3.py gwen.hapsolo.contigs.primary.lengths thille.hifiv3.asm.p_ctg.masked.lengths UoQ_sorted_contigLengths2.txt UoQ_sorted_scaffold_lengths.txt sorted_pame_lengths.txt
