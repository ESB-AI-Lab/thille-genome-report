#!/bin/bash
#SBATCH -J 3_sbatch_coverage		# jobname
#SBATCH -o step_3.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e step_3.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
###SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
###SBATCH --mem-per-cpu=4000
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL

title=$1
bedtools_module=$2
module load ${bedtools_module}
temp="${title}/temp"

bedtools genomecov -ibam ${temp}.sorted.bam -d > ${title}/coverage.txt
#output: chromosome, position, coverage

cut -f1 ${title}/coverage.txt | sort | uniq > ${title}/features.txt

if [[ $? -eq 0 ]]
then
    touch success.txt
else
    touch failure.txt
fi
