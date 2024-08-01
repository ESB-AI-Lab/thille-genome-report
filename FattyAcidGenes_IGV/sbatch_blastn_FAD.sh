#!/bin/bash
#SBATCH -J blastn_FAD		# jobname
#SBATCH -o blastn_FAD.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e blastn_FAD.e%A.%a	# error file name A is the jobid and a is the arraytaskid
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

MYDB=hass.hifi.pooled.asm.p_ctg.masked.fasta
QUERY=fatty_acid_genes.fasta
MYLABEL=HASS_BLAST_DB
SLURM_CPUS_PER_TASK=1

# BLAST
module load BLAST

makeblastdb -in ${MYDB} -input_type fasta -dbtype nucl -parse_seqids -out $(basename ${MYDB} .fasta) -title ${MYLABEL}

blastn -query ${QUERY} -db $(basename ${MYDB} .fasta)  -max_target_seqs 10 -outfmt 6 -evalue 1e-5 -num_threads ${SLURM_CPUS_PER_TASK} > $(basename ${QUERY} .fasta)_$(basename ${MYDB} .fasta)_blastn_5.outfmt6
