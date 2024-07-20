!/bin/bash
#SBATCH -J minimap2		# jobname
#SBATCH -o minimap2.o%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e minimap2.e%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 48:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=gasimmons@ucsd.edu
#SBATCH --mail-type=begin	# email me when the job starts
#SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL


UoQ_fasta="GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
gwen_fasta="gwen.scaffolds.20kNs.fasta"
gwen_gtf="gwen.scaffolds.20kNs.gene.annotation.genes.gtf"
temp_name="gwen_genes_UoQ"

PATH=/jet/home/simmonsg/shared/avocado/gavin/bin/minimap2/minimap2-2.24_x64-linux:$PATH

minimap2 -ax splice $UoQ_fasta $temp_name.fasta > $temp_name.sam
