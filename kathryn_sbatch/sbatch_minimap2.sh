#!/bin/bash
#SBATCH -J mm2.hass2thille                # jobname
#SBATCH -o mm2.hass2thille.o%A.%a         # jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e mm2.hass2thille.e%A.%a         # error file name A is the jobid and a is the arraytaskid
###SBATCH -c 64                    # start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH -N 1
#SBATCH -p RM	#RM-512  #p1,gcluster           # queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH --ntasks-per-node=128
###SBATCH --mem-per-cpu=20000
#SBATCH -t 18:00:00             # run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=b1sanchez@ucsd.edu
#SBATCH --mail-type=begin       # email me when the job starts
#SBATCH --mail-type=end         # email me when the job finishes
#SBATCH --export=ALL

# Set the number of threads per task(Default=1)
#export OMP_NUM_THREADS=1
#export MV2_SHOW_CPU_BINDING=1
#CPUS=$(($SLURM_CPUS_ON_NODE/2))
CPUS=$SLURM_CPUS_ON_NODE
#echo  $SLURM_CPUS_PER_TASK

echo ${CPUS}

##miniasm2.24 x64
PATH=/jet/home/bsanchez/bin/minimap2-2.24_x64-linux:$PATH

#REF=hass.hifi.pooled.asm.p_ctg.fasta
#REF=hass.hifi.pooled.asm.a_ctg.fasta
REF=thille.hifiv3.asm.p_ctg.fasta
#REF=thille.hifiv3.asm.a_ctg.fasta
variety=a048
variety=Eugenin
variety=Flavia
variety=mendez
variety=motherhass
variety=hass
variety=thille
#FASTQ=${variety}.pooled.hifi.fastq.gz
FASTQ=${variety}.hifi.fastq.gz

minimap2 -t ${CPUS} -ax map-hifi ${REF} ${FASTQ} > $(basename ${REF} .fasta)_$(basename ${FASTQ} .hifi.fastq.gz).hifi.aln.sam # PacBio HiFi/CCS genomic reads (v2.19 or later)
pigz -p ${CPUS} $(basename ${REF} .fasta)_$(basename ${FASTQ} .hifi.fastq.gz).hifi.aln.sam
#minimap2 -t ${CPUS} -ax asm20 ${REF} ${FASTQ} > $(basename ${REF} .fasta)_$(basename ${FASTQ} .hifi.fastq.gz).ccs.aln.sam    # PacBio HiFi/CCS genomic reads (v2.18 or earlier)
#pigz -p ${CPUS} $(basename ${REF} .fasta)_$(basename ${FASTQ} .hifi.fastq.gz).ccs.aln.sam
