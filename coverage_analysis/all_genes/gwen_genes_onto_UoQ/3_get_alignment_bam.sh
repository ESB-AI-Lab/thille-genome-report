#!/bin/bash

PATH=$(cat ../path.txt)$PATH

mdoule load samtools/1.13.0

uoq_fasta="../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
gwen_fasta="../symlinks/gwen.scaffolds.20kNs.fasta"
gwen_gtf="../symlinks/gwen.scaffolds.20kNs.gene.annotation.genes.gtf"
                  
samtools view -b temp.sam -o temp.bam
samtools sort temp.bam > gwen_onto_UoQ.sorted.bam
samtools index gwen_onto_UoQ.sorted.bam
#bam to bed
if [ $? = 0 ]
then
    rm minimap2.*
fi