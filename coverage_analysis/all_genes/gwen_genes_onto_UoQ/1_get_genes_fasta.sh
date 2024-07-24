#!/bin/bash


PATH=$(cat ../path.txt)$PATH

module load samtools/1.13.0

uoq_fasta="../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
gwen_fasta="../symlinks/gwen.scaffolds.20kNs.fasta"
gwen_gtf="../symlinks/gwen.scaffolds.20kNs.gene.annotation.genes.gtf"
                    
awk -F "\t" '$3 == "gene"' $gwen_gtf > temp.gtf
wc -l $gwen_gtf | awk '{ print $1 }' | xargs -I {} echo "there are" {} "annotations"
wc -l temp.gtf | awk '{ print $1 }' | xargs -I {} echo "there are" {} "genes"
bedtools getfasta -fi $gwen_fasta -bed temp.gtf -fo temp.fasta