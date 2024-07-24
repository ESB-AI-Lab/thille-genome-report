outdated

Goal: Find Locations of Cherry-Picked Gwen flowering genes mapped onto the hass genome from the University of Queensland (UoQ).

Input: 
    fasta of genome 
    gtf of gene
    fasta of the other genome
    
Output: 
    bed file of alignment
    
Procedure:
    sbatch 1_sbatch_genes.sh   <5m
    sbatch 2_sbatch_map.sh     <5m
    sbatch 3_sbatch_tools.sh   <5m

Comments:
    jobfile_genes.txt
        [genome file] [name of gene]
    jobfile_map.txt
        [name of gene] [other genome file] [name of other genome]

Genes:
    VOZ1 done
    CCR1
    EOBI
    SPA1 done
    KCS11
    AP2
    MTERF2
    NPR5
    CSU2 

Genomes I tried but tig not found:
    gwen.hapsolo.contigs.primary_new.fasta
    pamer.contigs.c21.consensus.consensus_pilon_pilon.gene.annotation.genes.avo_typeA_v_typeB.top1pct.20k.windowed.weir.fasta
    gwen.scaffolds.20kNs.fasta 
    
Last section of VOZ1
CM056820.1	46940787	46941030	tig00002445:264472-264715	60	+
Coincides with big drop in graph
Looking at coverage analysis, others are in regions of .2 to .5 median coverage