Goal: Visualize Gwen genes mapped onto the hass genome from the University of Queensland (UoQ).

Input: 
    UoQ genome:         ../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta 
    Gwen assembly:      ../symlinks/gwen.scaffolds.20kNs.fasta
    Gwen annotation:    ../symlinks/gwen.scaffolds.20kNs.gene.annotation.genes.gtf
    
Output: 
    alignment:          gwen_genes_onto_UoQ.sorted.bam
    alignment index:    gwen_genes_onto_UoQ.sorted.bam.bai
    
Procedure:
    bash 1_get_genes_fasta.sh
    sbatch 2_sbatch_get_alignment_sam.sh (check status with squeue)
    bash 3_get_alignment_bam.sh

Comments:
Check temp files if things get broken.
Backup folder has working output.