Goal: Map Thille ASM onto UoQ ASM. Gather Titlting Paths and Structural Variants.

Input:
    Reference genome:             ../GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta 
    Query genome:                ../thille.hifiv3.asm.p_ctg.fasta.out
    
Procedures: 

Using minimap2:
    sbatch sbatch_minimap2_dot_plot.sh
        minimap2_alignment_dot_plot.py
        
Using nucmer to map and return delta file:
    1: sbatch_mummer_THILLE.sh 
        UoQ_thille.hifiv3.asm.p_ctg.fasta.out.delta
        UoQ_thille.hifiv3.asm.p_ctg.fasta.outMM.mgaps
        UoQ_thille.hifiv3.asm.p_ctg.fasta.outMM.ntref
    
Getting Tilting Path:
    2: sbatch_show-tiling_UoQTHILLE.sh
        UoQ_thille.hifiv3.asm.p_ctg.fasta.out.delta.txt
Getting Bed file from Tilting Path:
    3: sbatch_tiling_to_bed.sh
        thille_hass_tilting_output.bed
Getting Structural Vvariants:
    4: sbatch_show-coords_UoQTHILLE.sh
        UoQ_thille.hifiv3.asm.p_ctg.fasta.out.delta.coords
    5: sbatch_show-snps_UoQTHILLE.sh
        UoQ_thille.hifiv3.asm.p_ctg.fasta.out.delta.snps
