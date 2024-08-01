Goal: Map Gwen reads onto the Thille genome. Compute and graph coverage.

Input:
    Thille genome:             ../thille.hifiv3.asm.p_ctg.fasta 
    Gwen reads:                ../gwen_pacb_subreads.fastq.gz
    Thille Lengths:            ../thille.hifiv3.asm.p_ctg.masked.lengths
    Python Graphing:           ../graphCoverage.py
    Python Dips:               ../findCoverageDips.py
Output:
    Coverage by Position:   gwen_thilleCoverage_output.cov.txt
    And Coverage Graphs:    ./graphs
   
Procedure:
    sbatch 1_sbatch_Gwen_ThilleCoverage.sh
        gwen_pacb_subreads.fastq.gz_thille.hifiv3.asm.p_ctg.sam
        gwen_pacb_subreads.fastq.gz_thille.hifiv3.asm.p_ctg.bam
        gwen_pacb_subreads.fastq.gz_thille.hifiv3.asm.p_ctg.sorted.bam
        gwen_pacb_subreads.fastq.gz_thille.hifiv3.asm.p_ctg.sorted.bam.bai
        thille.hifiv3.asm.p_ctg.fasta.fai
    sbatch 2_sbatch_Gwen_ThilleCoverage.sh
        gwen_pacb_subreads.fastq.gz_thille.hifiv3.asm.p_ctg.bed
        gwen_thilleCoverage_output.cov.txt
    sbatch 3_sbatch_graphCoverage.sh
        30k_graph_data
        graph_gwen_onto_thille_30k_median_ex-10
    sbatch 4_sbatch_findCoverageDrops.sh
        coverage_dips_1k
