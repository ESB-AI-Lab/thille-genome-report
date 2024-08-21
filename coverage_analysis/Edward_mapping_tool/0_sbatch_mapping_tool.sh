#!/bin/bash
#SBATCH -J 0_sbatch_mapping_tool		# jobname
#SBATCH -o step_0.%A.%a	# jobname.o%j for single (non array) jobs jobname.o%A.%a for array jobs
#SBATCH -e step_0.%A.%a	# error file name A is the jobid and a is the arraytaskid
#SBATCH -a 1		# start and stop of the array start-end
#SBATCH -A mcb180013p
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -p RM-shared		# queue (partition) -- compute, shared, large-shared, debug (30mins)
#SBATCH -t 6:00:00		# run time (dd:hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=e1tang@ucsd.edu
###SBATCH --mail-type=begin	# email me when the job starts
###SBATCH --mail-type=end	# email me when the job finishes
#SBATCH --export=ALL
#--------------------------------------------------------------------------------
#INPUTS

# Gwen -> UoQ
: '
title="gwen_reads_onto_UoQ"
reads_path="../symlinks/gwen_pacb_subreads.fastq.gz"
reference_path="../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
minimap2_path="/jet/home/etang1/bin/minimap2"
minimap2_options="map-pb"
samtools_module="samtools/1.13.0"
bedtools_module="bedtools/2.30.0"
window_size="30000"
gene_bed="selected_genes_onto_UoQ.bed"
gene_txt="selected_genes_onto_UoQ.txt"
graph_title="Gwen Raw Pacbio Data Mapped Onto Hass Reference, Median of Window Size 30k"

# Thille -> UoQ
title="thille_reads_onto_UoQ"
reads_path="../symlinks/thille.hifiv3.fastq.gz"
reference_path="../symlinks/GCA_029852735.1_ASM2985273v1_genomic.UoQ.fasta"
minimap2_path="/jet/home/etang1/bin/minimap2"
minimap2_options="map-hifi"
samtools_module="samtools/1.13.0"
bedtools_module="bedtools/2.30.0"
window_size="30000"
gene_bed="selected_genes_onto_UoQ.bed"
gene_txt="selected_genes_onto_UoQ.txt"
graph_title="Thille Raw HiFi Data Mapped Onto Hass Reference, Median of Window Size 30k"

# Gwen -> Thille (coverage.txt symlinked from gavin)
title="gwen_reads_onto_thille"
reads_path=""
reference_path=""
minimap2_path=""
minimap2_options=""
samtools_module="samtools/1.13.0"
bedtools_module="bedtools/2.30.0"
window_size="30000"
gene_bed="selected_genes_onto_thille.bed"
gene_txt="selected_genes_onto_thille.txt"
graph_title="Gwen Raw Pacbio Data Mapped Onto Thille Reference, Median of Window Size 30k"

'
# Thille -> Thille (coverage.txt symlinked from gavin)
title="thille_reads_onto_thille"
reads_path=""
reference_path=""
minimap2_path=""
minimap2_options=""
samtools_module="samtools/1.13.0"
bedtools_module="bedtools/2.30.0"
window_size="30000"
gene_bed="selected_genes_onto_thille.bed"
gene_txt="selected_genes_onto_thille.txt"
graph_title="Thille Raw HiFi Data Mapped Onto Thille Reference, Median of Window Size 30k"

#--------------------------------------------------------------------------------
# Do not comment out functions

function sbatch_check {
    if [[ $? -ne 0 ]]
    then
        echo "failure"
        exit 1
    fi
}

function error_check {
    while true
    do
        if [ -f "success.txt" ]; then
            rm success.txt
            break
        fi
        if [ -f "failure.txt" ]; then
            echo "failure"
            rm failure.txt
            exit 1
        fi
        sleep 5
    done
}

# seq the number of items in the sbatch array
function completion_check {



for i in $(seq 1 30);
do
    while true
    do
        if [ -f "complete_${i}.txt" ]; then
            break
        fi
        sleep 5
    done
done
rm "complete"*".txt"
}
: '
#--------------------------------------------------------------------------------1
if [ -f "jobfile_${title}.txt" ]; then
    echo "jobfile_${title} found"
else
    echo "jobfile_${title} not found"
    exit 1
fi

echo "attempting to create ${title}"
mkdir ${title}
if [[ $? -ne 0 ]]
then
    exit 1
fi
echo "${title} created"
echo

echo "calling 1_sbatch_alignment.sh with parameters:"
echo "reads_path = ${reads_path}"
echo "reference_path = ${reference_path}"
echo "minimap2_path = ${minimap2_path}"
echo "minimap2_options = ${minimap2_options}"
start=$SECONDS
sbatch 1_sbatch_alignment.sh ${title} ${reads_path} ${reference_path} ${minimap2_path} ${minimap2_options}  
sbatch_check
error_check
echo "completed in $(( SECONDS - start )) seconds"
echo

#--------------------------------------------------------------------------------2
echo "calling 2_sbatch_sort.sh with parameters:"
echo "samtools_module = ${samtools_module}"
start=$SECONDS
sbatch 2_sbatch_sort.sh ${title} ${samtools_module}
sbatch_check
error_check
echo "completed in $(( SECONDS - start )) seconds"
echo

#--------------------------------------------------------------------------------3
echo "calling 3_sbatch_coverage.sh with parameters:"
echo "bedtools_module = ${bedtools_module}"
start=$SECONDS 
sbatch 3_sbatch_coverage.sh ${title} ${bedtools_module}
sbatch_check
error_check
echo "completed in $(( SECONDS - start )) seconds"
echo

#--------------------------------------------------------------------------------4
echo "calling 4_sbatch_graph_coverage.sh with parameters:"
echo "window_size = ${window_size}"
echo "gene_txt = ${gene_txt}"
start=$SECONDS
sbatch 4_sbatch_graph_coverage.sh ${title} ${window_size} ${gene_txt}
sbatch_check
completion_check
echo "completed in $(( SECONDS - start )) seconds"
echo
'
#--------------------------------------------------------------------------------5
echo "calling 5_sbatch_stack_graphs.sh with parameters:"
echo "graph_title = ${graph_title}"
start=$SECONDS
sbatch 5_sbatch_stack_graphs.sh ${title} "${graph_title}"
sbatch_check
error_check
echo "completed in $(( SECONDS - start )) seconds"
echo

#--------------------------------------------------------------------------------6
echo "calling 6_sbatch_coverage_analysis.sh with parameters:"
echo "bedtools_module = ${bedtools_module}"
echo "gene_bed = ${gene_bed}"
start=$SECONDS
sbatch 6_sbatch_coverage_analysis.sh ${title} ${bedtools_module} ${gene_bed}
sbatch_check
completion_check
rm -f ${title}/intersection.bed
touch ${title}/intersection.bed
ls ${title}/temp_intersection_*.bed | xargs -I {} cat {} >> ${title}/intersection.bed
echo "completed in $(( SECONDS - start )) seconds"
echo

#--------------------------------------------------------------------------------
