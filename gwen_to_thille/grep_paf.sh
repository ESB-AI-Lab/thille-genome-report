# run: bash grep_paf.sh > grep_genes_paf.txt
PAF=gwen_to_thille.paf
CHROMS="CM056814.1 CM056818.1"
for n in ${CHROMS}
do
    echo $n
    grep $n $PAF
done
