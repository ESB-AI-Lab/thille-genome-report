#!/bin/sh
# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start (0-based; BED-like; closed)
# 4	int	Query end (0-based; BED-like; open)
# 5	char	Relative strand: "+" or "-"
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start on original strand (0-based)
# 9	int	Target end on original strand (0-based)
# 10	int	Number of residue matches
# 11	int	Alignment block length
# 12	int	Mapping quality (0-255; 255 for missing)

# want QRY positions (cols 1, 3, 4), relative strand (col 5)
# get REF??
# chrom, start, end, orientation

PAF="thille_to_UoQ.paf"
HEADER="chrom\tstart\tend\tstrand\n"
printf $HEADER
# cat $PAF | awk '{print $1 "\t" $3 "\t" $4 "\t" $5}'
cat $PAF | awk '{print $6 "\t" $8 "\t" $9 "\t" $5}'
# LINE1="CM056814.1\t21131165\t21131684\n"
# LINE2="CM056814.1\t21132796\t21133279\n" 
# LINE3="CM056818.1\t3269725\t3274107\n"    
# LINE4="CM056818.1\t3266660\t3266824"

# for l in $LINE1 $LINE2 $LINE3 $LINE4
# do
#     printf ${l}
# done
