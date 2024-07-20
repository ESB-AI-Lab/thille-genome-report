#!/usr/bin/env/python3

import sys

from typing import List, Tuple

# return list of tuples
# coverage categories:
# 0 = below 0 of median
# 1 = below 0.2 of median
# 2 = below 0.5 of median
# 3 = above 0.5 of median
# (chromosome, start position, region length, coverage category, average coverage, median value used for comparison)

fields = Tuple[str, int, int, int, float, int]

def Category(coverage: int, median: int) -> int:
    if coverage == 0:
        return 0
    elif coverage < median * 0.2:
        return 1
    elif coverage < median * 0.5:
        return 2
    else:
        return 3  # for coverage equal to or above 0.5 of the median

def CoverageAnalysis(coverages: List[str], chromosome: str, median: int) -> List[fields]:
    regions = []
    curr = 0.0
    length = 1
    total = 0.0
    prevCategory = 0
    startPosition = 0

    for i, coverage_line in enumerate(coverages):
        coverage_values = coverage_line.split('\t')
        if len(coverage_values) >= 2:
            curr = float(coverage_values[1])
            category = Category(int(curr), median)

            if category == prevCategory:
                # extend the region
                total += curr
                length += 1
            else:
                # write the region and start a new one
                average = round(total / length, 2)
                regions.append((chromosome, startPosition, length, prevCategory, average, median))
                total = curr
                length = 1
                prevCategory = category
                startPosition = i

    # write the last region
    average = round(total / length, 2)
    regions.append((chromosome, startPosition, length, prevCategory, average, median))

    return regions

# arguments for array sbatching
chromosome = str(sys.argv[1])
median = int(sys.argv[2])

# read coverage file into python list
coverageFile = '30k_graph_data/' + chromosome + '_graph_data' + '.txt'
coverage = []
with open(coverageFile, 'r') as f:
    coverages = f.read().splitlines()

# kowalski, analysis
regions = CoverageAnalysis(coverages, chromosome, median)

# write to file
with open('coverage_analysis/coverage_analysis_thille_reads_onto_thille_' + chromosome + '.bed', 'w') as f:
    f.write('# coverage categories:\n')
    f.write('# 0 = below 0 of median\n')
    f.write('# 1 = below 0.2 of median\n')
    f.write('# 2 = below 0.5 of median\n')
    f.write('# 3 = above 0.5 of median\n')
    f.write('# (chromosome, start position, region length, coverage category, average coverage, median value used for comparison)\n')
    for r in regions:
        # bed format, 6 fields
        f.write(str(r[0]) + '\t' + str(r[1]) + '\t' + str(r[2]) + '\t' + str(r[3]) + '\t' + str(r[4]) + '\t' + str(r[5]) + '\n')
