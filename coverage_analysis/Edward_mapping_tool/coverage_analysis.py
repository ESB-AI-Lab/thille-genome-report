#!/usr/bin/env/python3
import sys
from typing import List, Tuple

# arguments for array sbatching
title = sys.argv[1]
chromosome = str(sys.argv[2])
median = int(float(sys.argv[3]))

# return list of tuples 
# coverage categories: 
# 0 = 0 coverage
# 1 = below 0.2 of median
# 2 = between 0.2 and 0.5 of median
# 3 = above 0.5 of median
# (chromosome, start position, stop position, coverage category, average coverage, median value used for comparison)
fields = Tuple[str, int, int, int, float, int]
lengthRequirement = 1000

def Category(coverage: int, median: int) -> int:
    categoryLow = median * 0.2
    categoryMed = median * 0.5
    if coverage == 0:
        return 0
    elif coverage <= categoryLow:
        return 1
    elif coverage <= categoryMed:
        return 2
    else:
        return 3

def CoverageAnalysis(coverages: List[str], chromosome: str, median: int) -> List[fields]:
    regions = []
    curr = int(coverages[0])
    length = 1
    total = curr
    prevCategory = Category(curr, median)
    lenCoverages = len(coverages)
    for i in range(1, lenCoverages):
        curr = int(coverages[i])
        category = Category(curr, median)
        if category == prevCategory:
            # extend the region
            total += curr
            length += 1
        else:
            # write the region and start a new one
            startPosition = i - length
            average = round(total / length, 2)
            if prevCategory != 3 and length >= lengthRequirement:
                regions.append((chromosome, startPosition, startPosition + length, prevCategory, average, median))
            total = curr
            length = 1
            prevCategory = category
    # write the last region
    i = lenCoverages
    startPosition = i - length
    average = round(total / length, 2)
    if prevCategory != 3 and length >= lengthRequirement:
        regions.append((chromosome, startPosition, startPosition + length, prevCategory, average, median))
    return regions

# read coverage file into python list
coverageFile = title + "/temp_coverage_" + chromosome + '.txt'
coverage = []
with open(coverageFile, 'r') as f:
    coverages = f.read().splitlines()

# kowalski, analysis
regions = CoverageAnalysis(coverages, chromosome, median)

# write to file
with open(title + '/' + chromosome + '_coverage_analysis.bed', 'w') as f:
    for r in regions:
        # bed format, 6 fields
        f.write(str(r[0]) + '\t' + str(r[1]) + '\t' + str(r[2]) + '\t' + str(r[3]) + '\t' + str(r[4]) + '\t' + str(r[5]) + '\n')
