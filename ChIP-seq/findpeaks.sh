#!/bin/bash
# tag directories are like sorted bam files. It transform the sequence alignment into platform independent data structure. It places all relevant information about the experiment into a "Tag Directory"

makeTagDirectory $1 $1.nodup.bam 

# peak finding using homer
# $1 -> tag directory
# -L -> fold enrichment over local tag count, default: 4.0
# -C ->  fold enrichment limit of expected unique tag positions, default: 2.0
# -size -> Peak size, default: auto
# -mindist ->  (minimum distance between peaks, default: peak size x2)
# -tbp -> Maximum tags per bp to count, 0 = no limit, default: auto
# -o -> output filename

findPeaks $1 -L 0 -C 3 -size 1000 -minDist 1000 -tbp 3 -o $1.peaks
