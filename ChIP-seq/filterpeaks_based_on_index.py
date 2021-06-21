#!/usr/bin/env python

############################################################################################################
# filter peaks that overlap between replicates from each parent file based on indices from mergePeak 
# filterpeaks_based_on_index.py DM.H3K4me1.R1.indices DM.H3K4me1.R1.peaks  > DM.H3K4me1.R1.peaks.overlap
# Indices are retrieve using following command
# grep -v \# mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks | awk '{print $9}'|sort >  DM.H3K4me1.R1.indices
# DM.H3K4me1.R1.peaks -> from mergePeaks
############################################################################################################

import sys
index_list=[line.strip("\n") for line in open(sys.argv[1])]

with open(sys.argv[2]) as infile:
	for lines in infile:
		if not lines.startswith("#"):
			line=lines.strip("\n").split("\t")
			if line[0] in index_list:
				print lines.strip("\n")
