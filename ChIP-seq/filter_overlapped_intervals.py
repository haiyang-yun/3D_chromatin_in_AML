#!/usr/bin/env python

##############################################################################
# Compare the distance between the centers of consecutive genomic intervals
# It takes the sorted file after concatenation of peaks from two replicates
# and find the distance between the two centers. If the distance between 
# centers is less than 500,it keep the one with higher signal 
# else keep both the intervals
##############################################################################
import sys

with open(sys.argv[1]) as infile:
	line=next(infile)
	line1=line.strip("\n").split("\t")
	for lines in infile:
		line2=lines.strip("\n").split("\t")
		#print line1
		#print line2
		if line1[0]==line2[0] and (int(line2[1])+500)-(int(line1[1])+500) <= 500:
			if float(line1[4]) >= float(line2[4]):
				pass
			else:   
				line1=line2
		else:
			print "\t".join(line1)
			line1=line2
	print "\t".join(line1)
"""
with open(sys.argv[1]) as infile:
	for lines in infile:
		line2=next(infile)
		line1=lines.strip("\n").split("\t")
		line2=line2.strip("\n").split("\t")
		if float(line1[4]) >= float(line2[4]):
			print "\t".join(line1)
		else:   
			print "\t".join(line2)
	
		#raw_input("Press Enter to continue...")
"""
