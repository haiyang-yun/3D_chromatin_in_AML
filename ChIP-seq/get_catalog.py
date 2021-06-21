#!/usr/bin/env python

##############################################################################
# Compare the consecutive genomic intervals and keep the one with higher score
# if the distance between the intervals is less than 500. If the
# It takes the sorted file after concatenation of peaks from all 
# the samples
##############################################################################
import sys
with open(sys.argv[1]) as infile:
	line=next(infile)
	line1=line.strip("\n").split("\t")
	for lines in infile:
		line2=lines.strip("\n").split("\t")
		#print line1
		#print line2
		if line1[0]==line2[0] and int(line2[1])<= int(line1[2]):
			if int(line2[1])-int(line1[2]) < 500:
				if float(line1[5]) >= float(line2[5]):
					pass
				else:
					line1=line2
			else:   
				print "\t".join(line1)
				line1=line2	
		else:
			print "\t".join(line1)
			line1=line2
		#raw_input("Press Enter to continue...")
