#!/bin/bash

## Lineage Negative     
./process_aligned_reads.sh  -g WT -m WT.H3K27ac.R1 -x mm10
./process_aligned_reads.sh  -g WT -m WT.H3K27ac.R2 -x mm10
./process_aligned_reads.sh  -g WT -m WT.H3K4me1.R1 -x mm10
./process_aligned_reads.sh  -g WT -m WT.H3K4me1.R2 -x mm10
./process_aligned_reads.sh  -g WT -m WT.H3K4me3.R1 -x mm10
./process_aligned_reads.sh  -g WT -m WT.H3K4me3.R2 -x mm10
./process_aligned_reads.sh  -g WT -m WT.input -x mm10

## Neutrophills
./process_aligned_reads.sh  -g NE -m NE.H3K27ac.R1 -x mm10
./process_aligned_reads.sh  -g NE -m NE.H3K27ac.R2 -x mm10
/process_aligned_reads.sh  -g NE -m NE.H3K4me1.R1 -x mm10
./process_aligned_reads.sh  -g NE -m NE.H3K4me1.R2 -x mm10
./process_aligned_reads.sh  -g NE -m NE.H3K4me3.R1 -x mm10
./process_aligned_reads.sh  -g NE -m NE.H3K4me3.R2 -x mm10
./process_aligned_reads.sh  -g NE -m NE.input -x mm10

## Double mutant     
./process_aligned_reads.sh  -g DM -m DM.H3K27ac.R1 -x mm10
./process_aligned_reads.sh  -g DM -m DM.H3K27ac.R2 -x mm10
./process_aligned_reads.sh  -g DM -m DM.H3K4me1.R1 -x mm10
./process_aligned_reads.sh  -g DM -m DM.H3K4me1.R2 -x mm10
./process_aligned_reads.sh  -g DM -m DM.H3K4me3.R1 -x mm10
./process_aligned_reads.sh  -g DM -m DM.H3K4me3.R2 -x mm10
./process_aligned_reads.sh  -g DM -m DM.input -x mm10

## FLT3      
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K27ac.R1 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K27ac.R2 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K4me1.R1 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K4me1.R2 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K4me3.R1 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.H3K4me3.R2 -x mm10
./process_aligned_reads.sh  -g FLT3 -m FLT3.input -x mm10

## NPM1      
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K27ac.R1 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K27ac.R2 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K4me1.R1 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K4me1.R2 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K4me3.R1 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.H3K4me3.R2 -x mm10
./process_aligned_reads.sh  -g NPM1 -m NPM1.input -x mm10
