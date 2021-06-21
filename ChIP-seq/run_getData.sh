#!/bin/bash

## Wildtype
get_data.sh  -g WT -m WT.H3K27ac.R1  -i WT.H3K27ac.R1.fastq  -x mm10
get_data.sh  -g WT -m WT.H3K27ac.R2  -i WT.H3K27ac.R2.fastq  -x mm10
get_data.sh  -g WT -m WT.H3K4me1.R1  -i WT.H3K4me1.R1.fastq -x mm10
get_data.sh  -g WT -m WT.H3K4me1.R2  -i WT.H3K4me1.R2.fastq -x mm10
get_data.sh  -g WT -m WT.H3K4me3.R1  -i WT.H3K4me3.R1.fastq -x mm10
get_data.sh  -g WT -m WT.H3K4me3.R2  -i WT.H3K4me3.R2.fastq -x mm10
get_data.sh  -g WT -m WT.input -i WT.input.fastq -x mm10

## Neutrophil
get_data.sh  -g NE -m NE.H3K27ac.R1  -i NE.H3K27ac.R1.fastq  -x mm10
get_data.sh  -g NE -m NE.H3K27ac.R2  -i NE.H3K27ac.R2.fastq  -x mm10
get_data.sh  -g NE -m NE.H3K4me1.R1  -i NE.H3K4me1.R1.fastq -x mm10
get_data.sh  -g NE -m NE.H3K4me1.R2  -i NE.H3K4me1.R2.fastq -x mm10
get_data.sh  -g NE -m NE.H3K4me3.R1  -i NE.H3K4me3.R1.fastq -x mm10
get_data.sh  -g NE -m NE.H3K4me3.R2  -i NE.H3K4me3.R2.fastq -x mm10
get_data.sh  -g NE -m NE.input -i NE.input.fastq -x mm10

## Double mutant
get_data.sh  -g DM -m DM.H3K27ac.R1 -i DM.H3K27ac.R1.fastq  -x mm10
get_data.sh  -g DM -m DM.H3K27ac.R2 -i DM.H3K27ac.R2.fastq  -x mm10
get_data.sh  -g DM -m DM.H3K4me1.R1 -i DM.H3K4me1.R1.fastq -x mm10
get_data.sh  -g DM -m DM.H3K4me1.R2 -i DM.H3K4me1.R2.fastq -x mm10
get_data.sh  -g DM -m DM.H3K4me3.R1 -i DM.H3K4me3.R1.fastq -x mm10
get_data.sh  -g DM -m DM.H3K4me3.R2 -i DM.H3K4me3.R2.fastq -x mm10
get_data.sh  -g DM -m DM.input -i DM.input.fastq -x mm10

## FLT3
get_data.sh  -g FLT3 -m FLT3.H3K27ac.R1 -i FLT3.H3K27ac.R1.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.H3K27ac.R2 -i FLT3.H3K27ac.R2.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.H3K4me1.R1 -i FLT3.H3K4me1.R1.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.H3K4me1.R2 -i FLT3.H3K4me1.R2.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.H3K4me3.R1 -i FLT3.H3K4me3.R1.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.H3K4me3.R2 -i FLT3.H3K4me3.R2.fastq  -x mm10
get_data.sh  -g FLT3 -m FLT3.input -i FLT3.input.fastq  -x mm10

## NPM1
get_data.sh  -g NPM1 -m NPM1.H3K27ac.R1 -i NPM1.H3K27ac.R1.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.H3K27ac.R2 -i NPM1.H3K27ac.R2.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.H3K4me1.R1 -i NPM1.H3K4me1.R1.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.H3K4me1.R2 -i NPM1.H3K4me1.R2.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.H3K4me3.R1 -i NPM1.H3K4me3.R1.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.H3K4me3.R2 -i NPM1.H3K4me3.R2.fastq  -x mm10
get_data.sh  -g NPM1 -m NPM1.input -i NPM1.input.fastq  -x mm10

