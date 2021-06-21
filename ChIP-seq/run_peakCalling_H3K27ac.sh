#!/bin/bash

## DM
./peakCalling.sh	DM	DM.H3K27ac.R1	1e-20	DM.input.BED
./peakCalling.sh	DM	DM.H3K27ac.R2	1e-20	DM.input.BED

## FLT3
./peakCalling.sh	FLT3	FLT3.H3K27ac.R1	1e-20	FLT3.input.BED
./peakCalling.sh	FLT3	FLT3.H3K27ac.R2	1e-20	FLT3.input.BED
## WT
./peakCalling.sh	WT	WT.H3K27ac.R1	1e-20	WT.input.BED
./peakCalling.sh	WT	WT.H3K27ac.R2	1e-20	WT.input.BED

## NPM1
./peakCalling.sh	NPM1	NPM1.H3K27ac.R1	1e-20	NPM1.input.BED
./peakCalling.sh	NPM1	NPM1.H3K27ac.R2	1e-20	NPM1.input.BED
