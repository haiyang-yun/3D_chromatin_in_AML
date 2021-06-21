#!/bin/bash
### merge replicates
header=WT.CHiC.R1.r_1_2.hicup.sam
files="WT.CHiC.R1.r_1_2.hicup.sam WT.CHiC.R2.r_1_2.hicup.sam"
output=WT.CHiC.R1-2.sam
(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $output

header=NPM1.CHiC.R1.r_1_2.hicup.sam
files="NPM1.CHiC.R1.r_1_2.hicup.sam NPM1.CHiC.R2.r_1_2.hicup.sam"
output=NPM1.CHiC.R1-2.sam
(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $output

header=FLT3.CHiC.R1.r_1_2.hicup.sam
files="FLT3.CHiC.R1.r_1_2.hicup.sam FLT3.CHiC.R2.r_1_2.hicup.sam"
output=FLT3.CHiC.R1-2.sam
(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $output

header=DM.CHiC.R1.r_1_2.hicup.sam
files="DM.CHiC.R1.r_1_2.hicup.sam DM.CHiC.R2.r_1_2.hicup.sam"
output=DM.CHiC.R1-2.sam
(grep ^@ $header; for f in $files; do grep -v ^@ $f; done) > $output


### HOMER HiC PCA calculation
#WT.CHiC.R1-2
hicup2homer WT.CHiC.R1-2.sam
makeTagDirectory WT.CHiC.R1-2 -format HiCsummary WT.CHiC.R1-2.sam.homer
runHiCpca.pl WT.CHiC.R1-2 WT.CHiC.R1-2/ -res 50000 -cpu 12 -genome mm10
bedGraphToBigWig WT.CHiC.R1-2.PC1.bedGraph ../../mm10.chrom.sizes.v2 WT.CHiC.R1-2.PC1.bw
mkdir -p WT.CHiC/WT.CHiC.R1-2
mv WT.CHiC.R1-2* WT.CHiC/WT.CHiC.R1-2

# NPM1.CHiC.R1-2
hicup2homer NPM1.CHiC.R1-2.sam
makeTagDirectory NPM1.CHiC.R1-2 -format HiCsummary NPM1.CHiC.R1-2.sam.homer
runHiCpca.pl NPM1.CHiC.R1-2 NPM1.CHiC.R1-2/ -res 50000 -cpu 12 -genome mm10
bedGraphToBigWig NPM1.CHiC.R1-2.PC1.bedGraph ../../mm10.chrom.sizes.v2 NPM1.CHiC.R1-2.PC1.bw
mkdir -p NPM1.CHiC/NPM1.CHiC.R1-2
mv NPM1.CHiC.R1-2* NPM1.CHiC/NPM1.CHiC.R1-2

# FLT3.CHiC.R1-2
hicup2homer FLT3.CHiC.R1-2.sam
makeTagDirectory FLT3.CHiC.R1-2 -format HiCsummary FLT3.CHiC.R1-2.sam.homer
runHiCpca.pl FLT3.CHiC.R1-2 FLT3.CHiC.R1-2/ -res 50000 -cpu 12 -genome mm10
bedGraphToBigWig FLT3.CHiC.R1-2.PC1.bedGraph ../../mm10.chrom.sizes.v2 FLT3.CHiC.R1-2.PC1.bw
mkdir -p FLT3.CHiC/FLT3.CHiC.R1-2
mv FLT3.CHiC.R1-2* FLT3.CHiC/FLT3.CHiC.R1-2

# DM.CHiC.R1-2
hicup2homer DM.CHiC.R1-2.sam
makeTagDirectory DM.CHiC.R1-2 -format HiCsummary DM.CHiC.R1-2.sam.homer
runHiCpca.pl DM.CHiC.R1-2 DM.CHiC.R1-2/ -res 50000 -cpu 12 -genome mm10
bedGraphToBigWig DM.CHiC.R1-2.PC1.bedGraph ../../mm10.chrom.sizes.v2 DM.CHiC.R1-2.PC1.bw
mkdir -p DM.CHiC/DM.CHiC.R1-2
mv DM.CHiC.R1-2* DM.CHiC/DM.CHiC.R1-2


### Commpartment Differential analysis
annotatePeaks.pl WT.CHiC.R1.PC1.txt mm10 -noblanks -bedGraph WT.CHiC.R1.PC1.bedGraph WT.CHiC.R2.PC1.bedGraph NPM1.CHiC.R1.PC1.bedGraph NPM1.CHiC.R2.PC1.bedGraph > NPM1vsWT_output.txt

annotatePeaks.pl WT.CHiC.R1.PC1.txt mm10 -noblanks -bedGraph WT.CHiC.R1.PC1.bedGraph WT.CHiC.R2.PC1.bedGraph FLT3.CHiC.R1.PC1.bedGraph FLT3.CHiC.R2.PC1.bedGraph > FLT3vsWT_output.txt

annotatePeaks.pl WT.CHiC.R1.PC1.txt mm10 -noblanks -bedGraph WT.CHiC.R1.PC1.bedGraph WT.CHiC.R2.PC1.bedGraph DM.CHiC.R1.PC1.bedGraph DM.CHiC.R2.PC1.bedGraph > DMvsWT_output.txt

getDiffExpression.pl NPM1vsWT_output.txt WT WT NPM1 NPM1 -pc1 -export regions > NPM1vsWT_diffOutput.txt

getDiffExpression.pl FLT3vsWT_output.txt WT WT FLT3 FLT3  -pc1 -export regions > FLT3vsWT_diffOutput.txt

getDiffExpression.pl DMvsWT_output.txt WT WT DM DM -pc1 -export regions > DMvsWT_diffOutput.txt


### Compartment correlation
getHiCcorrDiff.pl Compartment.corr.DMvsWT WT.CHiC.R1-2 DM.CHiC.R1-2 -cpu 8 -res 50000
getHiCcorrDiff.pl Compartment.corr.FLT3vsWT WT.CHiC.R1-2 FLT3.CHiC.R1-2 -cpu 8 -res 50000
getHiCcorrDiff.pl Compartment.corr.NPM1vsWT WT.CHiC.R1-2 NPM1.CHiC.R1-2 -cpu 8 -res 50000


