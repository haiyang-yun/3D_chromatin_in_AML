#!/bin/bash

# convert bam to chicago chinput files
bam2chicago.sh WT.CHiC.R1.bam CHiC.mm10.baitmap Digest.mm10.rmap WT.CHiC.R1 nodelete

bam2chicago.sh WT.CHiC.R2.bam CHiC.mm10.baitmap Digest.mm10.rmap WT.CHiC.R2 nodelete

bam2chicago.sh DM.CHiC.R1.bam CHiC.mm10.baitmap Digest.mm10.rmap DM.CHiC.R1 nodelete

bam2chicago.sh DM.CHiC.R2.bam CHiC.mm10.baitmap Digest.mm10.rmap DM.CHiC.R2 nodelete

bam2chicago.sh FLT3.CHiC.R1.bam CHiC.mm10.baitmap Digest.mm10.rmap FLT3.CHiC.R1 nodelete

bam2chicago.sh FLT3.CHiC.R2.bam CHiC.mm10.baitmap Digest.mm10.rmap FLT3.CHiC.R2 nodelete

bam2chicago.sh NPM1.CHiC.R1.bam CHiC.mm10.baitmap Digest.mm10.rmap NPM1.CHiC.R1 nodelete

bam2chicago.sh NPM1.CHiC.R2.bam CHiC.mm10.baitmap Digest.mm10.rmap NPM1.CHiC.R2 nodelete

# Call interactions on chinput files
Rscript runChicago.R --design-dir .-design-dir WT.CHiC.R1.chinput,WT.CHiC.R2.chinput WT.CHiC.R1-2

Rscript runChicago.R --design-dir .-design-dir DM.CHiC.R1.chinput,DM.CHiC.R2.chinput DM.CHiC.R1-2

Rscript runChicago.R --design-dir .-design-dir NPM1.CHiC.R1.chinput,NPM1.CHiC.R2.chinput NPM1.CHiC.R1-2

Rscript runChicago.R --design-dir .-design-dir FLT3.CHiC.R1.chinput,FLT3.CHiC.R2.chinput FLT3.CHiC.R1-2

