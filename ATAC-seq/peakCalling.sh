#!/bin/bash
REPORT=MACS.$3.log #$3 is pvalue
mkdir -p $1/$2/PeakCalling/MACS.NoModel.NoLambda
cd $1/$2/PeakCalling/MACS.NoModel.NoLambda
run_MACS.NoModel.NoLambda.pl $2.nodup.BED $3  2> $REPORT
