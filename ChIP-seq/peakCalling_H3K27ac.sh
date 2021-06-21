#!/bin/bash
GSE=$1
GSM=$2
pval=$3
control=$4
REPORT=$GSM.$pval.log #$4=pvalue
mkdir -p $GSE/$GSM/PeakCalling/MACS.NoModel
cd $GSE/$GSM/PeakCalling/MACS.NoModel
./run_MACS.NoModel.pl $GSM.nodup.BED $pval $control  2> $REPORT
