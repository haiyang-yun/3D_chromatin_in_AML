#!/bin/bash

GSE=$1
GSM=$2
fileName=$GSE/$GSM/$GSM
picard-tools SortSam INPUT=$fileName.sam OUTPUT=$fileName.sorted.bam SO=coordinate
picard-tools MarkDuplicates INPUT=$fileName.sorted.bam OUTPUT=$fileName.nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=$fileName.metricN.log  VALIDATION_STRINGENCY=LENIENT
samtools index $fileName.nodup.bam
samtools view -bhS $fileName.sam > $fileName.bam
rm $fileName.sam $fileName.sorted.bam
