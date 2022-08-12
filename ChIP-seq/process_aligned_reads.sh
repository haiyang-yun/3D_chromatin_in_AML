#!/bin/bash

# Example usage of running this script: process_aligned_reads.sh -g [GENOTYPE] -m [OUTPUT_FOLDER] -x mm10
# -g [GENOTYPE]: top folder named as genotype (e.g. WT, FLT3, NPM1, DM)for input sam files 
# -m [OUTPUT_FOLDER]: sub folder (e.g. H3K4me1, ATAC-seq) where sam files are stored and serve as output folder to store filtered reads 

GSE=$1
GSM=$2
genome=$3
fileName=$GSE/$GSM/$GSM
picard-tools SortSam INPUT=$fileName.sam OUTPUT=$fileName.sorted.bam SO=coordinate
picard-tools MarkDuplicates INPUT=$fileName.sorted.bam OUTPUT=$fileName.nodup.bam REMOVE_DUPLICATES=true METRICS_FILE=$fileName.metricN.log  VALIDATION_STRINGENCY=LENIENT
samtools index $fileName.nodup.bam
rm $fileName.sorted.bam
samtools view -h $fileName.nodup.bam > $fileName.nodup.sam
