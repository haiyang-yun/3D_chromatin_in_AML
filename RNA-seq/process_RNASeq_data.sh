#!/bin/bash
GSE=$1
GSM=$2
file=$3
genome=$4
genomeDir=$5
annotationFile=$6


mkdir -p $GSE/$GSM
mv $file".r_1.fastq" $GSE/$GSM/$GSM".r_1.fastq"	
mv $file".r_2.fastq" $GSE/$GSM/$GSM".r_2.fastq"
cd $GSE/$GSM
runRNA_STAR_paired.pl  $GSM".r_1.fastq" $GSM".r_2.fastq" $GSM $genome $genomeDir  $annotationFile  "y"
