#!/bin/bash
GSE=$1
GSM=$2
fastq=$3
cd $GSE/$GSM/
file1="$GSM.r_1.fastq"
file2="$GSM.r_2.fastq"
digestFile="/serenity/data/reference-genomes/Digest_mm10_HindIII_None_15-51-59_29-04-2016.txt"
path="/usr/bin/bowtie2"
idx="/serenity/data/reference-genomes/mm10"

hicup --bowtie2 $path  --index $idx --digest $digestFile --format Sanger --longest 800 --shortest 150 --threads 12 $file1 $file2
