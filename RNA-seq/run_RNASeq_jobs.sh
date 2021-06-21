#!/bin/bash
# Wildtype
process_RNASeq_data.sh WT.EXPR WT.EXPR.R1 WT.EXPR.R1 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
process_RNASeq_data.sh WT.EXPR WT.EXPR.R2 WT.EXPR.R2 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
# Double Mutant
process_RNASeq_data.sh DM.EXPR DM.EXPR.R1 DM.EXPR.R1 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
process_RNASeq_data.sh DM.EXPR DM.EXPR.R2 DM.EXPR.R2 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
# NPM1 mutation
process_RNASeq_data.sh NPM1.EXPR NPM1.EXPR.R1 NPM1.EXPR.R1 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
process_RNASeq_data.sh NPM1.EXPR NPM1.EXPR.R2 NPM1.EXPR.R2 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
# FLT3 mutation
process_RNASeq_data.sh FLT3.EXPR FLT3.EXPR.R1 FLT3.EXPR.R1 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10
process_RNASeq_data.sh FLT3.EXPR FLT3.EXPR.R2 FLT3.EXPR.R2 STAR-GENOMES-mm10.gencode.vM7.comprehensive gencode.vM7.comprehensive.annotation.gtf  mm10


