# RNA-seq data mapping anlaysis

RNA-seq reads were quality filtered and mapped using STAR against the mouse genome (mm10). Uniquely mapped reads were quantified with HTSeq vand protein-coding genes with non-zero read count in wildtype or mutant HSPC (n = 16,771) were included for downstream analysis. Reads Per Kilobase Million (RPKM) mapped reads for each protein-coding genes were calculated using Bioconductor package edgeR.

run_RNASeq_jobs.sh
process_RNASeq_data.sh
runRNA_STAR_paired.pl

