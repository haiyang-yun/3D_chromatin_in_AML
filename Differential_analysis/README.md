# Differential_analysis_NGSdata

Differential expression of protein-coding genes was analysed with these counts using Bioconductor package DESeq2. Significantly differential expression was considered by setting adjusted P value (adjP) < 0.05 and fold change (FC) ≥ 1.5 between mutants and wildtype HSPC. 

RNAseq_differential_analysis.R

Differential enrichment was analysed using edgeR, with significant changes being defined by FDR value < 0.05 and FC ≥ 2 (gain or loss) in the presence of any mutations. 

ATACseq_differential_analysis.R


Differential enrichment of chromatin marks (H3K4me1 and H3K27ac) was analysed using edgeR, with significant changes being defined by FDR value < 0.05 and FC ≥ 1.5 (gain or loss) in the presence of any mutations. 

H3K4me1_differential_analysis.R
H3K27ac_differential_analysis.R

