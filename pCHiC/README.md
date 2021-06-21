# Chromatin interaction anlaysis on pCHiC data

## Promoter-anchored interaction analysis
Paired-end sequences of pCHiC were processed using the HiCUP pipeline with default parameters for the following steps: quality control, identification of reads containing HiC junctions, mapping to reference genome mm10, and filtering duplicated HiC di-tags. Output bam files containing valid HiC di-tags were processed by Bioconductor package CHiCAGO, to call significant promoter-based interactions. CHiCAGO considers distance effect on interaction frequencies by virtue of a convolution background model and a distance-weighted p-value. CHiCAGO scores represent -log weighted p-values, the higher the more likelihood of interaction formed. CHiCAGO scores were calculated for each pCHiC sample, and significant interactions were called when CHiCAGO scores ≥ 5. Significant interactions were computed per HSPC by merging their replicates, and were combined to form a matrix of total unique interactions. While, interactions present in both replicates per HSPC were considered with high confidence (HC). Total pCHiC reads at individual promoters were summed to perform differential analysis in mutant vs WT HSPC using edgeR. Differential total interaction reads were defined by adjP < 0.05 and absolute FC > 1. Rewired interactions were identified by comparing and ranking their CHiCAGO scores in WT and mutant cells. Those HC interactions, absent in WT (score < 5) but present in mutant (score ≥ 5), with scores ranked in the bottom quartile in WT but in the top quartile in mutant, were considered mutation-associated gained interactions. And vice versa, HC interactions that were present in WT but absent in mutant, ranked in the top quartile in WT but in the bottom quartile in mutant, were considered as lost interactions by mutations. 

## Chromatin compartment analysis
Sub-nuclear compartmentation represented by self-associating chromatin domains were analysed by means of principal component analysis (PCA) on capture HiC data in a similar way as described on HiC data, mainly using HOMER software. To do so, uniquely mapped reads from HiCUP analysis were used to create "Tag Directory" using makeTagDirectory from HOMER. Principal component (PC1) values were calculated by running runHiCpca.pl with default setting (resolution at 50 kb). This led to the separation of chromatin into two compartments, with positive PC1 regions reflecting "active" chromatin and negative PC1 regions indicative of "inactive" chromatin. Regions of continuous positive or negative PC1 values were stitched to be identified as A or B compartments, respectively. Genome-wide correlation of compartment PC1 values between mutant and WT cells was performed by running getHiCcorrDiff.pl from HOMER. Flipped compartments were identified using HOMER findHiCCompartments.pl. 

## pCHiC data mapping:
hicup.sh
runHiCU_pCHiCmapping.sh

## Chromatin compartment analysis:
runHOMER_compartment.sh

## Interaction calling:
runChicago_interactions_calling.sh

## Interactions viewed as 4C data
plotBaits.R

## Differential anlaysis at individual promoters:
runEdgeR_pCHiC_promoter_reads.R

## Interactions rewiring analysis:
pCHiC_rewire_analysis.R


