# ATAC-seq data mapping and peak calling

Chromatin accessibility probed by ATAC-seq was analysed in like manner as H3K27ac ChIP-seq analysis, including the procedures of trimming, mapping, filtering and peak calling. Briefly, trimmed sequences were mapped against mm10 reference genome using Bowtie2 and only uniquely mapped reads were kept. Peaks were called using MACS2 with the setting “-nomodel -nolambda” and only those with p-values less than 1e-20 were considered significant. A list of ATAC-seq consensus peak set was made using DiffBind. 

## reads mapping
get_data.sh
run_getData.sh

## remove duplicates
remdup.sh
run_remdup.sh

## peakcalling
peakCalling.sh
run_peakCalling.sh
run_MACS.NoModel.Nolambda.pl






