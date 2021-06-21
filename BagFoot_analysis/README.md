# Bivariate-Genomic-Footprinting

Code for footprinting and motif-flanking accessibility using BaGFoot (Baek et al)

Analysis was performed on ATAC-seq profiles at all accessible regions represented by ATAC-seq consensus peaks. In brief, known TF motifs were first scanned at ± 100 bp from ATAC-seq peak summit and their occurrence was aggregated. Next, for each motif, ATAC-seq reads were counted at aggregated motif-flanking regions (±200 bp from motif center) and normalised to total sequencing counts in each sample. Differential accessibility was then determined as altered read counts flanking each motif in mutant vs WT. To probe TF footprinting depth, expected cuts were computed by mapping ATAC-seq reads to each motif occurred genome-wide and was set as baseline, and observed cuts only counted the reads within ± 100 bp from ATAC-seq peak summit. Footprinting depth was calculated as cut bias represented by a log ratio difference of observed cuts divided by expected cuts. Together, for each comparison between mutant and WT HSPC, the differences of footprinting depth and flanking accessibility for each TF motif were plotted in a Bagplot, where the inner polygon (“bag”) encompasses at most half of the TFs and the outer polygon “fence” is formed by inflating the bag geometrically by a default factor of 2.5. The outliers with a p-value less than 0.05 were considered significant.


aux.R
bagfoot.R

Baek S, et al (2017) ***Genomic Footprinting Detects Changes in Transcription Factor Activity*** ([Cell Rep 19, 1710-1722](https://pubmed.ncbi.nlm.nih.gov/28538187/)).
