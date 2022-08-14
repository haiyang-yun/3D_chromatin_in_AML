# 3D_chromatin_in_AML
Release of custom code for the manuscript "Mutational synergy during leukemia induction remodels chromatin accessibility, histone modifications and three-dimensional DNA topology to alter gene expressionL" (Yun et al. Nature Genetics 2021)

The scripts are provided for the integrative analysis of chromatin accessibility (ATAC-seq), chromatin states (ChIP-seq), DNA looping (pCHiC) and transcriptome (RNA-se) across four cellular states in a murine allelic series that models gene mutational synergy in the induction of acute myeloid leukemia. 

Analytical procedures using these scripts cover from raw data processing all the way through to integrated data analysis. The raw data processing includes quality control, read mapping, data filtering, normalization, and statistical calling. Subsequently, the processed data are first subjected to integrated analysis on dynamic chromatin states, to reveal differential clusters of cis-regulatory elements (CREs) that demonstrate similar dynamic chromatin modifications. Afterwards, the specific clusters of CREs with characteristic gain or loss of enhancer signatures are annotated to target genes, using either linear or spatial proximity information. Differential mRNA expression is further analysed for these genes along with their associated functional network. The relevant biological information of the data used and their functional interpretation are discussed in great detail in (Yun et al. Nature Genetics 2021).

## contents
11 directories, 60 files
.
├── ATAC-seq
│   ├── README.md
│   ├── get_data.sh
│   ├── peakCalling.sh
│   ├── remdup.sh
│   ├── run_MACS.NoModel.NoLambda.pl
│   ├── run_getData.sh
│   ├── run_peakCalling.sh
│   └── run_remdup.sh
├── BagFoot_analysis
│   ├── README.md
│   ├── aux.R
│   ├── bagfoot.R
│   └── run_bagfoot.sh
├── ChIP-seq
│   ├── README.md
│   ├── enhancer_calling.R
│   ├── filter_overlapped_intervals.py
│   ├── filter_peaks.sh
│   ├── filterpeaks_based_on_index.py
│   ├── findpeaks.sh
│   ├── get_catalog.py
│   ├── get_data.sh
│   ├── peakCalling_H3K27ac.sh
│   ├── process_aligned_reads.sh
│   ├── run_MACS.NoModel.pl
│   ├── run_findpeaks.sh
│   ├── run_getData.sh
│   ├── run_peakCalling_H3K27ac.sh
│   └── run_process_aligned_reads.sh
├── Differential_analysis
│   ├── ATACseq_differential_analysis.R
│   ├── H3K27ac_differential_analysis.R
│   ├── H3K4me1_differential_analysis.R
│   ├── README.md
│   └── RNAseq_differential_analysis.R
├── Homer_motif_anlaysis
│   ├── Homer_motif_analysis.R
│   └── README.md
├── PCA_correlation_analysis
│   ├── H3K4me1_correlation.R
│   ├── PCA_analysis.R
│   └── README.md
├── README.md
├── RNA-seq
│   ├── README.md
│   ├── process_RNASeq_data.sh
│   ├── runRNA_STAR_paired.pl
│   └── run_RNASeq_jobs.sh
├── Seurat_analysis
│   ├── Multiomics_Seurat_analysis.R
│   └── README.md
├── Volcanoplots_MAplots
│   ├── ATAC_MAplots.R
│   ├── Chromatin_modifications_MAplots.R
│   ├── README.md
│   └── RNAseq_volcanoplots.R
├── deeptools_analysis
│   ├── README.md
│   ├── deeptools_ChIPseq_ATACseq.R
│   ├── deeptools_leukemia_specific_changes.R
│   └── deeptools_tsne_clusters.R
└── pCHiC
    ├── README.md
    ├── hicup.sh
    ├── pCHiC_rewire_analysis.R
    ├── plotBaits.R
    ├── runChicago_interactions_calling.sh
    ├── runEdgeR_pCHiC_promoter_reads.R
    ├── runHOMER_compartment.sh
    └── runHiCUP_pCHiCmapping.sh
