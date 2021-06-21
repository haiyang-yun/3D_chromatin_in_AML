# 3D_chromatin_in_AML
Release of custom code for paper "3D chromatin in AML" Yun et al. 2021

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
