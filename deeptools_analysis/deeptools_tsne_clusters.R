### deeptools for tsne clusters
# H3K4me1
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC_consensus_summit2kb_tsne.cluster-1.bed ATAC_consensus_summit2kb_tsne.cluster-2.bed ATAC_consensus_summit2kb_tsne.cluster-3.bed ATAC_consensus_summit2kb_tsne.cluster-4.bed ATAC_consensus_summit2kb_tsne.cluster-5.bed ATAC_consensus_summit2kb_tsne.cluster-6.bed ATAC_consensus_summit2kb_tsne.cluster-7.bed ATAC_consensus_summit2kb_tsne.cluster-8.bed ATAC_consensus_summit2kb_tsne.cluster-9.bed ATAC_consensus_summit2kb_tsne.cluster-10.bed -S LN.H3K4me1.R1-2.cpm_peaks.bw NPM1.H3K4me1.R1-2.cpm_peaks.bw FLT3.H3K4me1.R1-2.cpm_peaks.bw DM.H3K4me1.R1-2.cpm_peaks.bw -o TSNE_10clusters_H3K4me1.gz --missingDataAsZero -p 18

plotProfile -m TSNE_10clusters_H3K4me1.gz -out TSNE_10clusters_H3K4me1_profile.pdf --perGroup --numPlotsPerRow 10 --plotWidth 3.59 --plotHeight 4.72 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' -y 'H3K4me1' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --regionsLabel '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'

# H3K4me3
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC_consensus_summit2kb_tsne.cluster-1.bed ATAC_consensus_summit2kb_tsne.cluster-2.bed ATAC_consensus_summit2kb_tsne.cluster-3.bed ATAC_consensus_summit2kb_tsne.cluster-4.bed ATAC_consensus_summit2kb_tsne.cluster-5.bed ATAC_consensus_summit2kb_tsne.cluster-6.bed ATAC_consensus_summit2kb_tsne.cluster-7.bed ATAC_consensus_summit2kb_tsne.cluster-8.bed ATAC_consensus_summit2kb_tsne.cluster-9.bed ATAC_consensus_summit2kb_tsne.cluster-10.bed -S LN.H3K4me3.R1-2.cpm_peaks.bw NPM1.H3K4me3.R1-2.cpm_peaks.bw FLT3.H3K4me3.R1-2.cpm_peaks.bw DM.H3K4me3.R1-2.cpm_peaks.bw -o TSNE_10clusters_H3K4me3.gz --missingDataAsZero -p 18

plotProfile -m TSNE_10clusters_H3K4me3.gz -out TSNE_10clusters_H3K4me3_profile.pdf --perGroup --numPlotsPerRow 10 --plotWidth 4 --plotHeight 5 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' -y 'H3K4me3' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --regionsLabel '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'

# H3K27ac
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC_consensus_summit2kb_tsne.cluster-1.bed ATAC_consensus_summit2kb_tsne.cluster-2.bed ATAC_consensus_summit2kb_tsne.cluster-3.bed ATAC_consensus_summit2kb_tsne.cluster-4.bed ATAC_consensus_summit2kb_tsne.cluster-5.bed ATAC_consensus_summit2kb_tsne.cluster-6.bed ATAC_consensus_summit2kb_tsne.cluster-7.bed ATAC_consensus_summit2kb_tsne.cluster-8.bed ATAC_consensus_summit2kb_tsne.cluster-9.bed ATAC_consensus_summit2kb_tsne.cluster-10.bed -S LN.H3K27ac.R1-2.cpm_peaks.bw NPM1.H3K27ac.R1-2.cpm_peaks.bw FLT3.H3K27ac.R1-2.cpm_peaks.bw DM.H3K27ac.R1-2.cpm_peaks.bw -o TSNE_10clusters_H3K27ac.gz --missingDataAsZero -p 18

plotProfile -m TSNE_10clusters_H3K27ac.gz -out TSNE_10clusters_H3K27ac_profile.pdf --perGroup --numPlotsPerRow 10 --plotWidth 4 --plotHeight 5 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' -y 'H3K27ac' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --regionsLabel '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'

# ATAC-seq
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC_consensus_summit2kb_tsne.cluster-1.bed ATAC_consensus_summit2kb_tsne.cluster-2.bed ATAC_consensus_summit2kb_tsne.cluster-3.bed ATAC_consensus_summit2kb_tsne.cluster-4.bed ATAC_consensus_summit2kb_tsne.cluster-5.bed ATAC_consensus_summit2kb_tsne.cluster-6.bed ATAC_consensus_summit2kb_tsne.cluster-7.bed ATAC_consensus_summit2kb_tsne.cluster-8.bed ATAC_consensus_summit2kb_tsne.cluster-9.bed ATAC_consensus_summit2kb_tsne.cluster-10.bed -S LN.ATAC.R1-2.cpm_peaks.bw NPM1.ATAC.R1-2.cpm_peaks.bw FLT3.ATAC.R1-2.cpm_peaks.bw DM.ATAC.R1-2.cpm_peaks.bw  -o TSNE_10clusters_ATAC-seq.gz --missingDataAsZero -p 18

plotProfile -m TSNE_10clusters_ATAC-seq.gz -out TSNE_10clusters_ATAC-seq_profile.pdf --perGroup --numPlotsPerRow 10 --plotWidth 4 --plotHeight 5 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' -y 'ATAC-seq' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --regionsLabel '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'






