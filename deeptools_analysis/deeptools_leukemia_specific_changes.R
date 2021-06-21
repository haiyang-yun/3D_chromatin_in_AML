### analyse leukemia specific changes with deeptools
## up
## creat 2-clusters with NE profile only
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R atac_co20.summit200bp.ov.WD.enh.H3K4me1.up.bed -S NE.H3K4me1.cpm.bw -o atac_summit.ov.WD.enh.up_Ne.H3K4me1.gz --missingDataAsZero -p 8

plotHeatmap -m atac_summit.ov.WD.enh.up_Ne.H3K4me1.gz -out atac_summit.ov.WD.enh.up_Ne.H3K4me1.2cluster_heatmap.pdf --colorMap coolwarm --samplesLabel 'NE' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel '1' '2' --refPointLabel 'center' --legendLocation 'none' --outFileSortedRegions atac_summit.ov.Ne.enh.up_H3K4me1.2cluster.bed --kmeans 2 --whatToShow 'heatmap and colorbar'

sort -k1,1 -k2,2n atac_summit.ov.Ne.enh.up_H3K4me1.2cluster.bed | awk '{if ($13 == "1") print $1"\t"$2"\t"$3}' > Ne.up_WDNG_k4me1.cluster1.bed
sort -k1,1 -k2,2n atac_summit.ov.Ne.enh.up_H3K4me1.2cluster.bed | awk '{if ($13 == "2") print $1"\t"$2"\t"$3}' > Ne.up_WDNG_k4me1.cluster2.bed

### down
## creat 2-clusters with NE profile only
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R atac_co20.summit200bp.ov.WD.enh.H3K4me1.down.bed -S NE.H3K4me1.cpm.bw -o atac_summit.ov.WD.enh.down_Ne.H3K4me1.gz --missingDataAsZero -p 18

plotHeatmap -m atac_summit.ov.WD.enh.down_Ne.H3K4me1.gz -out atac_summit.ov.WD.enh.down_Ne.H3K4me1.2cluster_heatmap.pdf --colorMap coolwarm --samplesLabel 'NE' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel '3' '4' --refPointLabel 'center' --legendLocation 'none' --outFileSortedRegions atac_summit.ov.Ne.enh.down_H3K4me1.2cluster.bed --kmeans 2 --whatToShow 'heatmap and colorbar'

sort -k1,1 -k2,2n atac_summit.ov.Ne.enh.down_H3K4me1.2cluster.bed | awk '{if ($13 == "3") print $1"\t"$2"\t"$3}' > Ne.down_WDNG_k4me1.cluster3.bed
sort -k1,1 -k2,2n atac_summit.ov.Ne.enh.down_H3K4me1.2cluster.bed | awk '{if ($13 == "4") print $1"\t"$2"\t"$3}' > Ne.down_WDNG_k4me1.cluster4.bed

### combined plots
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R Ne.up_WDNG_k4me1.cluster1.bed Ne.up_WDNG_k4me1.cluster2.bed Ne.down_WDNG_k4me1.cluster3.bed Ne.down_WDNG_k4me1.cluster4.bed -S LN.H3K4me1.cpm.bw DM.H3K4me1.cpm.bw NE.H3K4me1.cpm.bw GMP.K4ME1.re.cpm_peaks_norm.bw LN.H3K27ac.cpm.bw DM.H3K27ac.cpm.bw NE.H3K27ac.cpm.bw GMP.H3K27AC.re.cpm_peaks_norm.bw LN.ATAC.cpm_peaks.bw DM.ATAC.cpm_peaks.bw NE.ATAC.cpm_peaks.bw GMP.ATAC_NEW.cpm_peaks.bw -o WDNG_all.gz --missingDataAsZero -p 20

plotHeatmap -m WDNG_all.gz -out WDNG_all_sortbyall_heatmap.pdf --colorMap Oranges Oranges Oranges Oranges Purples Purples Purples Purples Greens Greens Greens Greens --samplesLabel 'WT' 'DM' 'NEU' 'GMP' 'WT' 'DM' 'NEU' 'GMP' 'WT' 'DM' 'NEU' 'GMP' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel 'Gain-1' 'Gain-2' 'Loss-1' 'Loss-2' --zMin 0 0 0 0 0 0 0 0 0 0 0 0 --zMax 2.2 2.2 2.2 2.2 2.5 2.5 2.5 2.5 12 12 12 12 --refPointLabel 'center' --legendLocation 'none' --whatToShow 'heatmap and colorbar'
