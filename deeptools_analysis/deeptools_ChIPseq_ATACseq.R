### deeptools heatmap and profile plots

## ATAC
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R WD.ATAC.up_summit500.bed WD.ATAC.down_summit500.bed -S LN.ATAC.R1-2.cpm_peaks.bw NPM1.ATAC.R1-2.cpm_peaks.bw FLT3.ATAC.R1-2.cpm_peaks.bw DM.ATAC.R1-2.cpm_peaks.bw -o WNFD_ATAC.gz --missingDataAsZero -p 12

plotHeatmap -m WNFD_ATAC.gz -out WNFD_ATAC_heatmap.pdf --colorMap Greens --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel 'Accessibility gain in DM' 'Accessibility loss in DM' --refPointLabel 'center' --legendLocation 'none' --whatToShow 'heatmap and colorbar'

plotProfile -m WNFD_ATAC.gz -out WNFD_ATAC_profile.pdf --perGroup --numPlotsPerRow 1 --plotWidth 7.25 --plotHeight 6.5 --dpi 300 --refPointLabel 'center' --legendLocation 'upper-left' --plotType 'se' -y 'ATAC-seq enrichment' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --samplesLabel 'WT' 'NPM1' 'FLT3' 'DM' --regionsLabel 'Accessibility gain in DM' 'Accessibility loss in DM'

### H3K4me1
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC.summit_ov.WD.H3K4me1.up.bed ATAC.summit_ov.WD.H3K4me1.down.bed -S LN.H3K4me1.R1-2.cpm_peaks.bw NPM1.H3K4me1.R1-2.cpm_peaks.bw FLT3.H3K4me1.R1-2.cpm_peaks.bw DM.H3K4me1.R1-2.cpm_peaks.bw -o WNFD_H3K4me1.gz --missingDataAsZero -p 16

plotHeatmap -m WNFD_H3K4me1.gz -out WNFD_H3K4me1_heatmap.pdf --colorMap Oranges --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel 'H3K4me1 gain in DM' 'H3K4me1 loss in DM' --refPointLabel 'center' --legendLocation 'none' --whatToShow 'heatmap and colorbar'

plotProfile -m WNFD_H3K4me1.gz -out WNFD_H3K4me1_profile.pdf --perGroup --numPlotsPerRow 1 --plotWidth 7.25 --plotHeight 6.5 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' --yMin 0 0 --yMax 4 4 -y 'H3K4me1 enrichment' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --regionsLabel 'H3K4me1 gain in DM' 'H3K4me1 loss in DM'

### H3K27ac
computeMatrix reference-point --referencePoint center -b 1000 -a 1000 -R ATAC.summit_ov.WD.H3K27ac.up.bed  ATAC.summit_ov.WD.H3K27ac.down.bed -S LN.H3K27ac.R1-2.cpm_peaks.bw NPM1.H3K27ac.R1-2.cpm_peaks.bw FLT3.H3K27ac.R1-2.cpm_peaks.bw DM.H3K27ac.R1-2.cpm_peaks.bw -o WNFD_H3K27ac.gz --missingDataAsZero -p 16

plotHeatmap -m WNFD_H3K27ac.gz -out WNFD_H3K27ac_heatmap.pdf --colorMap Purples --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --heatmapHeight 10 --heatmapWidth 1.5 --dpi 300 -x '' --regionsLabel 'H3K27ac gain in DM' 'H3K27ac loss in DM' --refPointLabel 'center' --legendLocation 'none' --whatToShow 'heatmap and colorbar'

plotProfile -m WNFD_H3K27ac.gz -out WNFD_H3K27ac_profile.pdf --perGroup --numPlotsPerRow 1 --plotWidth 7.25 --plotHeight 6.5 --dpi 300 --refPointLabel 'center' --plotType 'se' --legendLocation 'upper-left' --yMin 0 0 --yMax 14 14 -y 'H3K27ac enrichment' --colors '#AAA9AD' '#0078D7' '#3CAEA3' '#ED553B' --samplesLabel 'WT' 'Npm1c' 'Flt3-ITD' 'DM' --regionsLabel 'H3K27ac gain in DM' 'H3K27ac loss in DM'

