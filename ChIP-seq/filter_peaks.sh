
### DM ###

mergePeaks -d 500 DM.H3K4me1.R1.peaks DM.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks | awk '{print $9}'|sort >  DM.H3K4me1.R1.indices
./filterpeaks_based_on_index.py DM.H3K4me1.R1.indices DM.H3K4me1.R1.peaks  > DM.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks | awk '{print $10}'|sort >  DM.H3K4me1.R2.indices
./filterpeaks_based_on_index.py DM.H3K4me1.R2.indices DM.H3K4me1.R2.peaks  > DM.H3K4me1.R2.peaks.overlap
grep -v \# DM.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"DM.H3K4me1.R1.peaks"}' > DM.H3K4me1.R1.bed 
grep -v \# DM.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"DM.H3K4me1.R2.peaks"}' > DM.H3K4me1.R2.bed
cat DM.H3K4me1.R1.bed DM.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  DM.H3K4me1.R1_DM.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py DM.H3K4me1.R1_DM.H3K4me1.R2.sorted.bed > DM.H3K4me1.bed

### WT ###
mergePeaks -d 500 WT.H3K4me1.R1.peaks WT.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_WT.H3K4me1.R1.peaks_WT.H3K4me1.R2.peaks | awk '{print $9}'|sort >  WT.H3K4me1.R1.indices
./filterpeaks_based_on_index.py WT.H3K4me1.R1.indices WT.H3K4me1.R1.peaks  > WT.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_WT.H3K4me1.R1.peaks_WT.H3K4me1.R2.peaks | awk '{print $10}'|sort >  WT.H3K4me1.R2.indices
./filterpeaks_based_on_index.py WT.H3K4me1.R2.indices WT.H3K4me1.R2.peaks  > WT.H3K4me1.R2.peaks.overlap
grep -v \# WT.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"WT.H3K4me1.R1.peaks"}' > WT.H3K4me1.R1.bed 
grep -v \# WT.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"WT.H3K4me1.R2.peaks"}' > WT.H3K4me1.R2.bed
cat WT.H3K4me1.R1.bed WT.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  WT.H3K4me1.R1_WT.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py WT.H3K4me1.R1_WT.H3K4me1.R2.sorted.bed > WT.H3K4me1.bed

### NPM1 ###
mergePeaks -d 500 NPM1.H3K4me1.R1.peaks NPM1.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks | awk '{print $9}'|sort >  NPM1.H3K4me1.R1.indices
./filterpeaks_based_on_index.py NPM1.H3K4me1.R1.indices NPM1.H3K4me1.R1.peaks  > NPM1.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks | awk '{print $10}'|sort >  NPM1.H3K4me1.R2.indices
./filterpeaks_based_on_index.py NPM1.H3K4me1.R2.indices NPM1.H3K4me1.R2.peaks  > NPM1.H3K4me1.R2.peaks.overlap
grep -v \# NPM1.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"NPM1.H3K4me1.R1.peaks"}' > NPM1.H3K4me1.R1.bed 
grep -v \# NPM1.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"NPM1.H3K4me1.R2.peaks"}' > NPM1.H3K4me1.R2.bed
cat NPM1.H3K4me1.R1.bed NPM1.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  NPM1.H3K4me1.R1_NPM1.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py NPM1.H3K4me1.R1_NPM1.H3K4me1.R2.sorted.bed > NPM1.H3K4me1.bed

### FLT3 ###
mergePeaks -d 500 FLT3.H3K4me1.R1.peaks FLT3.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks | awk '{print $9}'|sort >  FLT3.H3K4me1.R1.indices
./filterpeaks_based_on_index.py FLT3.H3K4me1.R1.indices FLT3.H3K4me1.R1.peaks  > FLT3.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks | awk '{print $10}'|sort >  FLT3.H3K4me1.R2.indices
./filterpeaks_based_on_index.py FLT3.H3K4me1.R2.indices FLT3.H3K4me1.R2.peaks  > FLT3.H3K4me1.R2.peaks.overlap
grep -v \# FLT3.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"FLT3.H3K4me1.R1.peaks"}' > FLT3.H3K4me1.R1.bed 
grep -v \# FLT3.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"FLT3.H3K4me1.R2.peaks"}' > FLT3.H3K4me1.R2.bed
cat FLT3.H3K4me1.R1.bed FLT3.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  FLT3.H3K4me1.R1_FLT3.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py FLT3.H3K4me1.R1_FLT3.H3K4me1.R2.sorted.bed > FLT3.H3K4me1.bed

###concatenate all the bed files
cat DM.H3K4me1.bed  FLT3.H3K4me1.bed  WT.H3K4me1.bed  NPM1.H3K4me1.bed | sort -k1,1 -k2,2n >H3K4me1.allpeaks.bed

### get the catalog
./get_catalog.py H3K4me1.allpeaks.bed > H3K4me1.bed

