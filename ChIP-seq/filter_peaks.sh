
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

### LN ###
mergePeaks -d 500 LN.H3K4me1.R1.peaks LN.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks | awk '{print $9}'|sort >  LN.H3K4me1.R1.indices
./filterpeaks_based_on_index.py LN.H3K4me1.R1.indices LN.H3K4me1.R1.peaks  > LN.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks | awk '{print $10}'|sort >  LN.H3K4me1.R2.indices
./filterpeaks_based_on_index.py LN.H3K4me1.R2.indices LN.H3K4me1.R2.peaks  > LN.H3K4me1.R2.peaks.overlap
grep -v \# LN.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"LN.H3K4me1.R1.peaks"}' > LN.H3K4me1.R1.bed 
grep -v \# LN.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"LN.H3K4me1.R2.peaks"}' > LN.H3K4me1.R2.bed
cat LN.H3K4me1.R1.bed LN.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  LN.H3K4me1.R1_LN.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py LN.H3K4me1.R1_LN.H3K4me1.R2.sorted.bed > LN.H3K4me1.bed

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
cat DM.H3K4me1.bed  FLT3.H3K4me1.bed  LN.H3K4me1.bed  NPM1.H3K4me1.bed | sort -k1,1 -k2,2n >H3K4me1.allpeaks.bed

### get the catalog
./get_catalog.py H3K4me1.allpeaks.bed > H3K4me1.bed

### NE #
mergePeaks -d 500 NE.H3K4me1.R1.peaks NE.H3K4me1.R2.peaks -prefix mergePeaks
grep -v \# mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks | awk '{print $9}'|sort >  NE.H3K4me1.R1.indices
./filterpeaks_based_on_index.py NE.H3K4me1.R1.indices NE.H3K4me1.R1.peaks  > NE.H3K4me1.R1.peaks.overlap
grep -v \# mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks | awk '{print $10}'|sort >  NE.H3K4me1.R2.indices
./filterpeaks_based_on_index.py NE.H3K4me1.R2.indices NE.H3K4me1.R2.peaks  > NE.H3K4me1.R2.peaks.overlap
grep -v \# NE.H3K4me1.R1.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"NE.H3K4me1.R1.peaks"}' > NE.H3K4me1.R1.bed 
grep -v \# NE.H3K4me1.R2.peaks.overlap |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8,"NE.H3K4me1.R2.peaks"}' > NE.H3K4me1.R2.bed
cat NE.H3K4me1.R1.bed NE.H3K4me1.R2.bed | sort -k1,1 -k2,2n >  NE.H3K4me1.R1_NE.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py NE.H3K4me1.R1_NE.H3K4me1.R2.sorted.bed > NE.H3K4me1.bed

cat LN.H3K4me1.bed  NE.H3K4me1.bed | sort -k1,1 -k2,2n >H3K4me1.LN_NEpeaks.bed
./get_catalog.py H3K4me1.LN_NEpeaks.bed > H3K4me1.LN_NE.bed
#########################################
# OLD COMMANDS
#########################################
# Find overlapping peaks
mergePeaks -d 500 DM.H3K4me1.R1.peaks DM.H3K4me1.R2.peaks -prefix mergePeaks
mergePeaks -d 500 LN.H3K4me1.R1.peaks LN.H3K4me1.R2.peaks -prefix mergePeaks
mergePeaks -d 500 NPM1.H3K4me1.R1.peaks NPM1.H3K4me1.R2.peaks -prefix mergePeaks
mergePeaks -d 500 NE.H3K4me1.R1.peaks NE.H3K4me1.R2.peaks -prefix mergePeaks
mergePeaks -d 500 FLT3.H3K4me1.R1.peaks FLT3.H3K4me1.R2.peaks -prefix mergePeaks
# Extract the genomic ranges from the overlap file
awk '{OFS="\t"};{print $2,$3,$4}' DM.H3K4me1_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks >DM.H3K4me1.R1_DM.H3K4me1.R2_overlap.bed && sed -i '1,+0d' DM.H3K4me1.R1_DM.H3K4me1.R2_overlap.bed

# convert homer peak file into bed file
grep -v \# DM.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8}' > DM.H3K4me1.R1.bed 
grep -v \# DM.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,$1,$6,$8}' > DM.H3K4me1.R2.bed

# Filter the overlapped peaks from each replicates
bedtools intersect -a  DM.H3K4me1.R1.bed -b DM.H3K4me1.R1_DM.H3K4me1.R2_overlap.bed -wa > DM.H3K4me1.R1.overlapped.bed
bedtools intersect -a  DM.H3K4me1.R2.bed -b DM.H3K4me1.R1_DM.H3K4me1.R2_overlap.bed -wa > DM.H3K4me1.R2.overlapped.bed

# concatenate both the files with overlapped peak ranges and sort it
cat DM.H3K4me1.R1.overlapped.bed DM.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  DM.H3K4me1.R1.R2.overlap.bed

# Filter the oerlapped peaks to keep the one with stronger signal
./filter_overlapped_intervals.py DM.H3K4me1.R1.R2.overlap.bed > DM.H3K4me1.bed
#########
# DM
#########
awk '{OFS="\t"};{print $2,$3,$4}' mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks >mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks.bed && sed -i '1,+0d' mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks.bed
grep -v \# DM.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,"DM.H3K4me1.R1.peaks",$6,$8}' > DM.H3K4me1.R1.bed 
grep -v \# DM.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,"DM.H3K4me1.R2.peaks",$6,$8}' > DM.H3K4me1.R2.bed
bedtools intersect -a  DM.H3K4me1.R1.bed -b mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks.bed -wa > DM.H3K4me1.R1.overlapped.bed
bedtools intersect -a  DM.H3K4me1.R2.bed -b mergePeaks_DM.H3K4me1.R1.peaks_DM.H3K4me1.R2.peaks.bed -wa > DM.H3K4me1.R2.overlapped.bed
cat DM.H3K4me1.R1.overlapped.bed DM.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  DM.H3K4me1.R1_DM.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py DM.H3K4me1.R1_DM.H3K4me1.R2.sorted.bed > DM.H3K4me1.bed

#########
# LN
#########
awk '{OFS="\t"};{print $2,$3,$4}' mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks >mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks.bed && sed -i '1,+0d' mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks.bed
grep -v \# LN.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,"LN.H3K4me1.R1.peaks",$6,$8}' > LN.H3K4me1.R1.bed 
grep -v \# LN.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,"LN.H3K4me1.R2.peaks",$6,$8}' > LN.H3K4me1.R2.bed
bedtools intersect -a  LN.H3K4me1.R1.bed -b mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks.bed -wa > LN.H3K4me1.R1.overlapped.bed
bedtools intersect -a  LN.H3K4me1.R2.bed -b mergePeaks_LN.H3K4me1.R1.peaks_LN.H3K4me1.R2.peaks.bed -wa > LN.H3K4me1.R2.overlapped.bed
cat LN.H3K4me1.R1.overlapped.bed LN.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  LN.H3K4me1.R1_LN.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py LN.H3K4me1.R1_LN.H3K4me1.R2.sorted.bed > LN.H3K4me1.bed

#########
# NPM1
#########
awk '{OFS="\t"};{print $2,$3,$4}' mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks >mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks.bed && sed -i '1,+0d' mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks.bed
grep -v \# NPM1.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,"NPM1.H3K4me1.R1.peaks",$6,$8}' > NPM1.H3K4me1.R1.bed 
grep -v \# NPM1.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,"NPM1.H3K4me1.R2.peaks",$6,$8}' > NPM1.H3K4me1.R2.bed
bedtools intersect -a  NPM1.H3K4me1.R1.bed -b mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks.bed -wa > NPM1.H3K4me1.R1.overlapped.bed
bedtools intersect -a  NPM1.H3K4me1.R2.bed -b mergePeaks_NPM1.H3K4me1.R1.peaks_NPM1.H3K4me1.R2.peaks.bed -wa > NPM1.H3K4me1.R2.overlapped.bed
cat NPM1.H3K4me1.R1.overlapped.bed NPM1.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  NPM1.H3K4me1.R1_NPM1.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py NPM1.H3K4me1.R1_NPM1.H3K4me1.R2.sorted.bed > NPM1.H3K4me1.bed

#########
# FLT3
#########
awk '{OFS="\t"};{print $2,$3,$4}' mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks >mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks.bed && sed -i '1,+0d' mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks.bed
grep -v \# FLT3.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,"FLT3.H3K4me1.R1.peaks",$6,$8}' > FLT3.H3K4me1.R1.bed 
grep -v \# FLT3.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,"FLT3.H3K4me1.R2.peaks",$6,$8}' > FLT3.H3K4me1.R2.bed
bedtools intersect -a  FLT3.H3K4me1.R1.bed -b mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks.bed -wa > FLT3.H3K4me1.R1.overlapped.bed
bedtools intersect -a  FLT3.H3K4me1.R2.bed -b mergePeaks_FLT3.H3K4me1.R1.peaks_FLT3.H3K4me1.R2.peaks.bed -wa > FLT3.H3K4me1.R2.overlapped.bed
cat FLT3.H3K4me1.R1.overlapped.bed FLT3.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  FLT3.H3K4me1.R1_FLT3.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py FLT3.H3K4me1.R1_FLT3.H3K4me1.R2.sorted.bed > FLT3.H3K4me1.bed

###concatenate all the bed files
cat DM.H3K4me1.bed  FLT3.H3K4me1.bed  LN.H3K4me1.bed  NPM1.H3K4me1.bed | sort -k1,1 -k2,2n >H3K4me1.allpeaks.bed

### get the catalog
./get_catalog.py H3K4me1.allpeaks.bed > H3K4me1.bed
#########
# NE
#########
awk '{OFS="\t"};{print $2,$3,$4}' mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks >mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks.bed && sed -i '1,+0d' mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks.bed
grep -v \# NE.H3K4me1.R1.peaks |awk '{OFS="\t"};{print $2,$3,$4,"NE.H3K4me1.R1.peaks",$6,$8}' > NE.H3K4me1.R1.bed 
grep -v \# NE.H3K4me1.R2.peaks |awk '{OFS="\t"};{print $2,$3,$4,"NE.H3K4me1.R2.peaks",$6,$8}' > NE.H3K4me1.R2.bed
bedtools intersect -a  NE.H3K4me1.R1.bed -b mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks.bed -wa > NE.H3K4me1.R1.overlapped.bed
bedtools intersect -a  NE.H3K4me1.R2.bed -b mergePeaks_NE.H3K4me1.R1.peaks_NE.H3K4me1.R2.peaks.bed -wa > NE.H3K4me1.R2.overlapped.bed
cat NE.H3K4me1.R1.overlapped.bed NE.H3K4me1.R2.overlapped.bed | sort -k1,1 -k2,2n >  NE.H3K4me1.R1_NE.H3K4me1.R2.sorted.bed
./filter_overlapped_intervals.py NE.H3K4me1.R1_NE.H3K4me1.R2.sorted.bed > NE.H3K4me1.bed

###concatenate all the bed files with NE
cat ../DM.H3K4me1.bed  ../FLT3.H3K4me1.bed  ../LN.H3K4me1.bed  ../NPM1.H3K4me1.bed  ../NE.H3K4me1.bed | sort -k1,1 -k2,2n >H3K4me1.allpeaks.bed

### get the catalog
../get_catalog.py H3K4me1.allpeaks.bed > H3K4me1.bed


