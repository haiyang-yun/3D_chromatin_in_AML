library(GenomicRanges)
library(GenomicAlignments)

#######################################################################
# H3K4me1 peakset
#######################################################################
peakFile=read.delim("./H3K4me1.bed", header = FALSE, sep ="\t")
peaks=GRanges(seqnames=peakFile$V1, IRanges(peakFile$V2,peakFile$V3),tagcount=peakFile$V5,peakName=peakFile$V7)
peakSet=keepStandardChromosomes(trim(resize(peaks, width=1000+width(peaks), fix="center"))) #146421 peaks

#######################################################################
#  Map H3K4me1 read counts on H3K4me1 peaks
#######################################################################
H3K4me1_bamFilePath="/serenity/data/HAYUN/ChIP-seq/catalog/bamfiles/H3K4me1"
H3K4me1_bamFiles = dir(H3K4me1_bamFilePath, pattern="*.bam$", full.name=T)
H3K4me1_aln <- lapply(H3K4me1_bamFiles, readGAlignments)
H3K4me1_counts <- lapply(H3K4me1_aln, GenomicRanges::countOverlaps, query=peakSet)
names(H3K4me1_counts) <- gsub(".nodup.bam","",basename(H3K4me1_bamFiles))
H3K4me1_df=data.frame(as.list(H3K4me1_counts),stringsAsFactors=FALSE)
H3K4me1_rowsum=rowSums (H3K4me1_df, na.rm = FALSE)
H3K4me1_maxlevel=apply(H3K4me1_df, 1, max)
H3K4me1_rawcount=cbind(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),H3K4me1_df,RowSum=H3K4me1_rowsum, Max=H3K4me1_maxlevel)
write.table(H3K4me1_rawcount,"H3K4me1.rawcount.csv",row.names=FALSE,quote=FALSE,sep="\t")

### Counts per million 
H3K4me1_fac=1000000/ lengths(H3K4me1_aln, use.names = TRUE)
H3K4me1_cpm=round(t(t(H3K4me1_df)*H3K4me1_fac),6)
H3K4me1_cpm_rowsum=rowSums (H3K4me1_cpm, na.rm = FALSE)
H3K4me1_maxlevel.cpm=apply(H3K4me1_cpm, 1, max)
H3K4me1_cpm_count=data.frame(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),H3K4me1_cpm,RowSum=H3K4me1_cpm_rowsum,Max=H3K4me1_maxlevel.cpm)
write.table(H3K4me1_cpm_count,"H3K4me1.cpm.csv",row.names=FALSE,quote=FALSE,sep="\t")


#######################################################################
# H3K4me3 peakset
#######################################################################
peakFile=read.delim("./H3K4me3.bed", header = FALSE, sep ="\t")
peaks=GRanges(seqnames=peakFile$V1, IRanges(peakFile$V2,peakFile$V3),tagcount=peakFile$V5,peakName=peakFile$V7)
peakSet=keepStandardChromosomes(trim(resize(peaks, width=1000+width(peaks), fix="center"))) #39819 peaks
#mergedPeakSet=reduce(peakSet) #37643 peaks

#######################################################################
#  Map H3K4me3 read counts on H3K4me3 peaks
#######################################################################
H3K4me3_bamFilePath="/serenity/data/HAYUN/ChIP-seq/catalog/bamfiles/H3K4me3"
H3K4me3_bamFiles = dir(H3K4me3_bamFilePath, pattern="*.bam$", full.name=T)
H3K4me3_aln <- lapply(H3K4me3_bamFiles, readGAlignments)
H3K4me3_counts <- lapply(H3K4me3_aln, GenomicRanges::countOverlaps, query=peakSet)
names(H3K4me3_counts) <- gsub(".nodup.bam","",basename(H3K4me3_bamFiles))
H3K4me3_df=data.frame(as.list(H3K4me3_counts),stringsAsFactors=FALSE)
H3K4me3_rowsum=rowSums (H3K4me3_df, na.rm = FALSE)
H3K4me3_maxlevel=apply(H3K4me3_df, 1, max)
H3K4me3_rawcount=cbind(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),H3K4me3_df,RowSum=H3K4me3_rowsum, Max=H3K4me3_maxlevel)
write.table(H3K4me3_rawcount,"H3K4me3.rawcount.csv",row.names=FALSE,quote=FALSE,sep="\t")

### Counts per million
H3K4me3_fac=1000000/ lengths(H3K4me3_aln, use.names = TRUE)
H3K4me3_cpm=round(t(t(H3K4me3_df)*H3K4me3_fac),6)
H3K4me3_cpm_rowsum=rowSums (H3K4me3_cpm, na.rm = FALSE)
H3K4me3_maxlevel.cpm=apply(H3K4me3_cpm, 1, max)
H3K4me3_cpm_count=data.frame(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),H3K4me3_cpm,RowSum=H3K4me3_cpm_rowsum,Max=H3K4me3_maxlevel.cpm)
write.table(H3K4me3_cpm_count,"H3K4me3.cpm.csv",row.names=FALSE,quote=FALSE,sep="\t")


############################################################################
# Defining enhancer based on the overlap between H3K4me3 and H3K4m1 peaks
############################################################################
H3K4me3peakFile=read.delim("../H3K4me3/based_on_HOMER/H3K4me3.cpm.csv", header = TRUE, sep ="\t")
H3K4me3peaks=GRanges(seqnames=H3K4me3peakFile$chr, IRanges(H3K4me3peakFile$start,H3K4me3peakFile$end),maxlevel=H3K4me3peakFile$Max)
H3K4me3peaks_gt16=H3K4me3peaks[H3K4me3peaks$maxlevel>16,]
enhancers=peakSet[!peakSet%over%H3K4me3peaks_gt16] #123231 enhancerPeaks
df.enhancer=data.frame(enhancers,stringsAsFactors=FALSE)
write.table(df.enhancer,"enhancerPeaks.bed",row.names=FALSE,quote=FALSE,sep="\t")
df.enhancer_merged=data.frame(reduce(enhancers),stringsAsFactors=FALSE) #98365 enhancerPeaks.merged
write.table(df.enhancer_merged,"enhancerPeaks.merged.bed",row.names=FALSE,quote=FALSE,sep="\t")