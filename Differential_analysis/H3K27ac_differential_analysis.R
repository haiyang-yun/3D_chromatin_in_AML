### Differential analysis by edgeR

library(GenomicRanges)
library(GenomicAlignments)
library("edgeR")
library("gplots")

getFit <- function(data, idx, lv, libsize){
  group <- factor(c(1,1,2,2))
  design <- model.matrix(~group)
  data <- data[,idx, drop=FALSE]
  y <- DGEList(counts=data, group=group,lib.size=libsize)
  keep <- rowSums(cpm(y)>1) >= 2
  y <- y[keep, , keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  y <- estimateGLMCommonDisp(y, design)
  y <- estimateGLMTrendedDisp(y, design)
  y <- estimateGLMTagwiseDisp(y, design)
  fit <- glmFit(y, design)
  lrt <- glmLRT(fit, coef=2)
  DEg <- as.data.frame(topTags(lrt, n=10^6))
  write.table(DEg,file=paste0(lv[1],".",lv[2],".diff_H3K27ac.edgeR.csv"), sep="\t", quote=F, row.names=TRUE)
}

peakFile=read.delim("./H3K27ac_consensus_peaks.bed", header = TRUE, sep ="\t")
peakSet=GRanges(seqnames=peakFile$chr, IRanges(peakFile$start,peakFile$end))
bamFilePath="/serenity/data/Haiyang/ChIP_HY/bamfiles/H3K27ac"
bamFiles = dir(bamFilePath, pattern="*.bam$", full.name=T)
aln <- lapply(bamFiles, readGAlignments)
counts <- lapply(aln, GenomicRanges::countOverlaps, query=peakSet)
names(counts) <- gsub(".nodup.bam","",basename(bamFiles))

df=data.frame(as.list(counts),stringsAsFactors=FALSE)
maxlevel=apply(df, 1, max)
rawcount=cbind(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),df,Width=width(peakSet), Max=maxlevel)

liblen=lengths(aln, use.names = TRUE)
data=df
peakId <- with(peakFile, paste(chr,"_",start,"_",end,sep=""))
row.names(data)=peakId

getFit(data, idx=c(5,6,7,8), lv=c("WT", "NPM1"), libsize=liblen[c(5,6,7,8)])
getFit(data, idx=c(5,6,3,4), lv=c("WT", "FLT3"), libsize=liblen[c(5,6,3,4)])
getFit(data, idx=c(5,6,1,2), lv=c("WT", "DM"), libsize=liblen[c(5,6,1,2)])