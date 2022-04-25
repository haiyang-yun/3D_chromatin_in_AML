########## Compute ATAC-seq consensus peak sets ##########
library(DiffBind)
library(GenomicRanges)
library(GenomicAlignments)

data = dba(sampleSheet="samplesheet_ATAC.csv")
data_consensus = dba.peakset(data, consensus=DBA_CONDITION, minOverlap=2)
data_consensus_all = dba(data_consensus, mask=data_consensus$masks$Consensus,minOverlap=1)
consensus_peaks_all = dba.peakset(data_consensus_all, bRetrieve=TRUE)
clist = data.frame(seqnames(consensus_peaks_all),start(consensus_peaks_all),end(consensus_peaks_all))
colnames(clist) = c("chr","start","end")
norchr = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chrX","chrY")
clist.filter = clist[is.element(clist$chr, norchr),]
write.table(clist.filter, "ATAC_consensus_peaks.bed", quote = FALSE, col.names = TRUE, row.names=FALSE, sep = "\t")

########## supplement info for peak maximal sample ##########
peakFile=read.delim("./ATAC_consensus_peaks.bed", header = TRUE, sep ="\t")
peakSet=GRanges(seqnames=peakFile$chr, IRanges(peakFile$start,peakFile$end))
bamFiles = dir([bamFilePath], pattern="*.bam$", full.name=T) ## set bamFilePath 
aln <- lapply(bamFiles, readGAlignments)
counts <- lapply(aln, GenomicRanges::countOverlaps, query=peakSet)
df=data.frame(as.list(counts),stringsAsFactors=FALSE)

### normalised read counts (RPKM)
fac=1000000/ lengths(aln, use.names = TRUE)
cpm=round(t(t(df)*fac),6)
rpkm=(cpm*1000)/width(peakSet)
atac.rpkm=data.frame(chr=seqnames(peakSet),start=start(ranges(peakSet)),end=end(ranges(peakSet)),round(rpkm,6))
atac.rpkm$max <- apply(atac.rpkm[,4:11], 1, max)    
atac.rpkm$max.ID <- colnames(atac.rpkm[,4:11])[apply(atac.rpkm[,4:11],1,which.max)]
atac.rpkm$peak.ID <- paste0(atac.rpkm$chr,":",atac.rpkm$start,"-",atac.rpkm$end) 
write.table(atac.rpkm[,c(1,2,3,14,13)], "ATAC_consensus_peakmax.bed", append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
