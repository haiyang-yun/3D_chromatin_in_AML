library(GenomicRanges)

######################## Differential chromatin marks
## define atac peaks and enhancer H3K4me1 peaks
atac.summit_co20 <- read.table("ATAC.consensus.peak.summit.bed", sep="\t",header=F)
atac <- read.delim("ATAC_consensus_peaks_filter.bed", sep="\t", header=TRUE)
atac$id <- paste0(atac$chr,"_",atac$start,"_",atac$end)
atac.peakset=GRanges(seqnames=atac$chr, IRanges(atac$start,atac$end))

atac.ov.atac.summit_co20 = atac[atac$id %in% atac.summit_co20$id,]

atac.summit_peak = GRanges(seqnames=atac.ov.atac.summit_co20$chr,IRanges(atac.ov.atac.summit_co20$start,atac.ov.atac.summit_co20$end),id=atac.ov.atac.summit_co20$id)

### H3K4me1
getList <- function(fileIn,fileOut1,group){
    tmp=read.delim(fileIn, sep="\t", header=TRUE, row.names=1)
    tmp1=strsplit(as.character(row.names(tmp)), "_") ### starsplit can only split character factors
    tmp2=GRanges(seqnames=sapply(tmp1, "[", 1), IRanges(as.numeric(sapply(tmp1, "[", 2)), as.numeric(sapply(tmp1, "[", 3))),
    logFC=tmp$logFC,logCPM=tmp$logCPM,FDR=tmp$FDR)
    tmp3=tmp2[tmp2 %over% atac.summit_peak,]
    tmp4=data.frame(seqnames(tmp3),start(tmp3),end(tmp3),tmp3$logFC,tmp3$logCPM,tmp3$FDR,
    stringsAsFactors = FALSE)
    colnames(tmp4)=c("chr","start","end","logFC","logCPM","FDR")
    ## MA plots
    pdf(fileOut1, width = 4.5, height = 4.5)
    par(mar=c(4.5,4.5,3,1),mgp=c(2.8,1,0),bty="l",las=1)
    with(tmp4, plot(logCPM, logFC, pch=20, col="grey",xlim=c(1,7),
    ylim = c(-4, 4),cex=1, cex.lab = 1.5, cex.main=1.5,font.main=2,
    cex.axis=1.2, font.axis=2,font.lab=2,
    xlab="log2 CPM H3K4me1",ylab="log2 FC",main=group,axes=FALSE))
    box(lwd=3)
    axis(1,at=c(2,4,6),labels=NULL,lwd=3,font=2,cex.axis=1.5)
    axis(2,at=c(-4,-2,0,2,4),labels=NULL,lwd=3,font=2,cex.axis=1.5)
    with(subset(tmp4, logFC >= log(1.5,2)), points(logCPM, logFC, pch=20, col="#db3236",cex=1))
    with(subset(tmp4, logFC <= -log(1.5,2)), points(logCPM, logFC, pch=20, col="#4885ed",cex=1))
    up <- subset(tmp4, logFC >= log(1.5,2))
    down <- subset(tmp4, logFC <= -log(1.5,2))
    Utext <- NROW(up)
    Dtext <- NROW(down)
    legend("topright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
    legend("bottomright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
    abline(h=c(-log(1.5,2),log(1.5,2)), lty=2,lwd=2)
    dev.off()
}
getList("WT.NPM1.diff_H3K4me1.edgeR.FDR0.05.csv","WN.H3K4me1.diff_MAplot.pdf","NPM1 vs WT")
getList("WT.FLT3.diff_H3K4me1.edgeR.FDR0.05.csv","WF.H3K4me1.diff_MAplot.pdf","FLT3 vs WT")
getList("WT.DM.diff_H3K4me1.edgeR.FDR0.05.csv","WD.H3K4me1.diff_MAplot.pdf","DM vs WT")



## define atac peaks and enhancer H3K27ac peaks
enh.all <- read.delim("enhancerPeaks.bed", sep="\t", header=T)
enh.all.peakset <- GRanges(seqnames = enh.all$chr, ranges = IRanges(enh.all$start,enh.all$end))

atac.summit.ov.enh = atac.summit_co20[atac.summit_co20$id %in% atac.summit_peak.ov.enh$id, ]
atac.summit.ov.enh_peak = GRanges(seqnames = atac.summit.ov.enh$chr,ranges = IRanges(atac.summit.ov.enh$start,atac.summit.ov.enh$end),id = atac.summit.ov.enh$id)

atac.summit500bp.ov.enh = resize(atac.summit.ov.enh_peak, width = 500, fix = "center")

enh.ov.atac.summit500 = enh.all.peakset[enh.all.peakset %over% atac.summit500bp.ov.enh, ]

#### H3K27ac
getList <- function(fileIn,fileOut1,group){
  tmp=read.delim(fileIn, sep="\t", header=TRUE, row.names=1)
  tmp1=strsplit(as.character(row.names(tmp)), "_") ### starsplit can only split character factors
  tmp2=GRanges(seqnames=sapply(tmp1, "[", 1), IRanges(as.numeric(sapply(tmp1, "[", 2)), as.numeric(sapply(tmp1, "[", 3))),
               logFC=tmp$logFC,logCPM=tmp$logCPM,FDR=tmp$FDR)
  tmp3=tmp2[tmp2 %over% enh.ov.atac.summit500,]
  tmp4=data.frame(seqnames(tmp3),start(tmp3),end(tmp3),tmp3$logFC,tmp3$logCPM,tmp3$FDR,
                  stringsAsFactors = FALSE)
  colnames(tmp4)=c("chr","start","end","logFC","logCPM","FDR")
  ## MA plots
  pdf(fileOut1, width = 4.5, height = 4.5)
  par(mar=c(4.5,4.5,3,1),mgp=c(2.8,1,0),bty="l",las=1)
  with(tmp4, plot(logCPM, logFC, pch=20, col="grey",xlim=c(1.5,8.5), 
                  ylim = c(-4.5, 4.5),cex=1, cex.lab = 1.5, cex.main=1.5,font.main=2,
                  cex.axis=1.2, font.axis=2,font.lab=2, 
                  xlab="log2 CPM H3K27ac Enh",ylab="log2 FC",main=group,axes=FALSE))
  box(lwd=3)
  axis(1,at=c(2,4,6,8),labels=NULL,lwd=3,font=2,cex.axis=1.5)
  axis(2,at=c(-4,-2,0,2,4),labels=NULL,lwd=3,font=2,cex.axis=1.5)
  with(subset(tmp4, logFC >= log(1.5,2)), points(logCPM, logFC, pch=20, col="#db3236",cex=1))
  with(subset(tmp4, logFC <= -log(1.5,2)), points(logCPM, logFC, pch=20, col="#4885ed",cex=1))
  up <- subset(tmp4, logFC >= log(1.5,2))
  down <- subset(tmp4, logFC <= -log(1.5,2))
  Utext <- NROW(up)
  Dtext <- NROW(down)
  legend("topright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
  legend("bottomright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
  abline(h=c(-log(1.5,2),log(1.5,2)), lty=2,lwd=2)
  dev.off()
}
getList("WT.NPM1.diff_H3K27ac.edgeR.FDR0.05.csv","WN.H3K27ac.diff_MAplot.pdf","NPM1 vs WT")
getList("WT.FLT3.diff_H3K27ac.edgeR.FDR0.05.csv","WF.H3K27ac.diff_MAplot.pdf","FLT3 vs WT")
getList("WT.DM.diff_H3K27ac.edgeR.FDR0.05.csv","WD.H3K27ac.diff_MAplot.pdf","DM vs WT")


