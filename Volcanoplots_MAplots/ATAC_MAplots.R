library(GenomicRanges)

### ATAC-seq MA plots

atac.new <- read.delim("ATAC_consensus_peaks.bed", sep="\t", header=FALSE,row.names=4)
getList <- function(fileIn,fileOut,plottitle){
    tmp0 <- read.delim(fileIn, sep="\t", header=TRUE, row.names=1)
    tmp1 = tmp0[row.names(tmp0) %in% row.names(atac.new),]
    pdf(fileOut, width = 4.5, height = 4.5)
    par(mar=c(4.5,4.8,2,1.7),mgp=c(2.8,1,0),bty="l",las=1)
    with(tmp1, plot(logCPM, logFC, pch=20, col="grey",xlim=c(0.5,8.5),
    ylim = c(-6, 6), cex=0.5, cex.lab = 1.5, cex.main=1.5,font.main=2,
    cex.axis=1.2, font.axis=2, font.lab=2, ylab="log2 FC",
    xlab="log2 CPM ATAC-seq",main=plottitle, axes=FALSE))
    box(lwd=3)
    axis(1,at=c(0,2,4,6,8),labels=NULL,lwd=3,font=2,cex.axis=1.5)
    axis(2,at=c(-6,-4,-2,0,2,4,6),labels=NULL,lwd=3,font=2,cex.axis=1.5)
    with(subset(tmp1, logFC >= 1), points(logCPM, logFC, pch=20, col="#db3236",cex=0.5))
    with(subset(tmp1, logFC <= -1), points(logCPM, logFC, pch=20, col="#4885ed",cex=0.5))
    up <- subset(tmp1, logFC >= 1)
    down <- subset(tmp1, logFC <= -1)
    Utext <- NROW(up)
    Dtext <- NROW(down)
    abline(h=c(1,-1),lty=2,lwd=2)
    legend("topright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
    legend("bottomright", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
    dev.off()
}
getList("WT.NPM1.diff_ATAC.edgeR.FDR0.05.csv","WN.ATAC.diff_MAplot.pdf","Npm1c vs WT")
getList("WT.FLT3.diff_ATAC.edgeR.FDR0.05.csv","WF.ATAC.diff_MAplot.pdf","Npm1c vs WT")
getList("WT.DM.diff_ATAC.edgeR.FDR0.05.csv","WD.ATAC.diff_MAplot.pdf","Npm1c vs WT")

