library("edgeR")
library("gplots")

getFit <- function(data, idx, lv, all=FALSE){
  g <- gsub("(.*)\\.CHiC.*", "\\1", colnames(data))
  group <- factor(g[idx],levels=c(lv[1],lv[2]))
  design <- model.matrix(~group)
  data <- data[,idx, drop=FALSE]
  y <- DGEList(counts=data, group=group)
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
  write.table(DEg,file=paste0(lv[1],".",lv[2],".totalreads",".diffExp.edgeR.csv"), sep="\t", quote=F, row.names=TRUE)
}
data=read.delim("chinput.totalreads.txt",sep="\t",header=TRUE,row.names=1)
getFit(data, idx=c(1,2,3,4), lv=c("WT", "NPM1"))
getFit(data, idx=c(1,2,5,6), lv=c("WT", "FLT3"))
getFit(data, idx=c(1,2,7,8), lv=c("WT", "DM"))


## volcano plots
getList <- function(fileIn, fileOut1, maintitle){
    tmp <- read.delim(fileIn, sep="\t",header=TRUE)
    tmp$FDR[tmp$FDR < 1e-20 ] <- 1e-20
    tmp.up=subset(tmp, FDR < 0.05 & logFC > 0)
    tmp.down=subset(tmp, FDR < 0.05 & logFC < 0)
    ### scatter plot
    pdf(fileOut1, width = 4.5, height = 4.5)
    #png(fileOut1, res=300, units = "in", width = 4.5, height = 4.5)
    par(mar=c(4.5,5,3,1),mgp=c(3.2,1,0),bty="l",las=1)
    with(tmp, plot(logFC, -log10(FDR), pch=20, col="grey", xlim=c(-7, 7), ylim=c(0,20),
    cex=1.2, cex.lab = 1.5, cex.axis=1.5, font.axis=2,font.lab=2,
    xlab="log2 FC pCHiC reads Pro", ylab="-log10(adjP)",axes=FALSE, main=maintitle, cex.main=1.6, font.main=2))
    box(lwd=3)
    axis(1,labels=TRUE,at=c(-6,-3,0,3,6),lwd=3,font=2,cex.axis=1.5)
    axis(2,labels=TRUE,lwd=3,font=2,cex.axis=1.5)
    with(tmp.up, points(logFC, -log10(FDR), pch=20, col="#db3236",cex=1.2))
    with(tmp.down, points(logFC, -log10(FDR), pch=20, col="#4885ed",cex=1.2))
    Utext <- nrow(tmp.up)
    Dtext <- nrow(tmp.down)
    legend("toprigh", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
    legend("topleft", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
    abline(v=0, h=-log10(0.05), lty=2,lwd=2)
    dev.off()
}
getList("WT.NPM1.totalreads.diffExp.edgeR.knowngenes.csv", "WN.diffpchic.pro.volcano.pdf",
"NPM1 vs WT")
getList("WT.FLT3.totalreads.diffExp.edgeR.knowngenes.csv", "WF.diffpchic.pro.volcano.pdf",
"FLT3 vs WT")
getList("WT.DM.totalreads.diffExp.edgeR.knowngenes.csv", "WD.diffpchic.pro.volcano.pdf",
"DM vs WT")

