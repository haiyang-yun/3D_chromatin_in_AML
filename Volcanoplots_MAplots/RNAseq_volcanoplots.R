setwd("~/Documents/20180920_manuscript/RNA-seq")
library(ggplot2)

getList <- function(fileIn, fileOut1, maintitle){
  tmp <- read.delim(fileIn, sep="\t",header=TRUE)
  tmp.up1=subset(tmp, padj < 0.05 & log2FC >= log(1.5,2))
  tmp.down1=subset(tmp, padj < 0.05 & log2FC <= -log(1.5,2))
  ### scatter plot
  tmp$padj[tmp$padj < 1e-100 ] <- 1e-100
  tmp[tmp == NA ] <- 1
  tmp.up=subset(tmp, padj < 0.05 & log2FC >= log(1.5,2))
  tmp.down=subset(tmp, padj < 0.05 & log2FC <= -log(1.5,2))
  pdf(fileOut1, width = 4.5, height = 4.5)
  #png(fileOut1, res=300, units = "in", width = 4.5, height = 4.5)
  par(mar=c(4.5,5,3,1),mgp=c(2.8,1,0),bty="l",las=1)
  with(tmp, plot(log2FC, -log10(padj), pch=20, col="grey", xlim=c(-10, 10), ylim=c(0,115),
                     cex=1, cex.lab = 1.5, cex.axis=1.5, font.axis=2,font.lab=2, 
                     xlab="log2 FC RNA-seq", ylab="-log10 adjP",axes=FALSE, main=maintitle, cex.main=1.6, font.main=2))
  box(lwd=3)
  axis(1,labels=TRUE,lwd=3,font=2,cex.axis=1.5)
  axis(2,labels=TRUE,lwd=3,font=2,cex.axis=1.5)
  with(tmp.up, points(log2FC, -log10(padj), pch=20, col="#db3236",cex=1))
  with(tmp.down, points(log2FC, -log10(padj), pch=20, col="#4885ed",cex=1))
  Utext <- nrow(tmp.up)
  Dtext <- nrow(tmp.down)
  legend("toprigh", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
  legend("topleft", x.intersp=-0.5, y.intersp=0.2, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
  abline(v=c(log(1.5,2),-log(1.5,2)), h=-log(0.05,10), lty=2,lwd=2)
  dev.off()
}
getList("WT_NPM1_de_forvolcano.txt", "WN.diffExp.volcano.pdf", "NPM1 vs WT")
getList("WT_FLT3_de_forvolcano.txt", "WF.diffExp.volcano.pdf", "FLT3 vs WT")
getList("WT_DM_de_forvolcano.txt", "WD.diffExp.volcano.pdf", "DM vs WT")