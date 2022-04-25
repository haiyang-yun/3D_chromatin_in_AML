library(GenomicRanges)
library(ggplot2)
library(dplyr) 

##################### Seurat analysis on merged cpm from replicates #####################
## all HindIII fragment
hf <- read.delim("Digest.mm10.rmap", sep="\t",header=FALSE)
colnames(hf) <- c("chr","start","end","FragID")
hf.peak=GRanges(seqnames=hf$chr, IRanges(hf$start,hf$end),names=hf$FragID)

## PCHiC annotation
PCHiC.anno <- read.delim("pCHiC_fragID_Gene.txt", sep="\t",header=FALSE)
colnames(PCHiC.anno) <- c("ID","Gene")

## Interaction profiles
int.data <- read.delim("pCHiC_matrix.txt", sep="\t",header=TRUE)
int.data.oeID.peak= GRanges(seqnames=int.data$oeChr, IRanges(int.data$oeStart,int.data$oeEnd),baitID=int.data$baitID,oeID=int.data$oeID)
int.data.baitID.peak= GRanges(seqnames=int.data$baitChr, IRanges(int.data$baitStart,int.data$baitEnd),baitID=int.data$baitID,oeID=int.data$oeID)

## Differential expression DM vs WT
rna <- read.csv("WT.DM.PC.diffExp.csv", sep=",",header=TRUE, row.names=1)
wd.exp <- rna[,c("Gene_name","log2FoldChange","padj")]
colnames(wd.exp) <- c("Gene_name","log2FC","padj")

getList <- function(clusterID){
  tmp1 <- read.delim(paste0("Multiomics_Cluster-",clusterID,"_summit200bp.bed"), sep="\t",header=FALSE)
  colnames(tmp1)=c("chr","start","end")
  tmp2=GRanges(seqnames=tmp1$chr, IRanges(tmp1$start,tmp1$end))
  ########## annotate by interaction baitID
  tmp3=(int.data.oeID.peak[int.data.oeID.peak%over%tmp2,])$baitID
  tmp5=PCHiC.anno[PCHiC.anno$ID%in%tmp3,]$Gene
  ########### annotate by BaitID
  tmp4=(hf.peak[hf.peak%over%tmp2,])$names
  tmp6=PCHiC.anno[PCHiC.anno$ID%in%tmp4,]$Gene
  ########### combined annotation
  tmp7=union(tmp5,tmp6)
  ## diff expression
  res=wd.exp[wd.exp$Gene_name%in%tmp7,]
  write.table(res, paste0("Cluster-",clusterID,"_genes_DMvsWT_diffexp.txt"),quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
  ## heatmap
  res[is.na(res)] <- 1
  res$padj[res$padj < 1e-100 ] <- 1e-100
  res.up=subset(res, padj < 0.05 & log2FC >= log(1.5,2))
  res.down=subset(res, padj < 0.05 & log2FC <= -log(1.5,2))
  pdf(paste0("Cluster-",clusterID,"_genes_DMvsWT_volcanoplot.pdf"), width = 4.5, height = 4.5)
  par(mar=c(4.5,5,3,1),mgp=c(2.8,1,0),bty="l",las=1)
  with(res, plot(log2FC, -log10(padj), pch=20, col="grey", xlim=c(-10, 10), ylim=c(0,110),
                 cex=0.5, cex.lab = 1.5, cex.axis=1.5, font.axis=2,font.lab=2, 
                 xlab="log2 FC RNA-seq", ylab="-log10 adjP",axes=FALSE, main="Cluster-6 genes\nDM vs WT", cex.main=1.5, font.main=2))
  box(lwd=3)
  axis(1,labels=TRUE,lwd=3,font=2,cex.axis=1.5)
  axis(2,labels=TRUE,lwd=3,font=2,cex.axis=1.5)
  with(res.up, points(log2FC, -log10(padj), pch=20, col="#db3236",cex=0.5))
  with(res.down, points(log2FC, -log10(padj), pch=20, col="#4885ed",cex=0.5))
  Utext <- nrow(res.up)
  Dtext <- nrow(res.down)
  legend("toprigh", x.intersp=-0.5, y.intersp=-0.5, inset=0.01,legend=Utext, text.col="#db3236", cex=1.5, text.font=2,box.lty=0)
  legend("topleft", x.intersp=-0.5, y.intersp=-0.5, inset=0.01,legend=Dtext, text.col="#4885ed", cex=1.5, text.font=2,box.lty=0)
  abline(v=c(log(1.5,2),-log(1.5,2)), h=-log(0.05,10), lty=2,lwd=2)
  dev.off()
  ###
  write.table(res.up$Gene_name, paste0("Cluster-",clusterID,"_DMvsWT_upgenes.txt"),quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}
getList(6)