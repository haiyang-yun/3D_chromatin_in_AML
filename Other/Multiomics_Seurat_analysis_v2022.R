library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(dplyr) 

##################### Seurat analysis on merged cpm from replicates #####################
atac.m.matrix <- read.delim("ATAC_consensus_summit2kb_adj_cpm_merge_transpose.txt", row.names = 1, head=TRUE)
rownames(atac.m.matrix) = c("WT-H3K4me1","Npm1c-H3K4me1","Flt3ITD-H3K4me1","DM-H3K4me1",
                            "WT-H3K4me3","Npm1c-H3K4me3","Flt3ITD-H3K4me3","DM-H3K4me3",
                            "WT-H3K27ac","Npm1c-H3K27ac","Flt3ITD-H3K27ac","DM-H3K27ac",
                            "WT-ATAC","Npm1c-ATAC","Flt3ITD-ATAC","DM-ATAC")
atac.m.object <- CreateSeuratObject(counts = atac.m.matrix, project = "atac.m.summit", min.cells = 0, min.features = 0)
atac.m.object <- NormalizeData(atac.m.object, normalization.method = "RC", scale.factor = 1000000)

atac.m.object <- FindVariableFeatures(atac.m.object, selection.method = "vst")
m.id <- rownames(atac.m.object)
atac.m.matrix.scale <- ScaleData(atac.m.object, features = m.id)
atac.m.matrix.pca <- RunPCA(atac.m.matrix.scale,features = m.id)

png("Multiomics_pca_cutoff.png", width = 4, height = 3.5, units = "in", res=300)
ElbowPlot(atac.m.matrix.pca, ndims = 15)
dev.off()
## suggesting a cutoff at 8 

m.tmp <- FindNeighbors(atac.m.matrix.pca, reduction = "pca", dims = 1:8, nn.eps = 0.5)
m.tmp <- FindClusters(m.tmp, resolution = 0.5)
## Number of communities: 10

m.tmp2 <- RunUMAP(m.tmp, dims = 1:10, min.dist = 0.75)
m.tmp3 <- RunTSNE(m.tmp, dims = 1:10)

p2 <- DimPlot(m.tmp2, reduction = "umap", pt.size = 0.1, label=TRUE, repel = TRUE, raster=FALSE) + ggtitle(label = "UMAP") + NoLegend()
p3 <- DimPlot(m.tmp3, reduction = "tsne", pt.size = 0.1, label=TRUE, repel = TRUE, raster=FALSE) + ggtitle(label = "tSNE") + NoLegend()

png("Multiomics_umap_tene.png", width = 7, height = 4, units = "in", res=300)
CombinePlots(plots = list(p2, p3))
dev.off()

## tsne clustering seperates better, choose tsne to make heatmap

png("Multiomics_heatmap.png", width = 9, height = 4, units = "in", res=300)
p4 <- DoHeatmap(m.tmp3, features = m.id, draw.lines = TRUE) + theme(legend.position="bottom") + 
  scale_color_manual(values=c('#F8766D','#D89000','#A3A500','#39B600','#00BF7D','#00BFC4',
                              '#00B0F6','#9590FF','#E76BF3','#FF62BC'))
ggsave("Multiomics_heatmap.png", width = 9, height = 4, units = "in", dpi=300) 

##################### Extract peak sets in Cluster-6 #####################
atac.summit <- read.table("ATAC_consensus_peak_summit.bed", sep="\t",header=F)
colnames(atac.summit) <- c("chr", "start", "end", "id")

atac_peak = GRanges(seqnames=atac.summit$chr, IRanges(atac.summit$start,atac.summit$end),id=atac.summit$id,strand = "+")

getList <- function(clusterID){
  tmp=subset(m.tmp3, idents = clusterID)
  tmp1=Cells(tmp)
  tmp3 = strsplit(tmp1, "_") 
  tmp4 = GRanges(seqnames=sapply(tmp3, "[", 1), IRanges(as.numeric(sapply(tmp3, "[", 2)),as.numeric(sapply(tmp3, "[", 3))))
  tmp5=atac_peak[atac_peak %over% tmp4,]
  tmp6=resize(tmp5, width = 200, fix = "center")
  tmp7=data.frame(seqnames(tmp6),start(tmp6),end(tmp6),stringsAsFactors = FALSE)
  write.table(tmp7, paste0("Multiomics_Cluster-",clusterID,"_summit200bp.bed"), append = FALSE, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
}
getList(6)
