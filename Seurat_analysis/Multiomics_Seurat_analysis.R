library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(dplyr) 

##################### Seurat analysis on merged cpm from replicates #####################
atac.m.matrix <- read.delim("ATAC_consensus_summit2kb_adj_cpm_merge_transpose.txt", row.names = 1, head=TRUE)
rownames(atac.m.matrix) = c("WT.H3K4me1","Npm1c.H3K4me1","Flt3-ITD.H3K4me1","DM.H3K4me1",
"WT.H3K4me3","Npm1c.H3K4me3","Flt3-ITD.H3K4me3","DM.H3K4me3",
"WT.H3K27ac","Npm1c.H3K27ac","Flt3-ITD.H3K27ac","DM.H3K27ac",
"WT.ATAC-seq","Npm1c.ATAC-seq","Flt3-ITD.ATAC-seq","DM.ATAC-seq")
atac.m.object <- CreateSeuratObject(counts = atac.m.matrix, project = "atac.m.summit", min.cells = 0, min.features = 0)
atac.m.object <- NormalizeData(atac.m.object, normalization.method = "RC", scale.factor = 1000000)

atac.m.object <- FindVariableFeatures(atac.m.object, selection.method = "vst")
m.id <- rownames(atac.m.object)
atac.m.matrix.scale <- ScaleData(atac.m.object, features = m.id)
atac.m.matrix.pca <- RunPCA(atac.m.matrix.scale,features = m.id,npcs = 50)

png("ATAC_consensus_summit2kb_adj_cpm_merg_pca_cutoff.png", width = 6, height = 4, units = "in", res=300)
ElbowPlot(atac.m.matrix.pca, ndims = 50)
dev.off()
## cutoff was set at 8 based on the graph

m.tmp <- FindNeighbors(atac.m.matrix.pca, reduction = "pca", dims = 1:8, nn.eps = 0.5)
m.tmp <- FindClusters(m.tmp, resolution = 0.5)
## Number of communities: 10

m.tmp3 <- RunTSNE(m.tmp, dims = 1:10)

##################### change the order of the clusters #####################
# first: to save old identity classes (the cluster labels) for reference.
m.tmp3[["old.ident"]] <- Idents(object = m.tmp3)

# Rename classes
m.tmp3 <- RenameIdents(object = m.tmp3, "0"="1","3"="2","7"="3","4"="4","9"="5","6"="6","5"="7","2"="8","1"="9","8"="10")
p1 <- DimPlot(m.tmp3, reduction = "tsne", pt.size = 0.1, label=TRUE, raster=FALSE,
order =  c("10","9","8","7","6","5","4","3","2","1")) +
ggtitle(label="Accessible chromatin regions") +
theme_bw() + theme_classic(base_size = 12, base_family = "") +
theme(axis.text = element_text(color="black", size=12, face="bold"),
axis.ticks = element_line(size = 1), axis.ticks.length = unit(.25, "cm"),
axis.title = element_text(color="black", size=14, face="bold"),
plot.title = element_text(color="black", size=14, face="bold",,hjust=0.5),
legend.title = element_text(color="black", size=14, face="bold",hjust=0.5),
legend.text = element_text(size=12, face="bold"))
ggsave("ATAC_consensus_summit2kb_adj_cpm_merg_tsne.2021.pdf", width = 5, height = 4.5)

pdf("ATAC_consensus_summit2kb_adj_cpm_merg_heatmap.pdf", width = 9.5, height = 5)
DoHeatmap(m.tmp3, features = m.id)
dev.off()

