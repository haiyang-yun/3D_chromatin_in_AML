library(ggplot2)
library(ggrepel)

##### PCA for RNA-seq, use normalised reads, filtered PC genes
data <- read.table("Normalized.filteredCounts.PC.csv", header = TRUE,sep=",")
df=data[,c(6,7,8,9,4,5,2,3)]
alldata <- df[ !rowSums(df[,colnames(df)[(1:ncol(df))]]==0)==ncol(df), ]
alldatapca <- prcomp(t(alldata), scale = TRUE)
alldatascores = as.data.frame(alldatapca$x) ### PC values
group = rownames(alldatascores)
summary(alldatapca) ### components frequencies
capture.output(summary(alldatapca), file ="RNA-seq_PCA_summary.txt", append = FALSE)
# Data design file
alldataDesign = data.frame(row.names = colnames(alldatascores),
condition = c( "WT", "WT", "NPM1", "NPM1", "FLT3", "FLT3", "DM", "DM"))
alldataDesign$condition = factor(alldataDesign$condition, levels=c("WT","NPM1","FLT3","DM"))
Sample = alldataDesign$condition
# PCA plot
bg <- ggplot(data = alldatascores, aes(x = PC1, y = PC2, label=group)) +
ggtitle("RNA-seq") + xlab(paste("PC1","(56.7%)")) + ylab(paste("PC2","(15.6%)")) +
geom_point(aes(color = Sample), size = 5) +
scale_color_manual(values=c("#AAA9AD","#0078D7","#3CAEA3","#ED553B")) +
theme_bw() + theme_classic(base_size = 12, base_family = "") +
theme(axis.text = element_text(color="black", size=12, face="bold"),
axis.ticks = element_line(size = 1), axis.ticks.length = unit(.25, "cm"),
axis.title = element_text(color="black", size=14, face="bold"),
plot.title = element_text(color="black", size=14, face="bold",,hjust=0.5),
legend.title = element_text(color="black", size=14, face="bold",hjust=0.5),
legend.text = element_text(size=12, face="bold"),
panel.border = element_rect(linetype = "solid",size = 1,fill = NA))
ggsave("RNA-seq_RPKM.PCA.pdf", width = 4.5, height = 3.5)


### PCA for ATAC-seq data
dataFile=read.delim("ATAC_consensus.cpm.csv", header = TRUE, sep =",",row.names=1)
data=dataFile[,c(1,2,3,4,5,6,7,8)]
secondlowest=min( data[data!=min(data)] ) ### find the 2nd lowest value
data[data==0]=secondlowest
colnames(data)=c("WT.R1","WT.R2","NPM1.R1","NPM1.R2","FLT3.R1","FLT3.R2","DM.R1","DM.R2")
alldata=log(data,2)
alldatapca <- prcomp(t(alldata), scale = TRUE)
alldatascores = as.data.frame(alldatapca$x) ### PC values
group = rownames(alldatascores)
summary(alldatapca) ### components frequencies
capture.output(summary(alldatapca), file ="ATAC_PCA_summary.txt", append = FALSE)
# Data design file
alldataDesign = data.frame(row.names = colnames(alldatascores),
condition = c( "WT", "WT", "NPM1", "NPM1", "FLT3", "FLT3", "DM", "DM"))
alldataDesign$condition = factor(alldataDesign$condition, levels=c("WT","NPM1","FLT3","DM"))
Sample = alldataDesign$condition
# PCA plot
bg <- ggplot(data = alldatascores, aes(x = PC1, y = PC2, label=group)) +
ggtitle("ATAC-seq") + xlab(paste("PC1","(54.6%)")) + ylab(paste("PC2","(29.7%)")) +
geom_point(aes(color = Sample), size = 5) +
scale_color_manual(values=c("#AAA9AD","#0078D7","#3CAEA3","#ED553B")) +
theme_bw() + theme_classic(base_size = 12, base_family = "") +
theme(axis.text = element_text(color="black", size=12, face="bold"),
axis.ticks = element_line(size = 1), axis.ticks.length = unit(.25, "cm"),
axis.title = element_text(color="black", size=14, face="bold"),
plot.title = element_text(color="black", size=14, face="bold",,hjust=0.5),
legend.title = element_text(color="black", size=14, face="bold",hjust=0.5),
legend.text = element_text(size=12, face="bold"),
panel.border = element_rect(linetype = "solid",size = 1,fill = NA))
ggsave("ATAC_PCA.pdf", width = 4.5, height = 3.5)


### pCHiC PCA
df <- read.table("PCHiC_scores_rep.overlap.all.noNE.txt", header = TRUE, row.names = 1)
alldata <- df[ !rowSums(df[,colnames(df)[(1:ncol(df))]]==0)==ncol(df), ]
alldatapca <- prcomp(t(alldata), scale = TRUE)
alldatascores = as.data.frame(alldatapca$x) ### PC values
group = rownames(alldatascores)
summary(alldatapca) ### components frequencies
capture.output(summary(alldatapca), file ="pCHiC_PCA_summary.txt", append = FALSE)
# Data design file
alldataDesign = data.frame(row.names = colnames(alldatascores),
condition = c( "WT", "WT", "NPM1", "NPM1", "FLT3", "FLT3", "DM", "DM"))
alldataDesign$condition = factor(alldataDesign$condition, levels=c("WT","NPM1","FLT3","DM"))
Sample = alldataDesign$condition
bg <- ggplot(data = alldatascores, aes(x = PC1, y = PC2, label=group)) +
ggtitle("pCHiC HC interaction scores") + xlab(paste("PC1","(25.4%)")) + ylab(paste("PC2","(19.1%)")) +
geom_point(aes(color = Sample), size = 5) +
scale_color_manual(values=c("#AAA9AD","#0078D7","#3CAEA3","#ED553B")) +
theme_bw() + theme_classic(base_size = 12, base_family = "") +
theme(axis.text = element_text(color="black", size=12, face="bold"),
axis.ticks = element_line(size = 1), axis.ticks.length = unit(.25, "cm"),
axis.title = element_text(color="black", size=14, face="bold"),
plot.title = element_text(color="black", size=14, face="bold",,hjust=0.5),
legend.title = element_text(color="black", size=14, face="bold",hjust=0.5),
legend.text = element_text(size=12, face="bold"),
panel.border = element_rect(linetype = "solid",size = 1,fill = NA))
ggsave("pCHiC_CHiCAGOscore_PCA.pdf", width = 4.5, height = 3.5)

