### H3K4me1 reads at new enhancers (ov ATAC-seq peak (cutoff20) summit +/-250bp)
library(flextable)
### replicates merged
k4me1_all <- read.delim("H3K4me1_CPM_newenhancers.csv", sep=",", header=TRUE)
nrow(k4me1_all) #41122
WT=k4me1_all$LN.H3K4me1.R1+k4me1_all$LN.H3K4me1.R2
NPM1=k4me1_all$NPM1.H3K4me1.R1+k4me1_all$NPM1.H3K4me1.R2
FLT3=k4me1_all$FLT3.H3K4me1.R1+k4me1_all$FLT3.H3K4me1.R2
DM=k4me1_all$DM.H3K4me1.R1+k4me1_all$DM.H3K4me1.R2
peaksID=paste(k4me1_all$chr,"_",k4me1_all$start,"_",k4me1_all$end,sep="")
k4me1_new = data.frame(WT,NPM1,FLT3,DM,row.names = peaksID)
colnames(k4me1_new) = c("WT","Npm1c","Flt3-ITD","DM")
#correlation heatmap
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex=1.8, font=2, col="#db3236")
}
panel.sc <-function(x, y,  ...)
{
    usr <- par("usr")
    points(x, y, col="#db3236", cex=.3, pch=20)
}

png("H3K4me1_newEnh_re.scatterplot.png", width = 4.5, height = 4.5, units = 'in', res = 300)
par(bg=NA)
pairs(log2(k4me1_new), lower.panel = panel.cor,upper.panel= panel.sc,
font.labels = 2,cex.labels=1.5,font.main=2,cex.main=1.5,
main="")
dev.off()
