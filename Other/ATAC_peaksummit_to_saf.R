########## convert peak summit to 2kb bins in .saf ##########
library(GenomicRanges)
library(GenomicAlignments)
atac.summit <- read.table("ATAC_consensus_peak_summit.bed", sep="\t",header=F)
colnames(atac.summit) <- c("chr", "start", "end", "id")

atac_peak = GRanges(seqnames=atac.summit$chr, IRanges(atac.summit$start,atac.summit$end),id=atac.summit$id,strand = "+")
atac_peak_2kb = resize(atac_peak, width = 2000, fix = "center")
atac_peak_adj = reduce(sort(atac_peak_2kb))
atac_peak_adj$id = c(paste0(seqnames(atac_peak_adj),"_",start(atac_peak_adj),"_",end(atac_peak_adj)))
atac_out=data.frame(atac_peak_adj$id,seqnames(atac_peak_adj),start(atac_peak_adj),end(atac_peak_adj),
                    strand(atac_peak_adj),stringsAsFactors = FALSE)
colnames(atac_out)=c("GeneID","Chr","Start","End","Strand")
write.table(atac_out, "ATAC_consensus_summit2kb_adj.saf", append = FALSE, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
