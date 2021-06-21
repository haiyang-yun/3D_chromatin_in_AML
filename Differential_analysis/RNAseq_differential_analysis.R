runDESeq<-function(s1,s2){
library(DESeq2)
############################################################
# List the files with read counts generated from HTSeq 
# This will include files for both conditions
# Each file contains read counts for all the transcripts
############################################################
sampleFiles =c(paste0(s1,".EXPR.R1.htseq.counts"),paste0(s1,".EXPR.R2.htseq.counts"),paste0(s2,".EXPR.R1.htseq.counts"),paste0(s2,".EXPR.R2.htseq.counts"))
print(sampleFiles)

###########################################################
# Get sample names from file Names
############################################################
sampleNames <- gsub(".htseq.counts","",sampleFiles)
sampleNames <- gsub(".EXPR","",sampleNames)
print(sampleNames)
###########################################################
# Define sample conditon for each sample in the order
# as the sample names are listed
############################################################
sampleCondition=c( s1,s1,s2,s2)
print(sampleCondition)
###########################################################
# Store all the information in a table
############################################################
sampleTable <- data.frame(sampleName = sampleNames,
                          fileName = sampleFiles,
                          condition = sampleCondition)
print(sampleTable)
#########################################################
# Read the counts for each transcript for all samples
#########################################################
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "./HTSeq-Counts",
                                       design= ~ condition)
print(ddsHTSeq)
########################################################
# Use only those transcripts which have a rowcount >1 
########################################################
geneFile=read.delim("gencode.vM7.proteinCoding_genes.list", header=F,row.name=1)
geneList=rownames(geneFile)
ddsHTSeq <- ddsHTSeq[rownames(ddsHTSeq)%in%geneList, ]
print(ddsHTSeq)
genesNameFile=read.csv("./gencode.vM7.geneNames.MOD.txt",sep="\t",header=F, row.name =1)
geneNames=genesNameFile[1]
rownames(ddsHTSeq)=geneNames[rownames(ddsHTSeq),]
########################################################
# Use only those transcripts which have a rowcount >1 
########################################################
dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]
print(dds)
dds$condition <- relevel(dds$condition, ref=s1)
print(dds$condition)
#######################################################
# Differential Expression analysis
#######################################################
ana<- DESeq(dds)
print(ana)
#######################################################
# Retrieve results of Differential Expression analysis
#######################################################
res <- results(ana)

#######################################################
# Summary of results of Differential Expression analysis
#######################################################
summary(res)

##########################################################
# Write output to a file
##########################################################
resOrdered <- res[order(res$log2FoldChange),]

write.table(as.data.frame(resOrdered),file=paste0(s1,".",s2,".PC.diffExp.csv"), sep="\t", quote=F, row.names=TRUE)

do.call(runDESeq, list(s1="WT",s2="DM"))
do.call(runDESeq, list(s1="WT",s2="NPM1"))
do.call(runDESeq, list(s1="WT",s2="FLT3"))


