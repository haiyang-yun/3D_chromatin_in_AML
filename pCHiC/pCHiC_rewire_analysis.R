setwd("~/Documents/Chicago_new/Differential.analysis/rewire_analysis")

All.score <- read.delim("R1-2_chicago.scores.matrix.txt", sep="\t",header=TRUE,row.names = 1)

### high means gained interations, low means lost interactions

############### WT vs NPM1
WT.pos.NPM1.neg <- read.delim("WT.pos.NPM1.neg.txt", sep="\t",header=FALSE)
WT.neg.NPM1.pos <- read.delim("WT.neg.NPM1.pos.txt", sep="\t",header=FALSE)
WT.NPM1.all <- rbind(WT.pos.NPM1.neg,WT.neg.NPM1.pos)
WT.NPM1.all.score <- All.score[row.names(All.score) %in% WT.NPM1.all[,1],1:2]

a = as.integer(nrow(WT.NPM1.all.score)*0.25) 
b = as.integer(nrow(WT.NPM1.all.score)*0.5) 
c = as.integer(nrow(WT.NPM1.all.score)*0.75)

WT.NPM1.tmp1 <- WT.NPM1.all.score[order(-WT.NPM1.all.score[,1]),] 
WT.NPM1.tmp1$WT.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.NPM1.tmp1)-c))

WT.NPM1.tmp2 <- WT.NPM1.tmp1[order(-WT.NPM1.tmp1[,2]),] 
WT.NPM1.tmp2$NPM1.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.NPM1.tmp2)-c))

WT.high.NPM1.low <- subset(WT.NPM1.tmp2,WT.group == "A" & NPM1.group == "D")
WT.low.NPM1.high <- subset(WT.NPM1.tmp2,WT.group == "D" & NPM1.group == "A")

write.table(WT.high.NPM1.low,file="WT.high.NPM1.low.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)
write.table(WT.low.NPM1.high,file="WT.low.NPM1.high.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)

############### WT vs FLT3
WT.pos.FLT3.neg <- read.delim("WT.pos.FLT3.neg.txt", sep="\t",header=FALSE)
WT.neg.FLT3.pos <- read.delim("WT.neg.FLT3.pos.txt", sep="\t",header=FALSE)
WT.FLT3.all <- rbind(WT.pos.FLT3.neg,WT.neg.FLT3.pos)
WT.FLT3.all.score <- All.score[row.names(All.score) %in% WT.FLT3.all[,1],c(1,3)] ## select specific columns

a = as.integer(nrow(WT.FLT3.all.score)*0.25) 
b = as.integer(nrow(WT.FLT3.all.score)*0.5) 
c = as.integer(nrow(WT.FLT3.all.score)*0.75)

WT.FLT3.tmp1 <- WT.FLT3.all.score[order(-WT.FLT3.all.score[,1]),] 
WT.FLT3.tmp1$WT.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.FLT3.tmp1)-c))

WT.FLT3.tmp2 <- WT.FLT3.tmp1[order(-WT.FLT3.tmp1[,2]),] 
WT.FLT3.tmp2$FLT3.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.FLT3.tmp2)-c))

WT.high.FLT3.low <- subset(WT.FLT3.tmp2,WT.group == "A" & FLT3.group == "D")
WT.low.FLT3.high <- subset(WT.FLT3.tmp2,WT.group == "D" & FLT3.group == "A")

write.table(WT.high.FLT3.low,file="WT.high.FLT3.low.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)
write.table(WT.low.FLT3.high,file="WT.low.FLT3.high.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)

############### WT vs DM
WT.pos.DM.neg <- read.delim("WT.pos.DM.neg.txt", sep="\t",header=FALSE)
WT.neg.DM.pos <- read.delim("WT.neg.DM.pos.txt", sep="\t",header=FALSE)
WT.DM.all <- rbind(WT.pos.DM.neg,WT.neg.DM.pos)
WT.DM.all.score <- All.score[row.names(All.score) %in% WT.DM.all[,1],c(1,4)] 

a = as.integer(nrow(WT.DM.all.score)*0.25) 
b = as.integer(nrow(WT.DM.all.score)*0.5) 
c = as.integer(nrow(WT.DM.all.score)*0.75)

WT.DM.tmp1 <- WT.DM.all.score[order(-WT.DM.all.score[,1]),] 
WT.DM.tmp1$WT.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.DM.tmp1)-c))

WT.DM.tmp2 <- WT.DM.tmp1[order(-WT.DM.tmp1[,2]),] 
WT.DM.tmp2$DM.group <- c(rep("A",a),rep("B",b-a),rep("C",c-b),rep("D",nrow(WT.DM.tmp2)-c))

WT.high.DM.low <- subset(WT.DM.tmp2,WT.group == "A" & DM.group == "D")
WT.low.DM.high <- subset(WT.DM.tmp2,WT.group == "D" & DM.group == "A")

write.table(WT.high.DM.low,file="WT.high.DM.low.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)
write.table(WT.low.DM.high,file="WT.low.DM.high.txt",row.names=TRUE,quote = FALSE, sep = "\t",col.names = TRUE)
