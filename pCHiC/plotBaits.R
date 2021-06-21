
library(devtools)
library(Chicago)
library(PCHiCdata)
########## plotBaits ##########
library(Chicago)
setwd("~/analyses/Chicago/enhancer_project")
WT.run <- readRDS("WT.run.rds")
NPM1.run <- readRDS("NPM1.run.rds")
FLT3.run <- readRDS("FLT3.run.rds")
DM.run <- readRDS("DM.run.rds")

## Irf8 729806
plotBaits(WT.run, n=1, Ncol = "N", baits = 729806,ylim=c(0,150),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Irf8_4C_norm_WT.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.1e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(NPM1.run, n=1, Ncol = "N", baits = 729806,ylim=c(0,170),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Irf8_4C_norm_NPM1.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.1e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(FLT3.run, n=1, Ncol = "N", baits = 729806,ylim=c(0,155),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Irf8_4C_norm_FLT3.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.1e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(DM.run, n=1, Ncol = "N", baits = 729806,ylim=c(0,120),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Irf8_4C_norm_DM.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.1e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

### Spi1 425092
plotBaits(WT.run, n=1, Ncol = "N", baits = 425092,ylim=c(0,160),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Spi1_4C_norm_WT.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.03e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(NPM1.run, n=1, Ncol = "N", baits = 425092,ylim=c(0,200),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Spi1_4C_norm_NPM1.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.03e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(FLT3.run, n=1, Ncol = "N", baits = 425092,ylim=c(0,200),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Spi1_4C_norm_FLT3.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.03e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(DM.run, n=1, Ncol = "N", baits = 425092,ylim=c(0,200),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Spi1_4C_norm_DM.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.03e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

### Hoxa9 617336
plotBaits(WT.run, n=1, Ncol = "N", baits = 617336, ylim=c(0,100), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa9_4C_norm_WT.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(NPM1.run, n=1, Ncol = "N", baits = 617336, ylim=c(0,100), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa9_4C_norm_NPM1.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(FLT3.run, n=1, Ncol = "N", baits = 617336, ylim=c(0,100),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa9_4C_norm_FLT3.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(DM.run, n=1, Ncol = "N", baits = 617336, ylim=c(0,100), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa9_4C_norm_DM.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

### Hoxa10 617342
plotBaits(WT.run, n=1, Ncol = "N", baits = 617342, ylim=c(0,100),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa10_4C_norm_WT.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(NPM1.run, n=1, Ncol = "N", baits = 617342, ylim=c(0,100),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa10_4C_norm_NPM1.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(FLT3.run, n=1, Ncol = "N", baits = 617342, ylim=c(0,100),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa10_4C_norm_FLT3.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(DM.run, n=1, Ncol = "N", baits = 617342, ylim=c(0,100),plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Hoxa10_4C_norm_DM.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 1.2e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

### Gfi1b 404376
plotBaits(WT.run, n=1, Ncol = "N", baits = 404376, ylim=c(0,115), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Gfi1b_4C_norm_WT.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.04e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(NPM1.run, n=1, Ncol = "N", baits = 404376, ylim=c(0,160), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Gfi1b_4C_norm_NPM1.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.04e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(FLT3.run, n=1, Ncol = "N", baits = 404376, ylim=c(0,160), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Gfi1b_4C_norm_FLT3.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.04e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

plotBaits(DM.run, n=1, Ncol = "N", baits = 404376, ylim=c(0,115), plotBaitNames = TRUE, plotBprof = TRUE, outfile = "Gfi1b_4C_norm_DM.pdf", width = 7, height = 3.5, removeBait2bait = FALSE, maxD = 0.04e6, bgCol = "black", plevel1 = 5, plevel2 = 3, lev1Col = "red", lev2Col = "blue",bgPch = 20, lev1Pch = 18, lev2Pch = 18)

