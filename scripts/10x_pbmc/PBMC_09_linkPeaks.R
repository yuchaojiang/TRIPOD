
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(presto)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)
library(gplots)
library(nbpMatching)
library(fields)
library(gridBase)
library(grid)

load('ATAC_RNA_WNN_integrated_high_var.rda') # Testing results
load('ATAC_RNA_WNN_footprint.rda') # Footprinting

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

### process the data
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RegionStats(object = pbmc, genome =BSgenome.Hsapiens.UCSC.hg38)
pbmc <- LinkPeaks(object = pbmc, peak.assay = "ATAC", expression.assay = "SCT",
                 pvalue_cutoff = 1, score_cutoff = 0,
                 distance = 2e+05, method = "pearson", genes.use = genes)
save(pbmc, file='linkPeak/pbmc.linkpeak.rda')

pearson.link <- Links(pbmc)
hist(pearson.link$pvalue)

pearson.link$pvalue.2 <- 2*pnorm(abs(pearson.link$zscore), lower.tail = FALSE) # two-sided
fdr.thresh <- 0.01
pearson.bh <- benjaminiHochsbergVectorAdjust(pearson.link$pvalue.2, fdr.thresh = fdr.thresh)
table(pearson.bh$which.reject)

identical(pearson.bh$pval.c.adj < fdr.thresh, pearson.bh$which.reject)
pearson.link$adj <- pearson.bh$pval.c.adj
hist(pearson.link$adj)

pearson.link.pos <- pearson.link[pearson.link$score > 0 & pearson.bh$which.reject]
pearson.link.neg <- pearson.link[pearson.link$score < 0 & pearson.bh$which.reject]

## output csv files
# gene, peak, coef, pval, adj
pos.df <- as.data.frame(mcols(pearson.link.pos))
pos.df <- pos.df[, c(2, 3, 1, 6, 7)]
colnames(pos.df) <- c("gene", "peak", "coef", "pval", "adj")

neg.df <- as.data.frame(mcols(pearson.link.neg))
neg.df <- neg.df[, c(2, 3, 1, 6, 7)]
colnames(neg.df) <- c("gene", "peak", "coef", "pval", "adj")

# get sets of peak-gene pair strings
# compute number of unique peak-gene pairs
convertPairToString <- function(x, col.1, col.2){
  gsub(" ", "", paste(x[col.1], x[col.2], sep = "_"))
}

linkpeaks.pos.strings <- unique(apply(pos.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
linkpeaks.neg.strings <- unique(apply(neg.df, 1, convertPairToString, col.1 = 1, col.2 = 2))
c(length(linkpeaks.pos.strings), length(linkpeaks.neg.strings))

write.table(pos.df, 'linkPeak/pos.df.csv', quote = F, sep = ",")
write.table(pos.df, 'linkPeak/neg.df.csv', quote = F, sep = ",")
save(linkpeaks.pos.strings, file='linkPeak/linkpeaks.pos.strings.rda')
save(linkpeaks.neg.strings, file='linkPeak/linkpeaks.neg.strings.rda')
