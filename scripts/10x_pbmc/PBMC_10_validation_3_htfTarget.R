
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

load('ATAC_RNA_WNN_integrated_high_var.rda') # Testing results
load('ATAC_RNA_WNN_footprint.rda') # Footprinting

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

# Highly variable genes
hTFtarget <- readRDS("data/hTFtarget.highly_variable_genes.rds")
hTFtarget=hTFtarget[names(hTFtarget) %in% names(xymats.all)]

for(i in 1:length(hTFtarget)){
  hTFtarget[[i]]=hTFtarget[[i]][hTFtarget[[i]]$tissue=='blood',]
}

get.precision.recall=function(goldCall, actualCall, total){
  TP=intersect(actualCall, goldCall)
  FP=setdiff(actualCall, TP)
  TN=intersect(setdiff(total, goldCall), setdiff(total, actualCall))
  FN=setdiff(setdiff(total, actualCall), TN)
  
  TP=length(TP)
  FP=length(FP)
  TN=length(TN)
  FN=length(FN)
  
  precision=TP/(TP+FP)
  recall=TP/(TP+FN)
  return(c(precision, recall))
}

gene.name='GNLY'

marginal=matrix(ncol=2, nrow=length(hTFtarget))
colnames(marginal)=c('precision','recall')
rownames(marginal)=names(hTFtarget)
TRIPOD=marginal

for(i in 1:length(hTFtarget)){
  gene.name =names(hTFtarget)[i]
  
  # gene range
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  peak.gr.g=xymats.all[[gene.name]]$xymats4$peak.gr.g
  TF.g=xymats.all[[gene.name]]$xymats4$TF.g
  hTFtarget.g=hTFtarget[[gene.name]]$TF
  
  beta.rej1=benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs, fdr.thresh = 0.05)$which.reject 
  
  gammas.rej1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats1$pvalgammas, do.plot = FALSE, fdr.thresh = 0.01)
  
  gammas.rej4=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, do.plot = FALSE, fdr.thresh = 0.05) |
                 benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5X=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5Y=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05) ) 
  
  TF.rej1=TF.g[sort(unique(which(beta.rej1)))]
  TF.rej5Xor5Y=TF.g[sort(unique(which(gammas.rej5X | gammas.rej5Y, arr.ind = T)[,2]))]

  marginal[i,]=get.precision.recall(goldCall = hTFtarget.g, actualCall = TF.rej1, total=TF.g)
  TRIPOD[i,]=get.precision.recall(goldCall = hTFtarget.g, actualCall = TF.rej5Xor5Y, total=TF.g)
}

pdf(file='hTFtarget.pdf', width=4, height=4)
plot(marginal[,1], marginal[,2], col='#2AB34B', pch=16, xlim=c(0,1), ylim=c(0,1),xlab='precision',ylab='recall', cex=0.8)
points(TRIPOD[,1], TRIPOD[,2], col='#7094CD', pch=16, cex=0.8)
dev.off()

