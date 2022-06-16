
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


#############################################
### Validating TF-gene: knockTFs
#############################################

# knockTF
# The names of the list elements represent target genes. In the database, there 
# are only 12 datasets for 12 TFs (BCL11A, ELK1, GATA3, JUN, MAF, MYB, NFATC3, 
# NFKB1, STAT3, STAT6, TAL1, and ZNF148) in the “peripheral blood” category and 
# One dataset in the “blood” category (HIF1A KD in primary human macrophage). 
# Other categories do not appear blood-related. 

knockTF.up=readRDS('data/knockTF.upregulated.0.05.list.rds')
knockTF.down=readRDS('data/knockTF.downregulated.0.05.list.rds')

knTF.up=vector(mode='list', length=11)
knTFs=c('BCL11A', 'ELK1', 'GATA3', 'JUN', 'MAF',
        'MYB', 'NFATC3', 'NFKB1', 'STAT3', 
        'TAL1', 'ZNF148')
names(knTF.up)=knTFs
# Removed STAT6, which does not have any up- or down-target gene.

rej1=rej4=rej5Xor5Y=rej5X=rej5Xand5Y=knTF.down=knTF.up

for(i in 1:length(knockTF.up)){
  temp=knockTF.up[[i]]
  for(t in 1:length(temp)){
    knTF.up[[temp[t]]]=c(knTF.up[[temp[t]]], names(knockTF.up)[i])
  }
}

for(i in 1:length(knockTF.down)){
  temp=knockTF.down[[i]]
  for(t in 1:length(temp)){
    knTF.down[[temp[t]]]=c(knTF.down[[temp[t]]], names(knockTF.down)[i])
  }
}

rm(knockTF.down); rm(knockTF.up)

length(knTF.up); length(knTF.down)
# 11 experiments for 11 TFs (STAT1 is removed because it does not have any up- or down-target genes)

for(i in 1:length(xymats.all)){
  if(i %%20 ==0) cat(i,'\t')
  gene.name=names(xymats.all)[i]
  
  # gene range
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  peak.gr.g=xymats.all[[gene.name]]$xymats4$peak.gr.g
  TF.g=xymats.all[[gene.name]]$xymats4$TF.g
  
  beta.rej1=benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs, fdr.thresh = 0.05)$which.reject &
                    xymats.all[[gene.name]]$xymats1$betaYijs>0
  gammas.rej4=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, do.plot = FALSE, fdr.thresh = 0.05) |
    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, do.plot = FALSE, fdr.thresh = 0.05)) &
    xymats.all[[gene.name]]$xymats1$betaYijs>0
    
  gammas.rej5X=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05)) &
    xymats.all[[gene.name]]$xymats5X$betahat.corrests>0
  
  gammas.rej5Y=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05) ) &
    xymats.all[[gene.name]]$xymats5Y$betahat.corrests>0
  
  TF.rej1=intersect(TF.g[sort(unique(which(beta.rej1)))], knTFs)
  TF.rej4=intersect(TF.g[sort(unique(which(gammas.rej4, arr.ind = T)[,2]))], knTFs)
  TF.rej5Xor5Y=intersect(TF.g[sort(unique(which(gammas.rej5X | gammas.rej5Y, arr.ind = T)[,2]))], knTFs)
  TF.rej5X=intersect(TF.g[sort(unique(which(gammas.rej5X, arr.ind = T)[,2]))], knTFs)
  TF.rej5Xand5Y=intersect(TF.g[sort(unique(which(gammas.rej5X & gammas.rej5Y, arr.ind = T)[,2]))], knTFs)
  
  
  if(length(TF.rej1)>0){
    for(t in 1:length(TF.rej1)){rej1[[TF.rej1[t]]]=c(rej1[[TF.rej1[t]]],gene.name)}
  }
  if(length(TF.rej4)>0){
    for(t in 1:length(TF.rej4)){rej4[[TF.rej4[t]]]=c(rej4[[TF.rej4[t]]],gene.name)}
  }
  if(length(TF.rej5Xor5Y)>0){
    for(t in 1:length(TF.rej5Xor5Y)){rej5Xor5Y[[TF.rej5Xor5Y[t]]]=c(rej5Xor5Y[[TF.rej5Xor5Y[t]]],gene.name)}
  }
  if(length(TF.rej5X)>0){
    for(t in 1:length(TF.rej5X)){rej5X[[TF.rej5X[t]]]=c(rej5X[[TF.rej5X[t]]],gene.name)}
  }
  if(length(TF.rej5Xand5Y)>0){
    for(t in 1:length(TF.rej5Xand5Y)){rej5Xand5Y[[TF.rej5Xand5Y[t]]]=c(rej5Xand5Y[[TF.rej5Xand5Y[t]]],gene.name)}
  }
  
}

precision=matrix(ncol=5, nrow=length(knTFs))
rownames(precision)=knTFs
colnames(precision)=c('rej1','rej4','rej5Xor5Y','rej5X','rej5Xand5Y')
recall=precision

for(i in 1:length(knTFs)){
  TF.name=knTFs[i]
  
  # Precision
  precision[i,1]=length(intersect(unique(knTF.up[[TF.name]]), rej1[[TF.name]]))/length( rej1[[TF.name]])
  precision[i,2]=length(intersect(unique(knTF.up[[TF.name]]), rej4[[TF.name]]))/length( rej4[[TF.name]])
  precision[i,3]=length(intersect(unique(knTF.up[[TF.name]]), rej5Xor5Y[[TF.name]]))/length( rej5Xor5Y[[TF.name]])
  precision[i,4]=length(intersect(unique(knTF.up[[TF.name]]), rej5X[[TF.name]]))/length( rej5X[[TF.name]])
  precision[i,5]=length(intersect(unique(knTF.up[[TF.name]]), rej5Xand5Y[[TF.name]]))/length( rej5Xand5Y[[TF.name]])
  
  # Recall
  recall[i,1]=length(intersect(unique(knTF.up[[TF.name]]), rej1[[TF.name]]))/length(unique(knTF.up[[TF.name]]))
  recall[i,2]=length(intersect(unique(knTF.up[[TF.name]]), rej4[[TF.name]]))/length(unique(knTF.up[[TF.name]]))
  recall[i,3]=length(intersect(unique(knTF.up[[TF.name]]), rej5Xor5Y[[TF.name]]))/length(unique(knTF.up[[TF.name]]))
  recall[i,4]=length(intersect(unique(knTF.up[[TF.name]]), rej5X[[TF.name]]))/length(unique(knTF.up[[TF.name]]))
  recall[i,5]=length(intersect(unique(knTF.up[[TF.name]]), rej5Xand5Y[[TF.name]]))/length(unique(knTF.up[[TF.name]]))
}

TF.filter=apply(precision, 1, function(x){any(is.nan(x)) | all(x==0)}) |
  apply(recall, 1, function(x){any(is.nan(x)) | all(x==0)})

precision=precision[!TF.filter,]
recall=recall[!TF.filter,]

precision=precision[,c(1,3)]
recall=recall[,c(1,3)]

precision
recall

library(reshape)
toplot=melt(precision)
colnames(toplot)=c('TF','Model','Precision')
toplot$Model=factor(toplot$Model, levels=c('rej1', 'rej5Xor5Y'))
p1 <- ggplot(data=toplot, aes(x=TF, y=Precision, fill=Model)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ylab('Precision = TP/(TP+FP)')+xlab('Transcription factor')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

toplot=melt(recall)
colnames(toplot)=c('TF','Model','Recall')
toplot$Model=factor(toplot$Model, levels=c('rej1', 'rej5Xor5Y'))
p2 <- ggplot(data=toplot, aes(x=TF, y=Recall, fill=Model)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ylab('Recall = TP/(TP+FN)')+xlab('Transcription factor') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2|p1

ggsave(filename = 'knockTF.precision_recall.pos.pdf', plot = p2|p1, width = 8, height = 3)

