
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

# LinkPeaks
linkpeak.pos=read.csv('linkPeak/pos.df.csv')
linkpeak.neg=read.csv('linkPeak/pos.df.csv')
linkpeak=rbind(linkpeak.pos, linkpeak.neg)

#############################################
### Validating peak - target genes using enhancer database
#############################################

# Enhancer database
fantom=readRDS('data/fantom.grl.rds')
fourDG=readRDS('data/four.d.genome.grl.rds')
enhancerAtlas=readRDS('data/enhanceratlas.combined.grl.rds')

# Look at a specific gene
gene.name='GNLY'
gene.name='CCR7'

pval.enhancerAtlas=matrix(ncol=5, nrow=length(xymats.all), data=NA)
rownames(pval.enhancerAtlas)=names(xymats.all)
colnames(pval.enhancerAtlas)=c('Model1','Model4','Model5XorModel5Y','Model5XandModel5Y', 'linkPeaks')
pval.fantom=pval.fourDG=pval.enhancerAtlas
for(i in 1:length(xymats.all)){
  if(i %%20 ==0) cat(i,'\t')
  gene.name=names(xymats.all)[i]
  
  # gene range
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  peak.gr.g=xymats.all[[gene.name]]$xymats4$peak.gr.g
  
  gammas.rej1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats1$pvalgammas, do.plot = FALSE, fdr.thresh = 0.01)
  
  gammas.rej4=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, do.plot = FALSE, fdr.thresh = 0.05) |
                 benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5X=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5Y=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05) ) 
  
  fantom.gene=fantom[[gene.name]]
  fourDG.gene=fourDG[[gene.name]]
  enhancerAtlas.gene=enhancerAtlas[[gene.name]]
  
  # Model 1:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$which.reject)))]
  pval.enhancerAtlas[i,1]=phyper.test(peak.gr.g, sig.peak.gr.g, enhancerAtlas.gene)
  pval.fantom[i,1]=phyper.test(peak.gr.g, sig.peak.gr.g, fantom.gene)
  pval.fourDG[i,1]=phyper.test(peak.gr.g, sig.peak.gr.g, fourDG.gene)
  
  # Model 4:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej4, arr.ind = T)[,1]))]
  pval.enhancerAtlas[i,2]=phyper.test(peak.gr.g, sig.peak.gr.g, enhancerAtlas.gene)
  pval.fantom[i,2]=phyper.test(peak.gr.g, sig.peak.gr.g, fantom.gene)
  pval.fourDG[i,2]=phyper.test(peak.gr.g, sig.peak.gr.g, fourDG.gene)
  
  # Model 5X or 5Y:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej5X | gammas.rej5Y, arr.ind = T)[,1]))]
  pval.enhancerAtlas[i,3]=phyper.test(peak.gr.g, sig.peak.gr.g, enhancerAtlas.gene)
  pval.fantom[i,3]=phyper.test(peak.gr.g, sig.peak.gr.g, fantom.gene)
  pval.fourDG[i,3]=phyper.test(peak.gr.g, sig.peak.gr.g, fourDG.gene)
  
  # Model 5X and 5Y:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej5X & gammas.rej5Y, arr.ind = T)[,1]))]
  pval.enhancerAtlas[i,4]=phyper.test(peak.gr.g, sig.peak.gr.g, enhancerAtlas.gene)
  pval.fantom[i,4]=phyper.test(peak.gr.g, sig.peak.gr.g, fantom.gene)
  pval.fourDG[i,4]=phyper.test(peak.gr.g, sig.peak.gr.g, fourDG.gene)
  
  # link peaks:
  temp=linkpeak[which(linkpeak$gene==gene.name),2]
  if(length(temp)>0){
    sig.peak.gr.g=StringToGRanges(temp)
    pval.enhancerAtlas[i,5]=phyper.test(peak.gr.g, sig.peak.gr.g, enhancerAtlas.gene)
    pval.fantom[i,5]=phyper.test(peak.gr.g, sig.peak.gr.g, fantom.gene)
    pval.fourDG[i,5]=phyper.test(peak.gr.g, sig.peak.gr.g, fourDG.gene)
  }
}

pval.fantom=pval.fantom[,c(5,1,3)] # linkpeak, marginal, model 5X or 5Y
pval.fourDG=pval.fourDG[,c(5,1,3)]
pval.enhancerAtlas=pval.enhancerAtlas[,c(5,1,3)]

par(mfrow=c(1,3))
par(mar=c(12, 4.1, 4.1, 2.1))
boxplot(pval.fantom, main='FANTOM', las=2, ylab='Hypergeometric test p-val')
boxplot(pval.fourDG, main='fourDG',las=2, ylab='Hypergeometric test p-val')
boxplot(pval.enhancerAtlas, main='enhancerAtlas',las=2, ylab='Hypergeometric test p-val')


#############################################
### Validating peak - target genes using GTEx
#############################################

# GTEx
load('gtex.rda') # Too large to be uploaded to GitHub
# The data is downloaded from the GTEx portal (whole blood), followed by processing.

pval.gtex=matrix(ncol=5, nrow=length(xymats.all), data=NA)
rownames(pval.gtex)=names(xymats.all)
colnames(pval.gtex)=c('Model1','Model4','Model5XorModel5Y','Model5XandModel5Y', 'linkPeaks')

for(i in 1:length(xymats.all)){
  if(i %%20 ==0) cat(i,'\t')
  gene.name=names(xymats.all)[i]
  # gene range
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  peak.gr.g=xymats.all[[gene.name]]$xymats4$peak.gr.g

  gammas.rej1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats1$pvalgammas, do.plot = FALSE, fdr.thresh = 0.01)
  
  gammas.rej4=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, do.plot = FALSE, fdr.thresh = 0.05) |
                 benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5X=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05)) 
  
  gammas.rej5Y=(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot = FALSE, fdr.thresh = 0.05) |
                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot = FALSE, fdr.thresh = 0.05) ) 

  gtex.gene=gtex[[gene.name]]
  
  if(length(gtex.gene)==0) next
  gtex.gene=gtex.gene[countOverlaps(gtex.gene,peak.gr.g)>0] # select only the eQTLs that overlap with peaks
  gtex.gene=gtex.gene[countOverlaps(gtex.gene, transcripts.gr)==0] # remove SNPs within gene bodies
  
  if(length(gtex.gene)==0) next

  # Model 1:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$which.reject)))]
  pval.gtex[i,1]=phyper.test(peak.gr.g, sig.peak.gr.g, gtex.gene)
  
  # Model 4:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej4, arr.ind = T)[,1]))]
  pval.gtex[i,2]=phyper.test(peak.gr.g, sig.peak.gr.g, gtex.gene)
  
  # Model 5X or 5Y:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej5X | gammas.rej5Y, arr.ind = T)[,1]))]
  pval.gtex[i,3]=phyper.test(peak.gr.g, sig.peak.gr.g, gtex.gene)
  
  # Model 5X and 5Y:
  sig.peak.gr.g=peak.gr.g[sort(unique(which(gammas.rej5X & gammas.rej5Y, arr.ind = T)[,1]))]
  pval.gtex[i,4]=phyper.test(peak.gr.g, sig.peak.gr.g, gtex.gene)
  
  # link peaks:
  temp=linkpeak[which(linkpeak$gene==gene.name),2]
  if(length(temp)>0){
    sig.peak.gr.g=StringToGRanges(temp)
    pval.gtex[i,5]=phyper.test(peak.gr.g, sig.peak.gr.g, gtex.gene)
  }
}

pval.gtex=pval.gtex[,c(5,1,3)]

par(mfrow=c(1,1))
par(mar=c(12, 4.1, 4.1, 2.1))
boxplot(pval.gtex, main='GTEx', las=2, ylab='Hypergeometric test p-val')

library(reshape)
pval.fantom.toplot=melt(pval.fantom)
colnames(pval.fantom.toplot)=c('Gene','Model','pval')
pval.fantom.toplot$Model=factor(pval.fantom.toplot$Model,
                                levels=c('linkPeaks','Model1','Model5XorModel5Y','Model5XandModel5Y'))
p1=ggplot(pval.fantom.toplot, aes(x=Model, y=pval, fill=Model)) +
  geom_boxplot(outlier.size = 0.5) + 
  ggtitle(paste('FANTOM5'))+
  ylab('Hypergeometric test p-value')

pval.fourDG.toplot=melt(pval.fourDG)
colnames(pval.fourDG.toplot)=c('Gene','Model','pval')
pval.fourDG.toplot$Model=factor(pval.fourDG.toplot$Model,
                                levels=c('linkPeaks','Model1','Model5XorModel5Y','Model5XandModel5Y'))
p2=ggplot(pval.fourDG.toplot, aes(x=Model, y=pval, fill=Model)) +
  geom_boxplot(outlier.size = 0.5) + 
  ggtitle(paste('4DGenome'))+
  ylab('Hypergeometric test p-value')

pval.enhancerAtlas.toplot=melt(pval.enhancerAtlas)
colnames(pval.enhancerAtlas.toplot)=c('Gene','Model','pval')
pval.enhancerAtlas.toplot$Model=factor(pval.enhancerAtlas.toplot$Model,
                                       levels=c('linkPeaks','Model1','Model5XorModel5Y','Model5XandModel5Y'))
p3=ggplot(pval.enhancerAtlas.toplot, aes(x=Model, y=pval, fill=Model)) +
  geom_boxplot(outlier.size = 0.5) + 
  ggtitle(paste('EnhancerAtlas'))+
  ylab('Hypergeometric test p-value')

pval.gtex.toplot=melt(pval.gtex)
colnames(pval.gtex.toplot)=c('Gene','Model','pval')
pval.gtex.toplot$Model=factor(pval.gtex.toplot$Model,
                                       levels=c('linkPeaks','Model1','Model5XorModel5Y','Model5XandModel5Y'))
p4=ggplot(pval.gtex.toplot, aes(x=Model, y=pval, fill=Model)) +
  geom_boxplot(outlier.size = 0.5) + 
  ggtitle(paste('GTEx'))+
  ylab('Hypergeometric test p-value')


ggsave(filename='validation_enhancer_databse.pdf', plot=p1|p2|p3|p4, width=13, height=4)

