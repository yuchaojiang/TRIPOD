
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

# LinkPeaks
linkpeak.pos=read.csv('linkPeak/pos.df.csv')
linkpeak.neg=read.csv('linkPeak/pos.df.csv')
linkpeak=rbind(linkpeak.pos, linkpeak.neg)

pbmc=SetIdent(pbmc,value = 'seurat.celltype')
DimPlot(pbmc, reduction = 'umap.rna', label=T)

gene.name='CCR7'
gene.name='GNLY'
gene.name='FCGR3A'
gene.name='MS4A1'

FeaturePlot(pbmc, features = paste("sct_",gene.name,sep=''), reduction = 'umap.rna', 
            label=T, label.size = 2, cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
VlnPlot(pbmc, features=paste("sct_",gene.name,sep=''))

# Peak Venn Diagram
alpha=0.05

for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){
  
  peak.gr.g=xymats.all[[gene.name]]$xymats1$peak.gr.g
  linkpeak.sig=unique(linkpeak[linkpeak$gene==gene.name & linkpeak$adj <=alpha,2])
  
  marginal.sig=unique(GRangesToString(peak.gr.g[benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$pval.c.adj<=alpha]))
  
  model4.sig= unique(GRangesToString(peak.gr.g[unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha, do.plot = FALSE) |
                                                       benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalXits, fdr.thresh=alpha, do.plot = FALSE), arr.ind = T)[,1])]))
  
  model5.sig= unique(GRangesToString(peak.gr.g[unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE)  |
                                                       benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, fdr.thresh=alpha, do.plot = FALSE) , arr.ind = T)[,1])]))
  x=list(linkpeak=linkpeak.sig,
         marginal=marginal.sig,
         nonparametric=model5.sig)
  library(ggvenn)
  p=ggvenn(x, 
           fill_color = c("#F8766D", "#00BA38", "#619CFF"),
           stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+ggtitle(gene.name)
  ggsave(filename = paste0('sig/sig_peak_venn_',gene.name,'.pdf'),p, height=2.5, width=2.5)
  
}

# TF venn diagram
for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){
  
  TF.g=xymats.all[[gene.name]]$xymats1$TF.g

  marginal.sig=unique(TF.g[benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs)$pval.c.adj<=alpha])
  
  model4.sig= unique(TF.g[unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha, do.plot = FALSE) |
                                                       benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, fdr.thresh=alpha, do.plot = FALSE), arr.ind = T)[,2])])
  
  model5.sig= unique(TF.g[unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE) |
                                                       benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, fdr.thresh=alpha, do.plot = FALSE) , arr.ind = T)[,2])])
  x=list(marginal=marginal.sig,
         nonparametric=model5.sig)
  library(ggvenn)
  p=ggvenn(x, 
           fill_color = c("#00BA38", "#619CFF"),
           stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+ggtitle(gene.name)
  ggsave(filename = paste0('sig/sig_TF_venn_',gene.name,'.pdf'),p, height=2.5, width=2.5)
  
}

# Trio overlap with conditional
for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){

  model4.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha, do.plot = FALSE) |
                                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalXits, fdr.thresh=alpha, do.plot = FALSE) |
                                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs, fdr.thresh=alpha, do.plot = FALSE)))
  
  model5X.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE) |
                                  benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest,fdr.thresh=alpha, do.plot = FALSE)))
  
  model5Y.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE)  |
                              benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, fdr.thresh=alpha, do.plot = FALSE)))
  
  x=list(model4.sig=model4.sig,
         nonparametricX=model5X.sig,
         nonparametricY=model5Y.sig)
  library(ggvenn)
  p=ggvenn(x, 
           fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
           stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+ggtitle(gene.name)
  ggsave(filename = paste0('sig/sig_trio_venn_',gene.name,'.pdf'),p, height=2.5, width=2.5)
  
}

# Trio overlap without conditional just interaction
for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){

  model4.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha, do.plot = FALSE)))
  
  model5X.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE)))
  
  model5Y.sig= unique(which(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha, do.plot = FALSE)))
  
  x=list(model4.sig=model4.sig,
         nonparametricX=model5X.sig,
         nonparametricY=model5Y.sig)
  library(ggvenn)
  p=ggvenn(x, 
           fill_color =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
           stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+ggtitle(gene.name)
  ggsave(filename = paste0('sig/sig_trio_venn_int_',gene.name,'.pdf'),p, height=2.5, width=2.5)
  
}
