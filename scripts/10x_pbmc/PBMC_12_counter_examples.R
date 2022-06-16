
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

# GNLY
gene.name='GNLY'
TF.name='TBX21'
peak.start=85693711
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==peak.start) # Peak access. chr2−85693711−85699520
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g==TF.name) # TF expr. TBX21 motif MA0690.1


plot.dir='./'
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4,
                  sel.t=peak.num, jj=TF.num,
                  plot.dir=plot.dir, to.plot =c('umap'), generate.pdf = TRUE)

pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
par(mfrow=c(3,4))
# Marginal: model 1
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 1 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 2 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 4 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
# Model 5: need to include this plot in the plot_gene_peak_TF function
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Yj",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Xt",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
dev.off()

# trios that are identified by interaction, not by marginal
# false negatives
for(gene.name in genes){
  alpha.int=0.0001 # FDR threshold for interactions
  alpha.mar=0.1 # FDR threshold for marginal
  dim(xymats.all[[gene.name]]$xymats4$Xit) # cells x peaks
  dim(xymats.all[[gene.name]]$xymats4$Y.TF) # cells x TFs
  
  # Marginal peaks
  length(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$pval.c.adj<=alpha.mar)
  
  # Marginal TFs
  length(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs)$pval.c.adj<=alpha.mar)
  
  # Trio
  dim(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE))
  dim(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE))
  
  
  trio.sig=which((benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE) |
                    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE)) &
                   (xymats.all[[gene.name]]$xymats5X$gammahat2s >0) &
                   (xymats.all[[gene.name]]$xymats5Y$gammahat2s >0), arr.ind = TRUE)
  
  peak.marginal.sig=which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$pval.c.adj<=alpha.mar)
  TF.marginal.sig=which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs)$pval.c.adj<=alpha.mar)
  
  trio.sig.non.marginal=trio.sig[(!trio.sig[,1] %in% peak.marginal.sig) &
                                   (!trio.sig[,2] %in% TF.marginal.sig) ,, drop=FALSE]
  
  if(length(trio.sig.non.marginal)>0){
    for(t in 1:nrow(trio.sig.non.marginal)){
      peak.num=trio.sig.non.marginal[t,1]
      TF.num=trio.sig.non.marginal[t,2]
      plot.dir='./counter_example_FN/'

      pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
      par(mfrow=c(3,4))
      # Marginal: model 1
      plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
      # Model 1 scatterplot
      plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
      # Model 2 scatterplot
      plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
      # Model 4 scatterplot
      plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
      # Model 5: need to include this plot in the plot_gene_peak_TF function
      res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                      Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                      Yg=xymats.all[[gene.name]]$xymats5X$Y,
                      match.by="Yj",
                      do.plot=TRUE,
                      col=wnn.celltype.col,
                      verbose=FALSE,
                      partition.screen=FALSE)
      res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                      Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                      Yg=xymats.all[[gene.name]]$xymats5X$Y,
                      match.by="Xt",
                      do.plot=TRUE,
                      col=wnn.celltype.col,
                      verbose=FALSE,
                      partition.screen=FALSE)
      dev.off()
      
    }
  }
}


# trios that are identified by marginal, not by interaction
# false positives
for(gene.name in genes){
  alpha.int=0.1 # FDR threshold for interactions (very insignificant for interaction)
  alpha.mar=0.01 # FDR threshold for marginal (very significant for marginal)
  dim(xymats.all[[gene.name]]$xymats4$Xit) # cells x peaks
  dim(xymats.all[[gene.name]]$xymats4$Y.TF) # cells x TFs
  
  # Marginal peaks
  length(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$pval.c.adj<=alpha.mar)
  
  # Marginal TFs
  length(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs)$pval.c.adj<=alpha.mar)
  
  # Trio
  dim(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE))
  dim(benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE))
  
  
  trio.sig=which((benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE) |
                    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, fdr.thresh=alpha.int, do.plot = FALSE) |
                    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha.int, do.plot = FALSE) |
                    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, fdr.thresh=alpha.int, do.plot = FALSE) |
                    benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, fdr.thresh=alpha.int, do.plot = FALSE)),
                 arr.ind = TRUE)
  
  peak.marginal.sig=which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalXits)$pval.c.adj<=alpha.mar &
                            xymats.all[[gene.name]]$xymats1$betaXits>0)
  TF.marginal.sig=which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalYijs)$pval.c.adj<=alpha.mar &
                          xymats.all[[gene.name]]$xymats1$betaYijs>0)
  peakxTF.marginal.sig=which(benjaminiHochsbergVectorAdjust(xymats.all[[gene.name]]$xymats1$pvalgammas)$pval.c.adj<=alpha.mar &
                               xymats.all[[gene.name]]$xymats1$gammas>0, arr.ind=TRUE)
  
  # sig in marginal, insig in interaction
  peak.marginal.sig.non.trio=peak.marginal.sig[(!peak.marginal.sig%in%trio.sig[,1])]
  TF.marginal.sig.non.trio=TF.marginal.sig[(!TF.marginal.sig%in%trio.sig[,2])]
  
  
  if(length(peak.marginal.sig.non.trio)>0 & length(TF.marginal.sig.non.trio)>0){
    for(peak.num in peak.marginal.sig.non.trio){
      for(TF.num in TF.marginal.sig.non.trio){
        if(xymats.all[[gene.name]]$xymats4$peakxmotif.g[peak.num, TF.num]==0) next
        plot.dir='./counter_example_FP/'

        pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
        par(mfrow=c(3,4))
        # Marginal: model 1
        plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
        # Model 1 scatterplot
        plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
        # Model 2 scatterplot
        plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
        # Model 4 scatterplot
        plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
        # Model 5: need to include this plot in the plot_gene_peak_TF function
        res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                        Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                        Yg=xymats.all[[gene.name]]$xymats5X$Y,
                        match.by="Yj",
                        do.plot=TRUE,
                        col=wnn.celltype.col,
                        verbose=FALSE,
                        partition.screen=FALSE)
        res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                        Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                        Yg=xymats.all[[gene.name]]$xymats5X$Y,
                        match.by="Xt",
                        do.plot=TRUE,
                        col=wnn.celltype.col,
                        verbose=FALSE,
                        partition.screen=FALSE)
        dev.off()
      }
    }
  }
}

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


# Plot the marginals with colors
gene.name='ZBTB16'
peak.num=25
TF.num=4

gene.name='RIN3'
peak.num=37
TF.num=12

gene.name='FGL2'
peak.num=12
TF.num=236

plot.dir='./'
pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot_recolored.pdf',sep=''), height=2.5, width=9)
par(mfrow=c(1,4))
source('plot_gene_peak_TF_col.R')
# Marginal: model 1
plot_gene_peak_TF_col(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)

source('matchedTest_col.R')
res=matchedTest_col(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Yj",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
# res=matchedTest_col(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
#                 Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
#                 Yg=xymats.all[[gene.name]]$xymats5X$Y,
#                 match.by="Xt",
#                 do.plot=TRUE,
#                 col=wnn.celltype.col,
#                 verbose=FALSE,
#                 partition.screen=FALSE)
dev.off()

par(mfrow=c(1,1))
a <- sample(1:100)
rbPal <- colorRampPalette(c('red','blue'))
b <- rbPal(10)[as.numeric(cut(a,breaks = 10))]
cuts<-levels(cut(a,breaks = 10))
cuts<-gsub(","," - ",cuts)
cuts<-gsub("\\(","[",cuts)
cuts=cuts[10:1]
pdf(file='counter_example_legend.pdf', width=4, height=4)
plot(a,col='white',pch=16)
legend("top",cuts,col=rbPal(10),pch=16)
dev.off()


