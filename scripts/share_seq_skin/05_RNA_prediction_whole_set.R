#######################################
###  Running MAESTRO
#######################################

# remotes::install_version("Seurat", version = "3.2.3")
# devtools::install_version(package = 'Signac', version = package_version('0.2.5'))
library(Seurat) # MAESTRO only works with Seurat V3!!
library(Signac)
# devtools::install_github("liulab-dfci/MAESTRO")
library(MAESTRO) 
library(dplyr)
library(ggplot2)
library(reticulate)

# python3.8 -m pip install tables
# python3.8 -m pip install h5py
# python3.8 -m pip install scipy
# python3.8 -m pip install pandas

# Need to load Python
use_python('/Users/yuchaojiang/pythonProject/bin/python', required = TRUE)

load('ATAC_RNA_WNN.rda')

for(dist in c(2000, 5000, 20000, 50000, 200000)){
  cat(dist,'\n')
  wnn.gene.maestro <- ATACCalculateGenescore(t(wnn.peak), organism = "GRCm38", decaydistance = dist, model = "Enhanced")
  wnn.gene.maestro <- t(as.matrix(wnn.gene.maestro))
  save(wnn.gene.maestro, file=paste0('prediction/wnn.gene.maestro.',dist,'.rda'))
}

install.packages('Seurat') # Re-install Seurat V4
install.packages('Signac')

#######################################
###  Look at prediction results across all models
#######################################

library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(presto)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

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

# Model 4 results
load('ATAC_RNA_WNN.rda')
save(peak.gr, file='prediction/peak.gr.rda') # Needed for testing on 3K cells

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

# Get highly variable genes
celltype_presto <- presto:::wilcoxauc.Seurat(X = skin, group_by = 'seurat.celltype', assay = 'data', seurat_assay = 'SCT')
celltype_genes=celltype_presto[order(celltype_presto[,'pval'])[1:1000],'feature'] # cell-type markers
highvar_genes=skin@assays$SCT@var.features[1:1000]
genes=unique(c(celltype_genes, highvar_genes))

# Visualize reduced dimensions and some gene markers
p1=DimPlot(skin, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE)+NoLegend()+ggtitle('ATAC reduced dim.')
p2=DimPlot(skin, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE)+NoLegend()+ggtitle('RNA reduced dim.')
p3=DimPlot(skin, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)+NoLegend()+ggtitle('WNN reduced dim.')
p4=DimPlot(skin, reduction = "wnn.umap", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+NoLegend()+ggtitle('WNN partition')
ggsave(filename ='skin_whole_set_reduced_dim.pdf',
       plot = p1+p2+p3+p4+ plot_layout(ncol=2),
       width = 8,
       height = 8)

p3=DimPlot(skin, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)+NoLegend()+ggtitle('WNN reduced dim.')
ggsave(filename='skin_whole_set_celltype.pdf',
       plot=p3,
       width=5, height=5)

p3=DimPlot(skin, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('WNN reduced dim.')
ggsave(filename='skin_whole_set_celltype_legend.pdf',
       plot=p3,
       width=12, height=12)

# Look at some gene markers
gene.name='Ano7'
gene.name='Padi3'
p1=FeaturePlot(skin, features = paste("sct_",gene.name,sep=''), reduction = 'wnn.umap', label=T, label.size = 3, cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
#RidgePlot(skin, features = paste("sct_",gene.name,sep=''), log = FALSE)
ggsave(filename=paste0('gene_rna_',gene.name,'_featureplot.pdf'),
       plot=p1, width=6, height=6)
p2=DotPlot(skin, features = paste("sct_",gene.name,sep=''))
ggsave(filename=paste0('gene_rna_',gene.name,'_dotplot.pdf'),
       plot=p2, width=6, height=6)

load(paste0('prediction/wnn.gene.maestro.',2000,'.rda'))
genes=intersect(genes, colnames(wnn.gene.maestro))

for(dist in c(2000, 5000, 20000, 50000, 200000)){
  cat('Distance:',dist,'\n\n')
  load(paste0('prediction//wnn.gene.maestro.',dist,'.rda'))
  ext.upstream=dist
  modelNames=c('-1a','0a','3a','maestro')
  pred.r=matrix(ncol=length(modelNames), nrow=length(genes))
  colnames(pred.r)=modelNames
  rownames(pred.r)=genes
  i=1
  for(gene.name in genes){
    if(i%%20==0) cat(i,'\t')
    for(m in 1:(length(modelNames)-1)){
      xymats=getXYMatrices(gene.name=gene.name,
                           ext.upstream=ext.upstream,
                           transcripts.gr = transcripts.gr,
                           peak.gr=peak.gr,
                           wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                           peakxmotif=peakxmotif, motifxTF=motifxTF,
                           wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
      if(ncol(xymats$Xit)==0 | ncol(xymats$X)==0) break # No peak overlap with the region
      xymats=fitModel(xymats, modelName=modelNames[m])
      pred.r[gene.name,modelNames[m]]=xymats$LOO.r
    }
    pred.r[gene.name,4]=round(cor(xymats$Y, wnn.gene.maestro[,gene.name], method = 'pearson'),3)
    if(all(!is.na(pred.r[gene.name,]))) i=i+1
    if(i >=400) break # Here I only looked at the first 200 genes within each dist.
  }
  # remove NAs
  pred.r=pred.r[apply(pred.r,1,function(x){all(!is.na(x))}),]
  colnames(pred.r)[1:3]=c('Gene.Activity','Peak.LASSO','Peak.TF.Lasso')
  pred.r=pred.r[,c(1,4,2,3)]
  save(pred.r, file=paste0('prediction/pred.r.',dist,'.rda'))
  cat('\n\n')
}

# load 20000 results
dist=20000
load(paste0('prediction/pred.r.',dist,'.rda'))
head(pred.r)

gene.name='Ano7'
gene.name='Padi3'
gene.name='Aldh1a3'

# This is to generate a plot for each gene with varying extended window size
for(dist in c(2000, 5000, 20000, 50000, 200000)){
  cat(dist,'\n')
  ext.upstream=dist
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr = transcripts.gr,
                       peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  if(ncol(xymats$Xit)==0 | ncol(xymats$Y.TF)==0) next
  pdf(file=paste0('prediction/RNA_prediction_',gene.name,'_',dist,'.pdf'), width=8, height=2)
  par(mfrow=c(1,4))
  # Model -1: sum of access. (alignment methods)
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr = transcripts.gr,
                       peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  xymats=fitModel(xymats, do.pred.plot = TRUE, modelName='-1a')
  
  # MAESTRO
  load(paste0('prediction/wnn.gene.maestro.',dist,'.rda'))
  plot(xymats$Y, wnn.gene.maestro[,gene.name], col=xymats$wnn.celltype.col, xlab='Gene expression', ylab='MAESTRO\'s regulatory potential', 
       pch=16, main=paste('MAESTRO\'s Regulatory Potential'))
  legend('bottomright',paste('r =',round(cor(xymats$Y, wnn.gene.maestro[,gene.name], method = 'pearson'),3)), bty='n')
  
  # Model 0: lasso on peak (sci-CAR)
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr = transcripts.gr,
                       peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  xymats=fitModel(xymats, do.pred.plot = TRUE, modelName='0a')
  
  # Model 3: lasso on peak and TF pair
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr = transcripts.gr,
                       peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  xymats=fitModel(xymats, do.pred.plot = TRUE, modelName='3a')
  dev.off() 
}

library(reshape)
library(ggplot2)
dist=2000
load(paste0('prediction/pred.r.',dist,'.rda')) 
pred.r=melt(pred.r)
pred.r=cbind(pred.r, dist=rep(dist, nrow(pred.r)))
pred.r.all=pred.r
for(dist in c(5000, 20000, 50000, 200000)){
  load(paste0('prediction/pred.r.',dist,'.rda'))  
  pred.r=melt(pred.r)
  pred.r=cbind(pred.r, dist=rep(dist, nrow(pred.r)))
  pred.r.all=rbind(pred.r.all, pred.r)
}

colnames(pred.r.all)=c('Gene','Model','Corr','Dist')
pred.r.all$Dist=as.factor(pred.r.all$Dist)

p=ggplot2::ggplot(pred.r.all, aes(x=Dist, y=Corr, fill=Model)) +
  geom_boxplot(outlier.size = 0.2, lwd=0.2) + ylim(-0.5, 1)+
  ggtitle(paste('RNA prediction using ATAC'))+
  ylab('Correlation coefficient between Y and Yhat')
ggsave(filename='skin_pred_cor.pdf',
       plot=p,
       width=5, height=3)


