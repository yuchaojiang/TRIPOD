
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
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
library(RColorBrewer)

load('ATAC_RNA_WNN_branch.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

# Reduced dimensions
p1=DimPlot(skin, reduction = "pseudo",group.by='celltype', label = TRUE, label.size = 3, repel = TRUE)+ggtitle('Annotated celltype') + NoLegend()
p2=DimPlot(skin, reduction = "pseudo",group.by='wsnn_res.0.8', label = TRUE, label.size = 3, repel = TRUE)+ggtitle('Seurat WNN clustering') + NoLegend()
p3=DimPlot(skin, reduction = "pseudo",group.by='ptime.clusters', label = TRUE, label.size = 3, repel = TRUE)+ggtitle('Kmeans clustering \n using Palantir\'s pseudotime and umaps') + NoLegend()
p4=DimPlot(skin, reduction = "pseudo",group.by='branch', label = TRUE, label.size = 3, repel = TRUE)+ggtitle('Branch from segmentation') + NoLegend()
# Pseudotime
p5=FeaturePlot(skin, features ='ptime', reduction = 'pseudo', label = TRUE, cols = c("lightgrey", "darkblue"))+ggtitle('Pseudotime')

ggsave(filename ='skin_hair_subset_reduced_dim.pdf',
       plot = p1+p2+p3+p4+p5+ plot_layout(ncol=3),
       width = 15,
       height = 10)

# Feature plots for some reported markers...
gene.name='Wnt3'
p1=FeaturePlot(skin, features = paste("sct_",gene.name,sep=''), reduction = 'pseudo', label = TRUE, cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
gene.name='Gpnmb'
p2=FeaturePlot(skin, features = paste("sct_",gene.name,sep=''), reduction = 'pseudo', label = TRUE, cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
p3=DotPlot(skin, features = paste("rna_",c('Wnt3','Gpnmb'), sep=''), group.by = 'seurat.celltype') + ggtitle('Gene Expression across cell types')+ylab('')
p4=DotPlot(skin, features = paste("rna_",c('Wnt3','Gpnmb'), sep=''), group.by = 'ptime.clusters') + ggtitle('Gene Expression across WNNs')+ylab('WNN')
ggsave(file='skin_hair_Wnt3_Gpnmb.pdf',
       plot=p1+p2+p3+p4+plot_layout(ncol=2),
       width=12, height=9)

# Trio regulatory model
# Remove peaks that are too wide
peak.filter=width(peak.gr)<6000
skin@assays$ATAC=subset(skin@assays$ATAC,rownames(skin@assays$ATAC)[width(skin@assays$ATAC@ranges)<=6000])
peak.gr=peak.gr[peak.filter]
wnn.peak=wnn.peak[,peak.filter]
peakxmotif=peakxmotif[peak.filter,]
rm(peak.filter)

ext.upstream=100000 # This is actually both up and downstream
transcripts <- transcripts.gr
TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
TSS <- GRanges(seqnames = seqnames(transcripts),
               ranges = IRanges(start = TSS.position, width = 1),
               strand = strand(transcripts),
               gene_name = transcripts$gene_name)
transcripts.ext.gr=getRegion(c(ext.upstream, ext.upstream), TSS)
rm(transcripts); rm(TSS.position); rm(TSS)

# DORC genes: these are what we want to focus on
genes = as.matrix(read.table('dorcs.txt'))[,1]
genes = genes[genes %in% rownames(skin@assays$RNA)]

# Below are the genes included in Sai's paper:
marker_genes=c('Gpnmb','Wnt3','Tubb6','Plpp3','Hoxc13','Foxp1',
               'Jag1','Krt14','Nfkbia','Wnt10b','Egr3',
               'Itga6','Id2','Ccnd2','Tnfaip3','Wnt5a')
marker_genes = marker_genes[marker_genes %in% rownames(skin@assays$RNA)] # marker genes in RNA
marker_genes %in% genes

genes=unique(c(marker_genes, genes))
gene.xymats4=vector(mode = "list", length = length(genes))
names(gene.xymats4)=genes
gene.xymats5X=gene.xymats5Y=gene.xymats2=gene.xymats1=gene.xymats4

for(i in 1:length(genes)){
  if(i%%50==0 | i<=5) cat(i,'\t')
  gene.name=genes[i]
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr=transcripts.gr, peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  if(ncol(xymats$Xit)==0 | ncol(xymats$Y.TF)==0) next
  
  # Model 1
  xymats1=fitModel(xymats, modelName='1a')
  gene.xymats1[[gene.name]]=xymats1
  
  # Model 2
  xymats2=fitModel(xymats, modelName='2a')
  gene.xymats2[[gene.name]]=xymats2
  
  # Model 4
  xymats4=fitModel(xymats, modelName='4a')
  gene.xymats4[[gene.name]]=xymats4
  
  # Model 5
  xymats5X=fitModel(xymats, modelName='5a', match.by='Xt', do.plot=FALSE, VERBOSE=FALSE)
  xymats5Y=fitModel(xymats, modelName='5a', match.by='Yj', do.plot=FALSE, VERBOSE=FALSE)
  gene.xymats5X[[gene.name]]=xymats5X
  gene.xymats5Y[[gene.name]]=xymats5Y
}

gene.filter=rep(TRUE, length(gene.xymats5X)) # Remove highly variable genes that dont have binding motifs for any TFs
for(i in 1:length(gene.xymats5X)){
  if(is.null(gene.xymats5X[[i]])) gene.filter[i]=FALSE
}
sum(gene.filter)
gene.xymats1=gene.xymats1[gene.filter]
gene.xymats2=gene.xymats2[gene.filter]
gene.xymats4=gene.xymats4[gene.filter]
gene.xymats5X=gene.xymats5X[gene.filter]
gene.xymats5Y=gene.xymats5Y[gene.filter]


# Combine all xymats
genes=names(gene.xymats1)
xymats.all=vector(mode = "list", length = length(genes))
names(xymats.all)=genes
for(i in 1:length(xymats.all)){
  if(i %%50==0) cat(i,'\t')
  gene.name=genes[i]
  xymats.all[[gene.name]]=list(xymats1=gene.xymats1[[gene.name]],
                               xymats2=gene.xymats2[[gene.name]],
                               xymats4=gene.xymats4[[gene.name]],
                               xymats5X=gene.xymats5X[[gene.name]],
                               xymats5Y=gene.xymats5Y[[gene.name]])
}
rm(gene.xymats1); rm(gene.xymats2); rm(gene.xymats4); rm(gene.xymats5X); rm(gene.xymats5Y)
rm(xymats1); rm(xymats2); rm(xymats4); rm(xymats5X); rm(xymats5Y)

save.image(file='ATAC_RNA_WNN_branch_integrated_high_var.rda')




library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
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
library(RColorBrewer)

load('ATAC_RNA_WNN_branch_integrated_high_var.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

# Gene Gpnmb
gene.name='Gpnmb'
#########################
# ADD Yj ~ Xt
#########################
names(xymats.all[[gene.name]]$xymats1) # Yg~Xt, Yg~Yj, Yg~Xt*Yj
names(xymats.all[[gene.name]]$xymats2) # Yg~Xt+Yj
names(xymats.all[[gene.name]]$xymats4) # Yg~Xt+Yj+Xt*Yj
names(xymats.all[[gene.name]]$xymats5X) # Match by Xt
names(xymats.all[[gene.name]]$xymats5Y) # Match by Yj

## PLOTTING  A TRIO
# Branch: TAC --> IRS; TAC --> inter --> Hair Shaft & Medulla

genes.selected=c('Gpnmb', # Medulla
                 'Krt71',# IRS
                 'Adamtsl1', # Hair Shaft
                 'Lgr6') # TAC

# specific genes
gene.name='Gpnmb'
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==49036381) # Peak access. 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='Fosl2') # TF expr. 

gene.name='Krt71'
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==101743171) # Peak access. 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='Gata3') # TF expr. 

plot.dir='./'
# umap plot for this trio
plot_gene_peak_TF(object=skin, group.by='celltype', reduction='pseudo', xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('umap'), generate.pdf = TRUE, featurePlot = 'Dot')

pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
par(mfrow=c(3,4))
# Marginal: model 1
plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 1 scatterplot
plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 2 scatterplot
plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 4 scatterplot
plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
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

for (gene.name in genes.selected){
  cat(gene.name,':\t\t')
  filter=which(xymats.all[[gene.name]]$xymats4$gammas>0 &
                 xymats.all[[gene.name]]$xymats1$betaXits.corr>0 &
                 xymats.all[[gene.name]]$xymats1$betaYijs.corr>0 &
                 xymats.all[[gene.name]]$xymats1$betaYijsOnXits.corr>0 &
                 xymats.all[[gene.name]]$xymats1$pvalBetaYijsOnXits<0.05, arr.ind=T)
  
  
  filter=filter[xymats.all[[gene.name]]$xymats1$pvalXits[filter[,1]]<0.05 &
                  xymats.all[[gene.name]]$xymats1$pvalYijs[filter[,2]]<0.05,, drop=FALSE]
  if(nrow(filter)==0) next
  gamma.pval=rep(NA, nrow(filter))
  for(i in 1:nrow(filter)){
    gamma.pval[i]=xymats.all[[gene.name]]$xymats4$pvalgammas[filter[i,1],filter[i,2]]
  }
  min.gamma.pval=min(gamma.pval)
  cat('Min p-value',min.gamma.pval,'\n')
  if(min.gamma.pval > 0.05) next
  peak.num=which(xymats.all[[gene.name]]$xymats4$pvalgammas==min.gamma.pval, arr.ind=T)[,1]
  TF.num=which(xymats.all[[gene.name]]$xymats4$pvalgammas==min.gamma.pval, arr.ind=T)[,2]
  
  plot.dir='./trio_selected_cases/'
  # umap plot for this trio
  plot_gene_peak_TF(object=skin, group.by='celltype', reduction='pseudo', xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('umap'), generate.pdf = TRUE, featurePlot = 'Dot')
  
  pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
  par(mfrow=c(3,4))
  # Marginal: model 1
  plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
  # Model 1 scatterplot
  plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
  # Model 2 scatterplot
  plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
  # Model 4 scatterplot
  plot_gene_peak_TF(object = skin, group.by='celltype', xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
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

