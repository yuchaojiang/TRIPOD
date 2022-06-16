
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
library(gridBase)
library(grid)

load('ATAC_RNA_WNN_footprint.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

# specific genes
gene.name='Krt71' # IRS
TF.name=toupper('Gata3')
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==101743171) # Peak access. 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='Gata3') # TF expr. 

gene.name='Gpnmb' # Medulla
TF.name=toupper('Fosl2')
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==49036381) # Peak access. 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='Fosl2') # TF expr. 

gene.name='Adamtsl1' # Hair Shaft
TF.name=toupper('Lef1')
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==86024870) # Peak access. 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='Lef1') # TF expr. 

# Get hierarchical clustering results for all trios
# should be using all the genes; gene-peak-TF cluster does not make much sense and would not be as stable
temp=aggregate(wnn.rna ~ wnn.celltype, FUN=sum)
rownames(temp)=temp[,1]
temp=temp[,-1]
temp=sweep(temp,1,rowSums(temp)/10^6,'/')
dend <- temp %>%scale %>%dist %>% hclust %>% as.dendrogram
rm(temp)

# Get WNN-specific reduction and find each WNN's nearest neighbors
reduction='pseudo'
skin.reduction=skin@reductions[[reduction]]@cell.embeddings
skin.pseudotime=as.matrix(skin$ptime)
wnn.reduction=matrix(nrow=length(levels(skin$seurat_clusters)), ncol=ncol(skin.reduction))
rownames(wnn.reduction)=paste('wnn_', levels(skin$seurat_clusters),sep='')
colnames(wnn.reduction)=colnames(skin.reduction)
wnn.pseudotime=matrix(nrow=length(levels(skin$seurat_clusters)), ncol=ncol(skin.pseudotime))
for(i in 1:nrow(wnn.rna)){
  cat(i,' ')
  wnn.reduction[i,]=apply(skin.reduction[skin$seurat_clusters==(i-1),], 2, mean)
  wnn.pseudotime[i,]=apply(skin.pseudotime[skin$seurat_clusters==(i-1),,drop=FALSE], 2, mean)
}

plot(skin.reduction[,1], skin.reduction[,2], col=skin$celltype, pch=16, cex=0.5)
points(wnn.reduction[,1],wnn.reduction[,2], cex=1, col='orange', pch=16)

library(FNN)
k=4 # How many wnn neighbors to include?
wnn.neighbors=get.knn(cbind(wnn.reduction, wnn.pseudotime), k=k)$nn.index
wnn.neighbors=cbind(1:nrow(wnn.neighbors), wnn.neighbors)

k=200 # How many cell neighbors to include?
cell.neighbors=get.knn(cbind(skin.reduction, skin.pseudotime), k=k)$nn.index
cell.neighbors=cbind(1:nrow(cell.neighbors), cell.neighbors)

# Let's remove one WNN at a time and plot it first
xymats.all=get.influence(xymats.all=xymats.all,
                         gene.name=gene.name, 
                         peak.num=peak.num, TF.num=TF.num, 
                         object=skin,
                         dend=dend,
                         wnn.neighbors = wnn.neighbors,
                         cell.neighbors = cell.neighbors,
                         reduction='pseudo')

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
