
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

# Get hierarchical clustering results for all trios
temp=aggregate(wnn.rna ~ wnn.celltype, FUN=sum)
rownames(temp)=temp[,1]
temp=temp[,-1]
temp=sweep(temp,1,rowSums(temp)/10^6,'/')
end <- temp %>%scale %>%dist %>% hclust %>% as.dendrogram
rm(temp)

# CCR7
gene.name='CCR7'
TF.name='LEF1'
peak.start=40596032

peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==peak.start) 
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g==TF.name) 

xymats.all=get.influence(xymats.all=xymats.all,
                         gene.name=gene.name, 
                         peak.num=peak.num, TF.num=TF.num, 
                         object=pbmc,
                         dend=dend,
                         to.test=c('WNN','CellType','Tree','Footprint'))

# Get WNN-specific reduction and find each WNN's nearest neighbors
reduction='wnn.umap'
pbmc.reduction=pbmc@reductions[[reduction]]@cell.embeddings
wnn.reduction=matrix(nrow=length(levels(pbmc$seurat_clusters)), ncol=ncol(pbmc.reduction))
rownames(wnn.reduction)=paste('wnn_', levels(pbmc$seurat_clusters),sep='')
colnames(wnn.reduction)=colnames(pbmc.reduction)
for(i in 1:nrow(wnn.rna)){
  cat(i,' ')
  wnn.reduction[i,]=apply(pbmc.reduction[pbmc$seurat_clusters==(i-1),], 2, mean)
}
plot(pbmc.reduction[,1], pbmc.reduction[,2], col=pbmc$celltype, pch=16, cex=0.5)
points(wnn.reduction[,1],wnn.reduction[,2], cex=1)
points(wnn.reduction[,1],wnn.reduction[,2], cex=1, col='orange', pch=16)


library(FNN)
k=4 # How many neighbors to include?
wnn.neighbors=get.knn(wnn.reduction, k=k)$nn.index
wnn.neighbors=cbind(1:nrow(wnn.neighbors), wnn.neighbors)

xymats.all=get.influence(xymats.all=xymats.all,
                         gene.name=gene.name, 
                         peak.num=peak.num, TF.num=TF.num, 
                         object=pbmc,
                         dend=dend,
                         wnn.neighbors=wnn.neighbors)

