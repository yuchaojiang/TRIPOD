
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
library(nbpMatching)
library(RColorBrewer)

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

load('skin.subset.qc.rda')

# k-means clustering: using the umap coordinates and pseudotime
set.seed(1)
k <- 40
ptime.kmeans=kmeans(x = cbind(skin$umap1, skin$umap2, skin$ptime), centers = k, nstart = 25)
ptime.clusters=as.factor(ptime.kmeans$cluster)
skin$ptime.clusters=ptime.clusters
table(ptime.clusters)

p1=DimPlot(skin, reduction = "pseudo",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p2=DimPlot(skin, reduction = "pseudo",group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Seurat WNN clustering')
p3=DimPlot(skin, reduction = "pseudo",group.by='ptime.clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Kmeans clustering \n using Palantir\'s pseudotime and umaps')
p1|p2|p3

skin=skin[,!skin$ptime.clusters %in% c(14,17,37)]

# Another round of re=normalization
# RNA analysis
DefaultAssay(skin) <- "RNA"
skin <- SCTransform(skin, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(skin) <- "ATAC"
skin <- RunTFIDF(skin)
skin <- FindTopFeatures(skin, min.cutoff = 'q0')
skin <- RunSVD(skin)
skin <- RunUMAP(skin, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

skin <- FindMultiModalNeighbors(skin, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
skin <- RunUMAP(skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 1)

p1=DimPlot(skin, reduction = "pseudo",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p2=DimPlot(skin, reduction = "pseudo",group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Seurat WNN clustering')
p3=DimPlot(skin, reduction = "pseudo",group.by='ptime.clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Kmeans clustering \n using Palantir\'s pseudotime and umaps')
p1|p2|p3


# Get the branches
# select clusters based on the cluster centers
centers <- as.data.frame(ptime.kmeans$centers)
colnames(centers)=c('umap1','umap2','pseudotime')
sel.rm <- centers$umap1 > 2 & centers$umap2 > -11.75 & centers$umap1 < 4.5 & centers$umap2 < -10.25 # 2, -11.75, 4.5, -10.25
sel.pre <- centers$umap1 > -2 & centers$umap2 > -11.5 & centers$umap1 < 0.5 & centers$umap2 < -9.75 # -2, -11.5, 0, -9.75
sel.branch1 <- (centers$umap1 > 0.5 & centers$umap2 > -10.75 & centers$umap1 < 2 & centers$umap2 < -8) | # 0, -10.75, 2, -8
  (centers$umap1 > 2 & centers$umap2 > -10 & centers$umap1 < 4 & centers$umap2 < -7.75) # 2, -10, 4, -7.75
sel.branch2 <- (centers$umap1 > 0.5 & centers$umap2 > -12.5 & centers$umap1 < 2 & centers$umap2 < -10.75) | # 0, -12.5, 2, -10.75
  (centers$umap1 > 2 & centers$umap2 > -14.25 & centers$umap1 < 5.5 & centers$umap2 < -11.75) # 2, -14.25, 6.5, -11.75
sel.branch2.1 <- centers$umap1 > 5.5 & centers$umap2 > -13.5 & centers$umap1 < 9 & centers$umap2 < -12.5 # 6.5, -13.25, 9, -12.5
sel.branch2.2 <- centers$umap1 > 5.5 & centers$umap2 > -14.25 & centers$umap1 < 9 & centers$umap2 < -13.5 # 6.5, -14.25, 9, -13.5
branches <- rep(NA, nrow(centers))
branches[sel.rm] <- "removed"
branches[sel.pre] <- "pre"
branches[sel.branch1] <- "branch1"
branches[sel.branch2] <- "branch2"
branches[sel.branch2.1] <- "branch2.1"
branches[sel.branch2.2] <- "branch2.2"

branch.df <- data.frame(cluster=1:k, branch=branches)
clusters=ptime.clusters
cluster.df <- data.frame(cluster=clusters, cell=names(clusters))
merged.df <- merge(cluster.df, branch.df, by=1, all.x=T)
rownames(merged.df) <- merged.df$cell
merged.df <- merged.df[colnames(skin), ]
skin$branch=merged.df$branch

p1=DimPlot(skin, reduction = "pseudo",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p2=DimPlot(skin, reduction = "pseudo",group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Seurat WNN clustering')
p3=DimPlot(skin, reduction = "pseudo",group.by='ptime.clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Kmeans clustering \n using Palantir\'s pseudotime and umaps')
p4=DimPlot(skin, reduction = "pseudo",group.by='branch', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Branch from segmentation')

table(skin$celltype, skin$branch) # Need to finetune this.
table(skin$seurat_clusters, skin$branch)

p1+p2+p3+p4+ plot_layout(ncol=2)

skin=skin[,(skin$branch=='branch1' & skin$celltype %in% c('IRS','TAC')) |
            (skin$branch=='branch2' & skin$celltype %in% c('Hair Shaft-cuticle.cortex','TAC')) |
            (skin$branch=='branch2.1' & skin$celltype %in% c('Hair Shaft-cuticle.cortex')) |
            (skin$branch=='branch2.2' & skin$celltype %in% c('Medulla')) |
            (skin$branch=='pre' & skin$celltype %in% c('TAC')) ]

# Another round of re=normalization
# RNA analysis
DefaultAssay(skin) <- "RNA"
skin <- SCTransform(skin, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(skin) <- "ATAC"
skin <- RunTFIDF(skin)
skin <- FindTopFeatures(skin, min.cutoff = 'q0')
skin <- RunSVD(skin)
skin <- RunUMAP(skin, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

skin <- FindMultiModalNeighbors(skin, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
skin <- RunUMAP(skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 1)

k <- 40
ptime.kmeans=kmeans(x = cbind(skin$umap1, skin$umap2, skin$ptime), centers = k, nstart = 25)
ptime.clusters=as.factor(ptime.kmeans$cluster)
skin$ptime.clusters=ptime.clusters
table(ptime.clusters)

skin$seurat.celltype=skin$celltype
save(skin, file='skin.subset.qc.branch.rda')

