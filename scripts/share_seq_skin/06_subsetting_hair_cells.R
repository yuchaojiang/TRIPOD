
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

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

load('skin.chromvar.qc.rda')
table(skin$seurat.celltype)

#############################################
# Read in Sai's coordiinates and cells (subset of cells)
#############################################

ptime=read.table('ptime.bed',sep='\t')

head(ptime)
ptime$bc=gsub('Trial60,skin.','',ptime$bc)
ptime$bc=gsub(',','.',ptime$bc)
ptime=ptime[ptime$celltype!='TAC-2',] # Remove TAC-2 group: cells are too noisy.
ptime$celltype[ptime$celltype=='TAC-1']='TAC' # Rename TAC-1 to TAC

ct=read.table('celltype_skin.txt', header = T, sep='\t')
ptime$rna.bc=ct[match(ptime$bc, ct[,2]),1]
any(is.na(ptime$rna.bc))

ptime=ptime[ptime$rna.bc %in% colnames(skin),]
skin=skin[,match(ptime$rna.bc, colnames(skin))] # Only select these subset of cells

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

dim(skin)
dim(ptime)
rownames(ptime)=ptime$rna.bc
all(colnames(skin)==rownames(ptime))
skin=AddMetaData(skin, ptime)
all(skin$celltype==ptime$celltype)
save(skin, file='skin.subset.rda')

#############################################
# Look at some reduced dimensions
#############################################

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

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

load('skin.subset.rda')
# Add Sai's umap coordinates as a reduction
skin[["pseudo"]] <- CreateDimReducObject(embeddings = cbind(umap1=skin$umap1, umap2=skin$umap2), key = "pseudoumap_", assay = DefaultAssay(skin))

# Ground truth umap coordinates: sent by Sai
DimPlot(skin, reduction = "pseudo", group.by = 'celltype', 
        label = TRUE, label.size = 3, repel = TRUE) + NoLegend()

# WNN clustering results 
p1=DimPlot(skin, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('WNN clustering')
p2=DimPlot(skin, reduction = "wnn.umap",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p1|p2

table(skin$celltype, skin$seurat_clusters)

skin=skin[,skin$seurat_clusters!=12]
skin$seurat_clusters=droplevels(skin$seurat_clusters)
table(skin$celltype, skin$seurat_clusters)

skin=skin[,(skin$celltype=='Hair Shaft-cuticle.cortex' & skin$seurat_clusters %in% c(2,3)) |
            (skin$celltype=='IRS' & skin$seurat_clusters %in% c(1,10))|
            (skin$celltype=='Medulla' & skin$seurat_clusters %in% c(6,11))|
            (skin$celltype=='TAC' & skin$seurat_clusters %in% c(0,4,5,7,8,9))]
skin$seurat_clusters=droplevels(skin$seurat_clusters)

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

p1=DimPlot(skin, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('WNN clustering')
p2=DimPlot(skin, reduction = "wnn.umap",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p1|p2
table(skin$celltype, skin$seurat_clusters)

p1=DimPlot(skin, reduction = "pseudo",group.by='celltype', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p2=DimPlot(skin, reduction = "pseudo",group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p1|p2

skin=skin[,!skin$seurat_clusters %in% c(9,11)]

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
p2=DimPlot(skin, reduction = "pseudo",group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE)+ggtitle('Annotated celltype')
p1|p2

save(skin, file='skin.subset.qc.rda')

