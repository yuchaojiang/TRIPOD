
library(Seurat)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)

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

load('pbmc_chromvar.rda')

pbmc=pbmc[,pbmc$celltype!='6_4'] # Remove a subcluster that failed to be annotated
pbmc$celltype=droplevels(pbmc$celltype)
DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

# Remove potential placelet cell: marker gene PPBP
table(pbmc@assays$RNA@counts['PPBP',])
pbmc=pbmc[,pbmc@assays$RNA@counts['PPBP',]==0]

dim(pbmc)
table(pbmc$celltype) # Annotation from Seurat's WNN tutorial; not well matched


# DC
VlnPlot(pbmc, features = c("rna_CST3"), slot = "counts", log = TRUE, group.by = 'celltype')
# NK 
VlnPlot(pbmc, features = c("rna_GNLY"), slot = "counts", log = TRUE, group.by = 'celltype')
# B cells
VlnPlot(pbmc, features = c("rna_MS4A1"), slot = "counts", log = TRUE, group.by = 'celltype')
# CD14+ mono
VlnPlot(pbmc, features = c("rna_CD14"), slot = "counts", log = TRUE, group.by = 'celltype')
# CD8 Naive, CD8 TEM_1, CD8 TEM_2
VlnPlot(pbmc, features = c("rna_CD8A"), slot = "counts", log = TRUE, group.by = 'celltype')


reference <- LoadH5Seurat("pbmc_multimodal_citeseq.h5seurat")

pdf("reference.celltype.pdf", width=12, height=6)
p1=DimPlot(object = reference, 
        reduction = "wnn.umap",
        group.by = "celltype.l1", 
        label = TRUE, 
        label.size = 3, 
        repel = TRUE) + NoLegend()
p2=DimPlot(object = reference, 
        reduction = "wnn.umap",
        group.by = "celltype.l2", 
        label = TRUE, 
        label.size = 3, 
        repel = TRUE) + NoLegend()
p1+p2
dev.off()

pbmc$seurat.celltype=pbmc$celltype

# find anchors between reference and query
# using a precomputed supervised PCA (spca) transformation for this example. 
DefaultAssay(pbmc) <- "SCT"
anchors <- FindTransferAnchors(
  reference = reference,
  query = pbmc,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# transfer cell type labels and protein data from the reference to the query
# project the query data onto the umap structure of the reference
pbmc <- MapQuery(
  anchorset = anchors,
  query = pbmc,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

pdf("pbmc_predicted_cell_typel1.pdf", width=18, height=6)
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("pbmc_predicted_cell_typel2.pdf", width=18, height=6)
p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# DC
VlnPlot(pbmc, features = c("rna_CST3"), slot = "counts", log = TRUE, group.by = 'predicted.celltype.l2')
# NK 
VlnPlot(pbmc, features = c("rna_GNLY"), slot = "counts", log = TRUE, group.by = 'predicted.celltype.l2')
# B cells
VlnPlot(pbmc, features = c("rna_MS4A1"), slot = "counts", log = TRUE, group.by = 'predicted.celltype.l2')
# CD14+ mono
VlnPlot(pbmc, features = c("rna_CD14"), slot = "counts", log = TRUE, group.by = 'predicted.celltype.l2')
# CD8 Naive, CD8 TEM_1, CD8 TEM_2
VlnPlot(pbmc, features = c("rna_CD8A"), slot = "counts", log = TRUE, group.by = 'predicted.celltype.l2')

pbmc$citeseq.celltype=factor(pbmc$predicted.celltype.l2)

p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat annotation")
p2<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "citeseq.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p1+p2

set.seed(1234)
pbmc=pbmc[,!pbmc$citeseq.celltype%in%c('dnT','Platelet','cDC1','ILC','CD4 Proliferating','ASDC')]
pbmc$citeseq.celltype=droplevels(pbmc$citeseq.celltype)
table(pbmc$citeseq.celltype)

pbmc=pbmc[,!pbmc$citeseq.celltype%in%c('CD4 CTL','CD8 TCM','NK Proliferating','NK_CD56bright')]
pbmc$citeseq.celltype=droplevels(pbmc$citeseq.celltype)
table(pbmc$citeseq.celltype)

table(pbmc$seurat.celltype, pbmc$citeseq.celltype)
levels(pbmc$seurat.celltype)=c(levels(pbmc$seurat.celltype), levels(pbmc$citeseq.celltype))
levels(pbmc$citeseq.celltype)=c(levels(pbmc$citeseq.celltype), levels(pbmc$seurat.celltype))


# CD27 is a marker for Memory B cells
VlnPlot(pbmc, features = c("rna_CD27"), slot = "counts", log = TRUE, group.by = 'seurat.celltype')
VlnPlot(pbmc, features = c("rna_CD27"), slot = "counts", log = TRUE, group.by = 'citeseq.celltype')
# Seurat mislabeled the naive and memory B cells
pbmc$seurat.celltype[pbmc$seurat.celltype=='Naive B']='B memory'
pbmc$seurat.celltype[pbmc$seurat.celltype=='Memory B']='B naive'
VlnPlot(pbmc, features = c("rna_CD27"), slot = "counts", log = TRUE, group.by = 'seurat.celltype')
VlnPlot(pbmc, features = c("rna_CD27"), slot = "counts", log = TRUE, group.by = 'citeseq.celltype')


pbmc$citeseq.celltype[pbmc$citeseq.celltype=='cDC2']='cDC'
pbmc$citeseq.celltype[pbmc$citeseq.celltype=='Plasmablast']='Plasma'
pbmc$seurat.celltype[pbmc$seurat.celltype=='Intermediate B']='B intermediate'
pbmc$seurat.celltype[pbmc$seurat.celltype=='CD8 TEM_1']='CD8 TEM' # combine Seurat CD8 TEM_1 and TEM_2
pbmc$seurat.celltype[pbmc$seurat.celltype=='CD8 TEM_2']='CD8 TEM'

pbmc$seurat.celltype=droplevels(pbmc$seurat.celltype)
pbmc$citeseq.celltype=droplevels(pbmc$citeseq.celltype)

levels(pbmc$seurat.celltype)%in%levels(pbmc$citeseq.celltype)

pbmc=pbmc[,as.character(pbmc$seurat.celltype)==as.character(pbmc$citeseq.celltype)]
pbmc$citeseq.celltype=factor(pbmc$citeseq.celltype, levels=levels(pbmc$seurat.celltype))

pbmc=pbmc[,!pbmc$seurat.celltype%in%c('Plasma','CD4 TEM','HSPC','gdT')]
pbmc$seurat.celltype=droplevels(pbmc$seurat.celltype)
pbmc$citeseq.celltype=droplevels(pbmc$citeseq.celltype)

table(pbmc$citeseq.celltype)
table(pbmc$seurat.celltype)



# Renormalize
pbmc <- subset(
  x = pbmc,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat annotation")
p2<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "citeseq.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p3<-DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat WNN clustering")
p1+p2+p3

table(pbmc$seurat.celltype, pbmc$seurat_clusters)


pbmc=pbmc[, (pbmc$seurat.celltype=='MAIT' & pbmc$seurat_clusters==15) |
            (pbmc$seurat.celltype=='CD8 Naive' & pbmc$seurat_clusters==2) |
            (pbmc$seurat.celltype=='CD4 TCM' & (pbmc$seurat_clusters==3 | pbmc$seurat_clusters==4)) |
            (pbmc$seurat.celltype=='Treg' & pbmc$seurat_clusters==13)   |
            (pbmc$seurat.celltype=='CD4 Naive' & pbmc$seurat_clusters==1)  |
            (pbmc$seurat.celltype=='NK' & pbmc$seurat_clusters==6)  |
            (pbmc$seurat.celltype=='CD14 Mono' & pbmc$seurat_clusters==0)  |
            (pbmc$seurat.celltype=='CD16 Mono' & pbmc$seurat_clusters==5)  |
            (pbmc$seurat.celltype=='pDC' & pbmc$seurat_clusters==16)  |
            (pbmc$seurat.celltype=='cDC' & pbmc$seurat_clusters==12) |
            (pbmc$seurat.celltype=='B intermediate' & pbmc$seurat_clusters==9) |
            (pbmc$seurat.celltype=='B memory' & pbmc$seurat_clusters==14) |
            (pbmc$seurat.celltype=='B naive' & pbmc$seurat_clusters==7) |
            (pbmc$seurat.celltype=='CD8 TEM' & (pbmc$seurat_clusters==8 | pbmc$seurat_clusters==10))]



# RNA analysis
DefaultAssay(pbmc) <- "RNA"
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

pbmc <- FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat annotation")
p2<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "citeseq.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p3<-DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat WNN clustering")
p1+p2+p3


p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "seurat.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat annotation")
p2<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "citeseq.celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p3<-DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat WNN clustering")
p1+p2+p3


# We will now get the gene activity matrix
# https://satijalab.org/signac/articles/pbmc_vignette.html
gene.activities <- GeneActivity(pbmc)
pbmc[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

# Perform normalization/scaling of the gene activity matrix
DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$ACTIVITY@scale.data[1:5,1:5]
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

save(pbmc, file='pbmc_chromvar_annotated.rda')




