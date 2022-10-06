library(Seurat)
library(SeuratDisk)
library(Signac)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(pheatmap)

# You can download the required data by running the following lines in a shell
# wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5
# wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz
# wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz.tbi

# The 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
pbmc[["ATAC"]] <- chrom_assay

VlnPlot(pbmc, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

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
pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% 
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

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

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(pbmc) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(pbmc), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)

pbmc <- SetAssayData(pbmc, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
# Note that this step can take 60 minutes 
pbmc <- RunChromVAR(
  object = pbmc,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat
reference <- LoadH5Seurat("pbmc_multimodal.h5seurat")

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

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.celltype.l1", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(pbmc, reduction = "wnn.umap", group.by = "predicted.celltype.l2", 
              label = TRUE, label.size = 3, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

pbmc$celltype=pbmc$predicted.celltype.l2
table(pbmc$celltype)

# Remove some minor cell types (with less than 100 cells)
pbmc=pbmc[,!pbmc$celltype%in%names(which(table(pbmc$celltype)<100))]
pbmc$celltype=droplevels(pbmc$celltype)
table(pbmc$celltype)

pbmc$celltype[pbmc$celltype=='cDC2']='cDC'


# Look at some cell-type markers
# DC
VlnPlot(pbmc, features = c("rna_CST3"), slot = "counts", log = TRUE, group.by = 'celltype')
# NK 
VlnPlot(pbmc, features = c("rna_GNLY"), slot = "counts", log = TRUE, group.by = 'celltype')
# B cells
VlnPlot(pbmc, features = c("rna_MS4A1"), slot = "counts", log = TRUE, group.by = 'celltype')
# CD14+ mono
VlnPlot(pbmc, features = c("rna_CD14"), slot = "counts", log = TRUE, group.by = 'celltype')
# CD8 Naive, CD8 TEM
VlnPlot(pbmc, features = c("rna_CD8A"), slot = "counts", log = TRUE, group.by = 'celltype')
# CD27 is a marker for Memory B cells (not Naive)
VlnPlot(pbmc, features = c("rna_CD27"), slot = "counts", log = TRUE, group.by = 'celltype')

# Renormalize after QC
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

p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p2<-DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat WNN clustering")
p1+p2

table(pbmc$seurat_clusters, pbmc$celltype)

# For each cluster, we only keep the main cell type
cell.keep=rep(FALSE, ncol(pbmc))
temp=table(pbmc$seurat_clusters, pbmc$celltype)

for(i in 1:nrow(temp)){
  clusteri=as.numeric(rownames(temp)[i])
  celltypei=colnames(temp)[which.max(temp[i,])]
  cell.keep[which(pbmc$seurat_clusters==clusteri & pbmc$celltype==celltypei)]=TRUE
}
pbmc=pbmc[,cell.keep]

# Renormalize after QC
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

p1<-DimPlot(pbmc, reduction = "wnn.umap", group.by = "celltype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("CITE-seq alignment")
p2<-DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("Seurat WNN clustering")
p1+p2

saveRDS(pbmc, file='pbmc.chromvar.annotated.rds')
