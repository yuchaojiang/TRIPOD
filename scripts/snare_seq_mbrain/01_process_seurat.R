
library(Signac)
library(Seurat)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
set.seed(1234)

# load processed data matrices for each assay
rna <- Read10X("./SNARE_seq/rna/", gene.column=1)
atac <- Read10X("./SNARE_seq/atac/", gene.column = 1)
fragments <- "./SNARE_seq/fragments.sort.bed.gz"

# create a Seurat object and add the assays
snare <- CreateSeuratObject(counts = rna)
snare[['ATAC']] <- CreateChromatinAssay(
  counts = atac,
  sep = c(":", "-"),
  genome = "mm10",
  fragments = fragments
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

# add the gene information to the object
Annotation(snare[["ATAC"]]) <- annotations

# Quality control
DefaultAssay(snare) <- "ATAC"
snare <- TSSEnrichment(snare)
snare <- NucleosomeSignal(snare)
snare$blacklist_fraction <- FractionCountsInRegion(
  object = snare,
  assay = 'ATAC',
  regions = blacklist_mm10
)
Idents(snare) <- "all"  # group all cells together, rather than by replicate
VlnPlot(
  snare,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment",
               "nucleosome_signal", "blacklist_fraction"),
  pt.size = 0.1,
  ncol = 5
)

snare <- subset(
  x = snare,
  subset = blacklist_fraction < 0.03 &
    TSS.enrichment < 20 &
    nCount_RNA > 800 &
    nCount_ATAC > 500
)
snare

# Gene expression data processing
DefaultAssay(snare) <- "RNA"

snare <- FindVariableFeatures(snare, nfeatures = 3000)
snare <- NormalizeData(snare)
snare <- ScaleData(snare)
snare <- RunPCA(snare, npcs = 30)
snare <- RunUMAP(snare, dims = 1:30, reduction.name = "umap.rna")
snare <- FindNeighbors(snare, dims = 1:30)
snare <- FindClusters(snare, resolution = 0.5, algorithm = 3)

p1 <- DimPlot(snare, label = TRUE) + NoLegend() + ggtitle("RNA UMAP")

# ATAC data processing
DefaultAssay(snare) <- 'ATAC'

snare <- FindTopFeatures(snare, min.cutoff = 10)
snare <- RunTFIDF(snare)
snare <- RunSVD(snare)
snare <- RunUMAP(snare, reduction = 'lsi', dims = 2:30, reduction.name = 'umap.atac')
p2 <- DimPlot(snare, reduction = 'umap.atac', label = TRUE) + NoLegend() + ggtitle("ATAC UMAP")

p1+p2

# Integrate with Allen scRNA-seq data for label transfer
# label transfer from Allen brain
allen <- readRDS("./allen_brain.rds")

# use the RNA assay in the SNARE-seq data for integration with scRNA-seq
DefaultAssay(snare) <- 'RNA'

transfer.anchors <- FindTransferAnchors(
  reference = allen,
  query = snare,
  dims = 1:30,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = allen$subclass,
  weight.reduction = snare[['pca']],
  dims = 1:30
)

snare <- AddMetaData(object = snare, metadata = predicted.labels)


# label clusters based on predicted ID
new.cluster.ids <- c(
  "L2/3 IT",
  "L4",
  "L6 IT",
  "L5 CT",
  "L4",
  "L5 PT",
  "Pvalb",
  "Sst",
  "Astro",
  "Oligo",
  "Vip/Lamp5",
  "L6 IT.2",
  "L6b",
  "NP"
)
names(x = new.cluster.ids) <- levels(x = snare)
snare <- RenameIdents(object = snare, new.cluster.ids)
snare$celltype <- Idents(snare)
DimPlot(snare, group.by = 'celltype', label = TRUE, reduction = 'umap.rna')

# jointly visualizing RNA and ATAC
DefaultAssay(snare) <- "ATAC"
CoveragePlot(snare, region = "chr2-22620000-22660000", features = "Gad2")


## calculate a wnn graph, representing a weighted combination of rna and atac-seq modalities
## use this graph for umap visualization and clustering
snare <- FindMultiModalNeighbors(snare, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
snare <- RunUMAP(snare, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
wnn.resolution=10
snare <- FindClusters(snare, resolution=wnn.resolution,graph.name = "wsnn", algorithm = 3, verbose = FALSE)
table(snare$seurat_clusters)

p1 <- DimPlot(snare, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(snare, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(snare, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))


saveRDS(snare, file='snare.rds')
