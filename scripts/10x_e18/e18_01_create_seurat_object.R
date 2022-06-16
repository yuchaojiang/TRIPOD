# load packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(patchwork)
library(biovizBase)

# set/create directories
dir.in <- "source_data"
if (!dir.exists(dir.in)) dir.create(data.in)
dir.out <- "derived_data"
if (!dir.exists(dir.out)) dir.create(data.out)
dir.fig <- "figures"
if (!dir.exists(dir.fig)) dir.create(data.fig)
dir.r <- "functions"

# download data into the data directory
# https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets
# e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5
# e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz
# e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz.tbi

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "e18_mouse_brain_fresh_5k_filtered_feature_bc_matrix.h5"
path <- file.path(dir.in, file)
inputdata.10x <- Read10X_h5(path)

# extract RNA and ATAC data
rna.counts <- inputdata.10x$`Gene Expression`
atac.counts <- inputdata.10x$Peaks

# create a Seurat object
e18 <- CreateSeuratObject(counts = rna.counts)
e18[["percent.mt"]] <- PercentageFeatureSet(e18, pattern = "^mt-")

# add in the ATAC-seq data
grange.counts <- StringToGRanges(rownames(atac.counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac.counts <- atac.counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
file <- "e18_mouse_brain_fresh_5k_atac_fragments.tsv.gz"
path <- file.path(dir.in, file)
chrom.assay <- CreateChromatinAssay(
  counts = atac.counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = path,
  min.cells = 10,
  annotation = annotations
)
e18[["ATAC"]] <- chrom.assay

# load cell type labels
# data(label.savercat)
file <- "label.savercat.rda"
path <- file.path(dir.in, file)
load(path); rm(path)
e18$label.savercat <- label.savercat

# perform RNA analysis
DefaultAssay(e18) <- "RNA"
e18 <- SCTransform(e18, verbose = FALSE)
e18 <- RunPCA(e18)
e18	<- RunUMAP(e18, reduction = "pca", dims = 1:50,
	reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

# perform ATAC analysis
DefaultAssay(e18) <- "ATAC"
e18 <- RunTFIDF(e18)
e18 <- FindTopFeatures(e18, min.cutoff = "q0")
e18 <- RunSVD(e18)
# exclude the first dimension as this is typically correlated with sequencing depth
e18 <- RunUMAP(e18, reduction = "lsi", dims = 2:50,
	reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# calculate a WNN graph
e18 <- FindMultiModalNeighbors(e18, reduction.list = list("pca", "lsi"),
	dims.list = list(1:50, 2:50))
e18 <- RunUMAP(e18, nn.name = "weighted.nn",
	reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# generate umap plots for sanity check
file <- "umap_rna_atac_before_qc.pdf"
path <- file.path(dir.fig, file)
pdf(path, width = 12, height = 6); rm(path)
p1 <- DimPlot(e18, reduction = "umap.rna",  group.by = "label.savercat",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("RNA")
p2 <- DimPlot(e18, reduction = "umap.atac",  group.by = "label.savercat",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("ATAC")
p3 <- DimPlot(e18, reduction = "wnn.umap", group.by = "label.savercat",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

# filter data based on cell types
celltypes.to.keep <- c("Radial glia", "Neuroblast",
	"Forebrain GABAergic", "Cortical or hippocampal glutamatergic",
	"Glioblast", "Oligodendrocyte", "Cajal-Retzius")

keep <- which(e18$label.savercat %in% celltypes.to.keep)
e18 <- subset(e18, cells = keep)
e18$celltype <- e18$label.savercat
e18$label.savercat <- NULL

# rename cell types for simplicity
e18$celltype[e18$celltype == "Cortical or hippocampal glutamatergic"] <- "Glutamatergic"
e18$celltype[e18$celltype == "Forebrain GABAergic"] <- "GABAergic"

# convert the cell type data to factors
levels.celltype  <- c(
	"Radial glia", "Neuroblast",
	"GABAergic", "Glutamatergic",
	"Glioblast", "Oligodendrocyte",
	"Cajal-Retzius"
)
e18$celltype <- factor(e18$celltype, levels = levels.celltype)

# filter data based on sequencing depths
e18 <- subset(
  x = e18,
  subset = nCount_ATAC < 1e5 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 1e5 &
    nCount_RNA > 1e3 &
    percent.mt < 20
)

# generate umap plots for sanity check
file <- "umap_rna_atac_after_qc.pdf"
path <- file.path(dir.fig, file)
pdf(path, width = 12, height = 6); rm(path)
p1 <- DimPlot(e18, reduction = "umap.rna", group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("RNA")
p2 <- DimPlot(e18, reduction = "umap.atac", group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("ATAC")
p3 <- DimPlot(e18, reduction = "wnn.umap", group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

file <- "e18.qc.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)
