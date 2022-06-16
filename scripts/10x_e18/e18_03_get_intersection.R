# load packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(patchwork)

# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "e18.chromvar.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

# get objects required for model fitting
DefaultAssay(e18) <- "ATAC"
transcripts.gr <- Signac:::CollapseToLongestTranscript(ranges = Annotation(e18))
transcripts.gr <- transcripts.gr[transcripts.gr$gene_biotype == "protein_coding"]
transcripts.gr <- transcripts.gr[seqnames(transcripts.gr) %in%
		paste0("chr", 1:19)]
transcripts.gr <- sort(transcripts.gr)

peaks.gr <- e18@assays$ATAC@ranges

motifxTF <- unlist(e18@assays$ATAC@motifs@motif.names)
motifxTF <- cbind(names(motifxTF), motifxTF)
colnames(motifxTF) <- c("motif", "TF")

peakxmotif <- e18@assays$ATAC@motifs@data

motifxTF[, 2] <- stringr::str_to_title(tolower(motifxTF[, 2]))

peakxmotif <- peakxmotif[, motifxTF[, 2] %in% rownames(e18@assays$RNA)]
motifxTF <- motifxTF[motifxTF[, 2] %in% rownames(e18@assays$RNA), ]

genes <- intersect(transcripts.gr$gene_name, rownames(e18@assays$SCT))
DefaultAssay(e18) <- "RNA"
e18@assays$RNA <- subset(e18@assays$RNA,
	features = match(genes, rownames(e18@assays$RNA)))
e18@assays$SCT <- subset(e18@assays$SCT,
	features = match(genes, rownames(e18@assays$SCT)))
transcripts.gr <- transcripts.gr[match(genes, transcripts.gr$gene_name)]
peakxmotif <- peakxmotif[, motifxTF[, 2] %in% genes]
motifxTF <- motifxTF[motifxTF[, 2] %in% genes, ]
e18@assays$chromvar <- subset(
	e18@assays$chromvar,
	features = match(motifxTF[, 1], rownames(e18@assays$chromvar)))

file <- "transcripts.gr.rds"
path <- file.path(dir.out, file)
saveRDS(transcripts.gr, path); rm(path)

file <- "peaks.gr.rds"
path <- file.path(dir.out, file)
saveRDS(peaks.gr, path); rm(path)

file <- "motifxTF.rds"
path <- file.path(dir.out, file)
saveRDS(motifxTF, path); rm(path)

file <- "peakxmotif.rds"
path <- file.path(dir.out, file)
saveRDS(peakxmotif, path); rm(path)

# perform RNA analysis again
DefaultAssay(e18) <- "RNA"
e18 <- SCTransform(e18, verbose = FALSE)
e18 <- RunPCA(e18)
e18	<- RunUMAP(e18, reduction = "pca", dims = 1:50,
	reduction.name = "umap.rna", reduction.key = "rnaUMAP_")

# perform ATAC analysis again
DefaultAssay(e18) <- "ATAC"
e18 <- RunTFIDF(e18)
e18 <- FindTopFeatures(e18, min.cutoff = "q0")
e18 <- RunSVD(e18)
e18 <- RunUMAP(e18, reduction = "lsi", dims = 2:50,
	reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# recalculate a WNN graph
e18 <- FindMultiModalNeighbors(e18, reduction.list = list("pca", "lsi"),
	dims.list = list(1:50, 2:50))
e18 <- RunUMAP(e18, nn.name = "weighted.nn",
	reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# generate umap plots for fig 1a and supplementary figure 1x
file <- "umap_rna_atac_after_re_norm.pdf"
path <- file.path(dir.fig, file)
pdf(path, width = 12, height = 6); rm(path)
p1 <- DimPlot(e18, reduction = "umap.rna",  group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("RNA")
p2 <- DimPlot(e18, reduction = "umap.atac",  group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("ATAC")
p3 <- DimPlot(e18, reduction = "wnn.umap", group.by = "celltype",
	label = TRUE, label.size = 2.5, repel = TRUE) +
	ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
dev.off()

file <- "e18.intersect.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)
