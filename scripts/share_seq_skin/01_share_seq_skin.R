
ct=read.table('celltype_skin.txt', header = T, sep='\t')
# first column: atac barcode
# second column: rna barcode
# third column: annotated cell type from the paper

rna.count=read.table('GSM4156608_skin.late.anagen.rna.counts.txt.gz', header = T, row.names = 1)
rna.count[1:5, 1:5]

# match the rna.count with the cell type labels
any(is.na(match(ct$rna.bc, colnames(rna.count))))
rna.count=rna.count[,match(ct$rna.bc, colnames(rna.count))]

# The barcode file is the same as the annotation ct; omitted
peak.barcodes <- scan('GSM4156597_skin.late.anagen.barcodes.txt.gz', what="")
all(ct$rna.bc==peak.barcodes)
rm(peak.barcodes)

# Read in the peak matrix: stored as Market matrix
# https://atlas.gs.washington.edu/mouse-atac/docs/
library(Matrix)
library(GenomicRanges)
peak.count= readMM("GSM4156597_skin.late.anagen.counts.txt.gz")
dim(peak.count)
peak.bed=read.delim('GSM4156597_skin.late.anagen.peaks.bed.gz', header=FALSE)
peak.granges=GRanges(seqnames=peak.bed$V1, ranges=IRanges(st=peak.bed$V2, end=peak.bed$V3))
rm(peak.bed)
rownames(peak.count)=paste(peak.granges)
colnames(peak.count)=ct$rna.bc # we will use the rna barcode from here and on

dim(ct) # barcode and cell type
dim(peak.count) # peak count
length(peak.granges) # GRanges for peaks
dim(rna.count) # rna count


library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
library(dplyr)
library(ggplot2)

# Create Seurat object
skin <- CreateSeuratObject(counts = rna.count)
skin[["percent.mt"]] <- PercentageFeatureSet(skin, pattern = "^MT-")
rm(rna.count)

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.use <- seqnames(peak.granges) %in% standardChromosomes(peak.granges)
peak.count <- peak.count[as.vector(grange.use), ]
peak.granges <- peak.granges[as.vector(grange.use)]
rm(grange.use)

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"

rm(peak.granges)
peak.count=peak.count[rowSums(peak.count)>10,]
# Leave the fragment file out for now.
# frag.file <- "GSM4156599_skin.atac.fragments.bed.gz"
chrom_assay <- CreateChromatinAssay(
  counts = peak.count,
  sep = c(":", "-"),
  genome = 'mm10',
  # fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
skin[["ATAC"]] <- chrom_assay
skin@assays$RNA
skin@assays$ATAC
rm(peak.count)
save(skin, file='skin_cluster.rda')


library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse skin
library(dplyr)
library(ggplot2)
load('skin_cluster.rda')

## perform basic qc based on the number of detected molecules for each modality
## as well as mitochondrial percentage
VlnPlot(skin, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()
skin <- subset(
  x = skin,
  subset = nCount_ATAC < 25000 &
    nCount_ATAC > 500 &
    nCount_RNA < 10000 &
    nCount_RNA > 500 &
    percent.mt < 20
)
VlnPlot(skin, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0) + NoLegend()

## perform pre-processing and dimensional reduction on both assays independently
## using standard approaches for rna and atac-seq data
DefaultAssay(skin) <- "RNA"
skin <- SCTransform(skin, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(skin) <- "ATAC"
skin <- RunTFIDF(skin)
skin <- FindTopFeatures(skin, min.cutoff = 'q0')
skin <- RunSVD(skin)
skin <- RunUMAP(skin, reduction = 'lsi', dims = 2:50,
                reduction.name = "umap.atac",
                reduction.key = "atacUMAP_")

## calculate a wnn graph, representing a weighted combination of rna and atac-seq modalities
## use this graph for umap visualization and clustering
skin <- FindMultiModalNeighbors(skin, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
skin <- RunUMAP(skin, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
table(skin$seurat_clusters)

p1 <- DimPlot(skin, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(skin, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(skin, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# We will now get the gene activity matrix
DefaultAssay(skin) <- "ATAC"
# gene.activities <- GeneActivity(skin) # This would not run because it requires fragment file to be loaded.

# Below is a short-cut but is not very precise: we only aggregate reads in the peak regions
# But there are reads off peak regions but within gene bodies that are defined
# Will load the fragment file in the second script.
annotation <- Annotation(object = skin)
transcripts <- Signac:::CollapseToLongestTranscript(ranges = annotation)
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"]
transcripts <- Extend(x = transcripts, upstream = 2000, 
                      downstream = 0)
gene.activities=matrix(nrow=length(transcripts), ncol=ncol(skin))
rownames(gene.activities)=transcripts$gene_name
colnames(gene.activities)=colnames(skin)
peak.count=skin@assays$ATAC@counts
peak.granges=skin@assays$ATAC@ranges
for(i in 1:nrow(gene.activities)){
  if(i %%200 ==0) cat(i,' ')
  peak.index=which(countOverlaps(peak.granges, transcripts[i])>0)
  gene.activities[i,]=apply(peak.count[peak.index,, drop=FALSE],2, sum)
}
gene.activities=gene.activities[apply(gene.activities,1,sum)>10,]
gene.activities=gene.activities[!duplicated(rownames(gene.activities)),] # remove duplicate genes
gene.activities=gene.activities[-which(rownames(gene.activities)==''),]

skin[['ACTIVITY']] <- CreateAssayObject(counts = gene.activities)

# Perform normalization/scaling of the gene activity matrix
DefaultAssay(skin) <- "ACTIVITY"
skin <- NormalizeData(skin)
skin <- FindVariableFeatures(skin)
all.genes <- rownames(skin)
skin <- ScaleData(skin, features = all.genes)
skin@assays$ACTIVITY@scale.data[1:5,1:5]
save(skin, file='skin.rda')


wnn.resolution=5
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p5=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
        label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p5

wnn.resolution=15
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p15=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
            label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p15

wnn.resolution=25
skin <- FindClusters(skin, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
p25=DimPlot(skin, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
            label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(skin@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()
p25


p5|p15|p25
