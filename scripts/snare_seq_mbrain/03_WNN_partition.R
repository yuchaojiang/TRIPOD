
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79) # For mouse snare
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

snare=readRDS("snare.chromvar.rds")

# SCtransform
DefaultAssay(snare) <- "RNA"
snare <- SCTransform(snare, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

####################################
##  Some general data structure
####################################

# This Seurat object has several assays: raw RNA, normalized/scaled RNA, ATAC peak,
# and ATAC gene activity
snare@assays$RNA # RNA 
snare@assays$SCT # RNA normalized by sctransform
snare@assays$ATAC # ATAC peak matrix
snare@assays$chromvar # Motif by cell matrix from chromVAR

snare@assays$RNA@counts[1:10, 1:10] # raw RNA count matrix
snare@assays$SCT@scale.data[1:10, 1:10] # RNA normalized by size factor, log, and scaled
snare@assays$ATAC@counts[1:10, 1:10] # raw ATAC peak matrix
snare@assays$chromvar@data[1:10,1:10] # chromVAR motif by cell matrix

# Dim. reduction plot
DimPlot(snare, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(snare, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(snare, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)

####################################
##  Get some additional annotations
####################################

DefaultAssay(snare) <- "ATAC"
transcripts <- Signac:::CollapseToLongestTranscript(ranges = Annotation(snare))
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"] # This is the gene coordinates
transcripts <- transcripts[seqnames(transcripts) %in% paste('chr', 1:19, sep='')]
transcripts <- sort(transcripts)
transcripts.gr <- transcripts; rm(transcripts)
transcripts.gr # GRanges obj for transcript/gene

peak.gr=snare@assays$ATAC@ranges # GRanges obj for ATAC peaks

motifxTF=unlist(snare@assays$ATAC@motifs@motif.names)
motifxTF=cbind(names(motifxTF), motifxTF)
colnames(motifxTF)=c('motif','TF')
head(motifxTF) # motif <-> TF/gene correspondence
peakxmotif=snare@assays$ATAC@motifs@data # peak by motif: overlap of ATAC peaks and motifs (binary)

# The mouse has gene symbol as Brca1 v.s. the gene/TF symbol as BRCA1 in human
motifxTF[,2]=stringr::str_to_title(tolower(motifxTF[,2]))

# Only keep the motifs/genes that we have in RNA
peakxmotif=peakxmotif[,motifxTF[,2]%in%rownames(snare@assays$RNA)]
motifxTF=motifxTF[motifxTF[,2]%in%rownames(snare@assays$RNA),] 

####################################
##  Generate some statistics/visualizations
####################################
ext.upstream=100000 # This is actually both up and downstream
transcripts <- transcripts.gr
TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
TSS <- GRanges(seqnames = seqnames(transcripts),
               ranges = IRanges(start = TSS.position, width = 1),
               strand = strand(transcripts),
               gene_name = transcripts$gene_name)
transcripts.ext.gr=getRegion(c(ext.upstream, ext.upstream), TSS)
transcripts.ext.gr # GRanges obj for extended gene bodies
transcripts.gr

# On average, how many peaks per extended gene region?
hist(countOverlaps(transcripts.ext.gr, peak.gr), xlab='Number of peaks', 
     main=paste('Number of peaks', ext.upstream/1000, 'Kb upstream of a gene'),
     xlim=c(0,100), breaks=seq(0,1000,2))

# On average, how many genes per extended gene region?
hist(countOverlaps(transcripts.ext.gr, transcripts.gr), xlab='Number of genes', 
     main=paste('Number of genes', ext.upstream/1000, 'Kb upstream of a gene'),
     breaks=seq(0,100,1), xlim=c(0,25))

# On average, how many motifs per peak region?
hist(apply(peakxmotif,1,sum), xlab='Number of motifs', 
     main=paste('Number of motifs per peak'),
     breaks=seq(0,1000,5), xlim=c(0,400))

# For peaks in the extended gene regions, how many fall into another gene?
hist(countOverlaps(subsetByOverlaps(peak.gr, transcripts.ext.gr), transcripts.gr),
     xlim=c(0,15), breaks=seq(-0.5, 100.5, 1),
     xlab='Number of other genes',
     main=paste('Number of additional genes peaks overlap in a', ext.upstream/1000, 'Kb window'))

# Take the intersection of the genes
genes=intersect(transcripts.gr$gene_name, rownames(snare@assays$SCT))
DefaultAssay(snare)='RNA'
snare@assays$RNA=subset(snare@assays$RNA,, match(genes, rownames(snare@assays$RNA)))
snare@assays$SCT=subset(snare@assays$SCT,, match(genes, rownames(snare@assays$SCT)))

transcripts.gr=transcripts.gr[match(genes,transcripts.gr$gene_name)]
transcripts.ext.gr=transcripts.ext.gr[match(genes,transcripts.ext.gr$gene_name)]

peakxmotif=peakxmotif[,motifxTF[,2]%in%genes]
motifxTF=motifxTF[motifxTF[,2]%in%genes,]
snare@assays$chromvar=subset(snare@assays$chromvar,, match(motifxTF[,1], rownames(snare@assays$chromvar)))

# Some re-normalization
# RNA analysis
DefaultAssay(snare) <- "RNA"
snare <- SCTransform(snare, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(snare) <- "ATAC"
snare <- RunTFIDF(snare)
snare <- FindTopFeatures(snare, min.cutoff = 'q0')
snare <- RunSVD(snare)
snare <- RunUMAP(snare, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

snare <- FindMultiModalNeighbors(snare, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
snare <- RunUMAP(snare, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
snare <- FindClusters(snare, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

snare$seurat.celltype=snare$celltype

DimPlot(snare, group.by='celltype', reduction = 'wnn.umap')
temp=table(snare$celltype, snare$seurat_clusters)
cell.filter=rep(FALSE, ncol(snare))

for(wnn in levels(snare$seurat_clusters)){
  wnn.celltype=names(which.max(temp[,wnn]))
  cell.filter[which(snare$seurat_clusters==wnn & snare$celltype==wnn.celltype)]=TRUE
}
snare=snare[,cell.filter]

# Some re-normalization
# RNA analysis
DefaultAssay(snare) <- "RNA"
snare <- SCTransform(snare, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(snare) <- "ATAC"
snare <- RunTFIDF(snare)
snare <- FindTopFeatures(snare, min.cutoff = 'q0')
snare <- RunSVD(snare)
snare <- RunUMAP(snare, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

snare <- FindMultiModalNeighbors(snare, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
snare <- RunUMAP(snare, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
snare <- FindClusters(snare, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
DimPlot(snare, group.by='celltype', reduction = 'wnn.umap')

save(snare, file='snare.chromvar.qc.rda')


####################################
##  Now let's get started: ATAC-RNA relationship
####################################
# Input 1: WNN resolution (i.e., cluster size)
wnn.resolution=10
snare <- FindClusters(snare, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
length(levels(snare@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])) # Number of clusters
table(snare@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])
table(snare$seurat_clusters, snare$celltype)

DimPlot(snare, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
        label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(snare@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()

wnn.numcells=table(snare$seurat_clusters) # Number of cells per wnn
names(wnn.numcells)=paste('wnn_', levels(snare$seurat_clusters),sep='')

wnn.rna=matrix(nrow=length(levels(snare$seurat_clusters)), ncol=nrow(snare@assays$RNA))
rownames(wnn.rna)=paste('wnn_', levels(snare$seurat_clusters),sep='')
colnames(wnn.rna)=rownames(snare@assays$RNA)
snare.rna.count=snare@assays$RNA@counts
for(i in 1:nrow(wnn.rna)){
  cat(i,' ')
  wnn.rna[i,]=apply(snare.rna.count[,snare$seurat_clusters==(i-1)], 1, sum)
  wnn.rna[i,]=wnn.rna[i,]/sum(wnn.rna[i,])*10^6 # adjust for lib. size for each wnn
}
rm(snare.rna.count); #rm(snare.rna.count.adj)

wnn.peak=matrix(nrow=length(levels(snare$seurat_clusters)), ncol=nrow(snare@assays$ATAC))
rownames(wnn.peak)=paste('wnn_', levels(snare$seurat_clusters),sep='')
colnames(wnn.peak)=rownames(snare@assays$ATAC)
snare.atac.count=snare@assays$ATAC@counts
for(i in 1:nrow(wnn.peak)){
  cat(i,' ')
  wnn.peak[i,]=apply(snare.atac.count[,snare$seurat_clusters==(i-1)], 1, sum)
  wnn.peak[i,]=wnn.peak[i,]/sum(wnn.peak[i,])*10^6 # adjust for lib. size for each wnn
}
rm(snare.atac.count); #rm(snare.atac.count.adj)

wnn.rna[1:5, 1:5] # wnn by gene matrix: sum of single-cell RNA read count
wnn.peak[1:5, 1:5] # wnn by peak matrix: sum of single-cell ATAC read count

# Get WNN cell types
temp=table(snare$seurat.celltype, snare$seurat_clusters)
wnn.celltype=rep(NA, nrow(wnn.rna))
for(i in 1:length(wnn.celltype)){
  temp.i_1=temp[,colnames(temp)==as.character(i-1)]
  wnn.celltype[i]=names(temp.i_1)[which.max(temp.i_1)]
}

# Get the corresponding color for each cell type from Seurat
p <- Seurat::DimPlot(snare, reduction = 'wnn.umap', label=T, group.by='seurat.celltype')
pbuild <- ggplot2::ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # this is to get the color palette by Seurat
pdata =cbind(snare$seurat.celltype, pdata)
wnn.celltype.col=rep(NA, length(wnn.celltype))
for(i in 1:length(wnn.celltype)){
  wnn.celltype.col[i]=pdata$colour[min(which(pdata$`snare$seurat.celltype`==wnn.celltype[i]))]
}

names(wnn.celltype) <- levels(snare)
snare <- RenameIdents(snare, wnn.celltype)
DimPlot(snare, reduction = "wnn.umap", label = TRUE, pt.size = 0.5)

# Remove unnecessary stored values/objects
rm(p); rm(pbuild); rm(pdata)
rm(temp)

save.image(file='ATAC_RNA_WNN.rda')
