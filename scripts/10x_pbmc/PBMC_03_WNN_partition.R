
library(Seurat)
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
library(DescTools)

load('pbmc_chromvar_annotated.rda')

####################################
##  Some general data structure
####################################

# This Seurat object has several assays: raw RNA, normalized/scaled RNA, ATAC peak,
# and ATAC gene activity
pbmc@assays$RNA # RNA 
pbmc@assays$SCT # RNA normalized by sctransform
pbmc@assays$ATAC # ATAC peak matrix
pbmc@assays$ACTIVITY # Gene activity from ATAC
pbmc@assays$chromvar # Motif by cell matrix from chromVAR

pbmc@assays$RNA@counts[1:10, 1:10] # raw RNA count matrix
pbmc@assays$SCT@scale.data[1:10, 1:10] # RNA normalized by size factor, log, and scaled
pbmc@assays$ATAC@counts[1:10, 1:10] # raw ATAC peak matrix
pbmc@assays$ACTIVITY@counts[1:10, 1:10] # ATAC gene activity matrix
pbmc@assays$ACTIVITY@scale.data[1:10, 1:10] # normalized/scaled ATAC activity matrix
pbmc@assays$chromvar@data[1:10,1:10] # chromVAR motif by cell matrix

# Dim. reduction plot
DimPlot(pbmc, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(pbmc, reduction = "umap.rna",  label = TRUE, label.size = 2.5, repel = TRUE)
DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE)

####################################
##  Get some additional annotations
####################################

DefaultAssay(pbmc) <- "ATAC"
transcripts <- Signac:::CollapseToLongestTranscript(ranges = Annotation(pbmc))
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"] # This is the gene coordinates
transcripts <- transcripts[seqnames(transcripts) %in% paste('chr', 1:22, sep='')]
transcripts <- sort(transcripts)
transcripts.gr <- transcripts; rm(transcripts)
transcripts.gr # GRanges obj for transcript/gene

peak.gr=pbmc@assays$ATAC@ranges # GRanges obj for ATAC peaks

motifxTF=unlist(pbmc@assays$ATAC@motifs@motif.names)
motifxTF=cbind(names(motifxTF), motifxTF)
colnames(motifxTF)=c('motif','TF')
head(motifxTF) # motif <-> TF/gene correspondence

peakxmotif=pbmc@assays$ATAC@motifs@data # peak by motif: overlap of ATAC peaks and motifs (binary)
# Only keep the motifs/genes that we have in RNA
peakxmotif=peakxmotif[,motifxTF[,2]%in%rownames(pbmc@assays$RNA)]
motifxTF=motifxTF[motifxTF[,2]%in%rownames(pbmc@assays$RNA),] 


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
genes=intersect(intersect(transcripts.gr$gene_name, rownames(pbmc@assays$SCT)),
                rownames(pbmc@assays$ACTIVITY))
DefaultAssay(pbmc)='RNA'
pbmc@assays$RNA=subset(pbmc@assays$RNA,, match(genes, rownames(pbmc@assays$RNA)))
pbmc@assays$SCT=subset(pbmc@assays$SCT,, match(genes, rownames(pbmc@assays$SCT)))
pbmc@assays$ACTIVITY=subset(pbmc@assays$ACTIVITY,, match(genes, rownames(pbmc@assays$ACTIVITY)))

transcripts.gr=transcripts.gr[match(genes,transcripts.gr$gene_name)]
transcripts.ext.gr=transcripts.ext.gr[match(genes,transcripts.ext.gr$gene_name)]

peakxmotif=peakxmotif[,motifxTF[,2]%in%genes]
motifxTF=motifxTF[motifxTF[,2]%in%genes,]
pbmc@assays$chromvar=subset(pbmc@assays$chromvar,, match(motifxTF[,1], rownames(pbmc@assays$chromvar)))

# Some re-normalization
# Perform normalization/scaling of the gene activity matrix
DefaultAssay(pbmc) <- "ACTIVITY"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$ACTIVITY@scale.data[1:5,1:5]
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- RunUMAP(pbmc, dims = 1:10)

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

pbmc@assays$ACTIVITY=NULL
pbmc@assays$prediction.score.celltype.l1=NULL
pbmc@assays$prediction.score.celltype.l2=NULL

levels(pbmc$citeseq.celltype)=c(levels(pbmc$citeseq.celltype),'CD4 Memory','CD8 Memory')
pbmc$citeseq.celltype[pbmc$citeseq.celltype=='CD4 TCM']='CD4 Memory'
pbmc$citeseq.celltype[pbmc$citeseq.celltype=='CD8 TEM']='CD8 Memory'
pbmc$citeseq.celltype=droplevels(pbmc$citeseq.celltype)
pbmc$citeseq.celltype=factor(pbmc$citeseq.celltype, levels=c('pDC','cDC','CD14 Mono','CD16 Mono','B naive','B intermediate','B memory',
                                                             'Treg','MAIT','NK','CD4 Naive','CD4 Memory','CD8 Naive','CD8 Memory'))
pbmc$seurat.celltype=pbmc$citeseq.celltype
pbmc$celltype=pbmc$citeseq.celltype
Idents(pbmc)='celltype'
table(pbmc$celltype)


####################################
##  Now let's get started: ATAC-RNA relationship
####################################
# Input 1: WNN resolution (i.e., cluster size)
wnn.resolution=15
pbmc <- FindClusters(pbmc, graph.name = "wsnn", resolution = wnn.resolution, algorithm = 3, verbose = FALSE)
length(levels(pbmc@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])) # Number of clusters
table(pbmc@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])

DimPlot(pbmc, reduction = "wnn.umap", group.by = paste('wsnn_res.',wnn.resolution,sep=''), 
              label = TRUE, label.size = 3, repel = TRUE) + 
  ggtitle(paste("WNN resolution", wnn.resolution, 'with', 
                length(levels(pbmc@meta.data[,paste('wsnn_res.',wnn.resolution,sep='')])), 'clusters'))+ NoLegend()

wnn.numcells=table(pbmc$seurat_clusters) # Number of cells per wnn
names(wnn.numcells)=paste('wnn_', levels(pbmc$seurat_clusters),sep='')

wnn.rna=matrix(nrow=length(levels(pbmc$seurat_clusters)), ncol=nrow(pbmc@assays$RNA))
rownames(wnn.rna)=paste('wnn_', levels(pbmc$seurat_clusters),sep='')
colnames(wnn.rna)=rownames(pbmc@assays$RNA)
pbmc.rna.count=pbmc@assays$RNA@counts
pbmc.rna.count.adj=sweep(pbmc.rna.count,2,colSums(pbmc.rna.count)/10^6,`/`) # adjust single-cell library size before computing the sd and mad
for(i in 1:nrow(wnn.rna)){
  cat(i,' ')
  wnn.rna[i,]=apply(pbmc.rna.count[,pbmc$seurat_clusters==(i-1)], 1, sum)
  wnn.rna[i,]=wnn.rna[i,]/sum(wnn.rna[i,])*10^6 # adjust for lib. size for each wnn
}
rm(pbmc.rna.count); rm(pbmc.rna.count.adj)

wnn.peak=matrix(nrow=length(levels(pbmc$seurat_clusters)), ncol=nrow(pbmc@assays$ATAC))
rownames(wnn.peak)=paste('wnn_', levels(pbmc$seurat_clusters),sep='')
colnames(wnn.peak)=rownames(pbmc@assays$ATAC)
pbmc.atac.count=pbmc@assays$ATAC@counts
pbmc.atac.count.adj=sweep(pbmc.atac.count,2,colSums(pbmc.atac.count)/10^6,`/`) # adjust single-cell library size before computing the sd and mad
for(i in 1:nrow(wnn.peak)){
  cat(i,' ')
  wnn.peak[i,]=apply(pbmc.atac.count[,pbmc$seurat_clusters==(i-1)], 1, sum)
  wnn.peak[i,]=wnn.peak[i,]/sum(wnn.peak[i,])*10^6 # adjust for lib. size for each wnn
}
rm(pbmc.atac.count); rm(pbmc.atac.count.adj)

wnn.rna[1:5, 1:5] # wnn by gene matrix: sum of single-cell RNA read count
wnn.peak[1:5, 1:5] # wnn by peak matrix: sum of single-cell ATAC read count

# Get WNN cell types
temp=table(pbmc$seurat.celltype, pbmc$seurat_clusters)
wnn.celltype=rep(NA, nrow(wnn.rna))
for(i in 1:length(wnn.celltype)){
  temp.i_1=temp[,colnames(temp)==as.character(i-1)]
  wnn.celltype[i]=names(temp.i_1)[which.max(temp.i_1)]
}

# Get the corresponding color for each cell type from Seurat
p <- Seurat::DimPlot(pbmc, reduction = 'wnn.umap', label=T, group.by='seurat.celltype')
pbuild <- ggplot2::ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # this is to get the color palette by Seurat
pdata =cbind(pbmc$seurat.celltype, pdata)
wnn.celltype.col=rep(NA, length(wnn.celltype))
for(i in 1:length(wnn.celltype)){
  wnn.celltype.col[i]=pdata$colour[min(which(pdata$`pbmc$seurat.celltype`==wnn.celltype[i]))]
}

names(wnn.celltype) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, wnn.celltype)
DimPlot(pbmc, reduction = "wnn.umap", label = TRUE, pt.size = 0.5)

# Remove unnecessary stored values/objects
rm(p); rm(pbuild); rm(pdata)
rm(temp)
rm(all.genes); rm(genes); rm(ext.upstream); rm(i); rm(temp.i_1)
save.image(file='ATAC_RNA_WNN.rda')

