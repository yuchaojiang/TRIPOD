
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
library(RColorBrewer)
library(plyr)

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

load('skin.subset.qc.branch.rda')
skin$celltype=as.factor(revalue(skin$celltype, c("Hair Shaft-cuticle.cortex"="Hair_Shaft")))
skin$seurat.celltype=skin$celltype
table(skin$celltype)

# Highly variable genes
celltype_presto <- presto:::wilcoxauc.Seurat(X = skin, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
celltype_genes=celltype_presto[order(celltype_presto[,'pval'])[1:1000],'feature'] # cell-type markers
highvar_genes=skin@assays$SCT@var.features[1:1000]
genes=unique(c(celltype_genes, highvar_genes)) # vector of highly variable genes
which(genes == 'Wnt3')

# Feature plots for some reported markers...
gene.name='Wnt3'
FeaturePlot(skin, features = paste("sct_",gene.name,sep=''), reduction = 'pseudo', cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
DotPlot(skin, features = paste("rna_",gene.name, sep=''), group.by = 'seurat.celltype') + ggtitle('Gene Expression')+ylab('Expression')
DotPlot(skin, features = paste("rna_",gene.name, sep=''), group.by = 'ptime.clusters') + ggtitle('Gene Expression')+ylab('Expression')
DotPlot(skin, features = paste("rna_",gene.name, sep=''), group.by = 'branch') + ggtitle('Gene Expression')+ylab('Expression')

####################################
##  Some general data structure
####################################

# This Seurat object has several assays: raw RNA, normalized/scaled RNA, ATAC peak,
# and ATAC gene activity
skin@assays$RNA # RNA 
skin@assays$SCT # RNA normalized by sctransform
skin@assays$ATAC # ATAC peak matrix
skin@assays$chromvar # Motif by cell matrix from chromVAR

skin@assays$RNA@counts[1:10, 1:10] # raw RNA count matrix
skin@assays$SCT@scale.data[1:10, 1:10] # RNA normalized by size factor, log, and scaled
skin@assays$ATAC@counts[1:10, 1:10] # raw ATAC peak matrix
skin@assays$chromvar@data[1:10,1:10] # chromVAR motif by cell matrix

####################################
##  Get some additional annotations
####################################

DefaultAssay(skin) <- "ATAC"
transcripts <- Signac:::CollapseToLongestTranscript(ranges = Annotation(skin))
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"] # This is the gene coordinates
transcripts <- transcripts[seqnames(transcripts) %in% paste('chr', 1:19, sep='')]
transcripts <- sort(transcripts)
transcripts.gr <- transcripts; rm(transcripts)
transcripts.gr # GRanges obj for transcript/gene

peak.gr=skin@assays$ATAC@ranges # GRanges obj for ATAC peaks

motifxTF=unlist(skin@assays$ATAC@motifs@motif.names)
motifxTF=cbind(names(motifxTF), motifxTF)
colnames(motifxTF)=c('motif','TF')
head(motifxTF) # motif <-> TF/gene correspondence
peakxmotif=skin@assays$ATAC@motifs@data # peak by motif: overlap of ATAC peaks and motifs (binary)

# The mouse has gene symbol as Brca1 v.s. the gene/TF symbol as BRCA1 in human
motifxTF[,2]=stringr::str_to_title(tolower(motifxTF[,2]))

# Only keep the motifs/genes that we have in RNA
peakxmotif=peakxmotif[,motifxTF[,2]%in%rownames(skin@assays$RNA)]
motifxTF=motifxTF[motifxTF[,2]%in%rownames(skin@assays$RNA),] 


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
genes=intersect(transcripts.gr$gene_name, rownames(skin@assays$SCT))
DefaultAssay(skin)='RNA'
skin@assays$RNA=subset(skin@assays$RNA,, match(genes, rownames(skin@assays$RNA)))
skin@assays$SCT=subset(skin@assays$SCT,, match(genes, rownames(skin@assays$SCT)))

transcripts.gr=transcripts.gr[match(genes,transcripts.gr$gene_name)]
transcripts.ext.gr=transcripts.ext.gr[match(genes,transcripts.ext.gr$gene_name)]

peakxmotif=peakxmotif[,motifxTF[,2]%in%genes]
motifxTF=motifxTF[motifxTF[,2]%in%genes,]
skin@assays$chromvar=subset(skin@assays$chromvar,, match(motifxTF[,1], rownames(skin@assays$chromvar)))

skin@assays$ATAC; length(peak.gr)
skin@assays$RNA; length(transcripts.gr); length(transcripts.ext.gr)
skin@assays$chromvar; dim(motifxTF); dim(peakxmotif)
nrow(skin@assays$RNA)
nrow(skin@assays$SCT)

# Some re-normalization
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
skin <- FindClusters(skin, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

# Generate WNN input for skin hair data only.
####################################
##  Now let's get started: ATAC-RNA relationship
####################################
skin$seurat_clusters=as.factor(as.numeric(skin$ptime.clusters)-1)
wnn.numcells=table(skin$seurat_clusters) # Number of cells per wnn
names(wnn.numcells)=paste('wnn_', levels(skin$seurat_clusters),sep='')

wnn.rna=matrix(nrow=length(levels(skin$seurat_clusters)), ncol=nrow(skin@assays$RNA))
rownames(wnn.rna)=paste('wnn_', levels(skin$seurat_clusters),sep='')
colnames(wnn.rna)=rownames(skin@assays$RNA)
skin.rna.count=skin@assays$RNA@counts
skin.rna.count.adj=sweep(skin.rna.count,2,colSums(skin.rna.count)/10^6,`/`) # adjust single-cell library size before computing the sd and mad
for(i in 1:nrow(wnn.rna)){
  cat(i,' ')
  wnn.rna[i,]=apply(skin.rna.count[,skin$seurat_clusters==(i-1)], 1, sum)
  wnn.rna[i,]=wnn.rna[i,]/sum(wnn.rna[i,])*10^6 # adjust for lib. size for each wnn
}
rm(skin.rna.count); rm(skin.rna.count.adj)

wnn.peak=matrix(nrow=length(levels(skin$seurat_clusters)), ncol=nrow(skin@assays$ATAC))
rownames(wnn.peak)=paste('wnn_', levels(skin$seurat_clusters),sep='')
colnames(wnn.peak)=rownames(skin@assays$ATAC)
skin.atac.count=skin@assays$ATAC@counts
skin.atac.count.adj=sweep(skin.atac.count,2,colSums(skin.atac.count)/10^6,`/`) # adjust single-cell library size before computing the sd and mad
for(i in 1:nrow(wnn.peak)){
  cat(i,' ')
  wnn.peak[i,]=apply(skin.atac.count[,skin$seurat_clusters==(i-1)], 1, sum)
  wnn.peak[i,]=wnn.peak[i,]/sum(wnn.peak[i,])*10^6 # adjust for lib. size for each wnn
}
rm(skin.atac.count); rm(skin.atac.count.adj)

wnn.rna[1:5, 1:5] # wnn by gene matrix: sum of single-cell RNA read count
wnn.peak[1:5, 1:5] # wnn by peak matrix: sum of single-cell ATAC read count

# Get WNN cell types: here I'm using the cell type (i.e., Medulla, TAC-1 etc.)
temp=table(skin$seurat.celltype, skin$seurat_clusters)
wnn.celltype=rep(NA, nrow(wnn.rna))
for(i in 1:length(wnn.celltype)){
  temp.i_1=temp[,colnames(temp)==as.character(i-1)]
  wnn.celltype[i]=names(temp.i_1)[which.max(temp.i_1)]
}

wnn.pseudotime=rep(NA, nrow(wnn.rna))
for(i in 1:length(wnn.celltype)){
  cell.i.index=which(skin$seurat_clusters==(i-1))
  wnn.pseudotime[i]=mean((skin$ptime)[cell.i.index])
}

boxplot(wnn.pseudotime~wnn.celltype)

# Get the corresponding color for each cell type from Seurat
p <- Seurat::DimPlot(skin, reduction = 'wnn.umap', label=T, group.by='seurat.celltype')
pbuild <- ggplot2::ggplot_build(p) # Use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # this is to get the color palette by Seurat
pdata =cbind(skin$seurat.celltype, pdata)
wnn.celltype.col=rep(NA, length(wnn.celltype))
for(i in 1:length(wnn.celltype)){
  wnn.celltype.col[i]=pdata$colour[min(which(pdata$`skin$seurat.celltype`==wnn.celltype[i]))]
}


skin=SetIdent(skin, value = 'seurat_clusters')
names(wnn.celltype) =names(wnn.celltype.col)= levels(skin$seurat_clusters)
skin <- RenameIdents(skin, wnn.celltype)
DimPlot(skin, reduction = "wnn.umap", label = TRUE, pt.size = 0.5)

# Remove unnecessary stored values/objects
rm(p); rm(pbuild); rm(pdata)
rm(temp)

save.image(file='ATAC_RNA_WNN_branch.rda')
