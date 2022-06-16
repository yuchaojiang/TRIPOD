#################################
### We will first process an independent PBMC 3K dataset
### and use it as a testing dataset
#################################

library(Seurat)
library(SeuratDisk)
library(patchwork)
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

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5")

# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# First, we need to calculate the peak coverage using the 10k cell peaks
frags.3k <- CreateFragmentObject(
  path = "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz",
  cells = colnames(atac_counts)
)

load('prediction/peak.gr.rda')
pbmc3k.counts <- FeatureMatrix(
  fragments = frags.3k,
  features = peak.gr,
  cells = colnames(atac_counts)
)

# Create Seurat object
pbmc <- CreateSeuratObject(counts = rna_counts)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
all(seqnames(peak.gr) %in% standardChromosomes(peak.gr))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

length(peak.gr)
dim(pbmc3k.counts)
dim(rna_counts)

frag.file <- "pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = pbmc3k.counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 0,
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

p1 <- DimPlot(pbmc, reduction = "umap.rna", group.by='seurat_clusters', label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X PBMC RNA") + NoLegend()
p2 <- DimPlot(pbmc, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("10X PBMC ATAC")+ NoLegend()
p1+p2+plot_layout(ncol=1)
p1|p2

wnn.resolution=10
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
wnn.rna.sd=wnn.rna.mad=wnn.rna.cv=wnn.rna.gini=wnn.rna
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
wnn.peak.sd=wnn.peak.mad=wnn.peak.cv=wnn.peak.gini=wnn.peak
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

pbmc3k=pbmc
wnn.rna3k=wnn.rna
wnn.peak3k=wnn.peak
save(pbmc3k, file='prediction/pbmc3k.rda')
save(wnn.rna3k, file='prediction/wnn.rna3k.rda')
save(wnn.peak3k, file='prediction/wnn.peak3k.rda')

#################################
### Train prediction model PBMC 10K
### Test PBMC 3K
#################################

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(presto)

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
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)
library(gplots)
library(nbpMatching)
library(fields)

load('ATAC_RNA_WNN.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

getXYMatrices_lasso<-function(gene.name, ext.upstream, ext.downstream=NULL, modelName, ...){
  if(is.null(ext.downstream)) ext.downstream=ext.upstream
  transcripts <- transcripts.gr
  TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
  TSS <- GRanges(seqnames = seqnames(transcripts),
                 ranges = IRanges(start = TSS.position, width = 1),
                 strand = strand(transcripts),
                 gene_name = transcripts$gene_name)
  transcripts.ext.gr=getRegion(c(ext.upstream, ext.downstream), TSS)
  
  Y=wnn.rna[,gene.name]
  Y.test=wnn.rna3k[,gene.name]
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  transcripts.gr.g # gene Granges
  transcripts.ext.gr.g # extended upstream gene Granges
  
  peak.gr.g=subsetByOverlaps(peak.gr, transcripts.ext.gr.g) # peaks within upstream gene regions
  Xit=wnn.peak[,countOverlaps(peak.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name])>0,drop=FALSE] # counts within these peaks.
  Xit.test=wnn.peak3k[,countOverlaps(peak.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name])>0,drop=FALSE]
  
  # Motifs within vicinity region of gene g
  peakxmotif.g=peakxmotif[countOverlaps(peak.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name])>0,, drop=FALSE]
  peakxmotif.g=peakxmotif.g[,apply(peakxmotif.g,2,sum)>0, drop=FALSE] # remove motifs that do not have binding sites in this gene region
  TF.g=motifxTF[match(colnames(peakxmotif.g),motifxTF[,1]),2]
  Y.TF=wnn.rna[,TF.g,drop=FALSE] # Y_ij: gene expression of TF j
  Y.TF.test=wnn.rna3k[,TF.g, drop=FALSE]
  
  TF.filter=apply(Y.TF,2,sum)>0 # Remove a TF if it is not expressed in any cells
  Y.TF=Y.TF[,TF.filter, drop=FALSE]
  Y.TF.test=Y.TF.test[,TF.filter, drop=FALSE]
  TF.g=TF.g[TF.filter, drop=FALSE]
  peakxmotif.g=peakxmotif.g[,TF.filter, drop=FALSE]
  rm(TF.filter)
  
  nonzero.peakxmotif.g=which(peakxmotif.g!=0, arr.ind=T)
  X=Xit[,nonzero.peakxmotif.g[,1], drop=FALSE]*Y.TF[,nonzero.peakxmotif.g[,2],drop=FALSE]
  X.test=Xit.test[,nonzero.peakxmotif.g[,1], drop=FALSE]*Y.TF.test[,nonzero.peakxmotif.g[,2],drop=FALSE]
  peakTF.filter=apply(X, 2, sum)>0 # additional QC: X_it * Y_ij can be zero across all i even though they themselves are non-zeros 
  nonzero.peakxmotif.g=nonzero.peakxmotif.g[peakTF.filter,,drop=FALSE]
  X=X[,peakTF.filter, drop=FALSE]
  X.test=X.test[,peakTF.filter, drop=FALSE]
  return(list(modelName=modelName, gene.name=gene.name, ext.upstream=ext.upstream,
              Y=Y, Y.test=Y.test, X=X,X.test=X.test, 
              Xit=Xit, Xit.test=Xit.test, Y.TF=Y.TF, Y.TF.test=Y.TF.test, 
              peakxmotif.g=peakxmotif.g, nonzero.peakxmotif.g=nonzero.peakxmotif.g, 
              TF.g=TF.g, peak.gr.g=peak.gr.g, 
              wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col))
}

celltype_presto <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
celltype_genes=celltype_presto[order(celltype_presto[,'pval'])[1:1000],'feature'] # cell-type markers
highvar_genes=pbmc@assays$SCT@var.features[1:1000]
genes=unique(c(celltype_genes, highvar_genes))

load('prediction/wnn.rna3k.rda')
load('prediction/wnn.peak3k.rda')
all(genes %in% colnames(wnn.rna3k))
dim(wnn.peak3k)

for(dist in c(2000, 5000, 20000, 50000, 200000)){
  cat('Dist',dist,'\n\n')
  
  gene.xymats.lasso=vector(mode = "list", length = length(genes))
  
  ext.upstream=dist
  for(i in 1:length(genes)){
    if(i%%50==0 | i<=5) cat(i,'\t')
    gene.name=genes[i]
    xymats=getXYMatrices_lasso(gene.name, ext.upstream=ext.upstream, modelName='3a')
    if(ncol(xymats$X) <= 2 | ncol(xymats$Xit) <= 2 | ncol(xymats$Y.TF) <=2 |
       sum(xymats$Y >0)<=10) next
    
    # Remove outliers
    xymats$X=apply(xymats$X, 2, capValues)
    xymats$Y=capValues(xymats$Y)
    xymats$Xit=apply(xymats$Xit, 2, capValues)
    xymats$Y.TF=apply(xymats$Y.TF, 2, capValues)
    
    # LASSO using peak x TF
    lasso.mod <- glmnet(xymats$X, xymats$Y, alpha = 1)  # 10-fold cross-validation to find the optimal lambda
    set.seed(1) 
    cv.out <- cv.glmnet(xymats$X, xymats$Y, alpha = 1)
    # plot(cv.out)
    bestlam <- cv.out$lambda.min
    if(max(lasso.mod$df)>1){
      bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df>1])) # Less penalty
    }
    fit <- glmnet(xymats$X, xymats$Y, alpha=1, lambda=bestlam)
    Y.test.pred <- predict(fit, newx=xymats$X.test)
    xymats$lasso.pred.peaktf=Y.test.pred
    
    # LASSO using peak only
    temp=xymats$X
    temp.test=xymats$X.test
    xymats$X=xymats$Xit
    xymats$X.test=xymats$Xit.test
    lasso.mod <- glmnet(xymats$X, xymats$Y, alpha = 1)  # 10-fold cross-validation to find the optimal lambda
    set.seed(1) 
    cv.out <- cv.glmnet(xymats$X, xymats$Y, alpha = 1)
    # plot(cv.out)
    bestlam <- cv.out$lambda.min
    if(max(lasso.mod$df)>1){
      bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df>1])) # Less penalty
    }
    fit <- glmnet(xymats$X, xymats$Y, alpha=1, lambda=bestlam)
    Y.test.pred <- predict(fit, newx=xymats$X.test)
    xymats$lasso.pred.peak=Y.test.pred
    
    # LASSO using TF only
    xymats$X=xymats$Y.TF
    xymats$X.test=xymats$Y.TF.test
    lasso.mod <- glmnet(xymats$X, xymats$Y, alpha = 1)  # 10-fold cross-validation to find the optimal lambda
    set.seed(1) 
    cv.out <- cv.glmnet(xymats$X, xymats$Y, alpha = 1)
    # plot(cv.out)
    bestlam <- cv.out$lambda.min
    if(max(lasso.mod$df)>1){
      bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df>1])) # Less penalty
    }
    fit <- glmnet(xymats$X, xymats$Y, alpha=1, lambda=bestlam)
    Y.test.pred <- predict(fit, newx=xymats$X.test)
    xymats$lasso.pred.tf=Y.test.pred
    
    xymats$X=temp
    xymats$X.test=temp.test
    gene.xymats.lasso[[gene.name]]=xymats
  }
  gene.filter=rep(TRUE, length(gene.xymats.lasso)) # Remove highly variable genes that dont have binding motifs for any TFs
  for(i in 1:length(gene.xymats.lasso)){
    if(is.null(gene.xymats.lasso[[i]])) gene.filter[i]=FALSE
  }
  sum(gene.filter)
  gene.xymats.lasso=gene.xymats.lasso[gene.filter]
  
  save(gene.xymats.lasso, file=paste0('prediction/gene.xymats.lasso_test3k.',dist,'.rda'))
  cat('\n\n')
}


dist=2000
load(paste0('prediction/gene.xymats.lasso_test3k.',dist,'.rda'))
genes=names(gene.xymats.lasso)
pred.cor=matrix(nrow=length(genes), ncol=3)
for(i in 1:nrow(pred.cor)){
  gene.name=genes[i]
  pred.cor[i,]=c(peak.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.peak),
                 tf.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.tf),
                 peaktf.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.peaktf))
}
colnames(pred.cor)=c('peak.cor','tf.cor','peaktf.cor')

library(reshape)
temp=melt(pred.cor)
temp=cbind(temp, dist)
temp=temp[,2:4]
colnames(temp)=c('model','corr','dist')

toplot=temp
for(dist in c(5000, 20000, 50000, 200000)){
  load(paste0('prediction/gene.xymats.lasso_test3k.',dist,'.rda'))
  genes=names(gene.xymats.lasso)
  pred.cor=matrix(nrow=length(genes), ncol=3)
  for(i in 1:nrow(pred.cor)){
    gene.name=genes[i]
    pred.cor[i,]=c(peak.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.peak),
                   tf.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.tf),
                   peaktf.cor=cor(gene.xymats.lasso[[gene.name]]$Y.test, gene.xymats.lasso[[gene.name]]$lasso.pred.peaktf))
  }
  colnames(pred.cor)=c('peak.cor','tf.cor','peaktf.cor')
  temp=melt(pred.cor)
  temp=cbind(temp, dist)
  temp=temp[,2:4]
  colnames(temp)=c('model','corr','dist')
  
  toplot=rbind(toplot,temp)
}
colnames(toplot)=c('Model','Corr','Dist')
toplot$Model=factor(toplot$Model, levels=c('peak.cor','tf.cor','peaktf.cor'))

toplot$Dist=as.factor(toplot$Dist)
toplot=toplot[toplot$Model!='tf.cor',] # Remove the TF prediction (only using TF exp. as features)
# Since it is not really looking at ATAC-RNA relationship.
p=ggplot(toplot, aes(x=Dist, y=Corr, fill=Model)) +
  geom_boxplot(outlier.size = 0.5) + ylim(-0.5, 1)+
  ggtitle(paste('RNA prediction using ATAC'))+
  ylab('Correlation coefficient between Y and Yhat')
ggsave(filename = 'pbmc_test_train_pred.pdf',plot = p, width=6, height = 4)

