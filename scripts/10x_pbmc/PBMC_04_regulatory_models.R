
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
library(plotrix)

load('ATAC_RNA_WNN.rda') 

# Remove peaks that are too wide
peak.filter=width(peak.gr)<6000
pbmc@assays$ATAC=subset(pbmc@assays$ATAC,rownames(pbmc@assays$ATAC)[width(pbmc@assays$ATAC@ranges)<=6000])
peak.gr=peak.gr[peak.filter]
wnn.peak=wnn.peak[,peak.filter]
peakxmotif=peakxmotif[peak.filter,]
rm(peak.filter)

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

ext.upstream=100000 # This is actually both up and downstream
transcripts <- transcripts.gr
TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
TSS <- GRanges(seqnames = seqnames(transcripts),
               ranges = IRanges(start = TSS.position, width = 1),
               strand = strand(transcripts),
               gene_name = transcripts$gene_name)
transcripts.ext.gr=getRegion(c(ext.upstream, ext.upstream), TSS)
rm(transcripts); rm(TSS.position); rm(TSS)

length(wnn.celltype); length(wnn.celltype.col)
dim(wnn.peak); pbmc@assays$ATAC
dim(wnn.rna); pbmc@assays$RNA; pbmc@assays$SCT
length(transcripts.gr); length(transcripts.ext.gr)
dim(peakxmotif); dim(motifxTF); pbmc@assays$chromvar

# Get highly variable genes
# Below are the cell-type markers by Seurat's tutorial
marker_genes=c('CEBPB', 'TBX21', 'IRF4', 'SOX4', 'EBF1', 'PAX5', 'IRF8', 'TCF4', 
              'IL7R', 'CCR7', 'IL7R', 'S100A4', 'CD14', 'LYZ', 'MS4A1', 'CD8A', 
              'FCGR3A', 'MS4A7', 'GNLY', 'NKG7', 'FCER1A', 'CST3')
celltype_presto <- presto:::wilcoxauc.Seurat(X = pbmc, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
celltype_genes=celltype_presto[order(celltype_presto[,'pval'])[1:1000],'feature'] # cell-type markers
highvar_genes=pbmc@assays$SCT@var.features[1:1000]
genes=unique(c(celltype_genes, highvar_genes))
rm(celltype_presto)
all(marker_genes %in% genes)

gene.xymats4=vector(mode = "list", length = length(genes))
names(gene.xymats4)=genes
gene.xymats5X=gene.xymats5Y=gene.xymats2=gene.xymats1=gene.xymats4

for(i in 1:length(genes)){
  if(i%%50==0 | i<=5) cat(i,'\t')
  gene.name=genes[i]
  xymats=getXYMatrices(gene.name=gene.name,
                       ext.upstream=ext.upstream,
                       transcripts.gr=transcripts.gr, peak.gr=peak.gr,
                       wnn.rna=wnn.rna, wnn.peak=wnn.peak,
                       peakxmotif=peakxmotif, motifxTF=motifxTF,
                       wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col)
  if(ncol(xymats$Xit)==0 | ncol(xymats$Y.TF)==0) next
  
  # Model 1
  xymats1=fitModel(xymats, modelName='1a')
  gene.xymats1[[gene.name]]=xymats1
  
  # Model 2
  xymats2=fitModel(xymats, modelName='2a')
  gene.xymats2[[gene.name]]=xymats2
  
  # Model 4
  xymats4=fitModel(xymats, modelName='4a')
  gene.xymats4[[gene.name]]=xymats4
  
  # Model 5
  xymats5X=fitModel(xymats, modelName='5a', match.by='Xt', do.plot=FALSE, VERBOSE=FALSE)
  xymats5Y=fitModel(xymats, modelName='5a', match.by='Yj', do.plot=FALSE, VERBOSE=FALSE)
  gene.xymats5X[[gene.name]]=xymats5X
  gene.xymats5Y[[gene.name]]=xymats5Y
}

gene.filter=rep(TRUE, length(gene.xymats5X)) # Remove highly variable genes that dont have binding motifs for any TFs
for(i in 1:length(gene.xymats5X)){
  if(is.null(gene.xymats5X[[i]])) gene.filter[i]=FALSE
}
sum(gene.filter)
gene.xymats1=gene.xymats1[gene.filter]
gene.xymats2=gene.xymats2[gene.filter]
gene.xymats4=gene.xymats4[gene.filter]
gene.xymats5X=gene.xymats5X[gene.filter]
gene.xymats5Y=gene.xymats5Y[gene.filter]


# Combine all xymats
genes=names(gene.xymats1)
xymats.all=vector(mode = "list", length = length(genes))
names(xymats.all)=genes
for(i in 1:length(xymats.all)){
  if(i %%50==0) cat(i,'\t')
  gene.name=genes[i]
  xymats.all[[gene.name]]=list(xymats1=gene.xymats1[[gene.name]],
                               xymats2=gene.xymats2[[gene.name]],
                               xymats4=gene.xymats4[[gene.name]],
                               xymats5X=gene.xymats5X[[gene.name]],
                               xymats5Y=gene.xymats5Y[[gene.name]])
}
rm(gene.xymats1); rm(gene.xymats2); rm(gene.xymats4); rm(gene.xymats5X); rm(gene.xymats5Y)
rm(xymats1); rm(xymats2); rm(xymats4); rm(xymats5X); rm(xymats5Y)

save.image(file='ATAC_RNA_WNN_integrated_high_var.rda')

# Gene CEBPB
gene.name='CEBPB'
#########################
# ADD Yj ~ Xt
#########################
names(xymats.all[[gene.name]]$xymats1) # Yg~Xt, Yg~Yj, Yg~Xt*Yj
names(xymats.all[[gene.name]]$xymats2) # Yg~Xt+Yj
names(xymats.all[[gene.name]]$xymats4) # Yg~Xt+Yj+Xt*Yj
names(xymats.all[[gene.name]]$xymats5X) # Match by Xt
names(xymats.all[[gene.name]]$xymats5Y) # Match by Yj

## PLOTTING  A TRIO
# CCR7
gene.name='CCR7'
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==40596032) # Peak access. chr17−40596032−40600718
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='LEF1') # TF expr. LEF1 motif MA0768.1

plot.dir='./'
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4,
                  sel.t=peak.num, jj=TF.num,
                  plot.dir=plot.dir, to.plot =c('umap'), generate.pdf = TRUE)

pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
par(mfrow=c(3,4))
# Marginal: model 1
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 1 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 2 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 4 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
# Model 5: need to include this plot in the plot_gene_peak_TF function
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Yj",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Xt",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
dev.off()

# GNLY
gene.name='GNLY'
peak.num=which( start(xymats.all[[gene.name]]$xymats4$peak.gr.g)==85693711) # Peak access. chr2−85693711−85699520
TF.num=which(xymats.all[[gene.name]]$xymats4$TF.g=='TBX21') # TF expr. TBX21 motif MA0690.1

plot.dir='./'
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4,
                  sel.t=peak.num, jj=TF.num,
                  plot.dir=plot.dir, to.plot =c('umap'), generate.pdf = TRUE)

pdf(file = paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_models1245_scatterplot.pdf',sep=''), height=8, width=12)
par(mfrow=c(3,4))
# Marginal: model 1
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('marginal'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 1 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats1, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model1.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 2 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats2, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model2.scatter'), generate.pdf = FALSE, partition.screen = FALSE)
# Model 4 scatterplot
plot_gene_peak_TF(object = pbmc, group.by='celltype',xymats = xymats.all[[gene.name]]$xymats4, sel.t=peak.num, jj=TF.num, plot.dir=plot.dir, to.plot =c('model4.scatter'),  generate.pdf = FALSE, partition.screen = FALSE)
# Model 5: need to include this plot in the plot_gene_peak_TF function
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Yj",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
res=matchedTest(Xt=xymats.all[[gene.name]]$xymats5X$X[,peak.num],
                Yj=xymats.all[[gene.name]]$xymats5X$Y.TF[,TF.num],
                Yg=xymats.all[[gene.name]]$xymats5X$Y,
                match.by="Xt",
                do.plot=TRUE,
                col=wnn.celltype.col,
                verbose=FALSE,
                partition.screen=FALSE)
dev.off()

