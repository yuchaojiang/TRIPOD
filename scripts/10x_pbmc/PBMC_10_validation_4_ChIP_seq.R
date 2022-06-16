
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

load('ATAC_RNA_WNN_integrated_high_var.rda')

# Remove all previously loaded functions
rm(list=lsf.str())
# Source updated functions
sapply(list.files(path='./R/',pattern='.R', full.names = TRUE),source,.GlobalEnv)

#############################################
### Validating TF-peak: TF ChIP-seq
#############################################
library(openxlsx)
chip.data=read.xlsx('../data/TF_ChIPseq/human_pbmc_filtered_list_20210312.xlsx')

cell.type='B Lymphocyte'
#cell.type='Monocyte'
#cell.type='T Lymphocyte'
factors=chip.data[chip.data$Cell_type==cell.type,]$Factor
factors=factors[isUnique(factors)]

wnn.rm=which(wnn.celltype%in% c('B memory','B naive','B intermediate'))
#wnn.rm=which(wnn.celltype%in% c('CD14 Mono','CD16 Mono'))
#wnn.rm=which(wnn.celltype%in% c('CD4 Naive', 'CD8 Naive','CD4 TCM','CD8 TEM'))
wnn.celltype[wnn.rm]

peak.overlap.chip=matrix(ncol=2, nrow=length(factors))
rownames(peak.overlap.chip)=factors
colnames(peak.overlap.chip)=c('AllPeaks','SigTFPeaks')

for(factor in factors){
  cat(factor,'\n')
  chip.id=chip.data[which(chip.data$Cell_type==cell.type & chip.data$Factor==factor),]$DCid
  chip.peak=read.table(paste0('../data/TF_ChIPseq/human_pbmc/',chip.id,'_peaks.bed'))
  chip.peak=GRanges(seqnames=chip.peak$V1, ranges = IRanges(start=chip.peak$V2, end = chip.peak$V3), intensity=chip.peak$V5)
  
# Add influential test for B cells
for(g in 1:length(xymats.all)){
  if (g%%20 ==0) cat(g," ")
  xymats=xymats.all[[g]]$xymats4
  pval_Bcells=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g),data=NA)

  j=which(xymats$TF.g==factor) # Which TF to look at ?
  Xit=xymats$X
  Yig=xymats$Y
  Yij=xymats$Y.TF[,j]
  if(sum(Yij>0)>4){
    for(t in which(xymats$peakxmotif.g[,j]==1)){
      if(all(Xit[,t]*Yij==0)) next
      pval_Bcells[t,j]=test.influential(Yig, Xit[,t], Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)['Yig']
    }
  }
  xymats.all[[g]]$xymats4$pval_Bcells=pval_Bcells
}

  peak.sig.bcell=peak.binding=peak.gr[0]
  for(i in 1:length(xymats.all)){
    if(i %%20 ==0) cat(i,'\t')
    gene.name=names(xymats.all)[i]
    xymats4=xymats.all[[gene.name]]$xymats4
    j=which(xymats4$TF.g==factor)
    if(length(j)==0) next
    peak.binding=c(peak.binding, xymats4$peak.gr.g[xymats4$peakxmotif.g[,j]>0])
    
    gammas.rej4=benjaminiHochsbergMatrix(xymats4$pvalgammas, do.plot = FALSE)
    beta.rej5X1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5X2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot=FALSE)
    beta.rej5Y1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5Y2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot=FALSE)
    
    peak.tf.sig=which(xymats4$pval_Bcells <0.05 & (gammas.rej4| beta.rej5X1 | beta.rej5Y1 | gammas.rej5X2 | gammas.rej5Y2), arr.ind = T)
    if(length(peak.tf.sig)>0){
      peak.tf.sig=peak.tf.sig[peak.tf.sig[,2]==j,, drop=FALSE]
      peak.sig.bcell=c(peak.sig.bcell, xymats4$peak.gr.g[peak.tf.sig[,1]])
    }
  }
  peak.sig.bcell=sort(peak.sig.bcell)
  peak.binding=sort(peak.binding)
  peak.insig.bcell=peak.binding[countOverlaps(peak.binding, peak.sig.bcell)==0]
  
  chip.peak=chip.peak[countOverlaps(chip.peak, peak.gr)>0]
  
  peak.overlap.chip[factor,]=c(sum(countOverlaps(peak.gr, chip.peak)>0)/length(peak.gr),
                               sum(countOverlaps(peak.sig.bcell, chip.peak)>0)/length(peak.sig.bcell))
  cat('\n\n')
}

peak.overlap.chip

save.image(file=paste0('peak.overlap.chip.bcell.rda'))



library(openxlsx)
chip.data=read.xlsx('../data/TF_ChIPseq/human_pbmc_filtered_list_20210312.xlsx')

#cell.type='B Lymphocyte'
cell.type='Monocyte'
#cell.type='T Lymphocyte'
factors=chip.data[chip.data$Cell_type==cell.type,]$Factor
factors=factors[isUnique(factors)]

#wnn.rm=which(wnn.celltype%in% c('B memory','B naive','B intermediate'))
wnn.rm=which(wnn.celltype%in% c('CD14 Mono','CD16 Mono'))
#wnn.rm=which(wnn.celltype%in% c('CD4 Naive', 'CD8 Naive','CD4 TCM','CD8 TEM'))
wnn.celltype[wnn.rm]

peak.overlap.chip=matrix(ncol=2, nrow=length(factors))
rownames(peak.overlap.chip)=factors
colnames(peak.overlap.chip)=c('AllPeaks','SigTFPeaks')

for(factor in factors){
  cat(factor,'\n')
  chip.id=chip.data[which(chip.data$Cell_type==cell.type & chip.data$Factor==factor),]$DCid
  chip.peak=read.table(paste0('../data/TF_ChIPseq/human_pbmc/',chip.id,'_peaks.bed'))
  chip.peak=GRanges(seqnames=chip.peak$V1, ranges = IRanges(start=chip.peak$V2, end = chip.peak$V3), intensity=chip.peak$V5)
  
  # Add influential test for B cells
  for(g in 1:length(xymats.all)){
    if (g%%20 ==0) cat(g," ")
    xymats=xymats.all[[g]]$xymats4
    pval_Bcells=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g),data=NA)

    j=which(xymats$TF.g==factor) # Which TF to look at ?
    Xit=xymats$X
    Yig=xymats$Y
    Yij=xymats$Y.TF[,j]
    if(sum(Yij>0)>4){
      for(t in which(xymats$peakxmotif.g[,j]==1)){
        if(all(Xit[,t]*Yij==0)) next
        pval_Bcells[t,j]=test.influential(Yig, Xit[,t], Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)['Yig']
      }
    }
    xymats.all[[g]]$xymats4$pval_Bcells=pval_Bcells
  }

  peak.sig.bcell=peak.binding=peak.gr[0]
  for(i in 1:length(xymats.all)){
    if(i %%20 ==0) cat(i,'\t')
    gene.name=names(xymats.all)[i]
    xymats4=xymats.all[[gene.name]]$xymats4
    j=which(xymats4$TF.g==factor)
    if(length(j)==0) next
    peak.binding=c(peak.binding, xymats4$peak.gr.g[xymats4$peakxmotif.g[,j]>0])
    
    gammas.rej4=benjaminiHochsbergMatrix(xymats4$pvalgammas, do.plot = FALSE)
    beta.rej5X1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5X2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot=FALSE)
    beta.rej5Y1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5Y2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot=FALSE)
    
    peak.tf.sig=which(xymats4$pval_Bcells <0.05 & (gammas.rej4| beta.rej5X1 | beta.rej5Y1 | gammas.rej5X2 | gammas.rej5Y2), arr.ind = T)
    if(length(peak.tf.sig)>0){
      peak.tf.sig=peak.tf.sig[peak.tf.sig[,2]==j,, drop=FALSE]
      peak.sig.bcell=c(peak.sig.bcell, xymats4$peak.gr.g[peak.tf.sig[,1]])
    }
  }
  peak.sig.bcell=sort(peak.sig.bcell)
  peak.binding=sort(peak.binding)
  peak.insig.bcell=peak.binding[countOverlaps(peak.binding, peak.sig.bcell)==0]
  
  chip.peak=chip.peak[countOverlaps(chip.peak, peak.gr)>0]
  
  peak.overlap.chip[factor,]=c(sum(countOverlaps(peak.gr, chip.peak)>0)/length(peak.gr),
                               sum(countOverlaps(peak.sig.bcell, chip.peak)>0)/length(peak.sig.bcell))
  cat('\n\n')
}

peak.overlap.chip

save.image(file='peak.overlap.chip.mono.rda')

library(openxlsx)
chip.data=read.xlsx('../data/TF_ChIPseq/human_pbmc_filtered_list_20210312.xlsx')

#cell.type='B Lymphocyte'
#cell.type='Monocyte'
cell.type='T Lymphocyte'
factors=chip.data[chip.data$Cell_type==cell.type,]$Factor
factors=factors[isUnique(factors)]

#wnn.rm=which(wnn.celltype%in% c('B memory','B naive','B intermediate'))
#wnn.rm=which(wnn.celltype%in% c('CD14 Mono','CD16 Mono'))
wnn.rm=which(wnn.celltype%in% c('CD4 Naive', 'CD8 Naive','CD4 TCM','CD8 TEM'))
wnn.celltype[wnn.rm]

peak.overlap.chip=matrix(ncol=2, nrow=length(factors))
rownames(peak.overlap.chip)=factors
colnames(peak.overlap.chip)=c('AllPeaks','SigTFPeaks')

for(factor in factors){
  cat(factor,'\n')
  chip.id=chip.data[which(chip.data$Cell_type==cell.type & chip.data$Factor==factor),]$DCid
  chip.peak=read.table(paste0('../data/TF_ChIPseq/human_pbmc/',chip.id,'_peaks.bed'))
  chip.peak=GRanges(seqnames=chip.peak$V1, ranges = IRanges(start=chip.peak$V2, end = chip.peak$V3), intensity=chip.peak$V5)
  
  # Add influential test for B cells
  for(g in 1:length(xymats.all)){
    if (g%%20 ==0) cat(g," ")
    xymats=xymats.all[[g]]$xymats4
    pval_Bcells=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g),data=NA)

    j=which(xymats$TF.g==factor) # Which TF to look at ?
    Xit=xymats$X
    Yig=xymats$Y
    Yij=xymats$Y.TF[,j]
    if(sum(Yij>0)>4){
      for(t in which(xymats$peakxmotif.g[,j]==1)){
        if(all(Xit[,t]*Yij==0)) next
        pval_Bcells[t,j]=test.influential(Yig, Xit[,t], Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)['Yig']
      }
    }
    xymats.all[[g]]$xymats4$pval_Bcells=pval_Bcells
  }
  
  peak.sig.bcell=peak.binding=peak.gr[0]
  for(i in 1:length(xymats.all)){
    if(i %%20 ==0) cat(i,'\t')
    gene.name=names(xymats.all)[i]
    xymats4=xymats.all[[gene.name]]$xymats4
    j=which(xymats4$TF.g==factor)
    if(length(j)==0) next
    peak.binding=c(peak.binding, xymats4$peak.gr.g[xymats4$peakxmotif.g[,j]>0])
    
    gammas.rej4=benjaminiHochsbergMatrix(xymats4$pvalgammas, do.plot = FALSE)
    beta.rej5X1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5X2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals, do.plot=FALSE)
    beta.rej5Y1=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest, do.plot=FALSE)
    gammas.rej5Y2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals, do.plot=FALSE)
    
    peak.tf.sig=which(xymats4$pval_Bcells <0.05 & (gammas.rej4| beta.rej5X1 | beta.rej5Y1 | gammas.rej5X2 | gammas.rej5Y2), arr.ind = T)
    if(length(peak.tf.sig)>0){
      peak.tf.sig=peak.tf.sig[peak.tf.sig[,2]==j,, drop=FALSE]
      peak.sig.bcell=c(peak.sig.bcell, xymats4$peak.gr.g[peak.tf.sig[,1]])
    }
  }
  peak.sig.bcell=sort(peak.sig.bcell)
  peak.binding=sort(peak.binding)
  peak.insig.bcell=peak.binding[countOverlaps(peak.binding, peak.sig.bcell)==0]
  
  chip.peak=chip.peak[countOverlaps(chip.peak, peak.gr)>0]
  
  peak.overlap.chip[factor,]=c(sum(countOverlaps(peak.gr, chip.peak)>0)/length(peak.gr),
                               sum(countOverlaps(peak.sig.bcell, chip.peak)>0)/length(peak.sig.bcell))
  cat('\n\n')
}

peak.overlap.chip

save.image(file='peak.overlap.chip.tcell.rda')


load('peak.overlap.chip.bcell.rda')
temp=peak.overlap.chip
temp=temp[apply(temp, 1, function(x){all(!is.nan(x))}),]
library(shape)
temp=melt(temp)
toplot=cbind(temp, celltype=rep('Bcell',nrow(temp)))

load('peak.overlap.chip.mono.rda')
temp=peak.overlap.chip
temp=temp[apply(temp, 1, function(x){all(!is.nan(x))}),]
temp=melt(temp)
toplot=rbind(toplot,cbind(temp, celltype=rep('Mono',nrow(temp))))

load('peak.overlap.chip.tcell.rda')
temp=peak.overlap.chip
temp=temp[apply(temp, 1, function(x){all(!is.nan(x))}),]
temp=melt(temp)
toplot=rbind(toplot,cbind(temp, celltype=rep('Tcell',nrow(temp))))

colnames(toplot)=c('TF','Peak','Percent','Celltype')

temp=toplot
toplot$TF=paste(toplot$Celltype,toplot$TF)

p <- ggplot(data=toplot, aes(x=TF, y=Percent, fill=Peak)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  theme_minimal()+ylab('Percentage of overlap')+xlab('Transcription factor')
p
ggsave(filename = 'chip_seq_percentage_overlap.pdf', plot = p, width = 7, height = 4)

