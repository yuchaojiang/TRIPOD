# For given gene, make the X, Y matrices, and other associated variables
# used for interpreting these matrices.  With the necessary filtering:
#   - Remove TF that do not have binding sites upstream of gene.
#   - Remove TF that are not expressed in any of the cells.
#   - For some models, remove peaks that overlap with genes.
# Currently this function is written so that it accesses
# global variables defined in this script above.  It should be rewritten
# to not rely on these global variables.
getXYMatrices<-function(gene.name, ext.upstream, ext.downstream=NULL,
	transcripts.gr, peak.gr,
	wnn.rna, wnn.peak,
	peakxmotif, motifxTF,
	wnn.celltype, wnn.celltype.col,
	...){
  if(is.null(ext.downstream)) ext.downstream=ext.upstream
  transcripts <- transcripts.gr
  TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
  TSS <- GRanges(seqnames = seqnames(transcripts),
                 ranges = IRanges(start = TSS.position, width = 1),
                 strand = strand(transcripts),
                 gene_name = transcripts$gene_name)
  transcripts.ext.gr=getRegion(c(ext.upstream, ext.downstream), TSS)
  
  Y=wnn.rna[,gene.name]
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  transcripts.gr.g # gene Granges
  transcripts.ext.gr.g # extended upstream gene Granges
  
  peak.gr.g=subsetByOverlaps(peak.gr, transcripts.ext.gr.g) # peaks within upstream gene regions
  Xit=wnn.peak[,countOverlaps(peak.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name])>0,drop=FALSE] # counts within these peaks.
  
  # Motifs within vicinity region of gene g
  peakxmotif.g=peakxmotif[countOverlaps(peak.gr, transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name])>0,, drop=FALSE]
  peakxmotif.g=peakxmotif.g[,apply(peakxmotif.g,2,sum)>0, drop=FALSE] # remove motifs that do not have binding sites in this gene region
  TF.g=motifxTF[match(colnames(peakxmotif.g),motifxTF[,1]),2]
  Y.TF=wnn.rna[,TF.g,drop=FALSE] # Y_ij: gene expression of TF j
  TF.filter=apply(Y.TF,2,sum)>0 # Remove a TF if it is not expressed in any cells
  Y.TF=Y.TF[,TF.filter, drop=FALSE]
  TF.g=TF.g[TF.filter, drop=FALSE]
  peakxmotif.g=peakxmotif.g[,TF.filter, drop=FALSE]
  nonzero.peakxmotif.g=which(peakxmotif.g!=0, arr.ind=T)
  rm(TF.filter)
  
  # Indicator for which peaks DOES NOT overlap with transcribed regions.
  peak.filter= !((countOverlaps(peak.gr.g, transcripts.gr)>0) & (countOverlaps(peak.gr.g, transcripts.gr.g) ==0))
  # Peak expression
  Y.peak=matrix(nrow=nrow(Xit), ncol=ncol(Xit), data=0) 
  colnames(Y.peak)=colnames(Xit); rownames(Y.peak)=rownames(Xit)
  for(j in which(!peak.filter)){
    genes.overlapped=transcripts.gr$gene_name[which(countOverlaps(transcripts.gr, peak.gr.g[j])>0)]
    genes.overlapped=intersect(unique(genes.overlapped), colnames(wnn.rna))
    Y.peak[,j]=apply(wnn.rna[,genes.overlapped, drop=FALSE],1,sum)
  }
  Y.peak[Y.peak<1]=1# epsilon=1 to avoid a zero in the denominator
  return(list(gene.name=gene.name, ext.upstream=ext.upstream,
              Y=Y, Xit=Xit, Y.TF=Y.TF, Y.peak=Y.peak, 
              peakxmotif.g=peakxmotif.g, nonzero.peakxmotif.g=nonzero.peakxmotif.g, 
              TF.g=TF.g, peak.gr.g=peak.gr.g, 
              wnn.celltype=wnn.celltype, wnn.celltype.col=wnn.celltype.col))
}
