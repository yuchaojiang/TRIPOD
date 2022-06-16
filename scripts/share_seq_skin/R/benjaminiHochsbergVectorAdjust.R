benjaminiHochsbergVectorAdjust<-function(pval.c, fdr.thresh=0.01){
  pval.c.rmNA=pval.c[!is.na(pval.c)]
  pvalvec.sorted=sort(pval.c.rmNA)
  pvalvec.ord=order(pval.c.rmNA)
  pval.c.adj.rmNA=pval.c.rmNA
  pval.c.adj.rmNA[pvalvec.ord]=pvalvec.sorted*length(pval.c.rmNA)/c(1:length(pval.c.rmNA))
  pval.c.adj=pval.c
  pval.c.adj[!is.na(pval.c)]=pval.c.adj.rmNA
  which.reject=pval.c.adj<fdr.thresh
  pval.c.adj=pmin(1, pval.c.adj)
  return(list(pval.c.adj=pval.c.adj,which.reject=which.reject))
}
