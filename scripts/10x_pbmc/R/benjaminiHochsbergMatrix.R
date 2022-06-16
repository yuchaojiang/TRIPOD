benjaminiHochsbergMatrix<-function(pvals, fdr.thresh=0.01,do.plot=NULL){
  if(is.null(do.plot)) do.plot=FALSE
  arr.ind=which(!is.na(pvals))
  pvalvec=pvals[arr.ind]
  pvalvec.sorted=sort(pvalvec)
  pvalvec.ord=order(pvalvec)
  kstar=max(which(pvalvec.sorted<= fdr.thresh*c(1:length(pvalvec))/length(pvalvec)))
  if(do.plot){
    plot(pvalvec.sorted,ylab="Sorted pvalues", main=paste(kstar,"passed", fdr.thresh, "B-H threshold"))
    abline(0, fdr.thresh/length(pvalvec), col="red")
    grid()
  }
  ## Plot 4:  How many significant peaks are there per TF?
  which.reject=matrix(nrow=nrow(pvals), ncol=ncol(pvals), data=FALSE)
  # cat(kstar,"rejections.\n")
  if(kstar>=1 & kstar<=length(pvals)) which.reject[arr.ind[pvalvec.ord[1:kstar]]]=TRUE
  which.reject  
}
