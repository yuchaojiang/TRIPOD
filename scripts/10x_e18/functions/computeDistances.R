computeDistances<-function(x){
  if(!is.matrix(x)){
    x=matrix(data=x,nrow=length(x),ncol=1)
  }
  
  dmat=matrix(nrow=nrow(x), ncol=nrow(x), data=NA)
  for(i in 1:nrow(x)){
    xi=matrix(nrow=nrow(x), ncol=ncol(x), data=x[i,], byrow=TRUE)  
    dmat[i,]=sqrt(rowSums((xi-x)^2)/ncol(x))
  }
  dmat
  
}
