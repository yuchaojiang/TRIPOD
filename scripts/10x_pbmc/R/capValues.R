capValues<-function(x, CAP.AT.QUANTILE=NULL){
  if(is.null(CAP.AT.QUANTILE)) CAP.AT.QUANTILE=0.02
  upperThresh=quantile(x, 1-CAP.AT.QUANTILE, na.rm=TRUE)
  lowerThresh=quantile(x, CAP.AT.QUANTILE, na.rm=TRUE)
  x=pmin(pmax(x, lowerThresh), upperThresh)
  x
}
