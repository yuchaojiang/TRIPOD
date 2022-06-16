test.influential = function(Yig, Xit, Yij, wnn.rm, plot.histogram=NULL,  
                            nsamp=NULL, ...){
  if(is.null(plot.histogram)) plot.histogram=FALSE
  if(is.null(nsamp)) nsamp=10000
  model4=lm(Yig~1+Xit+Yij+Xit*Yij) # Full model with all data points / WNNs
  # Model without the WNNs selected (observed)
  model4.delta.obs=lm(Yig[-wnn.rm]~1+Xit[-wnn.rm]+Yij[-wnn.rm]+(Xit*Yij)[-wnn.rm])
  
  # Yig.pred=(model4$fitted.values)[wnn.rm]
  # Yig.wnn.rm.pred=t(as.matrix(model4.delta.obs$coefficients))%*%rbind(rep(1, length(wnn.rm)), Xit[wnn.rm], Yij[wnn.rm], (Xit*Yij)[wnn.rm])
  Yig.pred=(model4$fitted.values)
  Yig.wnn.rm.pred=t(as.matrix(model4.delta.obs$coefficients))%*%rbind(rep(1, length(Yig)), Xit, Yij, (Xit*Yij))
  
  delta.coeff=matrix(ncol=5,nrow=nsamp)
  delta.coeff[1,1:4]=model4.delta.obs$coefficients-model4$coefficients
  delta.coeff[1,5]=mean(abs(Yig.pred-Yig.wnn.rm.pred))
  colnames(delta.coeff)=c(names(model4$coefficients), 'Yig')
  
  for(ii in 2:nsamp){
    wnn.rm.samp=sample(1:length(wnn.celltype),length(wnn.rm))
    model4.delta.samp=lm(Yig[-wnn.rm.samp]~1+Xit[-wnn.rm.samp]+
                           Yij[-wnn.rm.samp]+(Xit*Yij)[-wnn.rm.samp])
    # Yig.wnn.rm.samp.pred=t(as.matrix(model4.delta.samp$coefficients))%*%rbind(rep(1, length(wnn.rm)), Xit[wnn.rm], Yij[wnn.rm], (Xit*Yij)[wnn.rm])
    Yig.wnn.rm.samp.pred=t(as.matrix(model4.delta.samp$coefficients))%*%rbind(rep(1, length(Yig)), Xit, Yij, (Xit*Yij))
    
    delta.coeff[ii,1:4]=model4.delta.samp$coefficients-model4$coefficients
    delta.coeff[ii,5]=mean(abs(Yig.pred-Yig.wnn.rm.samp.pred))
  }
  delta.coeff.pval=apply(delta.coeff,2, function(x){sum(abs(x)>=abs(x[1]))/length(x)})
  delta.coeff.pval
  if(plot.histogram){
    par(mfrow=c(2,3))
    for(i in 1:5){
      hist(delta.coeff[,i],100, xlab=paste('Delta for',colnames(delta.coeff)[i]),
           main=paste('Delta coefficient for',colnames(delta.coeff)[i],': pval =', delta.coeff.pval[i]));abline(v=delta.coeff[1,i], col=2,lty=2)
    }
    par(mfrow=c(1,1))
  }
  return(delta.coeff.pval)
}
