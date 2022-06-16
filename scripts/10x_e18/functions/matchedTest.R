# Performs matched test.
# Xt: vector of accessibility for a given peak t
# Yj: vector of TF expression for a given TF j
# Yg: vector of gene expression for a given gene g
# Z: cell-level covariates (optional).
#     Xt, Yj, Yg, need to have the same length equal to nrow(Z).
# match.by: are we matching by Xt or Yj?  
#     TODO:  only matching by Xt is implemented, make it work when match.by=Yj.
### modified by yuriko trying to make it work when match.by=Yj
# Delta: % of Xt above which to test betahat. (as explained in slides)
# MIN.AFTER.FILTER: minimum number of matched cells after filtering in order do a test
# CAP.AT.QUANTILE: for robustness, cap values below CAP.AT.QUANTILE, or above 1-CAP.AT.QUANTILE

matchedTest<-function(Xt, Yj, Yg, Z=NULL, match.by, 
                      do.plot=TRUE, col="black",plot.dir=NULL, file.name=NULL,
                      partition.screen=TRUE, verbose=TRUE, 
                      Delta=NULL, MIN.POINTS.FOR.BETAHATMEAN=10,
                      MIN.POINTS.FOR.BETAHATCOR=10, MIN.NONZERO.FOR.GAMMAHAT2=10, CAP.AT.QUANTILE=0.02){
  
  to.match=NULL
  match.str=NULL
  if(is.null(Delta)) Delta=0.2
  if(match.by=="Xt"){
    U=Xt # matched predictor
    V=Yj # varied predictor
    match.str="Xt"
    varied.str="Yj"
  }
  if(match.by=="Yj"){
    U=Yj # matched predictor
    V=Xt # varied predictor
    match.str="Yj"
    varied.str="Xt"
  }
  # U=log(U+1) # take the log of the matched predictor # log-transformation is done outside the function
  if(!is.null(Z)){
    to.match=cbind(U, Z)
    match.str=paste(match.str, ",Z", sep="")
  } else {
    to.match=matrix(data=U,nrow=length(U),ncol=1)
  }
  if(is.null(to.match)){ error("to.match value is not valid.")}
  
  dmat = computeDistances(to.match)  # compute distances between the variables in to.match
  options(warn=-1) # suppress the warning on odd number of rows
  dmat2 = distancematrix(dmat)       # convert to distance matrix
  options(warn=0)
  nbm=nonbimatch(dmat2)              # perform matching
  # epsilon=quantile(abs(nbm$halves$Distance), 0.95) # 95th quantile
  epsilon=quantile(abs(nbm$halves$Distance), 0.75)+3*IQR(abs(nbm$halves$Distance))
  which.keep=which(abs(nbm$halves$Distance)<=epsilon )# throw away matches with difference > epsilon.
  which.keep=setdiff(which.keep, c(which(nbm$halves$Group1.ID=='ghost'), which(nbm$halves$Group2.ID=='ghost')))  # modified: also throw away the ghost matching
  
  matchedCells=data.frame(matrix(nrow=length(which.keep), ncol=8))  # one row per matched pair.
  names(matchedCells)=c("ID1", "ID2", "Yg1", "Yg2", "U1", "U2", "V1", "V2")
  matchedCells[,1]=nbm$halves$Group1.Row[which.keep]
  matchedCells[,2]=nbm$halves$Group2.Row[which.keep]
  matchedCells=cbind(matchedCells, CellType1=wnn.celltype[matchedCells[,1]])
  matchedCells=cbind(matchedCells, CellType2=wnn.celltype[matchedCells[,2]])
  
  matchedCells[,3]=Yg[matchedCells[,1]]
  matchedCells[,4]=Yg[matchedCells[,2]]
  matchedCells[,5]=U[matchedCells[,1]]
  matchedCells[,6]=U[matchedCells[,2]]
  matchedCells[,7]=V[matchedCells[,1]]
  matchedCells[,8]=V[matchedCells[,2]]
  matchedcol1=col[nbm$halves$Group1.Row[which.keep]]
  matchedcol2=col[nbm$halves$Group2.Row[which.keep]]
  
  # Compute differences and means
  dU=matchedCells$U2-matchedCells$U1
  meanU=(matchedCells$U1+matchedCells$U2)/2
  dYg=matchedCells$Yg2-matchedCells$Yg1
  dV=matchedCells$V2-matchedCells$V1
  # For robustness, cap values
  dYg=capValues(dYg, CAP.AT.QUANTILE)
  dV=capValues(dV, CAP.AT.QUANTILE)
  meanU=capValues(meanU, CAP.AT.QUANTILE)
  # throw away pairs for which dV==0 for computing pair-wise betahat.
  sel1= dV!=0
  sel1[is.na(sel1)]=FALSE
  Delta.val = max(quantile(meanU, Delta, na.rm=TRUE), 0.1)
  sel2= meanU>=Delta.val # only carry out level-one testing if matched variable U is greater than a threshold
  sel2[is.na(sel2)]=FALSE
  betahat2=rep(NA, length(dYg))
  betahat2[!sel2]=NA
  
  matchedCells=cbind(matchedCells, dU, meanU, dV, dYg, betahat2)
  num.matched.after.meanUfilter.dVfilter=sum(sel2 & sel1)
  
  # Initial values are NA, then these will be set to their appropriate estimates if
  # there are enough points for estimation.
  gammahat2=NA
  betahat.correst=NA
  gammahat.pval2=NA
  betahat.pval.rhotest=NA
  betahat.spearman.result.str <- ""
  gammahat2.result.str <- ""
  
  if(sum(sel2)>MIN.POINTS.FOR.BETAHATCOR){  # This is the more robust one
    res4=cor.test(dYg[sel2], dV[sel2], method="spearman") # matched variable is larger than Delta, only these points are used for β 
    betahat.correst=res4$estimate
    betahat.pval.rhotest=res4$p.value  # betahat estimate based on spearman correlation
    betahat.spearman.result.str=paste("Rank corr = ", format(betahat.correst, digits=2), 
                                      ", pval=", format(betahat.pval.rhotest, digits=2), sep="")
  }
  if(sum(sel1&sel2)>MIN.NONZERO.FOR.GAMMAHAT2){ # need enough nonzero V*meanU
    lmfit=lm(dYg~dV+dV:meanU) # All points are used for gamma: changed from dV*meanU
    gammahat2=summary(lmfit)$coefficients[3,1]
    gammahat.pval2=summary(lmfit)$coefficients[3,4]  # gammahat2 estimate based on regression
    gammahat2.result.str=paste("Partial regr. coeff. = ", format(gammahat2, digits=2), 
                               ", pval=", format(gammahat.pval2, digits=2), sep="")
  }
  
  if(do.plot){
    if(!is.null(plot.dir) && !is.null(file.name)){
      path=paste(plot.dir, file.name, sep="/")
      png(path, height=350, width=800) 
    }
    if(partition.screen) par(mfrow=c(1,2)) 
    # plot(nbm$halves$Distance, xlab="Pair", 
    #      ylab=paste("diff ",match.str," in pair",sep=""), main=paste("Model 5a d",match.str," across matched pairs",sep=""))
    # grid()
    # abline(h=0, col="gray")
    # abline(h=epsilon,col="red")
    
    xlab=ifelse(match.by=="Xt", "Diff in TF expression Yj", "Diff in peak levels Xt")
    plot(dV, dYg, pch=16,col='white',xlab=xlab, ylab="Diff in Gene Expression Yg", 
         main=paste('Model5 match by',match.by,'level 1\n',betahat.spearman.result.str), cex.main=1)
    # radius=min((max(dV)-min(dV))/100, (max(dYg)-min(dYg))/100)
    # for(t in 1:length(dV)){
    #   floating.pie(dV[t],dYg[t],x=c(1,1),radius=radius,
    #                col=c(matchedcol1[t],matchedcol2[t]), startpos=pi/2, border=FALSE)
    # }
    points(dV[sel2], dYg[sel2], pch=1, cex=2, lwd=1.2)
    points(dV, dYg, col=matchedcol1, pch=16, cex=2)
    points(dV, dYg, col='white', pch=1, cex=1, lwd=0.7)
    points(dV, dYg, col=matchedcol2, pch=16, cex=1)
    # points(dV, dYg,  col=matchedcol1,pch=16) # This is to only plot one color
    grid()
    if(sum(sel2)>=MIN.POINTS.FOR.BETAHATCOR){
      slope=lm(dYg[sel2]~dV[sel2])$coefficients[2]
      #if(!is.na(slope)) abline(0, slope, col="blue")
    }
    if(sum(sel1&sel2)>MIN.NONZERO.FOR.GAMMAHAT2){
      remainder1=residuals(lm(dYg~dV))
      remainder2=residuals(lm(meanU*dV~dV))
      ylab=ifelse(match.by=="Xt", "Partial residuals dYg on dYj", "Partial residuals dYg on dXt") # "Partial residuals Yg on V"
      xlab=ifelse(match.by=="Xt", "Partial residuals Xt*dYj on dYj", "Partial residuals Yj*dXt on dXt") # "Partial residuals U*V on V"
      plot(remainder2,remainder1, col='white',pch=16, ylab=ylab, xlab=xlab,
           main=paste('Model5 match by',match.by,'level 2\n', gammahat2.result.str), cex.main=1)
      # radius=min((max(remainder1)-min(remainder1))/100, (max(remainder2)-min(remainder2))/100)
      # for(t in 1:length(remainder2)){
      #   floating.pie(remainder2[t],remainder1[t],x=c(1,1),radius=radius,
      #                col=c(matchedcol1[t],matchedcol2[t]), startpos=pi/2, border=FALSE)
      # }
      points(remainder2, remainder1, col=matchedcol1, pch=16, cex=2)
      points(remainder2, remainder1, col='white', pch=1, cex=1, lwd=0.7)
      points(remainder2, remainder1, col=matchedcol2, pch=16, cex=1)
      grid()
      #abline(0, gammahat2, col="blue")
    } else{
      plot(0,0, cex=0.001, main=gammahat2.result.str)
    }
  }
  if(!is.null(plot.dir) && !is.null(file.name)){
    dev.off()
  }
  
  list(matchedCells=matchedCells, match.by=match.by, 
       betahat.correst=betahat.correst, betahat.pval.rhotest=betahat.pval.rhotest,
       gammahat2=gammahat2, gammahat.pval2=gammahat.pval2)
  
}
