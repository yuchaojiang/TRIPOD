fitModel<-function(xymats, modelName, do.pred.plot=NULL, do.plot=NULL, # fdr.thresh=NULL, 
                   match.by=NULL, plot.dir=NULL, VERBOSE=NULL, do.capValues=NULL,
                   log=NULL, Delta=NULL){
  
  if(is.null(VERBOSE)){VERBOSE=TRUE}
  if(is.null(Delta)){Delta=0.2}
  ###############################################################
  ### Let's for model variants a, b, and c first
  ### The difference is how Xit is specified
  ### And this applies to models 1-4 the same way 
  ###############################################################
  if(grepl('a',modelName)){ # Model a: keep all peaks as are
    xymats$Xit=xymats$Xit
  } else if (grepl('b',modelName)){ # Model b: remove peaks in other genes
    xymats$peak.gr.g=xymats$peak.gr.g[xymats$peak.filter]
    xymats$Xit=xymats$Xit[,xymats$peak.filter, drop=FALSE]
    xymats$peakxmotif.g=xymats$peakxmotif.g[xymats$peak.filter,, drop=FALSE]
    TF.filter=apply(xymats$peakxmotif.g,2,sum)>0
    xymats$Y.TF=xymats$Y.TF[,TF.filter, drop=FALSE]
    xymats$TF.g=xymats$TF.g[TF.filter]
    xymats$peakxmotif.g=xymats$peakxmotif.g[,TF.filter, drop=FALSE]
  } else{ # Model c: peaks overlapping genes not removed but weighted based on expression of the overlapped genes
    xymats$Xit=xymats$Xit/xymats$Y.peak # # Model c: remove peaks in other genes
  }
	if(!is.null(log)){
    if(all(log %in% c("Y", "Xit", "Y.peak", "Y.TF"))){
      if(any(log=="Y")){
        xymats$Y=log(xymats$Y+1)
      }
      if(any(log=="Xit")){
        xymats$Xit=log(xymats$Xit+1)
      }
      if(any(log=="Y.peak")){
        xymats$Y.peak=log(xymats$Y.peak+1)
      }
      if(any(log=="Y.TF")){
        xymats$Y.TF=log(xymats$Y.TF+1)
      }
    } else {
      stop('The log argument has to be "Y", "Xit", "Y.peak", or "Y.TF".')
    }
  }
  if(is.null(do.capValues)) {do.capValues=TRUE}
  if(do.capValues){ # Cap values
    xymats$Y=capValues(xymats$Y)
    xymats$Xit=apply(xymats$Xit, 2, capValues)
    xymats$Y.peak=apply(xymats$Y.peak, 2, capValues)
    xymats$Y.TF=apply(xymats$Y.TF, 2, capValues)
  }
  # if(is.null(fdr.thresh)){fdr.thresh=0.01}
  if(is.null(do.pred.plot)){do.pred.plot=FALSE}
  if(is.null(do.plot)){do.plot=FALSE}

  ###############################################################
  ### We include model -1 and model 0 as comparisons
  ### Model -1: take the sum of accessibility of all peaks (one vector)
  ### Model 0: LASSO on peaks
  ### Model 1: gene ~ peak, gene ~ TF, TF ~ peak, 
  ###          and gene ~ peak * TF without main terms; marginal model
  ### Model 2: gene ~ peak + TF without interaction, one trio at a time
  ### Model 3: LASSO on peak x TF
  ### Model 4: Fixed TF and fixed peak, regression with peak, TF, and peakxTF
  ### Model 5: potential outcomes model
  ###############################################################
  
  if(grepl('^-1',modelName)){
    ######## MODEL -1: Y on sum of accessibility
    X=as.matrix(apply(xymats$Xit, 1, sum))
    if(do.pred.plot){
      plot(xymats$Y, X, col=xymats$wnn.celltype.col, xlab='Gene expression', 
           ylab=paste('Sum of ATAC peaks'),
           pch=16, main=paste('Gene activity:',gene.name, xymats$ext.upstream/1000, 'Kb up and downstream'))
      legend('bottomright',paste('r =',round(cor(xymats$Y, X, method = 'pearson'),3)), bty = 'n')
      # text(xymats$Y, X, xymats$wnn.celltype, col='gray', cex=0.6)
      # points(xymats$Y, X, col=xymats$wnn.celltype.col, pch=16)
    }
    xymats$LOO.r=round(cor(xymats$Y, X, method = 'pearson'),3)
  }
  if(grepl('^0',modelName)){
    if(ncol(xymats$Xit)>1){ # at least two columns for LASSO
      lasso.mod <- glmnet(xymats$Xit, xymats$Y, alpha = 1)  # 10-fold cross-validation to find the optimal lambda
      set.seed(1) 
      cv.out <- cv.glmnet(xymats$Xit, xymats$Y, alpha = 1)
      # plot(cv.out)
      bestlam <- cv.out$lambda.min
      if(max(lasso.mod$df)>1){
        bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df>1])) # Less penalty
      }
      # With the selected lambda, leave-one-out prediction
      Y.pred=rep(NA, length(xymats$Y))
      for(i in seq(length(xymats$Y))){
        fit <- glmnet(xymats$Xit[-i,], xymats$Y[-i], alpha=1, lambda=bestlam)
        Y.pred[i] <- predict(fit, newx=t(xymats$Xit[i,]))
      }
      if(do.pred.plot){
        plot(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, xlab='Gene expression', ylab='Peak LASSO LOO predicted gene expression', 
            pch=16, main=paste('Peak LASSO'))
        legend('bottomright',paste('r =',round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)), bty = 'n')
        # text(xymats$Y, Y.pred, xymats$wnn.celltype, col='gray', cex=0.6)
        # points(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, pch=16)
      }
      xymats$LOO.r=round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)
      # To get the estimated beta's, we will use all the data points
      beta.pred.g=predict(lasso.mod, s = bestlam, type = "coefficients")[-1]
      betas=as.matrix(beta.pred.g)
      betas[betas==0]=NA
      xymats$betas=betas
    } else{
      Y.pred=rep(NA, length(xymats$Y))
      for(i in seq(length(xymats$Y))){
        fit <- lm(xymats$Y[-i]~xymats$Xit[-i,])
        Y.pred[i] <- xymats$Xit[i,]*fit$coefficients[2]+fit$coefficients[1]
      }
      if(do.pred.plot){
        plot(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, xlab='Gene expression', ylab='Peak LM LOO predicted gene expression', 
             pch=16, main=paste('Model', modelName,': only one peak; no LASSO.'))
        legend('bottomright',paste('r =',round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)), bty = 'n')
        # text(xymats$Y, Y.pred, xymats$wnn.celltype, col='gray', cex=0.6)
        # points(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, pch=16)
      }
      xymats$LOO.r=round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)
    }
  }
  
  if(grepl('^3',modelName)){
    X=xymats$Xit[,xymats$nonzero.peakxmotif.g[,1], drop=FALSE]*xymats$Y.TF[,xymats$nonzero.peakxmotif.g[,2],drop=FALSE]
    peakTF.filter=apply(X, 2, sum)>0 # additional QC: X_it * Y_ij can be zero across all i even though they themselves are non-zeros 
    nonzero.peakxmotif.g=xymats$nonzero.peakxmotif.g[peakTF.filter,,drop=FALSE]
    X=X[,peakTF.filter, drop=FALSE]
    
    if(ncol(X)>1){ # at least two columns for LASSO{
      lasso.mod <- glmnet(X, xymats$Y, alpha = 1)  # 10-fold cross-validation to find the optimal lambda
      set.seed(1) 
      cv.out <- cv.glmnet(X, xymats$Y, alpha = 1)
      # plot(cv.out)
      bestlam <- cv.out$lambda.min
      bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df>1]))# here is a threshold to result in not too much penalty: at least three peaks
      # With the selected lambda, leave-one-out prediction
      Y.pred=rep(NA, length(xymats$Y))
      for(i in seq(length(xymats$Y))){
        fit <- glmnet(X[-i,], xymats$Y[-i], alpha=1, lambda=bestlam)
        Y.pred[i] <- predict(fit, newx=t(X[i,]))
      }
      if(do.pred.plot){
        plot(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, xlab='Gene expression', ylab='Peak-TF LASSO LOO predicted gene expression',
            pch=16, main=paste('Peak-TF LASSO'))
        legend('bottomright',paste('r =',round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)), bty = 'n')
        # text(xymats$Y, Y.pred, xymats$wnn.celltype, col='gray', cex=0.6)
        # points(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, pch=16)
      }
      xymats$LOO.r=round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)
      # To get the estimated beta's, we will use all the data points
      beta.pred.g=predict(lasso.mod, s = bestlam, type = "coefficients")[-1]
      betas=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
      for(i in which(beta.pred.g!=0)){
        betas[xymats$nonzero.peakxmotif.g[i,1],xymats$nonzero.peakxmotif.g[i,2]]=beta.pred.g[i]
      }
      xymats$betas=betas
    } else if(ncol(X)==1){
      Y.pred=rep(NA, length(xymats$Y))
      for(i in seq(length(xymats$Y))){
        fit <- lm(xymats$Y[-i]~X[-i,])
        Y.pred[i] <- X[i,]*fit$coefficients[2]+fit$coefficients[1]
      }
      if(do.pred.plot){
        plot(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, xlab='Gene expression', ylab='Peak-TF LM LOO predicted gene expression', 
             pch=16, main=paste('Model', modelName,': only one peak-TF; no LASSO.'))
        legend('bottomright',paste('r =',round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)), bty = 'n')
        # text(xymats$Y, Y.pred, xymats$wnn.celltype, col='gray', cex=0.6)
        # points(xymats$Y, Y.pred, col=xymats$wnn.celltype.col, pch=16)
      }
      xymats$LOO.r=round(cor(xymats$Y, as.numeric(Y.pred), method = 'pearson'),3)
    } else{
      xymats$LOO.r=NA
    }
  }
  
  if(grepl('^1',modelName)){ # MODEL 1: marginal associations
    nTFs=ncol(xymats$Y.TF)
    npeaks=ncol(xymats$Xit)
    betaXits <- pvalXits <- betaXits.corr <- betaXits.spearman <- pvalXits.spearman <-
    	rep(NA, npeaks)
    betaYijs <- pvalYijs <- betaYijs.corr <- betaYijs.spearman <- pvalYijs.spearman <-
    	rep(NA, nTFs)
    for(t in 1:npeaks){
      Xit=xymats$Xit[,t]
      Yig=xymats$Y
      if(sum(Xit>0)>4){
        res=lm(Yig~1+Xit)
        if(any(is.na(res$coefficients))) next
        betaXits[t]=res$coefficients[2]
        pvalXits[t]=summary(res)$coefficients[2,4]
        betaXits.corr[t] <- cor(Xit, Yig)
        options(warn = -1)
        spearman <- cor.test(Xit, Yig, method = "spearman")
        options(warn = 0)
        betaXits.spearman[t] <- spearman$estimate
        pvalXits.spearman[t] <- spearman$p.value
      }
    }
    for(j in 1:nTFs){
      Yig=xymats$Y
      Yij=xymats$Y.TF[,j]
      if(sum(Yij>0)>4){
        res=lm(Yig~1+Yij)
        if(any(is.na(res$coefficients))) next
        betaYijs[j]=res$coefficients[2]
        pvalYijs[j]=summary(res)$coefficients[2,4]
        betaYijs.corr[j] <- cor(Yij, Yig)
        options(warn = -1)
        spearman <- cor.test(Yij, Yig, method = "spearman")
        options(warn = 0)
        betaYijs.spearman[j] <- spearman$estimate
        pvalYijs.spearman[j] <- spearman$p.value
      }
    }

    gammas <- pvalgammas <- betaYijsOnXits <- pvalBetaYijsOnXits <-
    	gammas.corr <- betaYijsOnXits.corr <- 
    	gammas.spearman <- pvalgammas.spearman <-
    	betaYijsOnXits.spearman <- pvalBetaYijsOnXits.spearman <- 
    	matrix(nrow = nrow(xymats$peakxmotif.g), ncol = ncol(xymats$peakxmotif.g), data = NA)
    
    for(j in 1:nTFs){
      # cat(j," ")
      Xit=xymats$Xit
      Yig=xymats$Y
      Yij=xymats$Y.TF[,j]
      if(sum(Yij>0)>4){  # arbitrary check that not too many zero Y's, so that there is enough degrees of freedom.
        for(t in which(xymats$peakxmotif.g[,j]==1)){
          if(all(Xit[,t]*Yij==0)) next  # Need to check: X_it * Y_ij can be zero across all i even though they themselves are non-zeros 
          if(sum(Xit[,t]>0)<=4) next
          res=lm(Yig~1+Xit[,t]:Yij) # Dropping main effects
          if(any(is.na(res$coefficients))) next
          gammas[t,j]=res$coefficients[2]
          pvalgammas[t,j]=summary(res)$coefficients[2,4]
          
          lmfit=summary(lm(Yij~Xit[,t]))
          betaYijsOnXits[t,j] = lmfit$coefficients[2,1]
          pvalBetaYijsOnXits[t,j]=lmfit$coefficients[2,4]
          
          # marginal terms
          # betaXits.corr[t,j]=cor(Xit[,t], Yig)
          # betaYijs.corr[t,j]=cor(Yij, Yig)
          gammas.corr[t,j]=cor(Xit[,t]*Yij, Yig)
          betaYijsOnXits.corr[t,j]=cor(Xit[,t], Yij)
          options(warn = -1)
          spearman <- cor.test(Xit[,t]*Yij, Yig, method = "spearman")
          options(warn = 0)
          gammas.spearman[t,j] <- spearman$estimate
          pvalgammas.spearman[t,j] <- spearman$p.value
          options(warn = -1)
          spearman <- cor.test(Xit[,t], Yij, method = "spearman")
          options(warn = 0)
          betaYijsOnXits.spearman[t,j] <- spearman$estimate
          pvalBetaYijsOnXits.spearman[t,j] <- spearman$p.value
        }
      }
    }
    xymats$betaXits=betaXits
    xymats$pvalXits=pvalXits
    xymats$betaYijs=betaYijs
    xymats$pvalYijs=pvalYijs
    xymats$gammas=gammas
    xymats$pvalgammas=pvalgammas
    xymats$betaYijsOnXits=betaYijsOnXits
    xymats$pvalBetaYijsOnXits=pvalBetaYijsOnXits
    xymats$betaXits.corr=betaXits.corr
    xymats$betaYijs.corr=betaYijs.corr
    xymats$gammas.corr=gammas.corr
    xymats$betaYijsOnXits.corr=betaYijsOnXits.corr
    xymats$betaXits.spearman=betaXits.spearman
    xymats$pvalXits.spearman=pvalXits.spearman
    xymats$betaYijs.spearman=betaYijs.spearman
    xymats$pvalYijs.spearman=pvalYijs.spearman
    xymats$gammas.spearman=gammas.spearman
    xymats$pvalgammas.spearman=pvalgammas.spearman
    xymats$betaYijsOnXits.spearman=betaYijsOnXits.spearman
    xymats$pvalBetaYijsOnXits.spearman=pvalBetaYijsOnXits.spearman
  }
  
  if(grepl('^2',modelName)){
    nTFs=ncol(xymats$Y.TF)
    betaXits=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    betaYijs=betaXits
    pvalXits=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    pvalYijs=pvalXits
    for(j in 1:nTFs){
      # cat(j," ")
      Xit=xymats$Xit
      Yig=xymats$Y
      Yij=xymats$Y.TF[,j]
      if(sum(Yij>0)>4){  # arbitrary check that not too many zero Y's, so that there is enough degrees of freedom.
        for(t in which(xymats$peakxmotif.g[,j]==1)){
          if(all(Xit[,t]*Yij==0)) next  # Need to check: X_it * Y_ij can be zero across all i even though they themselves are non-zeros 
          if(sum(Xit[,t]>0)<=4) next
          Delta.val = max(quantile(Xit[,t], Delta, na.rm=TRUE), 0.1)
          sel.Xit= Xit[,t]>=Delta.val # only carry out Yij testing if Xit is greater than a threshold
          res=lm(Yig[sel.Xit]~1+Xit[sel.Xit,t]+Yij[sel.Xit])
          if(any(is.na(res$coefficients))) next
          betaYijs[t,j]=res$coefficients[3]
          pvalYijs[t,j]=summary(res)$coefficients[3,4]
          
          Delta.val = max(quantile(Yij, Delta, na.rm=TRUE), 0.1)
          sel.Yij= Yij>=Delta.val # only carry out Xit testing if Yij is greater than a threshold
          res=lm(Yig[sel.Yij]~1+Xit[sel.Yij,t]+Yij[sel.Yij])
          if(any(is.na(res$coefficients))) next
          betaXits[t,j]=res$coefficients[2]
          pvalXits[t,j]=summary(res)$coefficients[2,4]
        }
      }
    }
    xymats$betaXits=betaXits
    xymats$betaYijs=betaYijs
    xymats$pvalXits=pvalXits
    xymats$pvalYijs=pvalYijs
    xymats$Delta=Delta
    #xymats$which.reject=benjaminiHochsbergMatrix(pvalXits, fdr.thresh, do.plot=FALSE) | benjaminiHochsbergMatrix(pvalYijs, fdr.thresh, do.plot=FALSE) 
  }
  
  if(grepl('^4',modelName)){
    nTFs=ncol(xymats$Y.TF)
    betaXits=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    gammas=betaYijs=betaXits
    pvalXits=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    pvalgammas=pvalYijs=pvalXits
    for(j in 1:nTFs){
      # cat(j," ")
      Xit=xymats$Xit
      Yig=xymats$Y
      Yij=xymats$Y.TF[,j]
      if(sum(Yij>0)>4){  # arbitrary check that not too many zero Y's, so that there is enough degrees of freedom.
        for(t in which(xymats$peakxmotif.g[,j]==1)){
          if(all(Xit[,t]*Yij==0)) next  # Need to check: X_it * Y_ij can be zero across all i even though they themselves are non-zeros 
          if(sum(Xit[,t]>0)<=4) next
          res=lm(Yig~1+Xit[,t]+Yij+Xit[,t]*Yij)
          if(any(is.na(res$coefficients))) next
          betaXits[t,j]=res$coefficients[2]
          betaYijs[t,j]=res$coefficients[3]
          gammas[t,j]=res$coefficients[4]
          pvalXits[t,j]=summary(res)$coefficients[2,4]
          pvalYijs[t,j]=summary(res)$coefficients[3,4]
          pvalgammas[t,j]=summary(res)$coefficients[4,4]
        }
      }
    }
    xymats$betaXits=betaXits
    xymats$betaYijs=betaYijs
    xymats$gammas=gammas
    xymats$pvalXits=pvalXits
    xymats$pvalYijs=pvalYijs
    xymats$pvalgammas=pvalgammas
  }
  
  if(grepl('^5',modelName)){  # matching-based potential outcomes model
    betahat.pvals.rhotest=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    betahat.corrests=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    gammahat2s=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    gammahat2.pvals=matrix(nrow=nrow(xymats$peakxmotif.g), ncol=ncol(xymats$peakxmotif.g), data=NA)
    if(VERBOSE) cat("Performing matched tests ... \n")
    for(t in 1:ncol(xymats$Xit)){
      if(VERBOSE) cat("Peak #", t,"out of", ncol(xymats$Xit),": ",sum(xymats$peakxmotif.g[t,]),"TF with motifs. Motif:")
      for(j in which(xymats$peakxmotif.g[t,]==1)){
        if(VERBOSE) cat(j," ")
        
        Xt = xymats$Xit[,t]   # peak accessibility * motif presence
        Yg = xymats$Y        # gene expression
        Yj = xymats$Y.TF[,j]  # transcription factor expression
        
        file.name=paste(c("matchedTest_peak", t, "_TF", j, ".png"), collapse="") ### modified
        res=matchedTest(Xt, Yj, Yg, Z=xymats$Z, match.by=match.by, ### modified 
                        CAP.AT.QUANTILE=0.02, Delta=Delta, 
                        col=wnn.celltype.col, do.plot=do.plot,
                        plot.dir=plot.dir, file.name=file.name,
                        partition.screen=TRUE, verbose=FALSE)
        
        betahat.corrests[t,j]=res$betahat.correst
        betahat.pvals.rhotest[t,j]=res$betahat.pval.rhotest
        gammahat2s[t,j]=res$gammahat2
        gammahat2.pvals[t,j]=res$gammahat.pval2
      }
      if(VERBOSE) cat("\n")
    }
    xymats$betahat.corrests=betahat.corrests
    xymats$betahat.pvals.rhotest=betahat.pvals.rhotest
    xymats$gammahat2s=gammahat2s
    xymats$gammahat2.pvals=gammahat2.pvals
    #xymats$reject.beta=benjaminiHochsbergMatrix(betahat.pvals.rhotest, fdr.thresh=fdr.thresh, do.plot=FALSE)
    #xymats$reject.gamma=benjaminiHochsbergMatrix(gammahat2.pvals, fdr.thresh=fdr.thresh, do.plot=FALSE)
    xymats$match.by=match.by
  }
  xymats
}
