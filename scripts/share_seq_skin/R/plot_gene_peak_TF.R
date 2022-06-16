# object: a Seurat object
# group.by: for UMAP plots
plot_gene_peak_TF=function(object, xymats, sel.t, jj,
	                         group.by=NULL,
                           plot.dir=NULL, to.plot=NULL, do.capValues=NULL,
                           partition.screen=NULL, generate.pdf=NULL,
                           plot.influential.points=NULL,
                           featurePlot=NULL, reduction=NULL){
  if(is.null(to.plot)){
    to.plot=c('umap','marginal', 'model1.scatter', 'model2.scatter','model4.scatter')
  }
  if(is.null(group.by)){
    object[["default.ident"]]=Idents(object)
    group.by='default.ident'
  }
  if(is.null(plot.influential.points)) plot.influential.points=FALSE
  if(plot.influential.points) wnn.sig=xymats$wnn.sig
  if(is.null(partition.screen)){partition.screen=TRUE}
  if(is.null(generate.pdf)){generate.pdf=TRUE}
  if(is.null(do.capValues)) {do.capValues=TRUE}
  if(do.capValues){ # Cap values
    xymats$Y=capValues(xymats$Y)
    xymats$Xit=apply(xymats$Xit, 2, capValues)
    xymats$Y.peak=apply(xymats$Y.peak, 2, capValues)
    xymats$Y.TF=apply(xymats$Y.TF, 2, capValues)
  }
  if(is.null(featurePlot)) {featurePlot='Violin'}
  if(is.null(reduction)) {reduction='wnn.umap'}
  if('umap'%in% to.plot){
    # UMAP + violin plot for the gene of interest
    p2=FeaturePlot(object, features = paste("sct_",gene.name,sep=''), reduction = reduction, label=T, label.size = 2, cols = c("lightgrey", "darkblue"))+ggtitle(paste('Gene expr.',gene.name))
    DefaultAssay(object)='ATAC'
    p3=FeaturePlot(object, features = rownames(xymats$peakxmotif.g)[sel.t], reduction = reduction, label=T, label.size = 2, cols = c("lightgrey", "darkgreen"))+ggtitle(paste('Peak access.',rownames(xymats$peakxmotif.g)[sel.t]))
    p4=FeaturePlot(object, features = paste("sct_",xymats$TF.g[jj],sep=''), reduction = reduction, label=T, label.size = 2, col=c('lightgrey','darkred'))+ggtitle(paste('TF expr.',xymats$TF.g[jj], 'motif', names(xymats$TF.g)[[jj]]))
    p5=FeaturePlot(object, features = motifxTF[motifxTF[,2]==xymats$TF.g[jj],1], min.cutoff = 0, cols = c("lightgrey", "darkorange"), reduction = reduction, label=T, label.size = 2)+ggtitle(paste('TF motif access.',motifxTF[motifxTF[,2]==xymats$TF.g[jj],1]))
    
    # violin plot!
    if(featurePlot=='Violin'){
      p6=VlnPlot(object, features = paste("rna_",gene.name, sep=''), slot = "counts", log = TRUE, group.by = group.by) + ggtitle('Gene Expression')+ylab('Expression')
      p7=VlnPlot(object, features = rownames(xymats$peakxmotif.g)[sel.t], slot = "counts", log = TRUE, group.by = group.by)+ggtitle('Peak Accessibility') + ylab('Peak Accessibility')
      p8=VlnPlot(object, features = paste("rna_",xymats$TF.g[jj], sep=''), slot = "counts", log = TRUE, group.by = group.by) + ggtitle('TF Expression')+ylab('Expression')
      p9=VlnPlot(object, features =paste("chromvar_",motifxTF[motifxTF[,2]==xymats$TF.g[jj],1],sep=''), group.by = group.by, pt.size = 0.1) + ggtitle('Motif Accessibility Deviation')+ylab('chromVAR Deviation')
    } else if(featurePlot=='Dot'){
      # p6=DotPlot(object, features = paste("rna_",gene.name, sep=''),  cols = c("lightgrey", "darkblue"), group.by = group.by) + ggtitle('Gene Expression')+ylab('Expression')
      p6=DotPlot(object, features = paste("sct_",gene.name, sep=''),  cols = c("lightgrey", "darkblue"), group.by = group.by) + ggtitle('Gene Expression')+ylab('Expression')
      p7=DotPlot(object, features = rownames(xymats$peakxmotif.g)[sel.t], cols = c("lightgrey", "darkgreen"), group.by = group.by)+ggtitle('Peak Accessibility') + ylab('Peak Accessibility')
      p8=DotPlot(object, features = paste("sct_",xymats$TF.g[jj], sep=''),  cols=c('lightgrey','darkred'), group.by = group.by) + ggtitle('TF Expression')+ylab('Expression')
      p9=VlnPlot(object, features =paste("chromvar_",motifxTF[motifxTF[,2]==xymats$TF.g[jj],1],sep=''), group.by = group.by, pt.size = 0.1) + ggtitle('Motif Accessibility Deviation')+ylab('chromVAR Deviation')
    }
    
    if(generate.pdf) {ggsave(filename = paste(plot.dir,'/',gene.name,'_',ext.upstream/1000,'kb_peak_',sel.t,'_TF_',jj,'_umap.pdf',sep=''),
           plot = p2+p3+p4+p5+p6+p7+p8+p9+ plot_layout(ncol=4),
           width = 20,
           height = 8)
      } else{
        plot = p2+p3+p4+p5+p6+p7+p8+p9+ plot_layout(ncol=4)
        print(plot)
      }
  }
  
  if('marginal' %in% to.plot){
    if(generate.pdf) {pdf(file= paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',sel.t,'_TF_',jj,'_marginal.pdf',sep=''),
        width=10, height=4)}
    # Now add the colors for the cell types
    if(partition.screen) {par(mfrow=c(1,3))}
    plot(xymats$Xit[,sel.t], xymats$Y, cex.main=1, xlab="X_it", ylab="Y_ig", 
         main=paste('Marginal Gene',gene.name,"expr. VS Peak", sel.t,"access.:\nr =", round(xymats$betaXits.corr[sel.t],2),
                    'p-val =', signif(xymats$pvalXits[sel.t],2)),
         col=wnn.celltype.col, pch=16); grid()
    text(xymats$Xit[,sel.t], xymats$Y, wnn.celltype, col='gray', cex=0.6)
    points(xymats$Xit[,sel.t], xymats$Y, cex.main=1, xlab="X_it", ylab="Y_ig",
           col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(xymats$Xit[wnn.sig,sel.t], xymats$Y[wnn.sig], cex.main=1, pch=0)}
    
    plot(xymats$Y.TF[, jj], xymats$Y, cex.main=1, col=wnn.celltype.col, pch=16,
         xlab="Y_ij", ylab="Y_ig", main=paste('Marginal Gene',gene.name,"expr. VS TF", jj, "expr.:\nr =", round(xymats$betaYijs.corr[jj],2), 
                                              'p-val =', signif(xymats$pvalYijs[jj],2))); grid()
    text(xymats$Y.TF[, jj], xymats$Y, wnn.celltype,col='gray', cex=0.6)
    points(xymats$Y.TF[, jj], xymats$Y, cex.main=1, col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(xymats$Y.TF[wnn.sig, jj], xymats$Y[wnn.sig], cex.main=1,pch=0)}
    
    plot(xymats$Xit[,sel.t], xymats$Y.TF[,jj], cex.main=1, col=wnn.celltype.col, pch=16,
         xlab="X_it", ylab="Y_ij", 
         main=paste('Marginal Peak access. VS TF  expr.:\nr =', round(xymats$betaYijsOnXits.corr[sel.t,jj],2), 
                    'p-val =', signif(xymats$pvalBetaYijsOnXits[sel.t,jj],2))); grid()
    text(xymats$Xit[,sel.t],xymats$Y.TF[,jj], wnn.celltype, col='gray', cex=0.6)
    points(xymats$Xit[,sel.t],xymats$Y.TF[,jj], cex.main=1, col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(xymats$Xit[wnn.sig,sel.t],xymats$Y.TF[wnn.sig,jj],cex.main=1, pch=0)}
    if(generate.pdf) {dev.off()}
  }
  
  if('model1.scatter' %in% to.plot){
    if(generate.pdf) {pdf(file= paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',sel.t,'_TF_',jj,'_model1_scatterplot.pdf',sep=''),
        width=4, height=4)}
    plot(xymats$Xit[,sel.t]*xymats$Y.TF[,jj], xymats$Y,cex.main=1, col=wnn.celltype.col, pch=16,
         xlab="X_it * Y_ij", ylab="Y_ig", 
         main=paste('Gene expr. VS (peak access. x TF expr.)\n',
                    'r =', round(cor(xymats$Xit[,sel.t]*xymats$Y.TF[,jj], xymats$Y),4),
                    'p-val =', signif(xymats$pvalgammas[sel.t,jj],2))); grid()
    text(xymats$Xit[,sel.t]*xymats$Y.TF[,jj], xymats$Y, wnn.celltype, col='gray', cex=0.6)
    points(xymats$Xit[,sel.t]*xymats$Y.TF[,jj], xymats$Y,cex.main=1, col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(xymats$Xit[wnn.sig,sel.t]*xymats$Y.TF[wnn.sig,jj], xymats$Y[wnn.sig],cex.main=1, pch=0)}
    if(generate.pdf) dev.off()
  }
  
  if('model2.scatter' %in% to.plot){
    if(generate.pdf) {pdf(file= paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',sel.t,'_TF_',jj,'_model2_scatterplot.pdf',sep=''),
        width=10, height=4)}
    # Now add the colors for the cell types
    if(partition.screen) {par(mfrow=c(1,3))}
    # Partial residual plot to visualize model 2
    # Compute partial residuals, and plot them against each other.
    remainY=residuals(lm(xymats$Y~xymats$Y.TF[,jj]))
    remainX=residuals(lm(xymats$Xit[,sel.t]~xymats$Y.TF[,jj]))
    temp=lm(remainY~remainX)
    partial.p = summary(temp)$coefficients[2,4]
    partial.beta=summary(temp)$coefficients[2,1]
    plot(remainX, remainY, cex.main=1,xlab="X_it partial resid", ylab="Y_ig partial resid", 
         col=wnn.celltype.col, pch=16,
         main=paste("Model 2a, alpha=", format(xymats$betaXits[sel.t,jj], digits=2), ", p-val=", format(xymats$pvalXits[sel.t,jj], digits=2))); grid()
    text(remainX, remainY, wnn.celltype,col='gray', cex=0.6)
    Delta.val = max(quantile(xymats$Y.TF[,jj], xymats$Delta, na.rm=TRUE), 0.1)
    sel.Yij= xymats$Y.TF[,jj]>=Delta.val # only carry out Xit testing if Yij is greater than a threshold
    points(remainX[sel.Yij], remainY[sel.Yij],pch=1)
    points(remainX, remainY, cex.main=1, 
           col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(remainX[wnn.sig], remainY[wnn.sig], cex.main=1, pch=0)}
    
    # Compute partial residuals, and plot them against each other.
    remainY=residuals(lm(xymats$Y~xymats$Xit[,sel.t]))
    remainX=residuals(lm(xymats$Y.TF[,jj]~xymats$Xit[,sel.t]))
    temp=lm(remainY~remainX)
    partial.p = summary(temp)$coefficients[2,4]
    partial.beta=summary(temp)$coefficients[2,1]
    plot(remainX, remainY, cex.main=1,xlab="Y_ij partial resid", ylab="Y_ig partial resid", 
         col=wnn.celltype.col, pch=16,
         main=paste("Model 2a, beta=", format(xymats$betaYijs[sel.t, jj], digits=2), ", p-val=", format(xymats$pvalYijs[sel.t, jj], digits=2))); grid()
    text(remainX, remainY, wnn.celltype,col='gray', cex=0.6)
    Delta.val = max(quantile(xymats$Xit[,sel.t], xymats$Delta, na.rm=TRUE), 0.1)
    sel.Xit= xymats$Xit[,sel.t]>=Delta.val # only carry out Yij testing if Xit is greater than a threshold
    points(remainX[sel.Xit], remainY[sel.Xit],pch=1)
    points(remainX, remainY, cex.main=1, 
           col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(remainX[wnn.sig], remainY[wnn.sig], cex.main=1, pch=0)}

    # Make a third empty plot to hold spot
    plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
    legend=unique(wnn.celltype)[match(levels(object$celltype), unique(wnn.celltype))]
    legend.col=unique(wnn.celltype.col)[match(levels(object$celltype), unique(wnn.celltype))]
    half.length=ceiling(length(legend)/2)
    legend('topleft', legend = legend[1:half.length],
           col=legend.col[1:half.length], pch=16, bty = 'n', cex=1)
    legend('topright', legend = legend[(half.length+1):length(legend)],
           col=legend.col[(half.length+1):length(legend)], pch=16, bty = 'n', cex=1)
    
    if(generate.pdf) dev.off()
  }
  
  if('model4.scatter' %in% to.plot){
    if(generate.pdf) {pdf(file= paste(plot.dir,gene.name,'_',ext.upstream/1000,'kb_peak_',sel.t,'_TF_',jj,'_model4_scatterplot.pdf',sep=''),
        width=4, height=4)}
    # Compute partial residuals, and plot them against each other.
    remainY=residuals(lm(xymats$Y~xymats$Xit[,sel.t]+xymats$Y.TF[,jj]))
    remainX=residuals(lm(xymats$Xit[,sel.t]*xymats$Y.TF[,jj]~xymats$Xit[,sel.t]+xymats$Y.TF[,jj]))
    temp=lm(remainY~remainX)
    partial.p = summary(temp)$coefficients[2,4]
    partial.beta=summary(temp)$coefficients[2,1]
    plot(remainX, remainY, cex.main=1,xlab="X_it * Y_ij partial resid", ylab="Y_ig partial resid", 
         col=wnn.celltype.col, pch=16,
         main=paste("Model 4a, gamma=", format(partial.beta, digits=2), ", p-val=", format(partial.p, digits=2))); grid()
    text(remainX, remainY, wnn.celltype,col='gray', cex=0.6)
    points(remainX, remainY, cex.main=1, 
           col=wnn.celltype.col, pch=16)
    if(plot.influential.points){
      points(remainX[wnn.sig], remainY[wnn.sig], cex.main=1, pch=0)}
    if(generate.pdf) dev.off()
  }
}
