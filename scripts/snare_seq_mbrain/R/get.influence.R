# object: a Seurat object
get.influence=function(xymats.all, gene.name, peak.num, TF.num, 
                       object, dend, reduction=NULL, model=NULL,
                       plot.dir=NULL, to.test=NULL, generate.pdf=NULL, 
                       partition.screen=NULL, wnn.neighbors=NULL, cell.neighbors=NULL){
  if(is.null(generate.pdf)) generate.pdf=TRUE
  if(is.null(partition.screen)) partition.screen=TRUE
  if(is.null(plot.dir)) plot.dir=getwd()
  if(is.null(to.test)) to.test=c('WNN','WNN.Sampling','CellType','Tree','Footprint')
  if(is.null(model)) model='4' # Can be chosen between model 4 and model 5; by default, model 4
  if('WNN.Sampling'%in%to.test & (is.null(wnn.neighbors) | is.null(cell.neighbors))) 
    stop('Need to specify WNN.neighbors or remove WNN.Sampling from to.test')
  Yig=xymats.all[[gene.name]]$xymats4$Y
  Xit=xymats.all[[gene.name]]$xymats4$Xit[,peak.num]
  Yij=xymats.all[[gene.name]]$xymats4$Y.TF[,TF.num]
  
  if(is.null(reduction)){reduction='wnn.umap'}  
  if(generate.pdf) {pdf(file =paste(plot.dir,'/', gene.name,'_',ext.upstream/1000,'kb_peak_',peak.num,'_TF_',TF.num,'_influence.pdf',sep=''), width=16, height=8)}
  if(partition.screen) {par(mfrow=c(2,4))}
  
  if('WNN' %in% to.test){
    # Remove one WNN at a time
    plot.new() 
    vps <- baseViewports()
    pushViewport(viewport(layout = grid.layout(2, 4)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    vp1 <-plotViewport(c(1,0.5,0,1))
    plotDFFIT=plot_cooks_DFFIT(Yig, Xit, Yij)
    p12=plotDFFIT
    print(p12[[1]], vp=vp1)
    popViewport()
    
    plot.new() 
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    print(p12[[2]], vp=vp1)
    popViewport()
  }
  
  if('WNN.Sampling' %in% to.test){
    delta.coeff.pval=matrix(nrow=length(wnn.celltype),ncol=5)
    colnames(delta.coeff.pval)=names(test.influential(Yig, Xit, Yij, 1, nsamp=10))
    rownames(delta.coeff.pval)=names(wnn.celltype)
    for(i in 1:length(wnn.celltype)){
      wnn.rm=wnn.neighbors[i,] # Remove k-nearest neighbors of the WNN. Can be inferred from umap coordinates + pseudotime
      delta.coeff.pval[i,]=test.influential(Yig, Xit, Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)
    }
    d = data.frame(obs=1:length(wnn.celltype),
                   cd=-log(delta.coeff.pval[,5],10),
                   txt=wnn.celltype)
    d$txt[d$cd<(-log(0.05,10))]=NA
    myplots=ggplot(d, aes(x = obs, y = cd, label = txt, ymin = min(cd), ymax = cd)) + 
      geom_linerange(colour = wnn.celltype.col) + 
      geom_point(shape = 16, colour = wnn.celltype.col) + geom_hline(yintercept = -log(0.05,10), colour = "black", linetype='dotted')+
      xlab("WNN") + ylab("-log(p-value)") + ggtitle(paste("WNN-specific p-values for ",colnames(delta.coeff.pval)[5],sep='')) + 
      geom_text(size = 3, family = "serif", fontface = "italic", colour =wnn.celltype.col, na.rm = TRUE) 
    xymats.all[[gene.name]]$xymats4$wnn.sampling.pval=delta.coeff.pval
    
    plot.new() 
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
    print(myplots, vp=vp1)
    popViewport()
    
    wnn.samp.logpval=-log(delta.coeff.pval[object$seurat_clusters,5],10)
    wnn.samp.logpval.avg=rep(NA, length(wnn.samp.logpval))
    for(tt in 1:length(wnn.samp.logpval.avg)){
      wnn.samp.logpval.avg[tt]=mean(wnn.samp.logpval[cell.neighbors[tt,]])
    }
    object$wnn.samp.logpval.avg=wnn.samp.logpval.avg
    p.feature=FeaturePlot(object, features='wnn.samp.logpval.avg', reduction=reduction, label=TRUE, label.size = 2)+
      ggtitle('WNN-specific p-values for Yig from sampling') +
      theme(plot.title = element_text(size = 10))
    plot.new() 
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
    print(p.feature, vp=vp1)
    popViewport()
  }
  
  if('CellType' %in% to.test){
    # Remove one cell type at a time
    delta.coeff.pval=matrix(nrow=length(unique(wnn.celltype)),ncol=5)
    colnames(delta.coeff.pval)=names(test.influential(Yig, Xit, Yij, 1, nsamp=10))
    rownames(delta.coeff.pval)=unique(wnn.celltype)
    for(wnn.rm.celltype in unique(wnn.celltype)){
      # cat(wnn.rm.celltype,'\t')
      wnn.rm=which(wnn.celltype==wnn.rm.celltype)
      delta.coeff.pval[wnn.rm.celltype,]=test.influential(Yig, Xit, Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)
    }
    delta.coeff.pval
    
    myplots <- list()
    for(i in 1:5){
      d = data.frame(obs=1:length(unique(wnn.celltype)),
                     cd=-log(delta.coeff.pval[,i],10),
                     txt=unique(wnn.celltype))
      d$txt[d$cd<(-log(0.05,10))]=NA
      d=d[levels(object$celltype),]
      d$obs=1:length(d$obs)
      myplots[[i]]=ggplot(d, aes(x = obs, y = cd, label = txt, ymin = min(cd), ymax = cd)) + 
        geom_linerange(colour = unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))]) + 
        geom_point(shape = 16, colour = unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))]) + geom_hline(yintercept = -log(0.05,10), colour = "black", linetype='dotted')+
        xlab("Cell types") + ylab("-log(p-value)") + ggtitle(paste("Cell-type-specific p-values for ",colnames(delta.coeff.pval)[i],sep='')) + 
        geom_text(size = 3, family = "serif", fontface = "italic", colour =unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))], na.rm = TRUE) 
    }
    p3=myplots[[5]] 
    
    plot.new() 
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
    print(p3, vp=vp1)
    popViewport()
    
    wnn.celltype.sig= rownames(delta.coeff.pval)[benjaminiHochsbergVectorAdjust(delta.coeff.pval[,5])$pval.c.adj<=0.05]
    wnn.sig=which(wnn.celltype%in%wnn.celltype.sig)
    xymats.all[[gene.name]]$xymats4$wnn.sig=wnn.sig
  }
  
  if('Tree' %in% to.test){
    subtrees <- partition_leaves(dend) # Subtrees, each corresponding to a split
    # Remove cell types at each split
    delta.coeff.pval=matrix(nrow=length(subtrees),ncol=5, data=NA)
    colnames(delta.coeff.pval)=names(test.influential(Yig, Xit, Yij, 1, nsamp=10))
    rownames(delta.coeff.pval)=1:length(subtrees)
    for(i in 2:length(subtrees)){
      # cat(i,'\t')
      wnn.rm=which(wnn.celltype%in%subtrees[[i]])
      delta.coeff.pval[i,]=test.influential(Yig, Xit, Yij, wnn.rm, plot.histogram = FALSE, nsamp=1000)
    }
    
    par(mar=c(6.1, 4.1, 2.1, 2.1))
    i=5
    nodes_col=rep('lightgray', length(subtrees))
    nodes_col[which(benjaminiHochsbergVectorAdjust(delta.coeff.pval[,i])$pval.c.adj<(0.05))]='red'
    nodes_col[which(is.na(delta.coeff.pval[,i]))]=1
    nodes_pch=rep(19, length(subtrees))
    nodes_pch[which(is.na(delta.coeff.pval[,i]))]=4
    wnn.celltype.col.unique=unique(wnn.celltype.col)
    names(wnn.celltype.col.unique)= unique(wnn.celltype)
    dend %>% set("nodes_pch", nodes_pch) %>% 
      set("nodes_col", nodes_col)%>% set("nodes_cex", 2) %>% set("labels_col", wnn.celltype.col.unique[dend %>% labels]) %>% # change color
      set("labels_cex", 1) %>% # Change size
      plot(main = paste("Hierarchical clustering with testing of",colnames(delta.coeff.pval)[i])) # plot
    legend('topright', col=c(1,2,'lightgray'), pch=c(4,19,19), legend=c('NA','Significant','Insignificant'), bty = 'n', cex=0.8)
    
    par(mar=c(5.1, 4.1, 4.1, 2.1))
  }
  
  if('Footprint' %in% to.test){
    vp1 <-plotViewport(c(0.1,0.1,0,1))
    # This is to plot the TF across all of its binding sites
    p5=PlotFootprint(object, features = TF.name, group.by='celltype', show.expected = FALSE) 
    p5=p5+ggtitle(paste('DNA footprinting for',TF.name))
    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3:4))
    print(p5, vp=vp1)
    popViewport()
  }
  if(generate.pdf) dev.off()
  return(xymats.all)
}
