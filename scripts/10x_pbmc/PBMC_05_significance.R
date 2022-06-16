
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)
library(ggplot2)
library(presto)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

library(glmnet)
library(pheatmap)
library(fields)
library(qvalue)
library(gplots)
library(patchwork)
library(olsrr)
library(pdftools)
library(ape)
library(dendextend)
library(gplots)
library(nbpMatching)
library(fields)
library(ggvenn)

load('ATAC_RNA_WNN_integrated_high_var.rda')

# How do the p-values compare between models
for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){
  pdf(file=paste0('linear_nonlinear_cor_',gene.name,'.pdf'), width=4, height=4)
  par(mfrow=c(1,1))
  # Gamma from model 4 and model 5X, model 5Y
  pairs(cbind(#model1=sign(c(xymats.all[[gene.name]]$xymats1$gammas))*transformPvalues(c(xymats.all[[gene.name]]$xymats1$pvalgammas)), # model 1
              model4=sign(c(xymats.all[[gene.name]]$xymats4$gammas))*transformPvalues(c(xymats.all[[gene.name]]$xymats4$pvalgammas)), # model 4
              model5X=sign(c(xymats.all[[gene.name]]$xymats5X$gammahat2s))*transformPvalues(c(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals)),  # model 5 match by X
              model5Y=sign(c(xymats.all[[gene.name]]$xymats5Y$gammahat2s))*transformPvalues(c(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals))), # model 5 match by Y
        upper.panel = panel.cor, lower.panel = function(x,y){smoothScatter(x,y,add=T)},
        main=paste(gene.name,'sign(gamma)*sqrt(-log(pval.gamma))'))
  # Alpha (for Xit) from model 2 and model 5Y
  pairs(cbind(model2=sign(c(xymats.all[[gene.name]]$xymats2$betaXits))*transformPvalues(c(xymats.all[[gene.name]]$xymats2$pvalXits)), # model 2
              #model4=sign(c(xymats.all[[gene.name]]$xymats4$betaXits))*transformPvalues(c(xymats.all[[gene.name]]$xymats4$pvalXits)), # model 4
              model5Y=sign(c(xymats.all[[gene.name]]$xymats5Y$betahat.corrests))*transformPvalues(c(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest))), # model 5 match by Y
        upper.panel = panel.cor, lower.panel = function(x,y){smoothScatter(x,y,add=T)},
        main=paste(gene.name,'sign(alpha)*sqrt(-log(pval.alpha))'))
  # Beta (for Yij) from model 2 and model 5X
  pairs(cbind(model2=sign(c(xymats.all[[gene.name]]$xymats2$betaYijs))*transformPvalues(c(xymats.all[[gene.name]]$xymats2$pvalYijs)), # model 2
              #model4=sign(c(xymats.all[[gene.name]]$xymats4$betaYijs))*transformPvalues(c(xymats.all[[gene.name]]$xymats4$pvalYijs)), # model 4
              model5X=sign(c(xymats.all[[gene.name]]$xymats5X$betahat.corrests))*transformPvalues(c(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest))), # model 5 match by X
        upper.panel = panel.cor, lower.panel = function(x,y){smoothScatter(x,y,add=T)},
        main=paste(gene.name,'sign(beta)*sqrt(-log(pval.beta))'))
  dev.off()
}

# Look at a specific gene
gene.name='CCR7'
gene.name='GNLY'
gene.name='FCGR3A'
gene.name='MS4A1'

alpha=0.05

for(gene.name in c('CCR7','GNLY','FCGR3A','MS4A1')){
  
  # gene range
  transcripts.gr.g=transcripts.gr[transcripts.gr$gene_name==gene.name]
  transcripts.ext.gr.g=transcripts.ext.gr[transcripts.ext.gr$gene_name==gene.name]
  peak.gr.g=xymats.all[[gene.name]]$xymats4$peak.gr.g
  
  # Gammas
  Linear.Int.=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats4$pvalgammas, fdr.thresh=alpha,do.plot = FALSE)
  NonLinear.Int.X=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$gammahat2.pvals,fdr.thresh=alpha, do.plot = FALSE)
  NonLinear.Int.Y=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$gammahat2.pvals,fdr.thresh=alpha, do.plot = FALSE)
  x=list(L.Int.=which(Linear.Int.),
         NL.Int.X=which(NonLinear.Int.X),
         NL.Int.Y=which(NonLinear.Int.Y))
  p1=ggvenn(x, 
            fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
            stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
  )+ggtitle(paste0(gene.name,': significant trios'))
  
  # Alphas
  beta.rej2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalXits,fdr.thresh=alpha, do.plot = FALSE) 
  beta.rej5Y=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5Y$betahat.pvals.rhotest,fdr.thresh=alpha, do.plot = FALSE)
  
  x=list(L.Cond.Y=which(beta.rej2),
         NL.Cond.Y=which(beta.rej5Y))
  p2=ggvenn(x, 
            fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
            stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
  )
  
  # Betas
  beta.rej2=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats2$pvalYijs,fdr.thresh=alpha, do.plot = FALSE)
  beta.rej5X=benjaminiHochsbergMatrix(xymats.all[[gene.name]]$xymats5X$betahat.pvals.rhotest,fdr.thresh=alpha, do.plot = FALSE)
  
  x=list(L.Cond.X=which(beta.rej2),
         NL.Cond.X=which(beta.rej5X))
  p3=ggvenn(x, 
            fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
            stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE
  )
  
  p=p1+p2+p3
  ggsave(filename=paste0('linear_nonlinear_venn_',gene.name,'.pdf'),plot = p,
         width = 12, height=3)
  
}
