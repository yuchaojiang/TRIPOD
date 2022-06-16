### Transform pvalues to make them easier to plot
transformPvalues<-function(pvals){
  sqrt(-log((pvals)))
}
