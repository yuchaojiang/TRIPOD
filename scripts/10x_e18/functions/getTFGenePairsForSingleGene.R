getTFGenePairsForSingleGene <- function(
  xymats,
	coef.name,
	pval.name,
	sign,
	fdr.thresh
) {
	gene.name <- xymats$gene.name
  Y.TF <- xymats$Y.TF
  coef <- xymats[[coef.name]]
  pval <- xymats[[pval.name]]
  adj <- benjaminiHochsbergVectorAdjust(pval, fdr.thresh = fdr.thresh)$pval.c.adj
  if (sign == "positive"){
    which.rej <- which(coef > 0 & adj < fdr.thresh)
  } else if (sign == "negative"){
    which.rej <- which(coef < 0 & adj < fdr.thresh)
  }
  results <- data.frame(
	  gene = rep(gene.name, length(which.rej)),
  	TF_num = which.rej,
    TF = colnames(Y.TF)[which.rej],
    coef = coef[which.rej],
    pval = pval[which.rej],
    adj = adj[which.rej])
  rownames(results) <- NULL
  results <- results[results$gene != results$TF, ]
  return(results)
}
