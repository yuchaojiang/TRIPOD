getPeakGenePairsForSingleGene <- function(
	xymats,
	coef.name,
	pval.name,
	sign,
	fdr.thresh
) {
	gene.name <- xymats$gene.name
  X <- xymats$X
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
  	peak_num = which.rej,
    peak = colnames(X)[which.rej],
    coef = coef[which.rej],
    pval = pval[which.rej],
    adj = adj[which.rej])
  rownames(results) <- NULL
  return(results)
}
