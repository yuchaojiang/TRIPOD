getTriosForSingleGene <- function(
	gene.name,
  xymats.list, 
	fdr.thresh, 
	sign, # "positive" or "negative"
  coef.name,
	pval.name
)	{
	xymats <- xymats.list[[gene.name]]
  peak.g <- rownames(xymats$peakxmotif.g)
  TF.g <- xymats$TF.g
  nonzero.peakxmotif.g <- xymats$nonzero.peakxmotif.g
  coef <- xymats[[coef.name]][nonzero.peakxmotif.g]
  pval <- xymats[[pval.name]][nonzero.peakxmotif.g]
  adj <- benjaminiHochsbergVectorAdjust(pval, fdr.thresh = fdr.thresh)$pval.c.adj
  if (sign == "positive") {
    which.rej <- which(coef > 0 & adj < fdr.thresh)
  } else if (sign == "negative") {
  		which.rej <- which(coef < 0 & adj < fdr.thresh)
  }
  df <- data.frame(gene = rep(gene.name, length(which.rej)),
    peak_num = nonzero.peakxmotif.g[which.rej, 1],
    TF_num = nonzero.peakxmotif.g[which.rej, 2], 
    peak = peak.g[nonzero.peakxmotif.g[which.rej, 1]],
    TF = TF.g[nonzero.peakxmotif.g[which.rej, 2]],
    coef = coef[which.rej],
  	pval = pval[which.rej],
  	adj = adj[which.rej])
  rownames(df) <- NULL
  # remove cases where the TF is the same as the target gene
  df <- df[df$gene != df$TF, ]
  return(df)
}
