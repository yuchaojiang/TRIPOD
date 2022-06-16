getPeakGenePairs <- function(
  xymats.list,
  fdr.thresh,
  sign, # "positive" or "negative"
  model.name
) {
  df <- data.frame()
  if (grepl("^1", model.name)){
    coef.name <- "betaXits"
    pval.name <- "pvalXits"
  } else {
    stop("This function is applicable to model 1 only.")
  }
  df.list <- bplapply(
    X = xymats.list,
    FUN = getPeakGenePairsForSingleGene,
    coef.name = coef.name,
    pval.name = pval.name,
    sign = sign,
    fdr.thresh = fdr.thresh
  )
  df <- do.call("rbind", df.list)
  rownames(df) <- NULL
  return(df)
}
