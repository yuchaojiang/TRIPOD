getEnrichmentScore <- function(
	data.frame.list # list of data frames
) {
	x <- data.frame(t(sapply(data.frame.list, colSums)))
  pval <- apply(x, 1, getHypergeometricPvalue)
  -log10(pval)
}
