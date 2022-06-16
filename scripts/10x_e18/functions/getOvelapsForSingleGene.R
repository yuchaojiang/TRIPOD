getOverlapsForSingleGene <- function(
	gene.name,
	TF.name,
	xymats.list,
	hits, # a data frame
	subject
) {
	xymats <- xymats.list[[gene.name]]
  peak.gr <- xymats$peak.gr.g # can select peaks containing the motif corresponding to the TF
  if (is.null(peak.gr)) {
  	print(paste0("peak.gr is NULL for ", gene.name))
  	return(NULL)
  }
  sig.peaks <- hits$peak[hits$gene == gene.name & hits$TF == TF.name]
  # if (length(sig.peaks) == 0) {return(NULL)}
  m <- sum(overlapsAny(peak.gr, subject))
  n <- length(peak.gr) - m
  if (length(sig.peaks) > 0) {
  	sig.peak.gr <- StringToGRanges(sig.peaks)
  	k <- length(sig.peak.gr)
  	q <- sum(overlapsAny(sig.peak.gr, subject))
  } else if (length(sig.peaks) == 0){
  	k <- q <- 0
  }
  df <- data.frame(q = q, m = m, n = n, k = k)
  return(df)
}

