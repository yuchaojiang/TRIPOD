getOverlapsAcrossGenes <- function(
	genes,
	TF.name,
	xymats.list,
	hits,
	subject
) {
	res.list <- bplapply(
		genes,
		getOverlapsForSingleGene,
		TF.name = TF.name,
  	xymats.list = xymats.list, # ...
		hits = hits,
  	subject = chip.gr
  )
  res <- do.call("rbind", res.list)
  return(res)
}
 