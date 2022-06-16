getSignificantPairs <- function(hits,
                                genes,
                                TSS,
                                resolution, # 1e4
                                gr.list) {
  res.list <- bplapply(
    genes,
    getSignificantPairsForSingleGene,
    TSS = TSS, # ...
    resolution = resolution,
    gr.list = gr.list,
    hits = hits
  )
  res <- unique(do.call("c", res.list))
  return(res)
}
