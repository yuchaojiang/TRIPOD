getInteractingPairs <- function(genes,
                                xymats.list,
                                TSS,
                                resolution, # 1e4
                                gr.list) {
  res.list <- bplapply(
    genes,
    getInteractingPairsForSingleGene,
    xymats.list = xymats.list, # ...
    TSS = TSS,
    resolution = resolution,
    gr.list = gr.list
  )
  res <- do.call("c", res.list)
  return(res)
}
