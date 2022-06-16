getAllPairsForSingleGene <- function(gene.name,
                                     xymats.list,
                                     TSS,
                                     resolution, # 1e4
                                     gr.list # interacting chromosome regions
) {
  xymats <- xymats.list[[gene.name]]
  peak.gr.g <- xymats$peak.gr.g
  TSS.g <- TSS[TSS$gene_name == gene.name]
  cnd.1 <- any(overlapsAny(TSS.g, gr.list[[1]]))
  cnd.2 <- any(overlapsAny(TSS.g, gr.list[[2]]))
  if (cnd.1 | cnd.2) {
    # remove peaks whose distance from TSS is below resolution
    # or remove peaks that overlap interacting chromosome regions containing TSS
    remove <- getRegion(rep(1, 2) * resolution, TSS.g)
    peak.all <- peak.gr.g[!overlapsAny(peak.gr.g, remove)]
    peak.all <- unique(Signac::GRangesToString(peak.all))
    peak.all <- paste0(gene.name, "_", peak.all)
    # return lists of
  } else {
    peak.all <- NULL
  }
  return(peak.all)
}
