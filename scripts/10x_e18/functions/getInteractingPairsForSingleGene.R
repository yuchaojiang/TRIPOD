getInteractingPairsForSingleGene <- function(gene.name,
                                             xymats.list,
                                             TSS,
                                             resolution, # 1e4
                                             gr.list) {
    xymats <- xymats.list[[gene.name]]
    peak.gr.g <- xymats$peak.gr.g
    if (is.null(peak.gr.g)) return(NULL)
  TSS.g <- TSS[TSS$gene_name == gene.name]
  # remove peaks whose distance from TSS is below resolution
    remove <- getRegion(rep(1, 2) * resolution, TSS.g)
  peak.gr.g <- peak.gr.g[!overlapsAny(peak.gr.g, remove)]
  # get the region 1 overlapping TSS if it exists
  fo <- as.data.frame(findOverlaps(TSS.g, gr.list[[1]]))
  if (!is.null(fo) && nrow(fo) > 0) {
      index <- unique(fo$subjectHits)
    peak.interact.1 <- peak.gr.g[overlapsAny(peak.gr.g, gr.list[[2]][index])]
    if (length(peak.interact.1) > 0) {
      peak.interact.1 <- unique(Signac::GRangesToString(peak.interact.1))
    } else {
      peak.interact.1 <- NULL
    }
  } else {
    peak.interact.1 <- NULL
  }
  # get the region 2 overlapping TSS if it exists
  fo <- as.data.frame(findOverlaps(TSS.g, gr.list[[2]]))
  if (!is.null(fo) && nrow(fo) > 0) {
      index <- unique(fo$subjectHits)
    peak.interact.2 <- peak.gr.g[overlapsAny(peak.gr.g, gr.list[[1]][index])]
    if (length(peak.interact.1) > 0) {
      peak.interact.2 <- unique(Signac::GRangesToString(peak.interact.2))
    } else {
      peak.interact.2 <- NULL
    }
  } else {
    peak.interact.2 <- NULL
  }
  # combine the results
  peak.interact <- c(peak.interact.1, peak.interact.2)
  if (length(peak.interact) > 0) {
    peak.interact <- paste0(gene.name, "_", peak.interact)
  } else {
    peak.interact <- NULL
  }
  return(peak.interact)
}
