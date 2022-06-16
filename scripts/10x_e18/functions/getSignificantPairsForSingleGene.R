getSignificantPairsForSingleGene <- function(gene.name,
                                             TSS,
                                             resolution, # 1e4
                                             gr.list,
                                             hits # data frame
) {
  TSS.g <- TSS[TSS$gene_name == gene.name]
  cnd.1 <- any(overlapsAny(TSS.g, gr.list[[1]]))
  cnd.2 <- any(overlapsAny(TSS.g, gr.list[[2]]))
  if (cnd.1 | cnd.2) {
    peak.sig.strings <- hits$peak[hits$gene == gene.name]
    if (length(peak.sig.strings) == 0) {
      return(NULL)
    }
    peak.sig <- StringToGRanges(peak.sig.strings)
    # remove peaks whose distance from TSS is below resolution
    remove <- getRegion(rep(1, 2) * resolution, TSS.g)
    peak.sig <- peak.sig[!overlapsAny(peak.sig, remove)]
    if (length(peak.sig) == 0) {
      return(NULL)
    }
    peak.sig <- unique(Signac::GRangesToString(peak.sig))
    peak.sig <- paste0(gene.name, "_", peak.sig)
    return(peak.sig)
  } else {
    return(NULL)
  }
}
