getRegion <- function(vicinity.size, TSS){
  upstream.size <- vicinity.size[1]
  downstream.size <- vicinity.size[2]
  if (upstream.size == 0) {
    upstream <- TSS
  } else if (upstream.size > 0) {
    upstream <- flank(TSS, width=upstream.size, start=T) # upstream excluding the TSS
  }
  if (downstream.size == 0) {
    downstream <- TSS
  } else if (downstream.size > 0) {
    downstream <- flank(TSS, width=downstream.size, start=F) # downstream excluding the TSS
  }
  vicinity <- punion(upstream, downstream, fill.gap=T)
  start(vicinity) <- pmax(0, start(vicinity))
  # end(vicinity) <- pmin(end(vicinity), choromosome.size))
  vicinity$gene_name <- TSS$gene_name
  return(vicinity)
}
