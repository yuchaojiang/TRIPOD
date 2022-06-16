convertTrioStringToDF <- function(x) {
  d <- data.frame(t(sapply(x, function(y) unlist(strsplit(y, "_")))))
  colnames(d) <- c("gene", "peak", "TF")
  rownames(d) <- NULL
  return(d)
}
