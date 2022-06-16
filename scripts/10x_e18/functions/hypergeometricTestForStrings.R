hypergeometricTestForStrings <- function(x, all, ground.truth) {
  q <- length(intersect(x, ground.truth)) # TP
  m <- length(ground.truth) # TP + FN
  n <- length(all) - m # TN + FP
  k <- length(x) # TP + FP
  if (m > 0) {
    pval <- phyper(q, m, n, k, lower.tail = FALSE)
  }
  return(c(q, m, n, k, pval))
}
