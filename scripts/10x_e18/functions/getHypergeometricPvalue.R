getHypergeometricPvalue <- function(x){
  ifelse(x[1] > 0, phyper(x[1], x[2], x[3], x[4], lower.tail = FALSE), NA)
}
