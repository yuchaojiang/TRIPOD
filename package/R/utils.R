#' Perform the Benjamini and Hochberg FDR procedure on a matrix
#'
#' @param pvals a matrix containing p-values.
#' @param fdr.thresh a numerical scalar specifying the FDR threshold.
#' @param do.plot a logical scalar indicating whether to
#' produce a diagnostic plot.
#'
#' @importFrom graphics grid
#' @return a logical matrix indicating whether the null hypothesis
#' has been rejected for each element.
performBHMatrix <- function(pvals, fdr.thresh = 0.01, do.plot = NULL) {
  if (is.null(do.plot)) do.plot <- FALSE
  arr.ind <- which(!is.na(pvals))
  pvalvec <- pvals[arr.ind]
  pvalvec.sorted <- sort(pvalvec)
  pvalvec.ord <- order(pvalvec)
  kstar <- max(which(pvalvec.sorted <= fdr.thresh * c(1:length(pvalvec)) / length(pvalvec)))
  if (do.plot) {
    plot(pvalvec.sorted, ylab = "Sorted pvalues",
    	main = paste(kstar, "passed", fdr.thresh, "B-H threshold"))
    abline(0, fdr.thresh / length(pvalvec), col = "red")
    grid()
  }
  which.reject <- matrix(nrow = nrow(pvals), ncol = ncol(pvals), data = FALSE)
  if (kstar >= 1 & kstar <= length(pvals)) which.reject[arr.ind[pvalvec.ord[1:kstar]]] <- TRUE
  return(which.reject)
}

#' Perform the Benjamini and Hochberg FDR procedure on a vector
#'
#' @param pval.c a numerical vector containing p-values.
#' @param fdr.thresh a numerical scalar specifying the FDR threshold.
#'
#' @return a list of two elements:
#' \describe{
#'   \item{pval.c.adj}{a numerical vector containing adjusted p-values}
#'   \item{which.reject}{a logical vector indicating whether the null hypothesis
#' has been rejected for each element.}
#' }
performBHVector <- function(pval.c, fdr.thresh = 0.01) {
  pval.c.rmNA <- pval.c[!is.na(pval.c)]
  pvalvec.sorted <- sort(pval.c.rmNA)
  pvalvec.ord <- order(pval.c.rmNA)
  pval.c.adj.rmNA <- pval.c.rmNA
  pval.c.adj.rmNA[pvalvec.ord] <- pvalvec.sorted * length(pval.c.rmNA) / c(1:length(pval.c.rmNA))
  pval.c.adj <- pval.c
  pval.c.adj[!is.na(pval.c)] <- pval.c.adj.rmNA
  which.reject <- pval.c.adj < fdr.thresh
  pval.c.adj <- pmin(1, pval.c.adj)
  return(list(pval.c.adj = pval.c.adj, which.reject = which.reject))
}

#' Cap values
#'
#' This function takes a numerical vector as an input and returns
#' values that are capped at a specific quantile.
#' This is useful for removing outliers.
#'
#' @param x a numerical vector
#' @param cap.at.quantile a quantile at which variables are capped
#' to remove outliers.
#'
#' @return a numerical vector containing capped values
#' @export
capValues <- function(x, cap.at.quantile = 0.02) {
  upper.thresh <- quantile(x, 1 - cap.at.quantile, na.rm = TRUE)
  lower.thresh <- quantile(x, cap.at.quantile, na.rm = TRUE)
  x <- pmin(pmax(x, lower.thresh), upper.thresh)
  return(x)
}

#' Compute distances
#'
#' @param x a numerical matrix or vector.
#'
#' @return a numerical matrix
computeDistances <- function(x) {
  if (!is.matrix(x)) {
    x <- matrix(data = x, nrow = length(x), ncol = 1)
  }

  dmat <- matrix(nrow = nrow(x), ncol = nrow(x), data = NA)
  for (i in 1:nrow(x)) {
    xi <- matrix(nrow = nrow(x), ncol = ncol(x), data = x[i, ], byrow = TRUE)
    dmat[i, ] <- sqrt(rowSums((xi - x)^2) / ncol(x))
  }
  return(dmat)
}

# performHypergeometricTest <- function(peaks.gr.g, sig.peaks.gr.g, enhancer.db) {
#   # peaks.gr.g: peak granges for gene g
#   # sig.peaks.gr.g: significant peak granges for gene g
#   # enhancer.db: enhancer database granges
#
#   # Hypergeomeetric test
#   # m+n: Total number of peaks (overlap or dont overlap with enhancers)
#   # m: Number peaks that overlap with enhancers
#   # k: Number of significant peaks
#   # q: Number of significant peaks that overlap with enhancers
#   if (length(enhancer.db) == 0) {
#     pval <- NA
#   } else {
#     m <- sum(countOverlaps(peaks.gr.g, enhancer.db) > 0)
#     n <- length(peaks.gr.g) - m
#     k <- length(sig.peaks.gr.g)
#     q <- sum(countOverlaps(sig.peaks.gr.g, enhancer.db) > 0)
#     # cat(m,'\t',m+n,'\t',q,'\t',k,'\n')
#     if (m == 0) { # no peak overlaps with enhancer databases
#       pval <- NA
#     } else {
#       pval <- phyper(q, m, n, k, lower.tail = F)
#     }
#   }
#   return(pval)
# }

#' Transform coefficients to make them easier to plot
#'
#' @param coefs a numerical matrix or vector.
#'
#' @return a numerical matrix or vector.
transformCoef <- function(coefs) {
  sign(coefs) * log(1 + abs(coefs))
}

#' Transform p-values to make them easier to plot
#'
#' @param pvals a numerical matrix or vector.
#'
#' @return a numerical matrix or vector.
transformPvalue <- function(pvals) {
  sqrt(-log((pvals)))
}
