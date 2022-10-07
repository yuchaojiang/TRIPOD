#' Get Cook's distance
#'
#' This function computes Cook's distance on the linear regression model
#' and prepares data for Cook's bar plot.
#'
#' @param ... further arguments to be passed.
#' @inheritParams testInfluence
#' @inheritParams getXYMatrices
#'
#' @import olsrr
#' @importFrom stats lm
#'
#' @export
getCooksD <- function(Yg, Xt, Yj, metacell.celltype, ...) {
  res <- lm(Yg ~ 1 + Xt + Yj + Xt * Yj) # res <- lm(Yg ~ Xt * Yj)
  k <- ols_prep_cdplot_data(res)
  d <- ols_prep_outlier_obs(k)
  d$txt[!is.na(d$txt)] <- metacell.celltype[!is.na(d$txt)]
  d.list <- list(d = d, k = k)
  return(d.list)
}

#' Get DFFIT
#'
#' This function computes DFFIT on the linear regression model
#' and prepares data for DFFIT plot.
#'
#' @inheritParams getCooksD
#'
#' @import olsrr
#' @importFrom stats dfbetas dffits lm
#'
#' @export
getDFFIT <- function(Yg, Xt, Yj, metacell.celltype, ...) {
	res <- lm(Yg ~ 1 + Xt + Yj + Xt * Yj) # res <- lm(Yg ~ Xt * Yj)
	dbetas <- NULL # why?
  obs <- NULL # why?
  txt <- NULL # why?
  dffitsm <- unlist(dffits(res))
  k <- olsrr:::model_n_coeffs(res)
  n <- olsrr:::model_rows(res)
  dffits.t <- sqrt(k / n) * 2
  title <- names(model.frame(res))[1]
  dffits.data <- data.frame(obs = seq_len(n), dbetas = dffitsm)
  d <- ols_prep_dfbeta_data(dffits.data, dffits.t)
  d$txt[!is.na(d$txt)] <- metacell.celltype[!is.na(d$txt)]
  d.list <- list(d = d, dffits.t = dffits.t, title = title)
  return(d.list)
}

#' Get p-values representing the influence of deleting a group of metacells
#'
#' This function computes sampling-based p-values representing the influence
#' of deleting a group of metacells on the linear interaction model.
#'
#' @param metacell.neighbors an integer vector containing indices of metacell
#' neighbors
#' @inheritParams getCooksD
#' @inheritParams testInfluence
#'
#' @return A m-by-4 or m-by-5 matrix containing p-values where m represents the
#' number of meta cells.
#'
#' @export
deleteNeighbor <- function(Yg, Xt, Yj, metacell.celltype, 
                           metacell.neighbors, seed = 1234, ...) {
  set.seed(seed)
  delta.coeff.pval <- matrix(nrow = length(metacell.celltype), ncol = 5)
  colnames(delta.coeff.pval) <- c("Intercept", "Xt", "Yj", "Xt:Yj", "Yg")
  for (i in 1:length(metacell.celltype)) {
    metacell.rm <- metacell.neighbors[i, ]
    delta.coeff.pval[i, ] <- testInfluence(Yg = Yg, Xt = Xt, Yj = Yj,
                                           metacell.celltype = metacell.celltype,
                                           metacell.rm = metacell.rm, 
                                           plot.histogram = FALSE,
                                           nsamp = 1000, seed = seed)
  }
  return(delta.coeff.pval)
}

#' Get p-values representing the influence of deleting a celltype
#'
#' This function computes sampling-based p-values representing the influence
#' of deleting metacells of a given cell type on the linear interaction model.
#'
#' @inheritParams getCooksD
#' @inheritParams testInfluence
#'
#' @return A data frame containing p-values.
#'
#' @export
deleteCellType <- function(Yg, Xt, Yj, metacell.celltype, seed = 1234, ...) {
  set.seed(seed)
  delta.coeff.pval <- matrix(nrow = length(unique(metacell.celltype)), ncol = 5)
  colnames(delta.coeff.pval) <- c("Intercept", "Xt", "Yj", "Xt:Yj", "Yg")
  rownames(delta.coeff.pval) <- unique(metacell.celltype)
  for (metacell.rm.celltype in unique(metacell.celltype)) {
    metacell.rm <- which(metacell.celltype == metacell.rm.celltype)
    delta.coeff.pval[metacell.rm.celltype, ] <- testInfluence(
      Yg = Yg, Xt = Xt, Yj = Yj, metacell.rm = metacell.rm, 
      metacell.celltype = metacell.celltype, 
      plot.histogram = FALSE, nsamp = 10000, seed = seed)
  }
  return(delta.coeff.pval)
}

#' Get p-values representing the influence of deleting branches of metacells
#'
#' This function computes sampling-based p-values representing influences
#' of deleting branches of metacells on the linear interaction model.
#' It requires a dendrogram of metacells.
#'
#' @param dend a dendrogram class object.
#' @inheritParams getCooksD
#' @inheritParams testInfluence
#'
#' @return A data frame containing p-values.
#'
#' @import dendextend
#' @export
deleteBranch <- function(Yg, Xt, Yj, dend, metacell.celltype, seed = 1234, ...) {
  set.seed(seed)
  # get subtrees, each of which corresponds to a split
  subtrees <- partition_leaves(dend)
  # remove cell types at each split
  delta.coeff.pval <- matrix(nrow = length(subtrees), ncol = 5, data = NA)
  colnames(delta.coeff.pval) <- c("Intercept", "Xt", "Yj", "Xt:Yj", "Yg")
  rownames(delta.coeff.pval) <- 1:length(subtrees)
  for (i in 2:length(subtrees)) {
    metacell.rm <- which(metacell.celltype %in% subtrees[[i]])
    delta.coeff.pval[i, ] <- testInfluence(
      Yg = Yg, Xt = Xt, Yj = Yj, metacell.rm = metacell.rm,
      metacell.celltype = metacell.celltype,
      plot.histogram = FALSE, nsamp = 1000, seed = seed)
  }
  return(delta.coeff.pval)
}

#' Get a dendrogram
#'
#' @param metacell.matrix a matrix containing RNA expression values for genes
#' or chromatin accessibility values for ATAC peaks.
#' Rows and columns represent cell clusters and feature names, respectively.
#' @inheritParams getXYMatrices
#'
#' @return A dendrogram.
#' @import dendextend
#' @importFrom stats aggregate as.dendrogram dist hclust
#' @export
getDendrogram <- function(metacell.matrix, metacell.celltype){
	tmp <- aggregate(metacell.matrix ~ metacell.celltype, FUN = sum)
  rownames(tmp) <- tmp[, 1]
  tmp <- tmp[, -1]
  # tmp <- tmp[ordered.celltype, ] # tmp <- tmp[col.map$celltype, ]
  tmp <- sweep(tmp, 1, rowSums(tmp)/10^6, "/")
  dend <- tmp %>% scale %>% dist %>% hclust %>% as.dendrogram
  rm(tmp)
  return(dend)
}

#' Get metacell neighbors
#'
#' @param reduction a character string specifying the reduced dimension data
#' @param k an integer value setting the number of metacell neighbors to include
#' @inheritParams getMetacellMatrix
#'
#' @return A matrix.
#' @import FNN Seurat
#' @export
getMetacellNeighbor <- function(object, reduction,
	cluster.name = "seurat_clusters", k = 4){
  reduction <- reduction
  cell.clusters <- object@meta.data[[cluster.name]]
  cell.reduction <- object@reductions[[reduction]]@cell.embeddings
  metacell.reduction <- matrix(
	  nrow = length(levels(cell.clusters)),
	  ncol = ncol(cell.reduction))
  rownames(metacell.reduction) <- paste0("metacell_", levels(cell.clusters))
  colnames(metacell.reduction) <- colnames(cell.reduction)
  for (i in 1:nrow(metacell.reduction)) {
    metacell.reduction[i, ] <- apply(cell.reduction[cell.clusters == (i - 1), ], 2, mean)
  }
  metacell.neighbors <- get.knn(metacell.reduction, k = k)$nn.index
  metacell.neighbors <- cbind(1:nrow(metacell.neighbors), metacell.neighbors)
  return(metacell.neighbors)
}

#' Get sampling-based p-values
#'
#' This function computes sampling-based p-values representing influences
#' of deleting metacells in a linear regression model.
#'
#' @param Yg a numeric vector containing RNA expression values.
#' @param Xt a numeric vector containing chromatin accessibility.
#' @param Yj a numeric vector containing TF expression.
#' @param alternative a character string specifying the alternative hypothesis.
#' This must be one of "two.sided", "greater", and "less".
#' @param metacell.rm an integer vector containing indices of metacells
#' to be removed.
#' @param metacell.celltype a character vector specifying cell types of the metacells.
#' @param plot.histogram an indicator as to whether to generate histograms
#' @param nsamp an integer indicating the number of sampling
#' @param seed an integer
#' @param ... further arguments passed to {\code{\link{hist}}}.
#'
#' @return A numerical vector containing p-values.
#' @importFrom stats lm model.frame
#' @importFrom graphics abline hist par
#'
#' @export
testInfluence <- function(Yg, Xt, Yj, metacell.rm, metacell.celltype,
                          plot.histogram = FALSE, nsamp = 10000, seed = 1234, ...) {
  set.seed(seed)
  # fit a model with all data points
  model <- lm(Yg ~ Xt * Yj)
  # fit a model without the selected metacells (observed)
  model.delta.obs <- lm(Yg[-metacell.rm] ~ Xt[-metacell.rm] * Yj[-metacell.rm])
  Yg.pred <- (model$fitted.values)
  Yg.metacell.rm.pred <- t(as.matrix(model.delta.obs$coefficients)) %*%
    rbind(rep(1, length(Yg)), Xt, Yj, (Xt * Yj))
  delta.coeff <- matrix(ncol = 5, nrow = nsamp)
  delta.coeff[1, 1:4] <- model.delta.obs$coefficients - model$coefficients
  delta.coeff[1, 5] <- mean(abs(Yg.pred - Yg.metacell.rm.pred))
  colnames(delta.coeff) <- c(names(model$coefficients), "Yg")
  for (ii in 2:nsamp) {
    metacell.rm.samp <- sample(1:length(metacell.celltype), length(metacell.rm))
    model.delta.samp <- lm(Yg[-metacell.rm.samp] ~ 1 + Xt[-metacell.rm.samp] +
                             Yj[-metacell.rm.samp] + (Xt * Yj)[-metacell.rm.samp])
    Yg.metacell.rm.samp.pred <- t(as.matrix(model.delta.samp$coefficients)) %*%
      rbind(rep(1, length(Yg)), Xt, Yj, (Xt * Yj))
    delta.coeff[ii, 1:4] <- model.delta.samp$coefficients - model$coefficients
    delta.coeff[ii, 5] <- mean(abs(Yg.pred - Yg.metacell.rm.samp.pred))
  }
  delta.coeff.pval <- apply(delta.coeff, 2, function(x) sum(abs(x) >= abs(x[1])) / length(x))
  if (plot.histogram) {
    par(mfrow = c(2, 3))
    for (i in 1:4) {
      hist(delta.coeff[, i], 100,
           xlab = paste("Delta coeff. for", colnames(delta.coeff)[i]),
           main = paste(metacell.celltype[metacell.rm][1],"\nDelta coeff. for",
                        colnames(delta.coeff)[i], "\npval =", delta.coeff.pval[i]), ...
      )
      abline(v = delta.coeff[1, i], col = 2, lty = 2)
    }
    i=5
    hist(delta.coeff[, i], 100,
         xlab = paste("Delta for", colnames(delta.coeff)[i]),
         main = paste(metacell.celltype[metacell.rm][1],"\nDelta for",
                      colnames(delta.coeff)[i], "\npval =", delta.coeff.pval[i]), ...
    )
    abline(v = delta.coeff[1, i], col = 2, lty = 2)
    par(mfrow = c(1, 1))
  }
  return(delta.coeff.pval)
}
