#' Create input matrices for model fitting for a given target gene
#'
#' getXYMatrices() creates input matrices, Yg, X, Yj, for model fitting
#' for a given target gene. Unlike {\code{\link{getXYMatrices}}},
#' the cases where Yg and Yj are the same are excluded.
#'
#' @param gene.name a character string containing the gene name. This must exist
#' in the "gene_name" column of the transcripts.gr object.
#' @param ext.upstream an integer representing the window size in bp.
#' @param ext.downstream an integer representing the window size in bp.
#' @param transcripts.gr a GRanges object containing annotated protein-coding transcripts.
#' @param peaks.gr a GRanges object containing ATAC peaks.
#' @param metacell.rna a matrix containing RNA expression values for genes.
#' Rows and columns represent cell clusters and gene names, respectively.
#' @param metacell.peak a matrix containing chromatin accessibility values for ATAC peaks.
#' Rows and columns represent cell clusters and ATAC peak positions, respectively.
#' @param peakxmotif a binary matrix containing indicators as to whether ATAC peaks contain motifs.
#' @param motifxTF a matrix containing names of motifs and TFs.
#' @param metacell.celltype a character vector specifying cell types of the metacells.
#' @param metacell.celltype.col a character vector specifying colors of the metacells
#' based on the cell types.
#'
#' @return A list with elements:
#' \item{gene.name}{a character string.}
#' \item{Yg}{ a vector containig RNA expression of the gene across metacells.}
#' \item{Xt}{a matrix containig chromatin accessibility of the ATAC peaks.}
#' \item{Yj}{a matrix containing TF RNA expression.
#' note that this is gene-specific since TFs that do no have motifs in the region
#' around TSS are excluded.}
#' \item{peakxmotif.g}{a gene-specific binary matrix containing indicators of
#' whether ATAC peaks contains motifs.
#' rows and columns represent peaks and motifs, respectively.}
#' \item{nonzero.peakxmotif.g}{a matrix with two columns.
#' The first and second columns contain row and column numbers
#' in the peakmotif.g matrix where the entry is non-zero.}
#'
#' @import GenomicRanges Matrix
#' @importFrom methods as
#'
#' @export
getXYMatrices <- function(gene.name, ext.upstream, ext.downstream = NULL,
                          transcripts.gr, peaks.gr,
                          metacell.rna, metacell.peak,
                          peakxmotif, motifxTF,
                          metacell.celltype = NULL, metacell.celltype.col = NULL) {
  if (is.null(ext.downstream)) ext.downstream <- ext.upstream
  transcripts.ext.gr <- promoters(transcripts.gr,
  	upstream = ext.upstream, downstream = ext.downstream + 1)

  # get a vector of the target gene expression
  Yg <- metacell.rna[, gene.name]
  # get ATAC peaks within the vicinity region of the target gene
  transcripts.gr.g <- transcripts.gr[transcripts.gr$gene_name == gene.name]
  transcripts.ext.gr.g <- transcripts.ext.gr[transcripts.ext.gr$gene_name == gene.name]
  peaks.gr.g <- subsetByOverlaps(peaks.gr, transcripts.ext.gr.g)
  Xt <- metacell.peak[, overlapsAny(peaks.gr,
  	transcripts.ext.gr[transcripts.ext.gr$gene_name == gene.name]), drop = FALSE]
  # get motifs within vicinity region of gene g
  peakxmotif.g <- peakxmotif[overlapsAny(peaks.gr,
  	transcripts.ext.gr[transcripts.ext.gr$gene_name == gene.name]), , drop = FALSE]
  # remove motifs that do not have binding sites in this gene region
  peakxmotif.g <- peakxmotif.g[, apply(peakxmotif.g, 2, sum) > 0, drop = FALSE]
  # get a vector containing TF names
  TF.g <- motifxTF[match(colnames(peakxmotif.g), motifxTF[, 1]), 2]
  # get gene expression of TF
  Yj <- metacell.rna[, TF.g, drop = FALSE]
  # get a logical vector for removing TFs if it is not expressed in any cells
  TF.filter1 <- apply(Yj, 2, sum) > 0
  # get a logical vector for removing a TF if it is the same as the target gene
  TF.filter2 <- TF.g != gene.name
  # get a logical vector
  TF.filter <- TF.filter1 & TF.filter2
  # remove the TFs
  Yj <- Yj[, TF.filter, drop = FALSE]
  TF.g <- TF.g[TF.filter, drop = FALSE]
  peakxmotif.g <- peakxmotif.g[, TF.filter, drop = FALSE]
  # get nonzero.peakxmotif.g
  peakxmotif.g <- as.matrix(peakxmotif.g)
  nonzero.peakxmotif.g <- which(peakxmotif.g != 0, arr.ind = T)
  rm(TF.filter)
  peakxmotif.g <- as(peakxmotif.g, "dgCMatrix")

  results <- list(
  	gene.name = gene.name,
  	ext.upstream = ext.upstream,
  	ext.downstream = ext.downstream,
    Yg = Yg,
  	Xt = Xt,
  	Yj = Yj,
    peakxmotif.g = peakxmotif.g,
  	nonzero.peakxmotif.g = nonzero.peakxmotif.g,
    TF.g = TF.g,
  	peaks.gr.g = peaks.gr.g,
    metacell.celltype = metacell.celltype,
  	metacell.celltype.col = metacell.celltype.col
  )
  return(results)
}

#' Perform model fitting for a given target gene
#'
#' fitModel() takes a list returned by {\code{\link{getXYMatrices}}}
#' as input, performs model fitting for a given target gene,
#' and returns a list object containing matrices of coefficients and p-values.
#' The (t, j) elements of the matrices
#' represent values for the t-th ATAC peak and the j-th TF.
#'
#' @param xymats a list returned by {\code{\link{getXYMatrices}}}.
#' @param model.name 	a character string specifying the model. This must be
#' "marginal", "conditional", "interaction", "TRIPOD".
#' @param match.by a character string specifying the variable to be matched.
#' This is necessary for "TRIPOD" only.
#' @param log a logical indicating whether to log-transform the variables.
#' @param cap.at.quantile a quantile at which variables are capped to remove outliers.
#' That is, values are capped if they are below cap.at.quantile
#' or above (1 - cap.at.quantile).
#' Set this to 0 if capping is unnecessary.
#' @param delta a numeric.
#' @param min.Xt an integer to set the minimum number of nonzero Xt.
#' @param min.Yj an integer to set the minimum number of nonzero Yj.
#' @param min.level1 an integer to set the minimum number of pairs.
#' @param min.level2 an integer to set the minimum number of pairs.
#' @param metacell.celltype.col a character vector specifying colors of the metacells
#' based on the cell types. This is necessary only for generating color-coded
#' scatter plots of pairs of metacells using {\code{\link{plotGenePeakTFScatter}}}.
#' @param verbose a logical indicating whether to print progress in "TRIPOD".
#'
#' @return A list containing some of the following elements. Note that the sets of
#' elements differ depending on the model.name input argument.
#' \item{gene.name}{a character string.}
#' \item{Yg}{ a vector containig capped RNA expression of the gene across metacells.}
#' \item{Xt}{a matrix containig capped chromatin accessibility of the ATAC peaks.}
#' \item{Yj}{a matrix containing capped TF RNA expression.
#' Note that TFs that do no have motifs in the region
#' around TSS are excluded.}
#' \item{peakxmotif.g}{a gene-specific binary matrix containing indicators of
#' whether ATAC peaks contains motifs.
#' rows and columns represent peaks and motifs, respectively.}
#' \item{nonzero.peakxmotif.g}{a matrix with two columns.
#' The first and second columns contain row and column numbers
#' in the peakmotif.g matrix where the entry is non-zero.}
#' \item{alpha.m}{a vector containing estimates of coefficient from the marginal model, Yg ~ Xt.}
#' \item{pval.alpha.m}{p-values associated with alpha.m.}
#' \item{beta.m}{a vector containing estimates of coefficient from the marginal model, Yg ~ Yj.}
#' \item{pval.beta.m}{p-values associated with beta.m.}
#' \item{gamma.m}{a matrix containing estimates of coefficient from the model, Yg ~ Xt:Yj.}
#' \item{pval.gamma.m}{p-values associated with gamma.m.}
#' \item{eta.m}{a matrix containing estimates of coefficient from the martinal model, Yj ~ Xt.}
#' \item{pval.eta.m}{p-values associated with eta.m.}
#' \item{cor.alpha}{Pearson correlation coefficient between Yg and Xt.}
#' \item{cor.beta}{Pearson correlation coefficient between Yg and Yj.}
#' \item{cor.gamma}{Pearson correlation coefficient between Yg and
#' the product of Xt and Yj.}
#' \item{cor.eta}{Pearson correlation coefficient between Yj and Xt.}
#' \item{alpha.c}{a matrix containing estimates of coefficient from the conditional model, Yg ~ Xt + Yj.}
#' \item{pval.alpha.c}{p-values associated with alpha.c.}
#' \item{beta.c}{a matrix containing estimates of coefficient from the conditional model, Yg ~ Xt + Yj.}
#' \item{pval.beta.c}{p-values associated with beta.c.}
#' \item{delta}{This applies to "conditional" and "TRIPOD" only.}
#' \item{alpha.i}{a matrix containing estimates of coefficient from the interation model.}
#' \item{pval.alpha.i}{p-values associated with alpha.i.}
#' \item{beta.i}{a matrix containing estimates of coefficient from the interation model.}
#' \item{pval.beta.i}{p-values associated with beta.i.}
#' \item{gamma.i}{a matrix containing estimates of coefficient from the interation model.}
#' \item{pval.gamma.i}{p-values associated with gamma.i.}
#' \item{coef.level1}{a matrix containing estimates of coefficient from TRIPOD level 1 test.}
#' \item{pval.level1}{p-values associated with coef.level1.}
#' \item{coef.level2}{a matrix containing estimates of coefficient from TRIPOD level 2 test.}
#' \item{pval.level2}{p-values associated with coef.level2.}
#' \item{match.by}{the matched variable. This applieds to "TRIPOD" only.}
#'
#' @import GenomicRanges Matrix nbpMatching
#' @importFrom stats IQR cor cor.test lm model.frame
#'
#' @export
fitModel <- function(xymats, model.name, match.by = NULL, log = NULL,
	                   cap.at.quantile = 0.02, delta = 0.2,
	                   min.Xt = 5, min.Yj = 5, min.level1 = 10, min.level2 = 10,
	                   metacell.celltype.col = NULL,
	                   # do.pred.plot = NULL, do.plot = NULL, plot.dir = NULL,
                     verbose = FALSE) {
	if (all(model.name != c("marginal", "conditional", "interaction", "TRIPOD"))) {
  	stop('The model.name argument must be one of "marginal", "conditional", "interaction", and "TRIPOD".')
  }
  if (!is.null(log)) {
    if (all(log %in% c("Yg", "Xt", "Yj"))) {
      if (any(log == "Yg")) {
        xymats$Yg <- log(xymats$Yg + 1)
      }
      if (any(log == "Xt")) {
        xymats$Xt <- log(xymats$Xt + 1)
      }
      if (any(log == "Yj")) {
        xymats$Yj <- log(xymats$Yj + 1)
      }
    } else {
      stop('The log argument must be one of "Yg", "Xt", and "Yj".')
    }
  }
	# cap values
	if (cap.at.quantile > 0) {
    xymats$Yg <- capValues(xymats$Yg, cap.at.quantile = cap.at.quantile)
    xymats$Xt <- apply(xymats$Xt, 2, capValues, cap.at.quantile = cap.at.quantile)
    xymats$Yj <- apply(xymats$Yj, 2, capValues, cap.at.quantile = cap.at.quantile)
  }

  if (model.name == "marginal") {
    nTFs <- ncol(xymats$Yj)
    npeaks <- ncol(xymats$Xt)
    alpha.m <- pval.alpha.m <- cor.alpha <- rep(NA, npeaks)
    beta.m <- pval.beta.m <- cor.beta <- rep(NA, nTFs)
    for (t in 1:npeaks) {
      Xt <- xymats$Xt[, t]
      Yg <- xymats$Yg
      # check the number of zeros in Xt to ensure that there is sufficient degrees of freedom
      if (sum(Xt > 0) >= min.Xt) {
        res <- lm(Yg ~ 1 + Xt)
        if (any(is.na(res$coefficients))) next
        alpha.m[t] <- res$coefficients[2]
        pval.alpha.m[t] <- summary(res)$coefficients[2, 4]
        cor.alpha[t] <- cor(Xt, Yg)
      }
    }
    for (j in 1:nTFs) {
      Yg <- xymats$Yg
      Yj <- xymats$Yj[, j]
      # check the number of zeros in Yj to ensure that there is sufficient degrees of freedom
      if (sum(Yj > 0) >= min.Yj) {
        res <- lm(Yg ~ 1 + Yj)
        if (any(is.na(res$coefficients))) next
        beta.m[j] <- res$coefficients[2]
        pval.beta.m[j] <- summary(res)$coefficients[2, 4]
        cor.beta[j] <- cor(Yj, Yg)
      }
    }

    gamma.m <- pval.gamma.m <- eta.m <- pval.eta.m <- cor.gamma <- cor.eta <-
      matrix(nrow = nrow(xymats$peakxmotif.g), ncol = ncol(xymats$peakxmotif.g), data = NA)

    for (j in 1:nTFs) {
      Xt <- xymats$Xt
      Yg <- xymats$Yg
      Yj <- xymats$Yj[, j]
      # check the number of zeros in Yj to ensure that there is sufficient degrees of freedom
      if (sum(Yj > 0) >= min.Yj) {
        for (t in which(xymats$peakxmotif.g[, j] == 1)) {
        	# check if Xt * Yj is zero across all metacells even though they themselves are non-zeros
          if (all(Xt[, t] * Yj == 0)) next
          if (sum(Xt[, t] > 0) < min.Xt) next
        	# drop the main effects
          res <- lm(Yg ~ 1 + Xt[, t]:Yj)
          if (any(is.na(res$coefficients))) next
          gamma.m[t, j] <- res$coefficients[2]
          pval.gamma.m[t, j] <- summary(res)$coefficients[2, 4]

          lmfit <- summary(lm(Yj ~ Xt[, t]))
          eta.m[t, j] <- lmfit$coefficients[2, 1]
          pval.eta.m[t, j] <- lmfit$coefficients[2, 4]

          cor.gamma[t, j] <- cor(Xt[, t] * Yj, Yg)
          cor.eta[t, j] <- cor(Xt[, t], Yj)
        }
      }
    }
    xymats$alpha.m <- alpha.m
    xymats$pval.alpha.m <- pval.alpha.m
    xymats$beta.m <- beta.m
    xymats$pval.beta.m <- pval.beta.m
    xymats$gamma.m <- gamma.m
    xymats$pval.gamma.m <- pval.gamma.m
    xymats$eta.m <- eta.m
    xymats$pval.eta.m <- pval.eta.m
    xymats$cor.alpha <- cor.alpha
    xymats$cor.beta <- cor.beta
    xymats$cor.gamma <- cor.gamma
    xymats$cor.eta <- cor.eta
  }

  if (model.name == "conditional") {
    nTFs <- ncol(xymats$Yj)
    alpha.c <- beta.c <- pval.alpha.c <- pval.beta.c <-
    	matrix(nrow = nrow(xymats$peakxmotif.g), ncol = ncol(xymats$peakxmotif.g), data = NA)

    for (j in 1:nTFs) {
      Xt <- xymats$Xt
      Yg <- xymats$Yg
      Yj <- xymats$Yj[, j]
      # check the number of zeros in Yj to ensure that there is sufficient degrees of freedom
      if (sum(Yj > 0) >= min.Yj) {
        for (t in which(xymats$peakxmotif.g[, j] == 1)) {
        	# check if Xt * Yj is zero across all metacells even though they themselves are non-zeros
          if (all(Xt[, t] * Yj == 0)) next
          if (sum(Xt[, t] > 0) < min.Yj) next
          delta.val <- max(quantile(Xt[, t], delta, na.rm = TRUE), 0.1)
          # carry out Yj testing only if Xt is greater than a threshold
          sel.Xt <- Xt[, t] >= delta.val
          res <- lm(Yg[sel.Xt] ~ 1 + Xt[sel.Xt, t] + Yj[sel.Xt])
          if (any(is.na(res$coefficients))) next
          beta.c[t, j] <- res$coefficients[3]
          pval.beta.c[t, j] <- summary(res)$coefficients[3, 4]

          delta.val <- max(quantile(Yj, delta, na.rm = TRUE), 0.1)
          # carry out Xt testing only if Yj is greater than a threshold
          sel.Yj <- Yj >= delta.val
          res <- lm(Yg[sel.Yj] ~ 1 + Xt[sel.Yj, t] + Yj[sel.Yj])
          if (any(is.na(res$coefficients))) next
          alpha.c[t, j] <- res$coefficients[2]
          pval.alpha.c[t, j] <- summary(res)$coefficients[2, 4]
        }
      }
    }
    xymats$alpha.c <- alpha.c
    xymats$beta.c <- beta.c
    xymats$pval.alpha.c <- pval.alpha.c
    xymats$pval.beta.c <- pval.beta.c
    xymats$delta <- delta
  }

  if (model.name == "interaction") {
    nTFs <- ncol(xymats$Yj)
    alpha.i <- beta.i <- gamma.i <- pval.alpha.i <- pval.gamma.i <- pval.beta.i <-
    	matrix(nrow = nrow(xymats$peakxmotif.g), ncol = ncol(xymats$peakxmotif.g), data = NA)

    for (j in 1:nTFs) {
      Xt <- xymats$Xt
      Yg <- xymats$Yg
      Yj <- xymats$Yj[, j]
      # check the number of zeros in Yj to ensure that there is sufficient degrees of freedom
      if (sum(Yj > 0) >= min.Yj) {
        for (t in which(xymats$peakxmotif.g[, j] == 1)) {
        	# check if Xt * Yj is zero across all metacells even though they themselves are non-zeros
          if (all(Xt[, t] * Yj == 0)) next
          if (sum(Xt[, t] > 0) < min.Yj) next
          res <- lm(Yg ~ 1 + Xt[, t] + Yj + Xt[, t] * Yj)
          if (any(is.na(res$coefficients))) next
          alpha.i[t, j] <- res$coefficients[2]
          beta.i[t, j] <- res$coefficients[3]
          gamma.i[t, j] <- res$coefficients[4]
          pval.alpha.i[t, j] <- summary(res)$coefficients[2, 4]
          pval.beta.i[t, j] <- summary(res)$coefficients[3, 4]
          pval.gamma.i[t, j] <- summary(res)$coefficients[4, 4]
        }
      }
    }
    xymats$alpha.i <- alpha.i
    xymats$beta.i <- beta.i
    xymats$gamma.i <- gamma.i
    xymats$pval.alpha.i <- pval.alpha.i
    xymats$pval.beta.i <- pval.beta.i
    xymats$pval.gamma.i <- pval.gamma.i
  }

  if (model.name == "TRIPOD") {
    pval.level1 <- coef.level1 <- coef.level2 <- pval.level2 <-
    	matrix(nrow = nrow(xymats$peakxmotif.g), ncol = ncol(xymats$peakxmotif.g), data = NA)
    if (verbose) cat("Performing matched tests ... \n")
    for (t in 1:ncol(xymats$Xt)) {
      if (verbose) cat("Peak #", t, "out of", ncol(xymats$Xt), ": ",
      	sum(xymats$peakxmotif.g[t, ]), "TF with motifs. Motif:")
      for (j in which(xymats$peakxmotif.g[t, ] == 1)) {
        if (verbose) cat(j, " ")
        Xt <- xymats$Xt[, t]
        Yg <- xymats$Yg
        Yj <- xymats$Yj[, j]
        # file.name <- paste(c("matchedTest_peak", t, "_TF", j, ".png"), collapse = "")
        res <- performMatchedTest(Xt = Xt, Yj = Yj, Yg = Yg,
          match.by = match.by,
          cap.at.quantile = cap.at.quantile, delta = delta,
        	min.level1 = min.level1, min.level2 = min.level2,
          # do.plot = do.plot, col = metacell.celltype.col,
          # plot.dir = plot.dir, file.name = file.name, partition.screen = TRUE,
          verbose = verbose
        )
        coef.level1[t, j] <- res$coef.level1 # <- res$betahat.correst
        pval.level1[t, j] <- res$pval.level1 # <- res$betahat.pval.rhotest
        coef.level2[t, j] <- res$coef.level2 # <- res$gammahat2
        pval.level2[t, j] <- res$pval.level2 # <- res$gammahat.pval2
      }
      if (verbose) cat("\n")
    }
    xymats$coef.level1 <- coef.level1
    xymats$pval.level1 <- pval.level1
    xymats$coef.level2 <- coef.level2
    xymats$pval.level2 <- pval.level2
    xymats$match.by <- match.by
    xymats$delta <- delta
  }

  return(xymats)
}

#' Perform matched test
#'
#' This function takes gene expression and chromatin accessibility vectors for
#' a given trio and performs matching-based nonparametric tests.
#'
#' @param Yg a numeric vector containing RNA expression values
#' for a given gene g.
#' @param Xt a numeric vector containing chromatin accessibility
#' for a given peak t.
#' @param Yj a numeric vector containing TF expression for a given TF j.
#' @param Z a matrix containing covariates.
#' @inheritParams fitModel
#'
#' @return A list containing the following elements.
#' \item{matched.cells}{a data frame with columns containing metacell indices
#' in matche pairs, differences in the outcome Yg,
#' differences in the varied predictor V,
#' and the average of the matched predictor U}
#' \item{match.by}{a character string specifying the matched variable.}
#' \item{coef.level1}{a numerical scalar representing an estimate of coefficient
#' from TRIPOD level 1 test.}
#' \item{pval.level1}{a numerical scalar representing a p-value associated with
#' coef.level1}
#' \item{coef.level2}{a numerical scalar representing an estimate of coefficient
#' from TRIPOD level 2 test.}
#' \item{pval.level2}{a numerical scalar representing a p-value associated with
#' coef.level2}
#' \item{matched.col1}{a character vector specifying colors of the first metacells
#' in the pairs based on the cell types.}
#' \item{matched.col2}{a character vector specifying colors of the second metacells
#' in the pairs based on the cell types.}
#'
#' @import nbpMatching
#' @importFrom stats IQR cor cor.test lm
#' @export
performMatchedTest <- function(Yg, Xt, Yj, Z = NULL, match.by,
	                      cap.at.quantile, delta,
                        min.level1, min.level2,
	                      metacell.celltype.col = NULL,
                        # do.plot = TRUE, col = "black", plot.dir = NULL, file.name = NULL,
                        # partition.screen = TRUE,
	                      verbose = FALSE
) {
  to.match <- NULL
  if (match.by == "Xt") {
  	# set the matched predictor
    U <- Xt
    # set the varied predictor
    V <- Yj
  }
  if (match.by == "Yj") {
  	# set the matched predictor
    U <- Yj
    # set the varied predictor
    V <- Xt
  }
  if (!is.null(Z)) {
    to.match <- cbind(U, Z)
  } else {
    to.match <- matrix(data = U, nrow = length(U), ncol = 1)
  }
  if (is.null(to.match)) {
    stop("to.match value is not valid.")
  }
  # compute distances between the variables in to.match
  dmat <- computeDistances(to.match)
  # suppress the warning on odd number of rows
  options(warn = -1)
  # convert to distance matrix
  dmat2 <- distancematrix(dmat)
  options(warn = 0)
  # perform matching
  nbm <- nonbimatch(dmat2)
  epsilon <- quantile(abs(nbm$halves$Distance), 0.75) + 3 * IQR(abs(nbm$halves$Distance))
  # exclude matches with difference greater than epsilon
  which.keep <- which(abs(nbm$halves$Distance) <= epsilon)
  # exclude ghost matching
  which.keep <- setdiff(which.keep, c(which(nbm$halves$Group1.ID == "ghost"),
  	which(nbm$halves$Group2.ID == "ghost")))
  # construct a matrix with each row containing a matched pair
  matched.cells <- data.frame(matrix(nrow = length(which.keep), ncol = 8))
  names(matched.cells) <- c("ID1", "ID2", "Yg1", "Yg2", "U1", "U2", "V1", "V2")
  matched.cells[, 1] <- nbm$halves$Group1.Row[which.keep]
  matched.cells[, 2] <- nbm$halves$Group2.Row[which.keep]
  matched.cells[, 3] <- Yg[matched.cells[, 1]]
  matched.cells[, 4] <- Yg[matched.cells[, 2]]
  matched.cells[, 5] <- U[matched.cells[, 1]]
  matched.cells[, 6] <- U[matched.cells[, 2]]
  matched.cells[, 7] <- V[matched.cells[, 1]]
  matched.cells[, 8] <- V[matched.cells[, 2]]
  if (!is.null(metacell.celltype.col)) {
    matched.col1 <- metacell.celltype.col[nbm$halves$Group1.Row[which.keep]]
    matched.col2 <- metacell.celltype.col[nbm$halves$Group2.Row[which.keep]]
  }

  # compute differences and means
  dU <- matched.cells$U2 - matched.cells$U1
  meanU <- (matched.cells$U1 + matched.cells$U2) / 2
  dYg <- matched.cells$Yg2 - matched.cells$Yg1
  dV <- matched.cells$V2 - matched.cells$V1
  # cap values for robustness
  dYg <- capValues(dYg, cap.at.quantile)
  dV <- capValues(dV, cap.at.quantile)
  meanU <- capValues(meanU, cap.at.quantile)
  # get a logical vector to select pairs for which dV does not equal 0
  sel1 <- dV != 0
  sel1[is.na(sel1)] <- FALSE
  delta.val <- max(quantile(meanU, delta, na.rm = TRUE), 0.1)
  # get a logical vector to select pairs for which the mean of U is greater than a threshold
  sel2 <- meanU >= delta.val
  sel2[is.na(sel2)] <- FALSE
  matched.cells <- cbind(matched.cells, dU, meanU, dV, dYg, sel1, sel2)
  # num.pairs <- sum(sel2 & sel1) # not used

  coef.level1 <- NA # gammahat2 <- NA
  pval.level1 <- NA # betahat.correst <- NA
  coef.level2 <- NA # gammahat.pval2 <- NA
  pval.level2 <- NA # betahat.pval.rhotest <- NA

  # carry out level 1 testing only if the mean of U is greater than a threshold
  if (sum(sel2) > min.level1) {
  	# use only points for which the mean of U is greater than a threshold
    res <- cor.test(dYg[sel2], dV[sel2], method = "spearman")
    coef.level1 <- res$estimate
    pval.level1 <- res$p.value
    rm(res)
  }
  # check if there are sufficient number of pairs that satisfy the conditions
  if (sum(sel1 & sel2) > min.level2) {
  	# use all points for gamma
    res <- lm(dYg ~ dV + dV:meanU)
    coef.level2 <- summary(res)$coefficients[3, 1]
    pval.level2 <- summary(res)$coefficients[3, 4]
    rm(res)
  }

  results <- list(
    matched.cells = matched.cells, match.by = match.by,
    coef.level1 = coef.level1, pval.level1 = pval.level1,
    coef.level2 = coef.level2, pval.level2 = pval.level2
  )
  if (!is.null(metacell.celltype.col)) {
    results$matched.col1 <- matched.col1
    results$matched.col2 <- matched.col2
  }
  return(results)
}

#' Get trios with significantly nonzero coefficients for a given target gene
#'
#' This function takes a list returned by {\code{\link{fitModel}}}
#' as input and generate a data frame containing trios with significantly
#' positive or negative coefficients for a given target gene.
#'
#' @param xymats a list returned by {\code{\link{fitModel}}}.
#' @param fdr.thresh a numerical scalar setting a FDR threshold.
#' @param sign a character string specifying the sign.
#' This should be either "positive" or "negative".
#' @param coef.name a character string specifying the name of the coefficient.
#' This must be one of the names of the elements in the xymats object.
#' @param pval.name a character string specifying the name of the element in
#' the xymats object that stores the P value.
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{peak_num}{an integer representing a peak index.
#' The mapping between the indices and the peak names are the same as that
#' in the xymats object.}
#' \item{TF_num}{an integer representing a TF index.
#' The mapping between the indices and the TF names are the same as that
#' in the xymats object.}
#' \item{peak}{a chracter string representingthe name of the peak.}
#' \item{TF}{a chracter string representing the name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#' @export
getTriosForSingleGene <- function(
#	gene.name,
# xymats.list,
	xymats,
	fdr.thresh,
	sign, # "positive" or "negative"
  coef.name,
	pval.name
)	{
	gene.name <- xymats$gene.name
  peak.g <- rownames(xymats$peakxmotif.g)
  TF.g <- xymats$TF.g
  nonzero.peakxmotif.g <- xymats$nonzero.peakxmotif.g
  coef <- xymats[[coef.name]][nonzero.peakxmotif.g]
  pval <- xymats[[pval.name]][nonzero.peakxmotif.g]
  adj <- performBHVector(pval, fdr.thresh = fdr.thresh)$pval.c.adj
  if (sign == "positive") {
    which.rej <- which(coef > 0 & adj < fdr.thresh)
  } else if (sign == "negative") {
  		which.rej <- which(coef < 0 & adj < fdr.thresh)
  }
  df <- data.frame(gene = rep(gene.name, length(which.rej)),
    peak_num = nonzero.peakxmotif.g[which.rej, 1],
    TF_num = nonzero.peakxmotif.g[which.rej, 2],
    peak = peak.g[nonzero.peakxmotif.g[which.rej, 1]],
    TF = TF.g[nonzero.peakxmotif.g[which.rej, 2]],
    coef = coef[which.rej],
  	pval = pval[which.rej],
  	adj = adj[which.rej])
  rownames(df) <- NULL
  # remove cases where the TF is the same as the target gene
  df <- df[df$gene != df$TF, ]
  return(df)
}

#' Get trios with significantly nonzero coefficients across genes
#'
#' This is a wrapper function for {\code{\link{getTriosForSingleGene}}}.
#'
#' @param xymats.list a list of lists (xymats objects) returned
#' by {\code{\link{fitModel}}}.
#' @param genes a character vector containing the symbols of target genes
#' to be analyzed.
#' @param model.name a character string specifying the model.
#' @param level an integer specifying whether level 1 or 2 test is used.
#' @inheritParams getTriosForSingleGene
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{peak_num}{an integer representing a peak index.
#' The mapping between the indices and the peak names are the same as that
#' in the xymats object.}
#' \item{TF_num}{an integer representing a TF index.
#' The mapping between the indices and the TF names are the same as that
#' in the xymats object.}
#' \item{peak}{a chracter string representingthe name of the peak.}
#' \item{TF}{a chracter string representing the name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#'
#' @import BiocParallel
#' @export
getTrios <- function(
	xymats.list,
	genes = NULL, # a character vector
	fdr.thresh,
	sign, # "positive" or "negative"
  model.name,
	level = NULL
	){
	df <- data.frame()
	if (length(model.name) > 1) {
		warning("The model.name argument takes only one element. The first element is used.")
		model.name <- model.name[1]
	}
	if (all(model.name != c("marginal", "conditional.on.Yj", "conditional.on.Xt",
		"interaction", "TRIPOD"))
	) {
		stop('The model.name argument should be set to "marginal", "conditional.on.Yj",
			"conditional.on.Xt", "interaction", or "TRIPOD".')
	}
  if (model.name != "TRIPOD" && !is.null(level)){
	  warning('The level argument is applicable to "TRIPOD" only.')
	}

	if (model.name == "marginal"){
		coef.name <- "gamma.m"
		pval.name <- "pval.gamma.m"
	} else if (model.name == "conditional.on.Yj") {
		coef.name <- "alpha.c"
		pval.name <- "pval.alpha.c"
	} else if (model.name == "conditional.on.Xt") {
		coef.name <- "beta.c"
		pval.name <- "pval.beta.c"
	} else if (model.name == "interaction") {
		coef.name <- "gamma.i"
		pval.name <- "pval.gamma.i"
	} else if (model.name == "TRIPOD") {
		if (level == 1) {
			coef.name <- "coef.level1"
		  pval.name <- "pval.level1"
		} else if (level == 2){
			coef.name <- "coef.level2"
		  pval.name <- "pval.level2"
		} else {
			stop("The level argument should be either 1 or 2.")
		}
	} # else {
	# 	stop("The model.name argument is not valid.")
	# }
	# if (is.null(genes)) {genes <- sapply(xymats.list, function(x) x$gene.name)}
  if (!is.null(genes)) {
  	xymats.list <- xymats.list[sapply(xymats.list, function(x) x$gene.name %in% genes)] # need to test this
  }

	df.list <- bplapply(
		X = xymats.list,
		FUN = getTriosForSingleGene, # try()
	  fdr.thresh = fdr.thresh,
	  sign = sign, # "positive" or "negative"
    coef.name = coef.name,
	  pval.name = pval.name
	)
  df <- do.call("rbind", df.list)
  rownames(df) <- NULL
  return(df)
}

#' Get peak-gene pairs with significantly nonzero coefficients for a given target gene
#'
#' This function takes a list returned by {\code{\link{fitModel}}}
#' as input and generate a data frame containing peak-gene pairs with significantly
#' positive or negative coefficients for a given target gene.
#'
#' @inheritParams getTriosForSingleGene
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{peak_num}{an integer representing a peak index.
#' The mapping between the indices and the peak names are the same as that
#' in the xymats object.}
#' \item{peak}{a chracter string representingthe name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#'
#' @export
getPeakGenePairsForSingleGene <- function(
	xymats,
	coef.name,
	pval.name,
	sign,
	fdr.thresh
) {
	gene.name <- xymats$gene.name
  Xt <- xymats$Xt
  coef <- xymats[[coef.name]]
  pval <- xymats[[pval.name]]
  adj <- performBHVector(pval, fdr.thresh = fdr.thresh)$pval.c.adj
  if (sign == "positive"){
    which.rej <- which(coef > 0 & adj < fdr.thresh)
  } else if (sign == "negative"){
    which.rej <- which(coef < 0 & adj < fdr.thresh)
  }
  results <- data.frame(
	  gene = rep(gene.name, length(which.rej)),
  	peak_num = which.rej,
    peak = colnames(Xt)[which.rej],
    coef = coef[which.rej],
    pval = pval[which.rej],
    adj = adj[which.rej])
  rownames(results) <- NULL
  return(results)
}

#' Get peak-gene pairs with significantly nonzero coefficients across genes
#'
#' This is a wrapper function for {\code{\link{getPeakGenePairsForSingleGene}}}.
#'
#' @inheritParams getTrios
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{peak_num}{an integer representing a peak index.
#' The mapping between the indices and the peak names are the same as that
#' in the xymats object.}
#' \item{peak}{a chracter string representingthe name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#'
#' @import BiocParallel
#' @export
getPeakGenePairs <- function(
	xymats.list,
	fdr.thresh,
	sign, # "positive" or "negative"
	model.name
) {
	df <- data.frame()
	if (model.name != "marginal"){
		stop('This function is applicable to "marginal" only.')
	}
  coef.name <- "alpha.m"
	pval.name <- "pval.alpha.m"
  df.list <- bplapply(
  	X = xymats.list,
  	FUN = getPeakGenePairsForSingleGene,
  	coef.name = coef.name,
	  pval.name = pval.name,
	  sign = sign,
	  fdr.thresh = fdr.thresh
  )
  df <- do.call("rbind", df.list)
  rownames(df) <- NULL
  return(df)
}

#' Get TF-gene pairs with significantly nonzero coefficients for a given target gene
#'
#' This function takes a list returned by {\code{\link{fitModel}}}
#' as input and generate a data frame containing TF-gene pairs with significantly
#' positive or negative coefficients for a given target gene.
#'
#' @inheritParams getTriosForSingleGene
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{TF_num}{an integer representing a TF index.
#' The mapping between the indices and the TF names are the same as that
#' in the xymats object.}
#' \item{TF}{a chracter string representing the name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#'
#' @export
getTFGenePairsForSingleGene <- function(
  xymats,
	coef.name,
	pval.name,
	sign,
	fdr.thresh
) {
	gene.name <- xymats$gene.name
  Yj <- xymats$Yj
  coef <- xymats[[coef.name]]
  pval <- xymats[[pval.name]]
  adj <- performBHVector(pval, fdr.thresh = fdr.thresh)$pval.c.adj
  if (sign == "positive"){
    which.rej <- which(coef > 0 & adj < fdr.thresh)
  } else if (sign == "negative"){
    which.rej <- which(coef < 0 & adj < fdr.thresh)
  }
  results <- data.frame(
	  gene = rep(gene.name, length(which.rej)),
  	TF_num = which.rej,
    TF = colnames(Yj)[which.rej],
    coef = coef[which.rej],
    pval = pval[which.rej],
    adj = adj[which.rej])
  rownames(results) <- NULL
  results <- results[results$gene != results$TF, ]
  return(results)
}

#' Get TF-gene pairs with significantly nonzero coefficients across genes
#'
#' This is a wrapper function for {\code{\link{getTFGenePairsForSingleGene}}}.
#'
#' @inheritParams getTrios
#'
#' @return A data frame with 8 columns containing the following:
#' \item{gene}{a character string containing the name of the target gene.}
#' \item{TF_num}{an integer representing a TF index.
#' The mapping between the indices and the TF names are the same as that
#' in the xymats object.}
#' \item{TF}{a chracter string representing the name of the peak.}
#' \item{coef}{a numerical scalar representing the value of the coefficient}
#' \item{pval}{a numerical scalar representing the nominal p value}
#' \item{adj}{a numerical scalar representing the BH-adjusted p value}
#'
#' @import BiocParallel
#' @export
getTFGenePairs <- function(
	xymats.list,
	fdr.thresh,
	sign, # "positive" or "negative"
  model.name
	){
	df <- data.frame()
	if ((model.name) != "marginal"){
	  stop('This function is applicable to "marginal" only.')
	}
	coef.name <- "beta.m"
	pval.name <- "pval.beta.m"
  df.list <- bplapply(
  	X = xymats.list,
  	FUN = getTFGenePairsForSingleGene,
  	coef.name = coef.name,
	  pval.name = pval.name,
	  sign = sign,
	  fdr.thresh = fdr.thresh
  )
  df <- do.call("rbind", df.list)
  rownames(df) <- NULL
  return(df)
}

#' Perform RNA prediction for a given target gene
#'
#' This function takes a list returned by {\code{\link{getXYMatrices}}}
#' as input and performs RNA prediction for a given target gene.
#'
#' @param xymats a list returned by {\code{\link{getXYMatrices}}}.
#' @param model.name 	a character string specifying the model. This must be
#' one of "gene.activity", "peak.LASSO", or "peak.TF.LASSO".
#' @param log a logical indicating whether to log-transform the variables.
#' This must be one of "Yg", "Xt", and "Yj".
#' @param cap.at.quantile a quantile at which variables are capped to remove outliers.
#' Set this to 0 if capping is unnecessary.
#' @param seed an integer.
#'
#' @return A list containing some of the following elements.
#' \item{gene.name}{a character string.}
#' \item{Yg}{ a vector containig capped RNA expression of the gene across metacells.}
#' \item{Xt}{a matrix containig capped chromatin accessibility of the ATAC peaks.}
#' \item{Yj}{a matrix containing capped TF RNA expression.
#' Note that TFs that do no have motifs in the region
#' around TSS are excluded.}
#' \item{peakxmotif.g}{a gene-specific binary matrix containing indicators of
#' whether ATAC peaks contains motifs.
#' rows and columns represent peaks and motifs, respectively.}
#' \item{nonzero.peakxmotif.g}{a matrix with two columns.
#' The first and second columns contain row and column numbers
#' in the peakmotif.g matrix where the entry is non-zero.}
#' \item{Y.pred}{A numerical vector containing predicted values of Yg.}
#' \item{cor}{A numerical scalar representing Pearson correlation between
#' the predicted and observed Yg values}
#'
#' @import GenomicRanges Matrix glmnet
#' @importFrom stats predict
#'
#' @export
performRNAPrediction <- function(xymats, model.name, log = NULL,
	cap.at.quantile = 0.02, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

	if (all(model.name != c("gene.activity", "peak.LASSO", "peak.TF.LASSO"))) {
  	stop('The model.name argument must be one of "gene.activity", "peak.LASSO", or "peak.TF.LASSO".')
  }
	xymats$model.name <- model.name

  if (all(log %in% c("Yg", "Xt", "Yj"))) {
    if (any(log == "Yg")) {
        xymats$Yg <- log(xymats$Yg + 1)
    }
    if (any(log == "Xt")) {
      xymats$Xt <- log(xymats$Xt + 1)
    }
    if (any(log == "Yj")) {
      xymats$Yj <- log(xymats$Yj + 1)
    }
  } else {
    stop('The log argument must be one of "Yg", "Xt", and "Yj".')
  }

	# cap values
	if (cap.at.quantile > 0) {
    xymats$Yg <- capValues(xymats$Yg, cap.at.quantile = cap.at.quantile)
    xymats$Xt <- apply(xymats$Xt, 2, capValues, cap.at.quantile = cap.at.quantile)
    xymats$Yj <- apply(xymats$Yj, 2, capValues, cap.at.quantile = cap.at.quantile)
	}

	if (model.name == "gene.activity") {
    Y.pred <- as.matrix(apply(xymats$Xt, 1, sum))
    xymats$Y.pred <- as.numeric(Y.pred)
    xymats$cor <- cor(xymats$Yg, as.numeric(Y.pred), method = "pearson")
  }

  if (model.name == "peak.LASSO") {
  	# check whether there are at least two columns to perform LASSO
    if (ncol(xymats$Xt) > 1) {
    	# perform 10-fold cross-validation to find the optimal lambda
      lasso.mod <- glmnet(xymats$Xt, xymats$Yg, alpha = 1)
      cv.out <- cv.glmnet(xymats$Xt, xymats$Yg, alpha = 1)
      bestlam <- cv.out$lambda.min
      if (max(lasso.mod$df) > 1) {
        bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df > 1]))
      }
      # perform leave-one-out prediction with the selected lambda
      Y.pred <- rep(NA, length(xymats$Yg))
      for (i in seq(length(xymats$Yg))) {
        fit <- glmnet(xymats$Xt[-i, ], xymats$Yg[-i], alpha = 1, lambda = bestlam)
        Y.pred[i] <- predict(fit, newx = t(xymats$Xt[i, ]))
      }
      xymats$Y.pred <- as.numeric(Y.pred)
      xymats$cor <- cor(xymats$Yg, as.numeric(Y.pred), method = "pearson")
    } else {
    	xymats$model.name <- "peak.linear"
      Y.pred <- rep(NA, length(xymats$Yg))
      for (i in seq(length(xymats$Yg))) {
        fit <- lm(xymats$Yg[-i] ~ xymats$Xt[-i, ])
        Y.pred[i] <- xymats$Xt[i, ] * fit$coefficients[2] + fit$coefficients[1]
      }
      xymats$Y.pred <- as.numeric(Y.pred)
      xymats$cor <- round(cor(xymats$Y, as.numeric(Y.pred), method = "pearson"), 3)
    }
  }

  if (model.name == "peak.TF.LASSO") {
    X <- xymats$Xt[, xymats$nonzero.peakxmotif.g[, 1], drop = FALSE] *
    	xymats$Yj[, xymats$nonzero.peakxmotif.g[, 2], drop = FALSE]
    # check if Xt * Yj can be zero across all metacells even though they themselves are non-zeros
    peakTF.filter <- apply(X, 2, sum) > 0
    nonzero.peakxmotif.g <- xymats$nonzero.peakxmotif.g[peakTF.filter, , drop = FALSE]
    X <- X[, peakTF.filter, drop = FALSE]

    # check whether there are at least two columns to perform LASSO
    if (ncol(X) > 1) {
    	# perform 10-fold cross-validation to find the optimal lambda
      lasso.mod <- glmnet(X, xymats$Yg, alpha = 1)
      cv.out <- cv.glmnet(X, xymats$Yg, alpha = 1)
      bestlam <- cv.out$lambda.min
      # set a threshold to lambda
      bestlam <- min(bestlam, max(lasso.mod$lambda[lasso.mod$df > 1]))
      # perform leave-one-out prediction with the selected lambda
      Y.pred <- rep(NA, length(xymats$Yg))
      for (i in seq(length(xymats$Yg))) {
        fit <- glmnet(X[-i, ], xymats$Yg[-i], alpha = 1, lambda = bestlam)
        Y.pred[i] <- predict(fit, newx = t(X[i, ]))
      }
      xymats$Y.pred <- Y.pred
      xymats$cor <- round(cor(xymats$Yg, as.numeric(Y.pred), method = "pearson"), 3)
    } else if (ncol(X) == 1) {
    	xymats$model.name <- "peak.TF.linear"
      Y.pred <- rep(NA, length(xymats$Yg))
      for (i in seq(length(xymats$Yg))) {
        fit <- lm(xymats$Yg[-i] ~ X[-i, ])
        Y.pred[i] <- X[i, ] * fit$coefficients[2] + fit$coefficients[1]
      }
      xymats$Y.pred <- as.numeric(Y.pred)
      xymats$cor <- cor(xymats$Yg, as.numeric(Y.pred), method = "pearson")
    } else {
    	xymats$Y.pred <- NA
      xymats$cor <- NA
    }
  }

  return(xymats)
}

