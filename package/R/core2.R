#' Create input matrices for model fitting for a given target gene
#'
#' This function creates input matrices, Yg, Xp, Yt, for model fitting
#' for a given target gene.
#' It allows for cases where a given TF corresponds to multiple motifs and
#' those where a given motif corresponds to multiple TFs.
#' Cases where Yg and Yt are the same are excluded.
#'
#' @param gene.name a character string containing the gene name. This must exist
#'   in the `gene_name` column of the `transcripts.gr` object.
#' @param ext.upstream an integer representing the window size in bp.
#' @param ext.downstream an integer representing the window size in bp.
#' @param tripod an object returned by {\code{\link{getObjectsForModelFit2}}}.
#' @param metacell an object returned by {\code{\link{getMetacellMatrices}}}.
#'
#' @return A list with elements:
#' \item{gene.name}{a character string.}
#' \item{Yg}{a vector containig RNA expression of the gene across metacells.}
#' \item{Xp}{a matrix containig chromatin accessibility of the ATAC peaks.}
#' \item{Yt}{a matrix containing TF RNA expression.
#' note that this is gene-specific since TFs that do no have motifs in the region
#' around TSS are excluded.}
#' \item{peak.motif.TF.g}{a gene-specific data frame containing mapping between
#' ATAC peaks, motif, and TFs. Each row contains a given ATAC peak region, a motif
#' that is present in the region, and the corresponding TF.}
#'
#' @import GenomicRanges Matrix
#' @importFrom methods as
#'
#' @export
getXYMatrices2 <- function(gene.name,
                           ext.upstream,
                           ext.downstream = NULL,
                           tripod,
                           metacell) {
    transcripts.gr <- tripod$transcripts.gr
    peaks.gr <- tripod$peaks.gr
    motif.TF.map <- tripod$motif.TF.map
    peakxmotif <- tripod$peakxmotif

    metacell.rna <- metacell$rna
    metacell.peak <- metacell$peak

    if (is.null(ext.downstream)) ext.downstream <- ext.upstream
    transcripts.ext.gr <- promoters(transcripts.gr,
                                    upstream = ext.upstream,
                                    downstream = ext.downstream + 1)

    # get a vector of the target gene expression
    Yg <- metacell.rna[, gene.name]

    # get ATAC peaks within the vicinity region of the target gene
    transcripts.gr.g <- transcripts.gr[transcripts.gr$gene_name == gene.name]
    transcripts.ext.gr.g <- transcripts.ext.gr[transcripts.ext.gr$gene_name == gene.name]
    peaks.gr.g <- subsetByOverlaps(peaks.gr, transcripts.ext.gr.g)
    Xp <- metacell.peak[, overlapsAny(peaks.gr, transcripts.ext.gr.g), drop = FALSE]

    # get motifs within the vicinity region
    peakxmotif.g <- peakxmotif[overlapsAny(peaks.gr, transcripts.ext.gr.g), , drop = FALSE]

    # remove motifs that do not have binding sites in this region
    peakxmotif.g <- peakxmotif.g[, colSums(peakxmotif.g) > 0, drop = FALSE]

    # get a vector containing unique TF names
    TF.g <- unique(motif.TF.map$TF[motif.TF.map$motif %in% colnames(peakxmotif.g)])

    # get gene expression of the TFs
    Yt <- metacell.rna[, TF.g, drop = FALSE]

    # get a logical vector for removing TFs if it is not expressed in any cells
    TF.filter1 <- colSums(Yt) > 0
    # get a logical vector for removing a TF if it is the same as the target gene
    TF.filter2 <- TF.g != gene.name
    # get a logical vector
    TF.filter <- TF.filter1 & TF.filter2

    # remove the TFs
    Yt <- Yt[, TF.filter, drop = FALSE]
    TF.g <- TF.g[TF.filter, drop = FALSE]

    # get nonzero.peakxmotif.g
    peakxmotif.g <- as.matrix(peakxmotif.g)
    nonzero.peakxmotif.g <- which(peakxmotif.g != 0, arr.ind = TRUE)
    peakxmotif.g <- as(peakxmotif.g, "dgCMatrix")

    # get a mapping between  motifs and corresponding TFs
    # note that a given motif may correspond to multiple TFs and
    # that a given TF may correspond to multiple motifs
    peak.names <- rownames(peakxmotif.g)[nonzero.peakxmotif.g[, 1]]
    motif.names <- colnames(peakxmotif.g)[nonzero.peakxmotif.g[, 2]]
    peak.motif.TF.g <- data.frame(peak = peak.names,
                                  motif = motif.names)
    peak.motif.TF.g <- merge(peak.motif.TF.g, motif.TF.map, by.x = "motif", by.y = "motif")

    # filter the data based on the list of TFs
    peak.motif.TF.g <- peak.motif.TF.g[peak.motif.TF.g$TF %in% TF.g, ]

    results <- list(
        gene.name = gene.name,
        ext.upstream = ext.upstream,
        ext.downstream = ext.downstream,
        Yg = Yg,
        Xp = Xp,
        Yt = Yt,
        peak.motif.TF.g = peak.motif.TF.g)

    results
}

#' Perform model fitting for a given target gene
#'
#' This function takes a list returned by {\code{\link{getXYMatrices}}}
#' as input, performs model fitting for a given target gene,
#' and returns a list or data fram object containing coefficients and p-values.
#'
#' @param xymats a list returned by {\code{\link{getXYMatrices}}}.
#' @param model.name a character string specifying the model. This must be
#'   "marginal", "conditional", "interaction", "TRIPOD".
#' @param match.by a character string specifying the variable to be matched.
#'   This is necessary for "TRIPOD" only.
#' @param log a logical indicating whether to log-transform the variables.
#' @param cap.at.quantile a quantile at which variables are capped to remove outliers.
#'   That is, values are capped if they are below `cap.at.quantile`
#'   or above (1 - `cap.at.quantile`).
#'   Set this to 0 if capping is unnecessary.
#' @param delta a numeric.
#' @param min.Xp an integer to set the minimum number of metacells
#'   with nonzero chromatin accessibility, `Xp`.
#' @param min.Yt an integer to set the minimum number of metacells
#'   with nonzero TF expression, `Yt`.
#' @param min.level1 an integer to set the minimum number of pairs.
#' @param min.level2 an integer to set the minimum number of pairs.
#' @param metacell.celltype.col a character vector specifying colors of the metacells
#'   based on the cell types. This is necessary only for generating color-coded
#'   scatter plots of pairs of metacells using {\code{\link{plotGenePeakTFScatter}}}.
#' @param verbose a logical indicating whether to print progress in "TRIPOD".
#'
#' @return This function returns a list containing some of the following elements. Note that the sets of
#' elements differ depending on the `model.name` input argument.
#' \item{gene.name}{a character string.}
#' \item{ext.upstream}{an integer representing the window size in bp.}
#' \item{ext.downstream}{an integer representing the window size in bp.}
#' \item{Yg}{a vector containig capped RNA expression of the gene across metacells.}
#' \item{Xp}{a matrix containig capped chromatin accessibility of the ATAC peaks.}
#' \item{Yt}{a matrix containing capped TF RNA expression.
#' Note that TFs that do no have motifs in the region around TSS are excluded.}
#' \item{peak.motif.TF.g}{a gene-specific data frame containing mapping between
#' ATAC peaks, motif, and TFs. Each row contains a given ATAC peak region, a motif
#' that is present in the region, and the corresponding TF.}
#' \item{coef}{a list or data fram object containing coefficients and p-values.}
#' \item{delta}{This applies to "conditional" and "TRIPOD" only.}
#' \item{match.by}{the matched variable. This applieds to "TRIPOD" only.}
#'
#' @import GenomicRanges Matrix nbpMatching
#' @importFrom stats IQR cor cor.test lm model.frame
#'
#' @export
fitModel2 <- function(xymats, model.name, match.by = NULL, log = NULL,
	              cap.at.quantile = 0.02, delta = 0.2,
	              min.Xp = 5, min.Yt = 5, min.level1 = 10, min.level2 = 10,
	              metacell.celltype.col = NULL,
	              verbose = FALSE) {
    if (all(model.name != c("marginal", "conditional", "interaction", "TRIPOD"))) {
  	stop('The model.name argument must be one of "marginal", "conditional", "interaction", and "TRIPOD".')
    }
    if (!is.null(log)) {
        if (all(log %in% c("Yg", "Xp", "Yt"))) {
            if (any(log == "Yg")) {
                xymats$Yg <- log(xymats$Yg + 1)
            }
            if (any(log == "Xp")) {
                xymats$Xp <- log(xymats$Xp + 1)
            }
            if (any(log == "Yt")) {
                xymats$Yt <- log(xymats$Yt + 1)
            }
        } else {
            stop('The log argument must be one of "Yg", "Xp", and "Yt".')
        }
    }
    # cap values
    if (cap.at.quantile > 0) {
        xymats$Yg <- capValues(xymats$Yg, cap.at.quantile = cap.at.quantile)
        xymats$Xp <- apply(xymats$Xp, 2, capValues, cap.at.quantile = cap.at.quantile)
        xymats$Yt <- apply(xymats$Yt, 2, capValues, cap.at.quantile = cap.at.quantile)
    }

    peak.motif.TF.g <- xymats$peak.motif.TF.g
    n <- nrow(xymats$peak.motif.TF.g)
    Yg <- xymats$Yg

    if (model.name == "marginal") {
        nTFs <- ncol(xymats$Yt)
        npeaks <- ncol(xymats$Xp)

        coef <- list(
            alpha.m = rep(NA, npeaks),
            pval.alpha.m = rep(NA, npeaks),
            beta.m = rep(NA, nTFs),
            pval.beta.m = rep(NA, nTFs),
            gamma.m = rep(NA, n),
            pval.gamma.m = rep(NA, n),
            eta.m = rep(NA, n),
            pval.eta.m = rep(NA, n)
        )

        for (p in 1:npeaks) {
            Xp <- xymats$Xp[, p]
            # check the number of zeros in Xp to ensure that there is sufficient degrees of freedom
            if (sum(Xp > 0) < min.Xp) next
            res <- lm(Yg ~ Xp)
            if (any(is.na(res$coefficients))) next
            coef$alpha.m[p] <- res$coefficients[2]
            coef$pval.alpha.m[p] <- summary(res)$coefficients[2, 4]
        }
        for (t in 1:nTFs) {
            Yt <- xymats$Yt[, t]
            # check the number of zeros in Yt to ensure that there is sufficient degrees of freedom
            if (sum(Yt > 0) < min.Yt) next
            res <- lm(Yg ~ Yt)
            if (any(is.na(res$coefficients))) next
            coef$beta.m[t] <- res$coefficients[2]
            coef$pval.beta.m[t] <- summary(res)$coefficients[2, 4]
        }
        for (i in 1:n) {
            peak.name <- peak.motif.TF.g$peak[i]
            TF.name <- peak.motif.TF.g$TF[i]
            Xp <- xymats$Xp[, peak.name]
            Yt <- xymats$Yt[, TF.name]
            # check that not too many zero Y's so that there is enough degrees of freedom
            if (sum(Yt > 0) < min.Yt) next
            # check if Xp * Yt is zero across all metacells even though they themselves are non-zeros
            if (all(Xp * Yt == 0)) next
            if (sum(Xp > 0) < min.Yt) next
            res <- lm(Yg ~ Xp:Yt)
            if (all(!is.na(res$coefficients))) {
                coef$gamma.m[i] <- res$coefficients[2]
                coef$pval.gamma.m[i] <- summary(res)$coefficients[2, 4]
            }
            summary <- summary(lm(Yt ~ Xp))
            if (all(!is.na(summary$coefficients))) {
                coef$eta.m[i] <- summary$coefficients[2, 1]
                coef$pval.eta.m[i] <- summary$coefficients[2, 4]
            }
        }
        xymats$coef <- coef
    }

    if (model.name == "conditional") {
        coef <- data.frame(
            alpha.c = rep(NA, n),
            beta.c = rep(NA, n),
            pval.alpha.c = rep(NA, n),
            pval.beta.c = rep(NA, n)
        )
        for (i in 1:n) {
            peak.name <- peak.motif.TF.g$peak[i]
            TF.name <- peak.motif.TF.g$TF[i]
            Xp <- xymats$Xp[, peak.name]
            Yt <- xymats$Yt[, TF.name]
            # check that not too many zero Y's so that there is enough degrees of freedom
            if (sum(Yt > 0) < min.Yt) next
            # check if Xp * Yt is zero across all metacells even though they themselves are non-zeros
            if (all(Xp * Yt == 0)) next
            if (sum(Xp > 0) < min.Yt) next

            delta.val <- max(quantile(Xp, delta, na.rm = TRUE), 0.1)
            # carry out Yt testing only if Xp is greater than a threshold
            sel.Xp <- Xp >= delta.val
            res <- lm(Yg[sel.Xp] ~ 1 + Xp[sel.Xp] + Yt[sel.Xp])
            if (all(!is.na(res$coefficients))) {
                coef$beta.c[i] <- res$coefficients[3]
                coef$pval.beta.c[i] <- summary(res)$coefficients[3, 4]
            }

            delta.val <- max(quantile(Yt, delta, na.rm = TRUE), 0.1)
            # carry out Xp testing only if Yt is greater than a threshold
            sel.Yt <- Yt >= delta.val
            res <- lm(Yg[sel.Yt] ~ 1 + Xp[sel.Yt] + Yt[sel.Yt])
            if (all(!is.na(res$coefficients))) {
                coef$alpha.c[i] <- res$coefficients[2]
                coef$pval.alpha.c[i] <- summary(res)$coefficients[2, 4]
            }
        }
        xymats$coef <- coef
        xymats$delta <- delta
    }

    if (model.name == "interaction") {
        coef <- data.frame(
            alpha.i = rep(NA, n),
            pval.alpha.i = rep(NA, n),
            beta.i = rep(NA, n),
            pval.beta.i = rep(NA, n),
            gamma.i = rep(NA, n),
            pval.gamma.i = rep(NA, n)
        )

        # iterate over the rows in peak.motif.TF.g (peak, motif, TF) (similar to nonzero.peakxmotif.g)
        for (i in 1:n) {
            peak.name <- peak.motif.TF.g$peak[i]
            TF.name <- peak.motif.TF.g$TF[i]
            Xp <- xymats$Xp[, peak.name]
            Yt <- xymats$Yt[, TF.name]
            # check that not too many zero Y's so that there is enough degrees of freedom
            if (sum(Yt > 0) < min.Yt) next
            # check if Xp * Yt is zero across all metacells even though they themselves are non-zeros
            if (all(Xp * Yt == 0)) next
            if (sum(Xp > 0) < min.Yt) next
            res <- lm(Yg ~ Xp * Yt)
            if (any(is.na(res$coefficients))) next

            coef$alpha.i[i]  <- res$coefficients[2]
            coef$beta.i[i] <- res$coefficients[3]
            coef$gamma.i[i] <- res$coefficients[4]
            coef$pval.alpha.i[i] <- summary(res)$coefficients[2, 4]
            coef$pval.beta.i[i] <- summary(res)$coefficients[3, 4]
            coef$pval.gamma.i[i] <- summary(res)$coefficients[4, 4]
        }
        xymats$coef <- coef
    }

    if (model.name == "TRIPOD") {
        coef <- data.frame(
            coef.level1 = rep(NA, n),
            pval.level1 = rep(NA, n),
            coef.level2 = rep(NA, n),
            pval.level2 = rep(NA, n)
        )

        # iterate over the rows in peak.motif.TF.g (peak, motif, TF) (similar to nonzero.peakxmotif.g)
        for (i in 1:n) {
            peak.name <- peak.motif.TF.g$peak[i]
            TF.name <- peak.motif.TF.g$TF[i]
            Xp <- xymats$Xp[, peak.name]
            Yt <- xymats$Yt[, TF.name]
            res <- performMatchedTest2(
                Xp = Xp, Yt = Yt, Yg = Yg,
                match.by = match.by,
                cap.at.quantile = cap.at.quantile, delta = delta,
                min.level1 = min.level1, min.level2 = min.level2,
                verbose = verbose
            )
            coef$coef.level1[i] <- res$coef.level1
            coef$pval.level1[i] <- res$pval.level1
            coef$coef.level2[i] <- res$coef.level2
            coef$pval.level2[i] <- res$pval.level2
        }
        xymats$coef <- coef
        xymats$match.by <- match.by
        xymats$delta <- delta
    }
    xymats
}

#' Perform matched test
#'
#' This function takes gene expression and chromatin accessibility vectors for
#' a given trio and performs matching-based nonparametric tests.
#'
#' @param Yg a numeric vector containing RNA expression values
#' for a given gene g.
#' @param Xp a numeric vector containing chromatin accessibility
#' for a given peak t.
#' @param Yt a numeric vector containing TF expression for a given TF j.
#' @param Z a matrix containing covariates.
#' @inheritParams fitModel2
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
performMatchedTest2 <- function(Yg, Xp, Yt,
                                Z = NULL,
                                match.by,
                                cap.at.quantile, delta,
                                min.level1, min.level2,
                                metacell.celltype.col = NULL,
	                        verbose = FALSE) {
  to.match <- NULL
  if (match.by == "Xp") {
  	# set the matched predictor
    U <- Xp
    # set the varied predictor
    V <- Yt
  }
  if (match.by == "Yt") {
  	# set the matched predictor
    U <- Yt
    # set the varied predictor
    V <- Xp
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
  results
}
