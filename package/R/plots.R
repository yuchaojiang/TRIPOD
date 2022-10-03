#' Plot Cook's distance
#'
#' @param d.list a list returned by {\code{\link{getCooksD}}}
#' @param ... further arguments to be passed.
#' @inheritParams getXYMatrices
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotCooksD <- function(d.list, metacell.celltype.col, ...) {
	d <- d.list$d
	k <- d.list$k
  p <- ggplot(d, aes(x = obs, y = cd, label = txt, ymin = 0, ymax = cd)) +
    geom_linerange(colour = metacell.celltype.col) +
    geom_point(shape = 16, colour = metacell.celltype.col) +
    geom_hline(yintercept = k$ts, colour = "red") +
    xlab("Metacell") +
    ylab("Cook's D") +
    ggtitle("Cook's D Chart for Yg") +
    geom_text(size = 3, family = "sans", colour = metacell.celltype.col, na.rm = TRUE) +
    annotate("text",
      x = Inf, y = Inf, hjust = 1.2, vjust = 1,
      family = "sans", colour = "darkred", label = paste("Threshold:", round(k$ts, 3))
    )
	return(p)
}

#' Plot DFFIT
#'
#' @param d.list a list returned by {\code{\link{getDFFIT}}}
#' @param ... further arguments to be passed.
#' @inheritParams getXYMatrices
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotDFFIT <- function(d.list, metacell.celltype.col, ...) {
	d <- d.list$d
	dffits.t <- d.list$dffits.t
  title <- d.list$title
  p <- ggplot(d, aes(x = obs, y = dbetas, label = txt, ymin = 0, ymax = dbetas)) +
    geom_linerange(colour = metacell.celltype.col) +
    geom_hline(yintercept = c(0, dffits.t, -dffits.t), colour = "red") +
    geom_point(colour = metacell.celltype.col, shape = 16) +
    xlab("Metacell") +
    ylab("DFFITS") +
    ggtitle(paste("DFFIT for", title)) +
    geom_text(
      hjust = -0.2, nudge_x = 0.15, size = 3, family = "sans",
      colour = metacell.celltype.col, na.rm = TRUE
    ) +
    annotate("text",
      x = Inf, y = Inf, hjust = 1.5, vjust = 2,
      family = "sans", colour = "darkred",
      label = paste("Threshold:", round(dffits.t, 2))
    )
	return(p)
}

#' Plot sampling-based p-values for deleting a group of metacells
#'
#' @param delta.coeff.pval a data frame returned by {\code{\link{deleteNeighbor}}}
#' @param coefficient a character string specifying the coefficient to be plotted.
#' This must be one of "Xt", "Yj", and "Xt:Yj".
#' If {\code{\link{deleteNeighbor}}} was run with an argument "two.sided",
#' this can be set to fitted values, "Yg".
#' @param ... further arguments to be passed.
#' @param threshold a value of nominal p-value at which a horizontal line is drawn.
#' @inheritParams getXYMatrices
#' @inheritParams testInfluence
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotInfluenceForNeighbor <- function(delta.coeff.pval, coefficient, metacell.celltype,
	metacell.celltype.col, threshold = 0.05, ...) {
	items <- c("Xt", "Yj", "Xt:Yj", "Yg")
	if (all(coefficient != items)) {
    stop('The coefficient argument must be one of "Xt", "Yj", "Xt:Yj", or "Yg".')
	}
	if (coefficient == "Yg" & ncol(delta.coeff.pval) == 4) {
		stop('The coefficient argument must be one of "Xt", "Yj", or "Xt:Yj" for one-sided tests.')
	}
  index <- which(items == coefficient) + 1
	d <- data.frame(
    obs = 1:length(metacell.celltype),
    cd = -log(delta.coeff.pval[, index], 10),
    txt = metacell.celltype
  )
	d$txt[d$cd < (-log(threshold, 10))] <- NA
	title <- paste("Metacell-specific p-values for ", colnames(delta.coeff.pval)[index], sep = "")
  p <- ggplot(d, aes(x = obs, y = cd, label = txt, ymin = 0, ymax = cd)) +
    geom_linerange(colour = metacell.celltype.col) +
    geom_point(shape = 16, colour = metacell.celltype.col) +
    geom_hline(yintercept = -log(threshold, 10), colour = "black", linetype = "dotted") +
    xlab("Metacell") +
    ylab("-log(p-value)") +
    ggtitle(title) +
    geom_text(size = 3, family = "sans",
      colour = metacell.celltype.col, na.rm = TRUE)
  return(p)
}

#' Generate a UMAP plot colored by sampling-based P values
#'
#' This function generates a UMAP plot colored by sampling-based P values
#' for deleting a group of metacells.
#'
#' @param object a Seurat object.
#' @param reduction a character string specifying the dimension reduction method.
#' The coordinates should be stored in the Seurat object.
#' @param ... further arguments to be passed.
#' @inheritParams plotInfluenceForNeighbor
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotFeatureForNeighbor <- function(delta.coeff.pval, object, reduction,
	coefficient, ...) {
	items <- c("Xt", "Yj", "Xt:Yj", "Yg")
	if (all(coefficient != items)) {
    stop('The coefficient argument must be one of "Xt", "Yj", "Xt:Yj", or "Yg".')
	}
	if (coefficient == "Yg" & ncol(delta.coeff.pval) == 4) {
		stop('The coefficient argument must be one of "Xt", "Yj", or "Xt:Yj" for one-sided tests.')
	}
  index <- which(items == coefficient) + 1
	object$metacell.samp.logpval <- -log(delta.coeff.pval[object$seurat_clusters, index], 10)
  p <- FeaturePlot(object, features = "metacell.samp.logpval",
  	reduction = reduction, label.size = 2, ...) +
    ggtitle(paste("Metacell-specific p-values for", colnames(delta.coeff.pval)[index], "from sampling")) +
    theme(plot.title = element_text(size = 10))
  return(p)
}

#' Plot sampling-based P-values for deleting a celltype
#'
#' @param ordered.celltype a chracter vector containing cell types. This
#' specifies the order of the cell types in the plot.
#' @inheritParams plotInfluenceForNeighbor
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plotInfluenceForCellType <- function(delta.coeff.pval, coefficient, threshold = 0.05,
	metacell.celltype, metacell.celltype.col, ordered.celltype) {
	items <- c("Xt", "Yj", "Xt:Yj", "Yg")
	if (all(coefficient != items)) {
    stop('The coefficient argument must be one of "Xt", "Yj", "Xt:Yj", or "Yg".')
	}
	if (coefficient == "Yg" & ncol(delta.coeff.pval) == 4) {
		stop('The coefficient argument must be one of "Xt", "Yj", or "Xt:Yj" for one-sided tests.')
	}
  index <- which(items == coefficient) + 1
  d <- data.frame(
    obs = 1:length(unique(metacell.celltype)),
    cd = -log(delta.coeff.pval[, index], 10),
    txt = unique(metacell.celltype))
  d$txt[d$cd < (-log(threshold, 10))] <- NA
  d <- d[ordered.celltype, ]
  d$obs <- 1:length(d$obs)
  p <- ggplot(d, aes(x = obs, y = cd, label = txt, ymin = 0, ymax = cd)) +
    geom_linerange(
      colour = unique(metacell.celltype.col)[
      match(rownames(d), unique(metacell.celltype))]) +
    geom_point(
      shape = 16,
      colour = unique(metacell.celltype.col)[
      match(rownames(d), unique(metacell.celltype))]) +
    geom_hline(yintercept = -log(threshold, 10), colour = "black", linetype = "dotted") +
    xlab("Cell types") + ylab("-log(p-value)") +
    ggtitle(paste0("Cell-type-specific p-values for ", colnames(delta.coeff.pval)[index])) +
    geom_text(size = 3, family = "sans",
      colour = unique(metacell.celltype.col)[
        match(rownames(d), unique(metacell.celltype))], na.rm = TRUE)
	return(p)
}

#' Plot sampling-based p-values for deleting a cell type
#'
#' @param dend a dendrogram.
#' @inheritParams plotInfluenceForNeighbor
#'
#' @return A ggplot object.
#'
#' @import dendextend
#' @importFrom graphics legend
#' @export
plotTree <- function(dend, delta.coeff.pval, coefficient,
	metacell.celltype, metacell.celltype.col, threshold = 0.05, ...) {
  par(mar = c(6.1, 4.1, 2.1, 2.1))
	items <- c("Xt", "Yj", "Xt:Yj", "Yg")
	if (all(coefficient != items)) {
    stop('The coefficient argument must be one of "Xt", "Yj", "Xt:Yj", or "Yg".')
	}
	if (coefficient == "Yg" & ncol(delta.coeff.pval) == 4) {
		stop('The coefficient argument must be one of "Xt", "Yj", or "Xt:Yj" for one-sided tests.')
	}
  index <- which(items == coefficient) + 1
  # get subtrees, each of which corresponds to a split
  subtrees <- partition_leaves(dend)
  nodes_col <- rep("lightgray", length(subtrees))
  nodes_col[which(performBHVector(delta.coeff.pval[, index])$pval.c.adj < threshold)] <- "red"
  nodes_col[which(is.na(delta.coeff.pval[, index]))] <- 1
  nodes_pch <- rep(19, length(subtrees))
  nodes_pch[which(is.na(delta.coeff.pval[, index]))] <- 4
  metacell.celltype.col.unique <- unique(metacell.celltype.col)
  names(metacell.celltype.col.unique) <- unique(metacell.celltype)
  dend %>%
    set("nodes_pch", nodes_pch) %>%
    set("nodes_col", nodes_col) %>%
    set("nodes_cex", 2) %>%
    set("labels_col", metacell.celltype.col.unique[dend %>% labels()]) %>%
    # change color
    set("labels_cex", 1) %>%
    # change size
    plot(main = paste("Hierarchical clustering with testing of", colnames(delta.coeff.pval)[index]))
  legend("topright", col = c(1, 2, "lightgray"), pch = c(4, 19, 19),
  legend = c("NA", "Significant", "Not significant"), bty = "n", cex = 0.8)
  par(mar = c(5.1, 4.1, 4.1, 2.1))
}

# plotFootprint <- function(object, features, assay = NULL, group.by = NULL, idents = NULL,
#                           label = TRUE, repel = TRUE, show.expected = TRUE, normalization = "subtract",
#                           label.top = 3, label.idents = NULL) {
#   plot.data <- GetFootprintData(
#     object = object, features = features,
#     assay = assay, group.by = group.by, idents = idents
#   )
#   motif.sizes <- Signac:::GetMotifSize(
#     object = object, features = features,
#     assay = assay
#   )
#   obs <- plot.data[plot.data$class == "Observed", ]
#   expect <- plot.data[plot.data$class == "Expected", ]
#   base <- ceiling(motif.sizes / 2)
#   obs$flanks <- sapply(
#     X = seq_len(length.out = nrow(x = obs)),
#     FUN = function(x) {
#       pos <- abs(obs[x, "position"])
#       size <- base[[obs[x, "feature"]]]
#       return((pos > size) & (pos < (size + 50)))
#     }
#   )
#   if (!is.null(normalization)) {
#     correction.vec <- expect$norm.value
#     names(correction.vec) <- paste(expect$position, expect$feature)
#     if (normalization == "subtract") {
#       obs$norm.value <- obs$norm.value - correction.vec[paste(
#         obs$position,
#         obs$feature
#       )]
#     } else if (normalization == "divide") {
#       obs$norm.value <- obs$norm.value / correction.vec[paste(
#         obs$position,
#         obs$feature
#       )]
#     } else {
#       stop("Unknown normalization method requested")
#     }
#   }
#   flanks <- obs[obs$flanks, ]
#   flanks <- group_by(.data = flanks, feature, group)
#   flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
#   topmean <- top_n(x = flankmeans, n = label.top, wt = mn)
#   ymax <- top_n(x = flankmeans, n = 1, wt = mn)
#   ymin <- top_n(x = flankmeans, n = 1, wt = -mn)
#   label.df <- data.frame()
#   sub <- obs[obs$position == 75, ]
#   for (i in seq_along(along.with = features)) {
#     if (is.null(x = label.idents)) {
#       groups.use <- topmean[topmean$feature == features[[i]], ]$group
#     } else {
#       groups.use <- label.idents
#     }
#     df.sub <- sub[(sub$feature == features[[i]]) & (sub$group %in%
#       groups.use), ]
#     label.df <- rbind(label.df, df.sub)
#   }
#   obs$label <- NA
#   label.df$label <- label.df$group
#   obs <- rbind(obs, label.df)
#   plotlist <- list()
#   for (i in seq_along(along.with = features)) {
#     df <- obs[obs$feature == features[[i]], ]
#     min.use <- ifelse(test = normalization == "subtract",
#       yes = -0.5, no = 0.5
#     )
#     axis.min <- min(min.use, ymin[ymin$feature == features[[i]], ]$mn)
#     axis.max <- ymax[ymax$feature == features[[i]], ]$mn +
#       0.5
#     p <- ggplot(data = df, mapping = aes(
#       x = position, y = norm.value,
#       color = group, label = label
#     ))
#     p <- p + geom_line(size = 0.2) + xlab("Distance from motif") +
#       ylab(label = "Tn5 insertion\nenrichment") + ylim(c(axis.min, axis.max)) #+ theme(legend.position = "none")#+ guides(color = guide_legend(override.aes = list(size = 1)))
#     if (label) {
#       if (repel) {
#         p <- p + ggrepel::geom_label_repel(
#           box.padding = 0.5,
#           show.legend = FALSE
#         )
#       } else {
#         p <- p + geom_label(show.legend = FALSE)
#       }
#     }
#     if (show.expected) {
#       df <- expect[expect$feature == features[[i]], ]
#       p1 <- ggplot(data = df, mapping = aes(
#         x = position,
#         y = norm.value
#       )) +
#         geom_line(size = 0.2) +
#         xlab("Distance from motif") +
#         ylab(label = "Expected\nTn5 enrichment") +
#         theme_classic()
#       p <- p + theme(
#         axis.title.x = element_blank(), axis.text.x = element_blank(),
#         axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank()
#       )
#       p <- p + p1 + plot_layout(ncol = 1, heights = c(3, 1))
#     }
#     plotlist[[i]] <- p
#   }
#   plots <- wrap_plots(plotlist)
#   return(plots)
# }

#' Generate UMAP plots for a trio
#'
#' This is a wrapper function of {\code{\link{FeaturePlot}}} from the Seurat
#' package for generating UMAP plots, in which colors represent degrees of
#' target gene expression,
#' chromatin accessibility of an ATAC peak region, and TF expression
#' across cells.
#'
#' @param label a logical scalar indicating whether to show cluster labels.
#' @param ... further arguments passed to {\code{\link{FeaturePlot}}}.
#' @inheritParams plotGenePeakTFDot
#' @inheritParams plotFeatureForNeighbor
#'
#' @export
plotGenePeakTFUMAP <- function(object, gene.name, peak.name, TF.name,
	reduction = "wnn.umap", label = FALSE,
	size.title = 15, size.axis = 15, ncol = 3, ...
) {
  p1 <- FeaturePlot(object, features = paste0("sct_", gene.name),
  	reduction = reduction, label = label, label.size = 2,
  	cols = c("lightgrey", "darkblue"), ...) +
  	ggtitle("Gene expression") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis))
  p2 <- FeaturePlot(object, features = peak.name,
  	reduction = reduction, label = label, label.size = 2,
  	cols = c("lightgrey", "darkgreen"), ...) +
  	ggtitle("Peak accessibility") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis))
  p3 <- FeaturePlot(object, features = paste0("sct_", TF.name),
  	reduction = reduction, label = label, label.size = 2,
  	cols = c("lightgrey", "darkred"), ...) +
  	ggtitle("TF expression") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis))
  plot <- p1 + p2 + p3 + plot_layout(ncol = ncol)
  print(plot)
}

#' Generate violin plots for a trio
#'
#' This is a wrapper function of {\code{\link{VlnPlot}}} from the Seurat
#' package for generating violin plots of target gene expression,
#' chromatin accessibility of an ATAC peak region, and TF expression
#' across specificed groups of cells.
#'
#' @param ... further arguments passed to {\code{\link{VlnPlot}}}.
#' @inheritParams plotGenePeakTFDot
#'
#' @import ggplot2 patchwork Seurat
#' @export
plotGenePeakTFVln <- function(object, gene.name, peak.name, TF.name,
	group.by, size.title = 15, size.axis = 15, ncol = 3, ...
) {
  p1 <- VlnPlot(object, features = paste0("rna_", gene.name),
      slot = "counts", log = TRUE, group.by = group.by, ...) +
      ggtitle("Gene expression") + ylab("Expression") +
      theme(plot.title = element_text(size = size.title, face = "plain"),
      	axis.title = element_text(size = size.axis),
      	legend.position = "none")
  p2 <- VlnPlot(object, features = peak.name,
      slot = "counts", log = TRUE, group.by = group.by, ...) +
      ggtitle("Peak accessibility") + ylab("Peak accessibility") +
  	  theme(plot.title = element_text(size = size.title, face = "plain"),
  	  	axis.title = element_text(size = size.axis),
  	  	legend.position = "none")
  p3 <- VlnPlot(object, features = paste0("rna_", TF.name),
    	slot = "counts", log = TRUE, group.by = group.by, ...) +
    	ggtitle("TF expression") + ylab("Expression") +
  	  theme(plot.title = element_text(size = size.title, face = "plain"),
  	  	axis.title = element_text(size = size.axis),
  	  	legend.position = "none")
  plot <- p1 + p2 + p3 + plot_layout(ncol = ncol)
  print(plot)
}

#' Generate dot plots for a trio
#'
#' This is a wrapper function of {\code{\link{DotPlot}}} from the Seurat
#' package for generating dot plots of target gene expression,
#' chromatin accessibility of an ATAC peak region, and TF expression
#' across specificed groups of cells.
#'
#' @param object a Seurat object.
#' @param gene.name a character string representing a target gene name.
#' @param group.by a character string specifying a meta data column
#' in the Seurat object containing groups on which plot stratification is based.
#' @param size.title a numerical scalar specifying the size of the titles
#' of the plots.
#' @param size.axis a numerical scalar specifying the size of the axis titles.
#' @param size.legend a numerical scalar specifying the size of the legend titles.
#' @param ncol an integer specifying the number of columns for arranging plots.
#' @param ... further arguments passed to {\code{\link{DotPlot}}}.
#' @inheritParams plotGenePeakTFScatter
#'
#' @import ggplot2 patchwork Seurat
#' @export
plotGenePeakTFDot <- function(object, gene.name, peak.name, TF.name,
	group.by, size.title = 15, size.axis = 15, size.legend = 15, ncol = 3, ...
) {
  p1 <- DotPlot(object, features = paste0("sct_", gene.name),
    cols = c("lightgrey", "darkblue"), group.by = group.by, ...) +
    ggtitle("Gene expression") + ylab("Expression") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis),
    	legend.title = element_text(size = size.axis))
  p2 <- DotPlot(object, features = peak.name,
    cols = c("lightgrey", "darkgreen"), group.by = group.by, ...) +
    ggtitle("Peak accessibility") + ylab("Peak accessibility") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis),
    	legend.title = element_text(size = size.axis))
  p3 <- DotPlot(object, features = paste0("sct_", TF.name),
    cols = c("lightgrey", "darkred"), group.by = group.by, ...) +
    ggtitle("TF Expression") + ylab("Expression") +
    theme(plot.title = element_text(size = size.title, face = "plain"),
    	axis.title = element_text(size = size.axis),
    	legend.title = element_text(size = size.axis))
  plot <- p1 + p2 + p3 + plot_layout(ncol = ncol)
  print(plot)
}

#' Generate scatter plots for analyzing regulatory relationships
#'
#' This function generated scatter plots of Yg, Xt, Yj, difference, or partial residuals.
#' Pearson correlation coefficient and associated P value are also printed.
#'
#' @param peak.name a character string containing a peak name,
#' such as "chr15-79201373-79203205".
#' @param TF.name a character string containing a TF name.
#' @param to.plot a character string specifying the model to be analyzed.
#' This must be one of "marginal.Yg.Xt", "marginal.Yg.Yj", "marginal.Yj.Xt",
#' "product", "conditional.on.Yj", "conditional.on.Xt", "interaction", or "TRIPOD".
#' @param level a numerical value specifying the level for TRIPOD. This must be
#' either 1 or 2.
#' @param metacell.celltype a character vector specifying cell types of the metacells.
#' @param metacell.celltype.col a character vector representing cell types of metacells.
#' @inheritParams fitModel
#'
#' @import nbpMatching
#' @importFrom graphics grid text points
#' @importFrom stats residuals
#' @export
plotGenePeakTFScatter <- function(
	xymats, peak.name, TF.name,
  to.plot, level = NULL, match.by = NULL,
	cap.at.quantile = 0.02, delta = 0.2,
	min.Xt = 5, min.Yj = 5, min.level1 = 10, min.level2 = 10,
  # generate.pdf = TRUE, plot.dir = getwd(), partition.screen = TRUE,
	metacell.celltype = NULL, metacell.celltype.col = "black"
) {
	if (all(to.plot != c("marginal.Yg.Xt", "marginal.Yg.Yj", "marginal.Yj.Xt",
		"product", "conditional.on.Yj", "conditional.on.Xt", "interaction", "TRIPOD"))
	) {
		stop('The to.plot argument must be one of "marginal.Yg.Xt", "marginal.Yg.Yj", "marginal.Yj.Xt",
		"product", "conditional.on.Yj", "conditional.on.Xt", "interaction", and "TRIPOD".')
	}
	if (to.plot != "TRIPOD" && !is.null(level)){
	  warning('The level argument is applicable to "TRIPOD" only.')
	}
	if (to.plot != "TRIPOD" && !is.null(match.by)){
	  warning('The match.by argument is applicable to "TRIPOD" only.')
	}
	peak.num <- which(rownames(xymats$peakxmotif) == peak.name) # need to check this
	TF.num <- which(xymats$TF.g == TF.name) # need to check this
	Yg <- xymats$Yg # may or may not need to be capped
	Xt <- xymats$Xt[, peak.num] # may or may not need to be capped
	Yj <- xymats$Yj[, TF.num] # may or may not need to be capped
	if (cap.at.quantile > 0) {
    Yg <- capValues(Yg, cap.at.quantile = cap.at.quantile)
    Xt <- capValues(Xt, cap.at.quantile = cap.at.quantile)
    Yj <- capValues(Yj, cap.at.quantile = cap.at.quantile)
	}

  if (to.plot == "marginal.Yg.Xt") {
    # plot Yg against Xt
  	tmp <- cor.test(Xt, Yg)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Marginal, Gene expr. VS Peak access.:\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(Xt, Yg,
      xlab = "Xt", ylab = "Yg",
      main = title, cex.main = 1, font.main = 1,
      col = metacell.celltype.col, pch = 16
    )
    grid()
    if (!is.null(metacell.celltype.col)) {
      text(Xt, Yg, metacell.celltype, col = "gray", cex = 0.6)
    	points(Xt, Yg,
    		col = metacell.celltype.col, pch = 16)
    }
    rm(title)
  }

  if (to.plot == "marginal.Yg.Yj") {
    # plot Yg against Yj
    tmp <- cor.test(Yj, Yg)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Marginal, Gene expr. VS TF expr.:\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(Yj, Yg,
      xlab = "Yj", ylab = "Yg",
    	main = title, cex.main = 1, font.main = 1,
    	col = metacell.celltype.col, pch = 16
    )
    grid()
    if (!is.null(metacell.celltype.col)) {
      text(Yj, Yg, metacell.celltype, col = "gray", cex = 0.6)
    	points(Yj, Yg,
    	  col = metacell.celltype.col, pch = 16)
    }
    rm(title)
  }

  if (to.plot == "marginal.Yj.Xt") {
    # plot Yj against Xt
    tmp <- cor.test(Xt, Yj)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Marginal, TF expr. VS Peak access.:\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(Xt, Yj,
      xlab = "Xt", ylab = "Yj",
      main = title, cex.main = 1, font.main = 1,
    	col = metacell.celltype.col, pch = 16
    )
    grid()
    if (!is.null(metacell.celltype.col)) {
    	text(Xt, Yj, metacell.celltype, col = "gray", cex = 0.6)
    	points(Xt, Yj,
    	  col = metacell.celltype.col, pch = 16)
    }
    rm(title)
  }

  if (to.plot == "product") {
    # plot Yg against Xt * Yj
  	tmp <- cor.test(Xt * Yj, Yg)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Gene expr. VS (peak access. x TF expr.):\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(Xt * Yj, Yg,
      xlab = "Xt * Yj", ylab = "Yg",
      main = title, cex.main = 1, font.main = 1,
    	col = metacell.celltype.col, pch = 16
    )
    grid()
    if (!is.null(metacell.celltype.col)) {
      text(Xt * Yj, Yg, metacell.celltype, col = "gray", cex = 0.6)
    	points(Xt * Yj, Yg,
        col = metacell.celltype.col, pch = 16)
    }
    rm(title)
  }

  if (to.plot == "conditional.on.Yj") {
    # compute partial residuals
    remainY <- residuals(lm(Yg ~ Yj))
    remainX <- residuals(lm(Xt ~ Yj))
    tmp <- cor.test(remainX, remainY)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Conditional on TF expr., Gene expr. VS Peak access.:\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(remainX, remainY,
      xlab = "Xt partial resid", ylab = "Yg partial resid",
      col = metacell.celltype.col, pch = 16,
      main = title, cex.main = 1, font.main = 1
    )
    grid()
    if (!is.null(metacell.celltype)) {
      text(remainX, remainY, metacell.celltype, col = "gray", cex = 0.6)
    	points(remainX, remainY,
    		col = metacell.celltype.col, pch = 16)
    }
    delta.val <- max(quantile(Yj, xymats$delta, na.rm = TRUE), 0.1)
    sel.Yj <- Yj >= delta.val
    points(remainX[sel.Yj], remainY[sel.Yj], pch = 1)
    rm(tmp); rm(title);  rm(delta.val)
  }

  if (to.plot == "conditional.on.Xt") {
    # compute partial residuals
    remainY <- residuals(lm(Yg ~ Xt))
    remainX <- residuals(lm(Yj ~ Xt))
    tmp <- cor.test(remainX, remainY)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Conditional on peak access., Gene expr. VS TF expr.:\n",
      "r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(remainX, remainY,
      xlab = "Yj partial resid", ylab = "Yg partial resid",
      col = metacell.celltype.col, pch = 16,
      main = title, cex.main = 1, font.main = 1
    )
    grid()
    if (!is.null(metacell.celltype)) {
      text(remainX, remainY, metacell.celltype, col = "gray", cex = 0.6)
    	points(remainX, remainY,
    		col = metacell.celltype.col, pch = 16)
    }
    delta.val <- max(quantile(Xt, xymats$delta, na.rm = TRUE), 0.1)
    sel.Xt <- Xt >= delta.val
    points(remainX[sel.Xt], remainY[sel.Xt], pch = 1)
    rm(tmp); rm(title); rm(delta.val)
  }

  if (to.plot == "interaction") {
    # compute partial residuals
    remainY <- residuals(lm(Yg ~ Xt + Yj))
    remainX <- residuals(lm(Xt * Yj ~ Xt + Yj))
    tmp <- cor.test(remainX, remainY)
    est <- tmp$estimate
    pval <- tmp$p.value
    title <- paste0("Interaction",
      ", r = ", format(est, digits = 2),
    	", p = ", format(pval, scientific = TRUE, digits = 2))
    plot(remainX, remainY,
      xlab = "Xt * Yj partial resid", ylab = "Yg partial resid",
    	main = title, cex.main = 1, font.main = 1,
      col = metacell.celltype.col, pch = 16,
    )
    grid()
    if (!is.null(metacell.celltype)){
      text(remainX, remainY, metacell.celltype, col = "gray", cex = 0.6)
    	points(remainX, remainY,
    		col = metacell.celltype.col, pch = 16)
    }
    rm(tmp); rm(est); rm(pval); rm(title)
  }

	# need to check if this runs when metacell.celltype.col is NULL
	if (to.plot == "TRIPOD") {
    res <- performMatchedTest(Xt = Xt, Yj = Yj, Yg = Yg,
      match.by = match.by,
      cap.at.quantile = cap.at.quantile, delta = delta,
      min.level1 = min.level1, min.level2 = min.level2,
    	metacell.celltype.col = xymats$metacell.celltype.col,
      verbose = FALSE)
    matched.cells <- res$matched.cells
    dYg <- matched.cells$dYg
    dV <- matched.cells$dV
    meanU <- matched.cells$meanU
    sel1 <- matched.cells$sel1
    sel2 <- matched.cells$sel2
    matched.col1 <- res$matched.col1
    matched.col2 <- res$matched.col2
    coef.level1 <- res$coef.level1
    pval.level1 <- res$pval.level1
    coef.level2 <- res$coef.level2
    pval.level2 <- res$pval.level2
		if (level == 1) {
			if (is.na(coef.level1)) stop("Level 1 test cannot be performed.")
			res <- cor.test(dYg[sel2], dV[sel2], method = "pearson")
      est <- res$estimate
      pval <- res$p.value
      ylab <- "Diff in Gene Expression Yg"
      xlab <- ifelse(match.by == "Xt", "Diff in TF expression Yj",
      	"Diff in peak levels Xt")
      title <- paste0("TRIPOD matching ", match.by, " level 1:\n",
      	"r = ", format(est, digits = 2),
      	", p = ", format(pval, scientific = TRUE, digits = 2))
      plot(dV, dYg, pch = 16, col = "white",
      	xlab = xlab, ylab = ylab,
        main = title, cex.main = 1, font.main = 1
      )
      points(dV[sel2], dYg[sel2], pch = 1, cex = 2, lwd = 1.2)
      points(dV, dYg, col = matched.col1, pch = 16, cex = 2)
      points(dV, dYg, col = "white", pch = 1, cex = 1, lwd = 0.7)
      points(dV, dYg, col = matched.col2, pch = 16, cex = 1)
      grid()
      rm(res); rm(est); rm(pval); rm(title)
		}
	  if (level == 2) {
      if (is.na(coef.level2)) stop("Level 2 test cannot be performed.")
	  	res <- lm(dYg ~ dV + dV:meanU)
      remainder1 <- residuals(lm(dYg ~ dV))
      remainder2 <- residuals(lm(meanU*dV ~ dV))
      res <- cor.test(remainder1, remainder2)
      est <- res$estimate
      pval <- res$p.value
      ylab <- ifelse(match.by == "Xt", "Partial residuals dYg on dYj",
      	"Partial residuals dYg on dXt") # "Partial residuals Yg on V"
      xlab <- ifelse(match.by == "Xt", "Partial residuals Xt*dYj on dYj",
      	"Partial residuals Yj*dXt on dXt") # "Partial residuals U*V on V"
      title <- paste0("TRIPOD matching ", match.by, " level 2:\n",
      	"r = ", format(est, digits = 2),
      	", p = ", format(pval, scientific = TRUE, digits = 2))
      plot(remainder2, remainder1, col = "white", pch = 16,
      	xlab = xlab, ylab = ylab,
        main = title, cex.main = 1, font.main = 1
      )
      points(remainder2, remainder1, col = matched.col1, pch = 16, cex = 2)
      points(remainder2, remainder1, col = "white", pch = 1, cex = 1, lwd = 0.7)
      points(remainder2, remainder1, col = matched.col2, pch = 16, cex = 1)
      grid()
      rm(res); rm(est); rm(pval); rm(title)
    }
  }
}

#' Generate scatter plots for analyzing RNA prediction
#'
#' This function takes an list object returned by
#' {\code{\link{performRNAPrediction}}} and generates scatter plots
#' comparing predicted and observed RNA expression values.
#' Pearson correlation coefficient is also printed.
#'
#' @param xymats a list object returned by {\code{\link{performRNAPrediction}}}.
#' @param metacell.celltype a character vector specifying cell types of the metacells.
#' @param metacell.celltype.col a character vector representing cell types of metacells.
#' @param ... further arguments passed to {\code{\link{plot}}}
#'
#' @importFrom graphics legend
#' @export
plotRNAPrediction <- function(
	xymats, metacell.celltype = NULL, metacell.celltype.col = "black", ...
) {
	model.name <- xymats$model.name
  if (model.name == "gene.activity") {
    title <- "Gene activity:"
    ylab <- "Sum of ATAC peaks"
  } else if (model.name == "peak.LASSO") {
    title <- "Peak LASSO:"
    ylab <- "Peak LASSO LOO predicted gene expression"
  } else if (model.name == "peak.linear") {
    title <- "Peak LM:"
    ylab <- "Peak LM LOO predicted gene expression"
  } else if (model.name == "peak.TF.LASSO") {
    title <- "Peak-TF LASSO:"
    ylab <- "Peak-TF LASSO LOO predicted gene expression"
  } else if (model.name == "peak.TF.linear") {
    title <- "Peak-TF LM:"
    ylab <- "Peak-TF LM LOO predicted gene expression"
  } else {
    stop("RNA prediction plot cannot be generated.")
  }
	gene.name <-  xymats$gene.name
	title <- paste(title, gene.name, xymats$ext.upstream / 1000, "Kb up and downstream")
  plot(xymats$Yg, xymats$Y.pred,
  	xlab = "Gene expression", ylab = ylab,
    col = xymats$metacell.celltype.col, pch = 16,
  	main = title, font.main = 1, ...
  )
  legend("bottomright", paste("r =", round(xymats$cor, 3)), bty = "n")
}

# addCellTypeLegend <- function(celltype.color.map, num.col) {
# 	# make anempty plot to hold spot
#   plot(NULL, xaxt = "n", yaxt = "n", bty = "n", ylab = "", xlab = "",
#   	xlim = 0:1, ylim = 0:1)
# 	legend <- color.celltype.map$celltype
#   legend.col <- color.celltype.map$color
#   length.per.col <- ceiling(length(legend) / num.col)
#   legend("topleft",
#     legend = legend[1:half.length],
#     col = legend.col[1:half.length], pch = 16, bty = "n", cex = 1
#   )
#   legend("topright",
#     legend = legend[(half.length + 1):length(legend)],
#     col = legend.col[(half.length + 1):length(legend)], pch = 16, bty = "n", cex = 1
#   )
# }
