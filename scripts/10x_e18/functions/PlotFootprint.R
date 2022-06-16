PlotFootprint <- function (object, features, assay = NULL, group.by = NULL, idents = NULL, 
                           label = TRUE, repel = TRUE, show.expected = TRUE, normalization = "subtract", 
                           label.top = 3, label.idents = NULL) 
{
  plot.data <- GetFootprintData(object = object, features = features, 
                                assay = assay, group.by = group.by, idents = idents)
  motif.sizes <- Signac:::GetMotifSize(object = object, features = features, 
                                       assay = assay)
  obs <- plot.data[plot.data$class == "Observed", ]
  expect <- plot.data[plot.data$class == "Expected", ]
  base <- ceiling(motif.sizes/2)
  obs$flanks <- sapply(X = seq_len(length.out = nrow(x = obs)), 
                       FUN = function(x) {
                         pos <- abs(obs[x, "position"])
                         size <- base[[obs[x, "feature"]]]
                         return((pos > size) & (pos < (size + 50)))
                       })
  if (!is.null(normalization)) {
    correction.vec <- expect$norm.value
    names(correction.vec) <- paste(expect$position, expect$feature)
    if (normalization == "subtract") {
      obs$norm.value <- obs$norm.value - correction.vec[paste(obs$position, 
                                                              obs$feature)]
    }
    else if (normalization == "divide") {
      obs$norm.value <- obs$norm.value/correction.vec[paste(obs$position, 
                                                            obs$feature)]
    }
    else {
      stop("Unknown normalization method requested")
    }
  }
  flanks <- obs[obs$flanks, ]
  flanks <- group_by(.data = flanks, feature, group)
  flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
  topmean <- top_n(x = flankmeans, n = label.top, wt = mn)
  ymax <- top_n(x = flankmeans, n = 1, wt = mn)
  ymin <- top_n(x = flankmeans, n = 1, wt = -mn)
  label.df <- data.frame()
  sub <- obs[obs$position == 75, ]
  for (i in seq_along(along.with = features)) {
    if (is.null(x = label.idents)) {
      groups.use <- topmean[topmean$feature == features[[i]], 
      ]$group
    }
    else {
      groups.use <- label.idents
    }
    df.sub <- sub[(sub$feature == features[[i]]) & (sub$group %in% 
                                                      groups.use), ]
    label.df <- rbind(label.df, df.sub)
  }
  obs$label <- NA
  label.df$label <- label.df$group
  obs <- rbind(obs, label.df)
  plotlist <- list()
  for (i in seq_along(along.with = features)) {
    df <- obs[obs$feature == features[[i]], ]
    min.use <- ifelse(test = normalization == "subtract", 
                      yes = -0.5, no = 0.5)
    axis.min <- min(min.use, ymin[ymin$feature == features[[i]], 
    ]$mn)
    axis.max <- ymax[ymax$feature == features[[i]], ]$mn + 
      0.5
    p <- ggplot(data = df, mapping = aes(x = position, y = norm.value, 
                                         color = group, label = label))
    p <- p + geom_line(size = 0.2) + xlab("Distance from motif") + 
      ylab(label = "Tn5 insertion\nenrichment") + ylim(c(axis.min, axis.max)) #+ theme(legend.position = "none")#+ guides(color = guide_legend(override.aes = list(size = 1)))
    if (label) {
      if (repel) {
        p <- p + ggrepel::geom_label_repel(box.padding = 0.5, 
                                           show.legend = FALSE)
      }
      else {
        p <- p + geom_label(show.legend = FALSE)
      }
    }
    if (show.expected) {
      df <- expect[expect$feature == features[[i]], ]
      p1 <- ggplot(data = df, mapping = aes(x = position, 
                                            y = norm.value)) + geom_line(size = 0.2) + xlab("Distance from motif") + 
        ylab(label = "Expected\nTn5 enrichment") + theme_classic()
      p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                     axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
      p <- p + p1 + plot_layout(ncol = 1, heights = c(3,  1))
    }
    plotlist[[i]] <- p
  }
  plots <- wrap_plots(plotlist)
  return(plots)
}
