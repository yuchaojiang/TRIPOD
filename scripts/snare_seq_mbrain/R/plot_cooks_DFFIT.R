plot_cooks_DFFIT=function(Yig, Xit, Yij,plot.all=NULL,...){
  if(is.null(plot.all)){plot.all=FALSE}
  model4=lm(Yig~1+Xit+Yij+Xit*Yij)
  #ols_plot_cooksd_bar(model4)
  #ols_plot_cooksd_chart(model4)
  
  k <- ols_prep_cdplot_data(model4)
  d <- ols_prep_outlier_obs(k)
  d$txt[!is.na(d$txt)]=wnn.celltype[!is.na(d$txt)]
  f <- ols_prep_cdplot_outliers(k)
  p <- ggplot(d, aes(x = obs, y = cd, label = txt, ymin = min(cd), ymax = cd)) + 
    geom_linerange(colour = wnn.celltype.col) + 
    geom_point(shape = 16, colour = wnn.celltype.col) + geom_hline(yintercept = k$ts, colour = "red") + 
    xlab("WNN") + ylab("Cook's D") + ggtitle("Cook's D Chart for Yig") + 
    geom_text(size = 3, family = "serif", fontface = "italic", colour = wnn.celltype.col, na.rm = TRUE) + 
    annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 1, 
             family = "serif", fontface = "italic", colour = "darkred", label = paste("Threshold:", round(k$ts, 3)))
  p1=p
  
  #ols_plot_dffits(model4)
  dbetas <- NULL
  obs <- NULL
  txt <- NULL
  dffitsm <- unlist(dffits(model4))
  k <- olsrr:::model_n_coeffs(model4)
  n <- olsrr:::model_rows(model4)
  dffits_t <- sqrt(k/n) * 2
  title <- names(model.frame(model4))[1]
  dfits_data <- data.frame(obs = seq_len(n), dbetas = dffitsm)
  d <- ols_prep_dfbeta_data(dfits_data, dffits_t)
  d$txt[!is.na(d$txt)]=wnn.celltype[!is.na(d$txt)]
  f <- ols_prep_dfbeta_outliers(d)
  p <- ggplot(d, aes(x = obs, y = dbetas, label = txt, ymin = 0, 
                     ymax = dffitsm)) + geom_linerange(colour = wnn.celltype.col) + 
    geom_hline(yintercept = c(0, dffits_t, -dffits_t), colour = "red") + 
    geom_point(colour = wnn.celltype.col, shape = 16) + xlab("WNN") + 
    ylab("DFFITS") + ggtitle(paste("DFFIT for",  title)) +
    geom_text(hjust = -0.2, nudge_x = 0.15, size = 3, family = "serif", 
              fontface = "italic", colour = wnn.celltype.col, na.rm = TRUE) + 
    annotate("text", x = Inf, y = Inf, hjust = 1.5, vjust = 2, 
             family = "serif", fontface = "italic", colour = "darkred", 
             label = paste("Threshold:", round(dffits_t, 2)))
  p2=p
  
  # ols_plot_dfbetas(model4)
  obs <- NULL
  txt <- NULL
  dfb <- dfbetas(model4)
  n <- nrow(dfb)
  np <- ncol(dfb)
  threshold <- 2/sqrt(n)
  myplots <- list()
  outliers <- list()
  for (i in seq_len(np)) {
    dbetas <- dfb[, i]
    df_data <- data.frame(obs = seq_len(n), dbetas = dbetas)
    d <- ols_prep_dfbeta_data(df_data, threshold)
    d$txt[!is.na(d$txt)]=wnn.celltype[!is.na(d$txt)]
    f <- ols_prep_dfbeta_outliers(d)
    p <- eval(substitute(ggplot(d, aes(x = obs, y = dbetas, label = txt, ymin = 0, ymax = dbetas)) + 
                           geom_linerange(colour = wnn.celltype.col) + 
                           geom_hline(yintercept = c(0, threshold, -threshold), colour = "red") + 
                           geom_point(colour = wnn.celltype.col, shape = 16) + xlab("WNN") + ylab("DFBETAS") + 
                           ggtitle(paste("DFBETA for", colnames(dfb)[i])) + 
                           geom_text(hjust = -0.2, nudge_x = 0.15, size = 2, 
                                     family = "serif", fontface = "italic", colour = wnn.celltype.col, 
                                     na.rm = TRUE) + 
                           annotate("text", x = Inf, y = Inf, hjust = 1.5, vjust = 2, 
                                    family = "serif", fontface = "italic", colour = "darkred", 
                                    label = paste("Threshold:", round(threshold, 2))), list(i = i)))
    myplots[[i]] <- p
    outliers[[i]] <- f
  }
  if(plot.all){
    plot = wrap_plots(myplots)+p2+p1+plot_layout(ncol=3)
  } else{
    plot = p1+p2+plot_layout(ncol=2)
  }
  return(plot)
}
