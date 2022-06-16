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
	if (grepl("^4", model.name)){
		coef.name <- "gammas"
		pval.name <- "pvalgammas"
		if (!is.null(level)){
			warning("The level argument is applicable to model 5 only.")
		}
	} else if (grepl("^5", model.name)) {
		if (level == 1) {
			coef.name <- "betahat.corrests"
		  pval.name <- "betahat.pvals.rhotest"
		} else if (level == 2){
			coef.name <- "gammahat2s"
		  pval.name <- "gammahat2.pvals"
		} else {
			stop("The level argument should be either 1 or 2.")
		}
	} else {
		stop("This function is applicable to model 4 and 5 only.")
	}
	if (is.null(genes)) {genes <- sapply(xymats.list, function(x) x$gene.name)}
	
	df.list <- bplapply(
		X = genes,
		FUN = getTriosForSingleGene, # try()
		xymats.list = xymats.list, 
	  fdr.thresh = fdr.thresh, 
	  sign = sign, # "positive" or "negative"
    coef.name = coef.name,
	  pval.name = pval.name
	)
  df <- do.call("rbind", df.list)
  rownames(df) <- NULL
  return(df)
}
