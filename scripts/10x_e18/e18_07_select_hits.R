# load packages
library(BiocParallel)

# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "xymats1.list.rds"
path <- file.path(dir.out, file)
xymats1.list <- readRDS(path); rm(path)

file <- "xymats4.list.rds"
path <- file.path(dir.out, file)
xymats4.list <- readRDS(path); rm(path)

file <- "xymats5X.list.rds"
path <- file.path(dir.out, file)
xymats5X.list <- readRDS(path); rm(path)

file <- "xymats5Y.list.rds"
path <- file.path(dir.out, file)
xymats5Y.list <- readRDS(path); rm(path)

# set FDR < 0.01
fdr.thresh <- 0.01

### get hits with significantly positive coefficient ###

# set the sign coefficient
sign <- "positive"; sign.str <- "pos"

# get peak-gene pairs with significantly non-zero alpha coefficients
# from model 1 (marginal Yg ~ Xt)
xymats1.alpha.df <- getPeakGenePairs(
  xymats.list = xymats1.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "1a")
file <- paste0("xymats1.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.alpha.df, path, sep = ",", quote = FALSE)

# get TF-gene pairs with significantly non-zero beta coefficients
# from model 1 (marginal Yg ~ Yj)
xymats1.beta.df <- getTFGenePairs(
  xymats.list = xymats1.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "1a")
file <- paste0("xymats1.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.beta.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 4 (interaction)
# xymats4.gamma.df <- getTrios(
#   xymats.list = xymats4.list,
#   fdr.thresh = fdr.thresh,
#   sign = sign,
#   model.name = "4a")

tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats4.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats4.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammas", pval.name = "pvalgammas")
}
xymats4.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats4.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats4.gamma.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero beta coefficients
# from model 5X level 1 (TRIPOD level 1 matching Xt)
# xymats5X.beta.df <- getTrios(
#  xymats.list = xymats5X.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 1)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5X.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5X.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "betahat.corrests", pval.name = "betahat.pvals.rhotest")
}
xymats5X.beta.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5X.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.beta.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 5X level 2 (TRIPOD level 2 matching Xt)
# xymats5X.gamma.df <- getTrios(
#  xymats.list = xymats5X.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 2)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5X.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5X.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammahat2s", pval.name = "gammahat2.pvals")
}
xymats5X.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5X.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.gamma.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero alpha coefficients
# from model 5Y level 1 (TRIPOD level 1 matching Yj)
# xymats5Y.alpha.df <- getTrios(
#  xymats.list = xymats5Y.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 1)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5Y.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5Y.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "betahat.corrests", pval.name = "betahat.pvals.rhotest")
}
xymats5Y.alpha.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5Y.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.alpha.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 5Y level 2 (TRIPOD level 2 matching Yj)
# xymats5Y.gamma.df <- getTrios(
#  xymats.list = xymats5Y.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 2)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5Y.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5Y.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammahat2s", pval.name = "gammahat2.pvals")
}
xymats5Y.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5Y.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.gamma.df, path, sep = ",", quote = FALSE)

### get hits with significantly negative coefficient ###

# set the sign coefficient
sign <- "negative"; sign.str <- "neg"

# get peak-gene pairs with significantly non-zero alpha coefficients
# from model 1 (marginal Yg ~ Xt)
xymats1.alpha.df <- getPeakGenePairs(
  xymats.list = xymats1.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "1a")
file <- paste0("xymats1.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.alpha.df, path, sep = ",", quote = FALSE)

# get TF-gene pairs with significantly non-zero beta coefficients
# from model 1 (marginal Yg ~ Yj)
xymats1.beta.df <- getTFGenePairs(
  xymats.list = xymats1.list,
  fdr.thresh = fdr.thresh,
  sign = sign,
  model.name = "1a")
file <- paste0("xymats1.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats1.beta.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 4 (interaction)
# xymats4.gamma.df <- getTrios(
#   xymats.list = xymats4.list,
#   fdr.thresh = fdr.thresh,
#   sign = sign,
#   model.name = "4a")

tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats4.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats4.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammas", pval.name = "pvalgammas")
}
xymats4.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats4.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats4.gamma.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero beta coefficients
# from model 5X level 1 (TRIPOD level 1 matching Xt)
# xymats5X.beta.df <- getTrios(
#  xymats.list = xymats5X.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 1)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5X.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5X.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "betahat.corrests", pval.name = "betahat.pvals.rhotest")
}
xymats5X.beta.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5X.beta.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.beta.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 5X level 2 (TRIPOD level 2 matching Xt)
# xymats5X.gamma.df <- getTrios(
#  xymats.list = xymats5X.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 2)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5X.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5X.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammahat2s", pval.name = "gammahat2.pvals")
}
xymats5X.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5X.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5X.gamma.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero alpha coefficients
# from model 5Y level 1 (TRIPOD level 1 matching Yj)
# xymats5Y.alpha.df <- getTrios(
#  xymats.list = xymats5Y.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 1)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5Y.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5Y.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "betahat.corrests", pval.name = "betahat.pvals.rhotest")
}
xymats5Y.alpha.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5Y.alpha.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.alpha.df, path, sep = ",", quote = FALSE)

# get trios with significantly non-zero gamma coefficients
# from model 5Y level 2 (TRIPOD level 2 matching Yj)
# xymats5Y.gamma.df <- getTrios(
#  xymats.list = xymats5Y.list,
#  fdr.thresh = fdr.thresh,
#  sign = sign,
#  model.name = "5a",
#  level = 2)
tmp.list <- list()
for (i in 1:length(xymats.list)) {
    xymats <- xymats5Y.list[[i]]
    gene.name <- xymats$gene.name
    tmp.list[[i]] <- getTriosForSingleGene(
        gene.name = gene.name, xymats.list = xymats5Y.list, fdr.thresh = fdr.thresh,
        sign = sign, coef.name = "gammahat2s", pval.name = "gammahat2.pvals")
}
xymats5Y.gamma.df <- do.call("rbind", tmp.list)
file <- paste0("xymats5Y.gamma.", sign.str, ".", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.table(xymats5Y.gamma.df, path, sep = ",", quote = FALSE)


