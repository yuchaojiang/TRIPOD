# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# load packages
library(dplyr)
library(ggplot2)
library(Seurat)
library(Signac)
library(olsrr)
library(dendextend)
library(gridBase)
library(grid)
library(FNN)

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
hit.list <- readRDS(path); rm(path)

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "metacell.rna.rds"
path <- file.path(dir.out, file)
wnn.rna <- readRDS(path); rm(path)

file <- "metacell.peak.rds"
path <- file.path(dir.out, file)
wnn.peak <- readRDS(path); rm(path)

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
wnn.celltype <- readRDS(path); rm(path)

file <- "metacell.celltype.col.rds"
path <- file.path(dir.out, file)
wnn.celltype.col <- readRDS(path); rm(path)

file <- "color.map.rds"
path <- file.path(dir.out, file)
col.map <- readRDS(path); rm(path)

file <- "e18.footprint.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

file <- "influential.sig.csv"
path <- file.path(dir.out, file)
sig <- read.csv(path); rm(path)

### select trios for plotting ###

# Pax6 -> Neurog2 via chr3-127639075-127643093 in Glioblast/Neuroblast
# Pax6 -> Kif23 via chr9-61945892-61948778 in Glioblast
# Neurog2 -> Mt1 via chr8-94097159-94099509 in Glioblast/Neuroblast
# Neurog2 -> Myo1h via chr5-114263473-114265217 in Neuroblast
# Sox9 -> Btg2 via chr1-134273297-134275102 in Neuroblast
# Sox9 -> Gli3 via chr13-15540753-15546568 in Glioblast
# Ascl1 -> Olig1 via chr16-91265825-91266538 in Glioblast
# Ascl1 -> Cbfa2t3 via chr8-122737655-122739809 in GABAergic

genes <- c("Neurog2", "Kif23", "Mt1", "Myo1h", "Btg2", "Gli3", "Olig1", "Cbfa2t3")
tfs <- rep(c("Pax6", "Neurog2", "Sox9", "Ascl1"), each = 2)
peaks <- c(
	"chr3-127639075-127643093",
	"chr9-61945892-61948778",
	"chr8-94097159-94099509",
	"chr5-114263473-114265217",
	"chr1-134273297-134275102",
	"chr13-15540753-15546568",
	"chr16-91265825-91266538",
	"chr8-122737655-122739809"
)

sig <- sig[!is.na(sig$X), ]
rownames(sig) <- sig$X
sig <- sig[, -1]
sig <-  sig[, 1:5]
sig.list <- split(sig, sig$TF)

### generate objects required for feature plots ###

# get meta-cell-specific reduction, and find each meta-cell's nearest neighbors
reduction <- "umap.rna"
cell.reduction <- e18@reductions[[reduction]]@cell.embeddings

metacell.reduction <- matrix(
	nrow = length(levels(e18$seurat_clusters)),
	ncol = ncol(cell.reduction)
)
rownames(metacell.reduction) <- paste0("metacell_", levels(e18$seurat_clusters))
colnames(metacell.reduction) <- colnames(cell.reduction)

for (i in 1:nrow(metacell.reduction)) {
  metacell.reduction[i, ] <- apply(
  	cell.reduction[e18$seurat_clusters == (i - 1), ], 2, mean)
}

file <- "cell.reduction.rds"
path <- file.path(dir.out, file)
saveRDS(cell.reduction, path); rm(path)

file <- "metacell.reduction.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.reduction, path); rm(path)

# set the number of metacell neighbors to include
k <- 4
metacell.neighbors <- get.knn(metacell.reduction, k = k)$nn.index
metacell.neighbors <- cbind(1:nrow(metacell.neighbors), metacell.neighbors)

file <- "metacell.neighbors.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.neighbors, path); rm(path)

# get a dendrogram based on RNA expression of all genes
metacell.celltype <- wnn.celltype
metacell.rna <- wnn.rna
tmp <- aggregate(metacell.rna ~ metacell.celltype, FUN = sum)
rownames(tmp) <- tmp[, 1]
tmp <- tmp[, -1]
tmp <- tmp[col.map$celltype, ]
tmp <- sweep(tmp, 1, rowSums(tmp)/10^6, "/")
dend <- tmp %>% scale %>% dist %>% hclust %>% as.dendrogram
rm(tmp)

file <- "dend.rds"
path <- file.path(dir.out, file)
saveRDS(dend, path); rm(path)

### determine axis ranges for scatter plots ###

range.matrix <- matrix(NA, nrow = length(genes), ncol = 4)
for (i in 1:length(genes)) {
	gene.name <- genes[i]
	TF.name <- tfs[i]
	peak.name <- peaks[i]
  X <- wnn.peak[, colnames(wnn.peak) == peak.name]
  Y.TF <- wnn.rna[, colnames(wnn.rna) == TF.name]
  Y <- wnn.rna[, colnames(wnn.rna) == gene.name]
  # cap values for models 1, 2, and 4
  X.capped <- capValues(X)
  Y.TF.capped <- capValues(Y.TF)
  Y.capped <- capValues(Y)
  remainX <- residuals(lm(X.capped * Y.TF.capped ~ X.capped + Y.TF.capped)) # regress out Xt and Yj
  remainY <- residuals(lm(Y.capped ~ X.capped + Y.TF.capped)) # regress out Xt and Yj
  range.matrix[i, ] <- c(range(remainX), range(remainY))
}
range.matrix
#            [,1]       [,2]        [,3]       [,4]
# [1,] -5786.8063  6367.6030 -115.194664 182.115696
# [2,] -7725.0819 10943.2342 -144.804552 333.612721
# [3,] -1798.2325  2372.5533  -30.352507  88.985140
# [4,] -1632.8376  1693.7822   -6.688546   9.253752
# [5,]  -691.3929  1417.6811  -88.301286 153.863458
# [6,] -4863.8450  2657.3170 -605.294592 768.423986
# [7,]  -275.7213   286.9157  -42.338282 146.685975
# [8,]  -298.9485   735.1378 -111.339770 207.912140

# manually set axis ranges
at.y.list <- list(
	c(-100, 0, 100),
	c(0, 150, 300),
  c(0, 50),
	c(-5, 0, 5),
	c(0, 100),
	c(-500, 0, 500),
	c(0, 100),
	c(-100, 0, 100, 200)
)

at.x.list <- list(
	c(-5000, 0, 5000),
	c(-5000, 0, 5000, 10000),
  c(-1000, 0, 1000, 2000),
	c(-1000, 0, 1000),
	c(-1000, 0, 1000),
	c(-4000, -2000, 0, 2000),
	c(-200, 0, 200),
	c(0, 500)
)

### generate plots ###

# set window size
ext.upstream <- 2e5
object <- e18; rm(e18)

# set plotting parameters
cex.axis <- 1.5
cex.lab <- 1.8
cex <- 2
scale.pt <- 2
cex.title <- 1.5
cex.title2 <- 0.5
line.title <- 0.5
line.title2 <- 2.5
line.xlab <- 3
line.ylab <- 3
lwd <- 1.5
mgp.x <- c(2, 0.8, 0) # mgp.x <- c(0.8, 0, 0)
mgp.y <- c(2, 0.8, 0) # mgp.y <- c(0.8, 0.2, 0)
xlab <- expression(italic(X[t])~italic(Y[j])~"partial residual")
ylab <- expression(italic(Y[g])~"partial residual")
reduction <- "umap.rna"
row.feature <- rep(1:4, each = 2)
col.feature <- rep(c(3, 6), 4)
row.pval <- rep(1:4, each = 2)
col.pval <- rep(c(2, 5), 4)
# object$seurat_clusters <- object$SCT_snn_res.15
pos.x <- c(0.3, 0.2, 0.15, 0.6, 0.4, 0.25, 0.2)
pos.y <- c(0.6, 0.75, 0.5, 0.65, 0.7, 0.85, 0.15)

set.seed(1)
file <- "scatter_feature_inflential.pdf"
path <- file.path(dir.fig, file)
pdf(file = path, width = 18, height = 12.5)
par(mfrow = c(4, 6))
par(mar = c(4, 4, 4, 0), pty = "s")
par(mgp = c(2, 0.8, 0))
for (i in 1:length(genes)) {
	gene.name <- genes[i]
	TF.name <- tfs[i]
	print(c(TF.name, gene.name))
	peak.name <- peaks[i]
  X <- wnn.peak[, colnames(wnn.peak) == peak.name]
  Y.TF <- wnn.rna[, colnames(wnn.rna) == TF.name]
  Y <- wnn.rna[, colnames(wnn.rna) == gene.name]
  X.capped <- capValues(X)
  Y.TF.capped <- capValues(Y.TF)
  Y.capped <- capValues(Y)
  remainX <- residuals(lm(X.capped * Y.TF.capped ~ X.capped + Y.TF.capped))
  remainY <- residuals(lm(Y.capped ~ X.capped + Y.TF.capped))
  # range.matrix[i, ] <- c(range(remainX), range(remainY))
  partial.cor <- cor.test(remainX, remainY)$estimate
  partial.cor.pval <- cor.test(remainX, remainY)$p.value
  plot(
    remainX,
    remainY,
	  xlab = "",
	  ylab = "",
    cex = cex,
	  cex.lab = cex.lab,
	  cex.axis = cex.axis,
	  main = "",
	  col = wnn.celltype.col,
	  pch = 16,
	  xaxt = "n",
	  yaxt = "n",
	  bty = "n")
  grid()
  box(lwd = lwd)
  par(mgp = mgp.x)
  axis(
	side = 1,
	at = at.x.list[[i]],
	las = 1,
	lwd = lwd,
	lwd.ticks = lwd,
	tck = -0.02,
	cex.axis = cex.axis)
  title(xlab = xlab, line = line.xlab, cex.lab = cex.lab)
  par(mgp = mgp.y)
  axis(
	  side = 2,
	  at = at.y.list[[i]],
	  las = 2,
	  lwd = lwd,
	  lwd.ticks = lwd,
	  tck = -0.01,
	  cex.axis = cex.axis)
  title(ylab = ylab, line = line.ylab, cex.lab = cex.lab)
  title <- paste0(
	  "r = ",
    format(partial.cor, digits = 2),
    " p = ", format(partial.cor.pval, scientific = TRUE, digits = 2))
  mtext(title, cex = cex.title, line = line.title)
  title.str <- paste0(TF.name, " -> ", gene.name, " via ", peak.name)
  mtext(title.str, cex = cex.title2, line = line.title2, adj = 0)
  rm(partial.cor)
  rm(partial.cor.pval)

  # cell-type specific p-values
  delta.coeff.pval <- matrix(nrow=length(wnn.celltype), ncol=5)
  Xit <- X.capped
  Yij <- Y.TF.capped
  Yig <- Y.capped
  rownames(delta.coeff.pval) <- names(wnn.celltype)
  delta.coeff.pval <- matrix(nrow = length(unique(wnn.celltype)), ncol = 5)
  colnames(delta.coeff.pval) <- c("Intercept", "Xit", "Yij", "Xit:Yij", "Yig")
  rownames(delta.coeff.pval) <- unique(wnn.celltype)
  for(wnn.rm.celltype in unique(wnn.celltype)){
    wnn.rm <- which(wnn.celltype == wnn.rm.celltype)
    delta.coeff.pval[wnn.rm.celltype, ] <- test.influential(Yig, Xit, Yij,
    	wnn.rm, plot.histogram = FALSE, nsamp = 1000)
  }

  myplots <- list()
  for (j in 1:5){
    d = data.frame(obs = 1:length(unique(wnn.celltype)),
                     cd = -log(delta.coeff.pval[, j], 10),
                     txt = unique(wnn.celltype))
    d$txt[d$cd<(-log(0.05,10))]=NA
    d=d[levels(object$celltype),]
    d$obs=1:length(d$obs)
    myplots[[j]] <- ggplot(d, aes(x = obs, y = cd, label = txt, ymin = 0, ymax = cd)) +
    geom_linerange(colour = unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))]) +
    geom_point(shape = 16, colour = unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))]) +
    geom_hline(yintercept = -log(0.01, 10), colour = "black", linetype='dotted')+
    xlab("Cell types") + ylab("-log(p-value)") +
    ggtitle("") +  #	ggtitle(paste0("Cell-type-specific p-values for ", colnames(delta.coeff.pval)[j])) +
    geom_text(size = 4, family = "Helvetica",
    	colour = unique(wnn.celltype.col)[match(rownames(d), unique(wnn.celltype))], na.rm = TRUE) +
    theme(
    	axis.text = element_text(size = 14),
      axis.title = element_text(size = 14, face = "plain"))
  }
  p3 <- myplots[[5]]

  plot.new()
  if (i == 1) {
    vps <- baseViewports()
    pushViewport(viewport(layout = grid.layout(4, 6)))
    vp1 <- plotViewport(c(0, 0, 0, 0)) # vp1 <- plotViewport(c(1, 0.5, 0, 1))
  }
  pushViewport(viewport(layout.pos.row = row.pval[i], layout.pos.col = col.pval[i]))
  print(p3, vp = vp1)
  popViewport()

  # meta-cell sampling
  delta.coeff.pval <- matrix(nrow=length(wnn.celltype), ncol=5)
  Xit <- X.capped
  Yij <- Y.TF.capped
  Yig<- Y.capped
  # colnames(delta.coeff.pval) <- c("(Intercept)", "Xit", "Yij", "Xit:Yij",  "Yig")
  rownames(delta.coeff.pval) <- names(wnn.celltype)
  for (j in 1:length(wnn.celltype)) {
    wnn.rm <- metacell.neighbors[j, ]
    delta.coeff.pval[j, ] <- test.influential(Yig, Xit, Yij, wnn.rm,
    	plot.histogram = FALSE, nsamp = 1000)
  }

  object$wnn.samp.logpval <- -log(delta.coeff.pval[object$seurat_clusters, 5], 10)
  p.feature <- FeaturePlot(
      object = object,
      features = "wnn.samp.logpval",
      reduction = reduction,
      label = FALSE,
  	  order = TRUE,
      label.size = 2) +
      ggtitle("") + # ggtitle('Meta-cell-specific p-values for Yig from sampling') +
      theme(plot.title = element_text(size = 10, face = "plain"))
  plot.new()
  pushViewport(viewport(layout.pos.row = row.feature[i], layout.pos.col = col.feature[i]))
  print(p.feature, vp = vp1)
  popViewport()
  for (j in 1:nrow(col.map)) {
    text(pos.x[j], pos.y[j], col.map$celltype[j], cex = 1)
  	if (col.map$celltype[j] == "Neuroblast") {
  		text(0.35, 0.4, col.map$celltype[j], cex = 1)
  	}
  }
}
dev.off()
