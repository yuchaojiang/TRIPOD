# set directories
dir.in <- "sourcedata"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# define functions
getPvalueforCellTypes <- function(
	x,
	xymats.list,
	wnn.celltype
) {
	gene.name <- x[1]
	peak.num <- as.integer(x[2])
	TF.num <- as.integer(x[3])
	print(c(peak.num, TF.num))
	TF.name <- x[5]
	xymats <- xymats.list[[gene.name]]
	Yg <- capValues(xymats$Y)
	Xt <- capValues(xymats$Xit[, peak.num])
	Yj <- capValues(xymats$Y.TF[, TF.num])
	celltypes <- levels(factor(wnn.celltype))
	delta.coeff.pval <- matrix(nrow = length(celltypes), ncol = 5)
  rownames(delta.coeff.pval) <- celltypes
  for (wnn.rm.celltype in unique(wnn.celltype)) {
  	wnn.rm <- which(wnn.celltype == wnn.rm.celltype)
  	delta.coeff.pval[wnn.rm.celltype, ] <- test.influential(
		Yig = Yg,
		Xit = Xt,
		Yij = Yj,
		wnn.rm = wnn.rm,
		nsamp = 1000)
  }
  return(delta.coeff.pval[, 5])
}

# read in data
file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
hit.list <- readRDS(path); rm(path)

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
wnn.celltype <- readRDS(path); rm(path)

file <- "metacell.celltype.col.rds"
path <- file.path(dir.out, file)
wnn.celltype.col <- readRDS(path); rm(path)

file <- "color.map.rds"
path <- file.path(dir.out, file)
col.map <- readRDS(path); rm(path)

# file <- "e18.metacell.rds"
# path <- file.path(dir.out, file)
# e18 <- readRDS(path); rm(path)

### perform influential analysis ###

# focus on neurogenesis and gliogenesis TFs
tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")
tfs <- c(tfs.neuro, tfs.glio)

names(hit.list)

# focus on model 4 (interaction) hits
df <- hit.list$m4.gamma.pos
df <- df[df$TF %in% tfs, ]

celltype.pvals <- t(apply(
  X = df,
  MARGIN = 1,
  FUN = getPvalueforCellTypes,
  xymats.list = xymats.list,
  wnn.celltype = wnn.celltype)
)

# create a col.map object
# celltypes <- levels(e18$celltype)
# col.map <- data.frame(celltype = wnn.celltype, color = wnn.celltype.col)
# col.map <- unique(col.map)
# rownames(col.map) <- col.map$celltype
# col.map <- col.map[celltypes, ]
# rownames(col.map) <- NULL
# file <- "color.map.rds"
# path <- file.path(dir.out, file)
# saveRDS(col.map, path); rm(path)

# order the columns
tmp <- celltype.pvals[, col.map$celltype]
celltype.pvals <- tmp; rm(tmp)
df <- cbind(df, celltype.pvals)
# file <- "influential.celltype.pvals.rds"
# path <- file.path(dir.out, file)
# saveRDS(df, path); rm(path)

# celltype.pvals <- df[, -(1:6)]
celltype.adj <- t(apply(celltype.pvals, 1, p.adjust, method = "bonferroni"))
colnames(celltype.adj) <- paste0("adj_", colnames(celltype.adj))
df <- cbind(df, celltype.adj)
file <- "influential.celltype.pvals.rds"
path <- file.path(dir.out, file)
saveRDS(df, path); rm(path)

fdr.thresh <- 0.01
# tmp <- apply(celltype.adj, 1, function(x) sum(any(x < fdr.thresh)))
tmp <- apply(df[, grep("^adj_", colnames(df))], 1, function(x) sum(any(x < fdr.thresh)))
df$sig <- tmp; rm(tmp)

# create a list of data frames for individual TFs
df$TF <- factor(df$TF, levels = tfs)
df$TF <- droplevels(df$TF)
celltype.pvals.list <- split(df, df$TF)

file <- "influential.celltype.pvals.list.rds"
path <- file.path(dir.out, file)
saveRDS(celltype.pvals.list, path); rm(path)

num.sig <- sapply(
	X = celltype.pvals.list,
	FUN = function(x) apply(x[, 16:22], 2, function(x) sum(x < fdr.thresh, na.rm = TRUE))
)
# num.sig[is.na(num.sig)] <- 0
rownames(num.sig) <- gsub("adj_", "", rownames(num.sig))

file <- "influential.num.sig.csv"
path <- file.path(dir.out, file)
write.csv(num.sig, path); rm(path)

sig <- df[df$sig == 1, ]
file <- "influential.sig.csv"
path <- file.path(dir.out, file)
write.csv(sig, path); rm(path)

### generate bar plots ###

titles <- colnames(num.sig)
titles <- gsub("[.]", "-", titles)

file <- "barplot_influential.pdf"
path <- file.path(dir.fig, file)
# pdf(path, width = 12, height = 2); rm(path)
# pdf(path, width = 12, height = 2.5); rm(path)
pdf(path, width = 12, height = 3); rm(path)
par(mgp = c(2, 0.8, 0))
par(mar = c(4, 0.5, 3, 0))
par(oma = c(0, 0, 0, 0.5))
# par(mfrow = c(1, 4))
mat <- matrix(1:12, nrow = 1)
widths <- c(1.1, rep(1, ncol(num.sig)))
layout(mat, widths = widths)
plot.new()
barplot(
	rev(num.sig[, 1]),
	beside = TRUE,
	horiz = TRUE,
	col = rev(col.map$color),
	xaxs = "r",
	yaxs = "r",
	# legend.text = c("Negative", "Positive"),
	# args.legend = list(bty = "n", xpd = TRUE),
	names.arg = rev(rownames(num.sig)),
  las = 2 # ,
	# xaxt = "n"
)
box()
mtext(titles[1], cex = 0.75, line = 0.3)
for (i in 2:ncol(num.sig)) {
	if (i != 8) {
	  barplot(
	    rev(num.sig[, i]),
	    beside = TRUE,
	    horiz = TRUE,
	    col = rev(col.map$color),
	    xaxs = "r",
	    yaxs = "r",
		  las = 2,
	    # legend.text = c("Negative", "Positive"),
	    # args.legend = list(bty = "n")
	    names.arg = "",
    )
	} else if (i == 8) {
		barplot(
	    rev(num.sig[, i]),
	    beside = TRUE,
	    horiz = TRUE,
	    col = rev(col.map$color),
	    xaxs = "r",
	    yaxs = "r",
	    # legend.text = c("Negative", "Positive"),
	    # args.legend = list(bty = "n", xpd = TRUE),
	    names.arg = "",
      las = 2,
	    xaxt = "n"
    )
		axis(1, at = 0:max(num.sig[, i]), las = 2)
	}
  box()
  mtext(titles[i], cex = 0.75, line = 0.3)
}
# legend(x = 10000, y = 34,
# 	legend = c("Positive", "Negative"),
# 	fill = rev(cols),
# 	bty = "n",
# 	xpd = TRUE)
dev.off()



