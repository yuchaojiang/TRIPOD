# load packages
library(BiocParallel)
library(GenomicRanges)

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

files <- list.files(dir.out, pattern = "[.]0[.]01[.]csv$", full.names = TRUE)
files <- c(files[grep("linkpeaks", files)], files[grep("xymats", files)])
names <- gsub(paste0(dir.out, "/"), "", gsub("[.]0[.]01[.]csv$", "", files))
names <- gsub("xymats", "m", names)
hit.list <- lapply(files, read.csv)
names(hit.list) <- names
hit.list <- hit.list[c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15)]
file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(hit.list, path); rm(path)

# get the total numbers of trios, peak-gene pairs, and TF-gene pairs
num.all.trio <- sapply(xymats.list, function(x) nrow(x$nonzero.peakxmotif.g))
tot.trio <- sum(num.all.trio)
num.all.peak.gene <- sapply(xymats.list, function(x) length(x$peak.gr.g))
tot.peak.gene <- sum(num.all.peak.gene)
num.all.TF.gene <- sapply(xymats.list, function(x) length(x$TF.g))
tot.TF.gene <- sum(num.all.TF.gene)

num.all.list <- list(num.all.trio, num.all.peak.gene, num.all.TF.gene)
file <- "num.all.list.rds"
path <- file.path(dir.out, file)
saveRDS(num.all.list, path); rm(path)

total <- data.frame(
	num = c(tot.trio, tot.peak.gene, tot.TF.gene)
)
rownames(total) <- c("trio", "peak_gene", "TF_gene")
colnames(total) <- "number"
file <- "total.hits.0.01.csv"
path <- file.path(dir.out, file)
write.csv(total, path); rm(path)

# generate bar plots
n <- 2
cols <- ggColorHue(n)

# compare the number of trios
trio.string.list <- bplapply(
	X = hit.list[-(1:6)],
	FUN = function(x) unique(apply(x, 1, convertTrioToString,
		col.1 = 1, col.2 = 4, col.3 = 5))
)

# take the union
trio.string.list$m5X.union.pos <- with(trio.string.list,
	union(m5X.beta.pos, m5X.gamma.pos))
trio.string.list$m5X.union.neg <- with(trio.string.list,
	union(m5X.beta.neg, m5X.gamma.neg))
trio.string.list$m5Y.union.pos <- with(trio.string.list,
	union(m5Y.alpha.pos, m5Y.gamma.pos))
trio.string.list$m5Y.union.neg <- with(trio.string.list,
	union(m5Y.alpha.neg, m5Y.gamma.neg))
trio.string.list$m5.union.pos <- with(trio.string.list,
	union(m5X.union.pos, m5Y.union.pos))
trio.string.list$m5.union.neg <- with(trio.string.list,
	union(m5X.union.neg, m5Y.union.neg))

file <- "trio.string.list.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(trio.string.list, path); rm(path)

num.trio <- sapply(trio.string.list, length)

# compare the number of unique peak-gene pairs
peak.string.list1 <- bplapply(
	X = hit.list[1:2],
	FUN = function(x) unique(apply(x, 1, convertPairToString, col.1 = 1, col.2 = 2))
)
peak.string.list2 <- bplapply(
	X = hit.list[3:4],
	FUN = function(x) unique(apply(x, 1, convertPairToString, col.1 = 1, col.2 = 3))
)
peak.string.list3 <- bplapply(
	X = hit.list[-(1:6)],
	FUN = function(x) unique(apply(x, 1, convertPairToString, col.1 = 1, col.2 = 4))
)
peak.string.list <- c(peak.string.list1, peak.string.list2, peak.string.list3)

# take the union
peak.string.list$m5X.union.pos <- with(peak.string.list,
	union(m5X.beta.pos, m5X.gamma.pos))
peak.string.list$m5X.union.neg <- with(peak.string.list,
	union(m5X.beta.neg, m5X.gamma.neg))
peak.string.list$m5Y.union.pos <- with(peak.string.list,
	union(m5Y.alpha.pos, m5Y.gamma.pos))
peak.string.list$m5Y.union.neg <- with(peak.string.list,
	union(m5Y.alpha.neg, m5Y.gamma.neg))
peak.string.list$m5.union.pos <- with(peak.string.list,
	union(m5X.union.pos, m5Y.union.pos))
peak.string.list$m5.union.neg <- with(peak.string.list,
	union(m5X.union.neg, m5Y.union.neg))

file <- "peak.string.list.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(peak.string.list, path); rm(path)

num.peak <- sapply(peak.string.list, length)

# compare the number of unique TF-gene pairs
TF.string.list1 <- bplapply(
	X = hit.list[5:6],
	FUN = function(x) unique(apply(x, 1, convertPairToString, col.1 = 1, col.2 = 3))
)
TF.string.list2 <- bplapply(
	X = hit.list[-(1:6)],
	FUN = function(x) unique(apply(x, 1, convertPairToString, col.1 = 1, col.2 = 5))
)
TF.string.list <- c(TF.string.list1, TF.string.list2)

# take the union
TF.string.list$m5X.union.pos <- with(TF.string.list,
	union(m5X.beta.pos, m5X.gamma.pos))
TF.string.list$m5X.union.neg <- with(TF.string.list,
	union(m5X.beta.neg, m5X.gamma.neg))
TF.string.list$m5Y.union.pos <- with(TF.string.list,
	union(m5Y.alpha.pos, m5Y.gamma.pos))
TF.string.list$m5Y.union.neg <- with(TF.string.list,
	union(m5Y.alpha.neg, m5Y.gamma.neg))
TF.string.list$m5.union.pos <- with(TF.string.list,
	union(m5X.union.pos, m5Y.union.pos))
TF.string.list$m5.union.neg <- with(TF.string.list,
	union(m5X.union.neg, m5Y.union.neg))

file <- "TF.string.list.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(TF.string.list, path); rm(path)

num.TF <- sapply(TF.string.list, length)

# generate a figure
combined <- merge(num.trio, num.peak, by = 0, all = TRUE)
combined <- merge(combined, num.TF, by.x = 1, by.y = 0, all = TRUE)
rownames(combined) <- combined[, 1]
combined <- combined[, -1]
models <- c(names(hit.list[1:6]), names(trio.string.list))
combined <- combined[models, ]
colnames(combined) <- c("trio", "peak", "TF")
file <- "num.hits.0.01.csv"
path <- file.path(dir.out, file)
write.csv(combined, path); rm(path)

file <- "num.hits.0.01.csv"
path <- file.path(dir.out, file)
combined <- read.csv(path, row.names = 1); rm(path)
# remove m5X.union.pos, m5X.union.neg, m5Y.union.pos, and m5Y.union.neg
combined <- combined[-c(17:20), ]

matrix.list <- list()
matrix.list$trio <- cbind(
	rev(combined[seq(2, nrow(combined), 2), 1]), # negative
	rev(combined[seq(1, nrow(combined), 2), 1]) # positive
)
matrix.list$peak <- cbind(
	rev(combined[seq(2, nrow(combined), 2), 2]), # negative
	rev(combined[seq(1, nrow(combined), 2), 2]) # positive
)
matrix.list$TF <- cbind(
	rev(combined[seq(2, nrow(combined), 2), 3]), # negative
	rev(combined[seq(1, nrow(combined), 2), 3]) # positive
)

model.names <- c(
	"LinkPeaks",
	"Marginal (Yg ~ Xt)",
	"Marginal (Yg ~ Yj)",
	"Interaction",
	"TRIPOD conditional matching Xt",
	"TRIPOD interaction matching Xt",
	"TRIPOD conditional matching Yj",
	"TRIPOD interaction matching Yj",
  "TRIPOD union"
)

titles <- c(
	paste0("Trio\n(n = ", tot.trio, ")"),
	paste0("Peak-gene pair\n(n = ", tot.peak.gene, ")"),
	paste0("TF-gene pair\n(n = ", tot.TF.gene, ")")
)

file <- "barplot_hits_fdr_0.01.pdf"
path <- file.path(dir.fig, file)
pdf(path, width = 6, height = 4); rm(path)
par(mgp = c(2, 0.8, 0))
par(mar = c(4, 0.5, 3, 0))
par(oma = c(0, 0, 0, 0.5))
par(lwd = 0.8)
mat <- matrix(1:4, nrow = 1)
widths <- c(2, 1, 1, 1)
layout(mat, widths = widths)
plot.new()
barplot(
	t(matrix.list[[1]]),
	beside = TRUE,
	horiz = TRUE,
	col = cols,
	xaxs = "r",
	yaxs = "r",
	# legend.text = c("Negative", "Positive"),
	# args.legend = list(bty = "n", xpd = TRUE),
	names.arg = rev(model.names),
  las = 2,
	lwd = 0.8
)
box()
mtext(titles[1], cex = 0.75, line = 0.3)
for (i in 2:length(matrix.list)) {
	barplot(
	  t(matrix.list[[i]]),
	  beside = TRUE,
	  horiz = TRUE,
	  col = cols,
	  xaxs = "r",
	  yaxs = "r",
		las = 2
	  # legend.text = c("Negative", "Positive"),
	  # args.legend = list(bty = "n")
	  # names.arg
  )
  box()
  mtext(titles[i], cex = 0.75, line = 0.3)
}
legend(x = 10000, y = 28,
	legend = c("Positive", "Negative"),
	fill = rev(cols),
	bty = "n",
	xpd = TRUE)
dev.off()
