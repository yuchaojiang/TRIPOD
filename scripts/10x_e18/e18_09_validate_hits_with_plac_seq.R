# load packages
library(Signac)
library(GenomicRanges)
library(ggplot2)
library(ggvenn)
library(cowplot)
library(patchwork)
library(BiocParallel)
library(BiocParallel)
library(GenomicRanges)
library(pheatmap)
library(RColorBrewer)

# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# download data into the dir.in directory
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3819641
# GSM3819641_RH_025_final.10k.2.peaks.bedpe.gz
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3819642
# GSM3819642_RH_026_final.10k.2.peaks.bedpe.gz

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "transcripts.gr.rds"
path <- file.path(dir.out, file)
transcripts.gr <- readRDS(path); rm(path)

file <- "peaks.gr.rds"
path <- file.path(dir.out, file)
peaks.gr <- readRDS(path); rm(path)

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "genes.rds"
path <- file.path(dir.out, file)
genes <- readRDS(path); rm(path)

file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
hit.list <- readRDS(path); rm(path)

file <- "GSM3819641_RH_025_final.10k.2.peaks.bedpe"
path <- file.path(dir.in, file)
plac.df.1 <- read.table(path, head = T); rm(path)

file <- "GSM3819642_RH_026_final.10k.2.peaks.bedpe"
path <- file.path(dir.in, file)
plac.df.2 <- read.table(path, head = T); rm(path)

### perform analysis for combined PLAC-seq replicates ###

# combine replicates
plac.df <- rbind(plac.df.1, plac.df.2)
plac.df <- plac.df[order(plac.df$start1), ]
plac.df <- plac.df[order(plac.df$chr1), ]

plac.1.gr <- GRanges(seqnames = plac.df$chr1,
	ranges = IRanges(start = plac.df$start1 + 1, end = plac.df$end1),
	index = 1:nrow(plac.df))
plac.2.gr <- GRanges(seqnames = plac.df$chr2,
	ranges = IRanges(start = plac.df$start2 + 1, end = plac.df$end2),
	index = 1:nrow(plac.df))

plac.gr.list <- list(plac.1.gr, plac.2.gr)
names(plac.gr.list) <- c("plac.1", "plac.2")
file <- "plac.merged.gr.list.rds"
path <- file.path(dir.out, file)
saveRDS(plac.gr.list, path); rm(path)

# set resolution
resolution <- 1e4

# get TSS
TSS <- promoters(transcripts.gr, upstream = 0, downstream = 1)

# get a list of pairs for all target genes whose TSS overlaps one of PLAC-seq regions
all.pairs <- getAllPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "all.pairs.rds"
path <- file.path(dir.out, file)
saveRDS(all.pairs, path); rm(path)

# obtain gene-peak pairs supported by PLAC-seq
plac.pairs <- getInteractingPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "plac.pairs.rds"
path <- file.path(dir.out, file)
saveRDS(plac.pairs, path); rm(path)

# select hits with positive coefficients
hit.list <- hit.list[c(1, 3, 7, 9, 11, 13, 15)]

sig.string.list <- lapply(
	X = hit.list,
	FUN = getSignificantPairs,
	genes,
	TSS,
	resolution,
	gr.list = plac.gr.list
)

names(sig.string.list) <- gsub("[.]pos", "", names(sig.string.list))
sig.string.list$m5X.union <- with(sig.string.list, union(m5X.beta, m5X.gamma))
sig.string.list$m5Y.union <- with(sig.string.list, union(m5Y.alpha, m5Y.gamma))
sig.string.list$m5.union <- with(sig.string.list, union(m5X.union, m5Y.union))
sig.string.list$union <- with(sig.string.list, union(m1.alpha, m5Y.union))

file <- "sig.string.list.pos.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(sig.string.list, path); rm(path)

res.matrix <- t(sapply(sig.string.list, hypergeometricTestForStrings,
	all = all.pairs, ground.truth = plac.pairs))
rownames(res.matrix) <- names(sig.string.list)
colnames(res.matrix) <- c("q", "m", "n", "k", "pval")
res <- as.data.frame(res.matrix)

# perform Bonferroni correction
res$adj <- p.adjust(res$pval, method = "bonferroni")

# take the negative log of p-values
res$log <- -log10(res$adj)

file <- "res.plac.seq.pos.0.01.csv"
path <- file.path(dir.out, file)
write.csv(res, path); rm(path)

# rename objects
res.merged <- res; rm(res)
sig.string.list.merged <- sig.string.list; rm(sig.string.list)

# draw three circle Venn diagram
int.string.list <- lapply(sig.string.list.merged, intersect, all.pairs)
set1 <- int.string.list$m1.alpha
set2 <- int.string.list$m5Y.union
set3 <- plac.pairs
set.list <- list(set1, set2, set3)
union <- Reduce("union", set.list)
factor.list <- lapply(set.list, factor, levels = union)
bool.list <- lapply(factor.list, function(x) as.logical(as.vector(table(x))))
t <- tibble(
  value = union,
	set1 = bool.list[[1]],
	set2 = bool.list[[2]],
  set3 = bool.list[[3]])
p2 <- ggplot(t, aes(A = set1, B = set2, C = set3)) +
  geom_venn(
    show_percentage = FALSE,
  	set_names = c(
		  "Marginal      ",
		  "     TRIPOD",
		  "PLAC-seq"),
    set_name_color = "black",
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    text_size = 1.8, stroke_size = 0.2, set_name_size = 2.5) +
	theme_void() +
  theme(
	  # plot.title = element_text(size = 6, hjust = 0.5),
	  # margin = margin(t = 10, b = -20)),
    plot.margin = margin(t = 0, r = 1, b = 0, l = 0, unit = "cm"),
    panel.spacing = unit(0, "lines")) +
  coord_fixed()

file <- "venn_plac_seq.pdf"
path <- file.path(dir.fig, file)
ggsave(filename = path, plot = p2, width = 2, height = 2); rm(path)

### perform analysis for PLAC-seq replicate 1 ###

# generate a Granges object for replicate 1
plac.df <- plac.df.1
plac.df <- plac.df[order(plac.df$start1), ]
plac.df <- plac.df[order(plac.df$chr1), ]
plac.1.gr <- GRanges(seqnames = plac.df$chr1,
	ranges = IRanges(start = plac.df$start1 + 1, end = plac.df$end1),
	index = 1:nrow(plac.df))
plac.2.gr <- GRanges(seqnames = plac.df$chr2,
	ranges = IRanges(start = plac.df$start2 + 1, end = plac.df$end2),
	index = 1:nrow(plac.df))
plac.gr.list <- list(plac.1.gr, plac.2.gr)
names(plac.gr.list) <- c("plac.1", "plac.2")
file <- "plac.rep1.gr.list.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(plac.gr.list, path)}; rm(path)

# analyze replicate 1
all.pairs <- getAllPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "all.pairs.rep1.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(all.pairs, path)}; rm(path)

plac.pairs <- getInteractingPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "plac.pairs.rep1.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(plac.pairs, path)}; rm(path)

sig.string.list <- lapply(
	X = hit.list,
	FUN = getSignificantPairs,
	genes,
	TSS,
	resolution,
	gr.list = plac.gr.list
)

names(sig.string.list) <- gsub("[.]pos", "", names(sig.string.list))
sig.string.list$m5X.union <- with(sig.string.list, union(m5X.beta, m5X.gamma))
sig.string.list$m5Y.union <- with(sig.string.list, union(m5Y.alpha, m5Y.gamma))
sig.string.list$m5.union <- with(sig.string.list, union(m5X.union, m5Y.union))
sig.string.list$union <- with(sig.string.list, union(m1.alpha, m5Y.union))

file <- "sig.string.list.pos.0.01.rep1.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(sig.string.list, path)}; rm(path)

res.matrix <- t(sapply(sig.string.list, hypergeometricTestForStrings,
	all = all.pairs, ground.truth = plac.pairs))
rownames(res.matrix) <- names(sig.string.list)
colnames(res.matrix) <- c("q", "m", "n", "k", "pval")
# create data frames
res <- as.data.frame(res.matrix)
# perform Bonferroni correction
res$adj <- p.adjust(res$pval, method = "bonferroni")
# take the negative log of p-values
res$log <- -log10(res$adj)
# save the result
file <- "res.plac.seq.pos.0.01.rep1.csv"
path <- file.path(dir.out, file)
write.csv(res, path); rm(path)
# rename objects
res1 <- res; rm(res)
sig.string.list1 <- sig.string.list; rm(sig.string.list)

### perform analysis for PLAC-seq replicate 2 ###

# generate a Granges object for replicate 2
plac.df <- plac.df.2
plac.df <- plac.df[order(plac.df$start1), ]
plac.df <- plac.df[order(plac.df$chr1), ]
plac.1.gr <- GRanges(seqnames = plac.df$chr1,
	ranges = IRanges(start = plac.df$start1 + 1, end = plac.df$end1),
	index = 1:nrow(plac.df))
plac.2.gr <- GRanges(seqnames = plac.df$chr2,
	ranges = IRanges(start = plac.df$start2 + 1, end = plac.df$end2),
	index = 1:nrow(plac.df))
plac.gr.list <- list(plac.1.gr, plac.2.gr)
names(plac.gr.list) <- c("plac.1", "plac.2")
file <- "plac.rep2.gr.list.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(plac.gr.list, path)}; rm(path)

# analyze replicate 2
all.pairs <- getAllPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "all.pairs.rep2.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(all.pairs, path)}; rm(path)

plac.pairs <- getInteractingPairs(
  genes = genes,
	xymats.list = xymats.list,
	TSS = TSS,
	resolution = resolution,
	gr.list = plac.gr.list)
file <- "plac.pairs.rep2.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(plac.pairs, path)}; rm(path)

sig.string.list <- lapply(
	X = hit.list,
	FUN = getSignificantPairs,
	genes,
	TSS,
	resolution,
	gr.list = plac.gr.list
)

names(sig.string.list) <- gsub("[.]pos", "", names(sig.string.list))
sig.string.list$m5X.union <- with(sig.string.list, union(m5X.beta, m5X.gamma))
sig.string.list$m5Y.union <- with(sig.string.list, union(m5Y.alpha, m5Y.gamma))
sig.string.list$m5.union <- with(sig.string.list, union(m5X.union, m5Y.union))
sig.string.list$union <- with(sig.string.list, union(m1.alpha, m5Y.union))

file <- "sig.string.list.pos.0.01.rep2.rds"
path <- file.path(dir.out, file)
if (!file.exists(path)) {saveRDS(sig.string.list, path)}; rm(path)

res.matrix <- t(sapply(sig.string.list, hypergeometricTestForStrings,
	all = all.pairs, ground.truth = plac.pairs))
rownames(res.matrix) <- names(sig.string.list)
colnames(res.matrix) <- c("q", "m", "n", "k", "pval")
# create data frames
res <- as.data.frame(res.matrix)
# perform Bonferroni correction
res$adj <- p.adjust(res$pval, method = "bonferroni")
# take the negative log of p-values
res$log <- -log10(res$adj)
# save the result
file <- "res.plac.seq.pos.0.01.rep2.csv"
path <- file.path(dir.out, file)
write.csv(res, path); rm(path)
# rename objects
res2 <- res; rm(res)
sig.string.list2 <- sig.string.list; rm(sig.string.list)

### combine and visualize the results ###

# read back in data
log.pval.matrix <- cbind(res1$log, res2$log, res.merged$log)
rownames(log.pval.matrix) <- rownames(res.merged)
colnames(log.pval.matrix) <- c("rep1", "rep2", "merged")
log.pval.matrix <- log.pval.matrix[c(1:7, 10:11), ]

file <- "log.pval.matrix.csv"
path <- file.path(dir.out, file)
if (!file.exists(path)) {write.csv(log.pval.matrix, path)}; rm(path)

# generate heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "Greens"))(100)
limit <- ceiling(max(log.pval.matrix, na.rm = TRUE)/5)*5
breaks <- seq(0, limit, limit/100)
labels.row <- c(
  "LinkPeaks",
	"Marginal",
	"Interaction",
	"TRIPOD conditional matching Xt",
  "TRIPOD interaction matching Xt",
  "TRIPOD conditional matching Yj",
  "TRIPOD interaction matching Yj",
  "TRIPOD union",
  "Union"
)
labels.col <- c(
	"Replicate 1",
	"Replicate 2",
	"Merged"
)
v <- c(log.pval.matrix)
number.color <- rep("gray10", length(v))
number.color[v > 20] <- "gray90"

file <- "heatmap_placseq_log_pval.pdf"
path <- file.path(dir.fig, file)
file.exists(path)
x <- pheatmap(
	log.pval.matrix,
  color = color,
	breaks = breaks,
  cluster_rows = FALSE,
	cluster_cols = FALSE,
  labels_row = labels.row,
	labels_col = labels.col,
  # show_rownames = FALSE,
	# show_colnames = FALSE,
  border_color = NA,
	fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
	display_numbers = TRUE,
	number_color = number.color, legend = FALSE
)
path <- file.path(dir.fig, file)
pdf(path, width = 2.5, height = 3.25)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
# dev.print(pdf, path, width = 2.5, height = 3.25); rm(path)
# graphics.off()

# draw three circle Venn diagram
int.string.list <- lapply(sig.string.list.merged, intersect, all.pairs)
marginal <- int.string.list[[2]]
string.list <- int.string.list[c(1, 3:7, 10)]

titles <- c(
  "LinkPeaks",
	"Interaction",
	"TRIPOD conditional\nmatching Xt",
  "TRIPOD interaction\nmatching Xt",
  "TRIPOD conditional\nmatching Yj",
  "TRIPOD interaction\nmatching Yj",
  "TRIPOD union"
)

p.list <- vector("list", length(string.list))
for (i in 1:length(string.list)) {
	set1 <- marginal
  set2 <- string.list[[i]]
  set3 <- plac.pairs
  set.list <- list(set1, set2, set3)
  union <- Reduce("union", set.list)
  factor.list <- lapply(set.list, factor, levels = union)
  bool.list <- lapply(factor.list, function(x) as.logical(as.vector(table(x))))
  t <- tibble(
    value = union,
	  set1 = bool.list[[1]],
	  set2 = bool.list[[2]],
  	set3 = bool.list[[3]])
	p.list[[i]] <- ggplot(t, aes(A = set1, B = set2, C = set3)) +
    geom_venn(
    	show_percentage = FALSE,
    	set_names = c(
		    "Marginal      ",
		    "Modeling",
		    "PLAC-seq"),
    	# set_name_color = "white",
    	fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF"),
    	# text_size = 1.8, stroke_size = 0.2, set_name_size = 1) +
    	text_size = 2.5, stroke_size = 0.2, set_name_size = 3) +
		theme_void() +
		theme(plot.title = element_text(size = 9.5, hjust = 0.5), # , margin = margin(t = 10, b = -20)),
			plot.margin = unit(c(0, 0, 0, 0), "cm"),
      panel.spacing = unit(0, "lines")) +
		coord_fixed() +
	  ggtitle(titles[i])
}

file <- "venn_plac_seq_all.pdf"
path <- file.path(dir.fig, file)
file.exists(path)
pdf(path)
wrap_plots(p.list) + plot_layout(ncol = 4)
dev.off(); rm(path)
