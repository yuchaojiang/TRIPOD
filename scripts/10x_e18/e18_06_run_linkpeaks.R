# load packages
library(Seurat)
library(Signac)
library(BSgenome.Mmusculus.UCSC.mm10)

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
file <- "e18.metacell.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path) ; rm(path)

file <- "genes.rds"
path <- file.path(dir.out, file)
genes <- readRDS(path) ; rm(path)

# run linkpeaks
DefaultAssay(e18) <- "ATAC"
e18 <- RegionStats(
	object = e18,
	genome = BSgenome.Mmusculus.UCSC.mm10)
e18 <- LinkPeaks(
  object = e18,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  pvalue_cutoff = 1,
  score_cutoff = 0,
  distance = 2e+05,
  method = "pearson",
  genes.use = genes)
links <- Links(e18)
# correct p-values for two-sided tests
links$pvalue <- 2*pnorm(abs(links$zscore), lower.tail = FALSE)
file <- "links.rds"
path <- file.path(dir.out, file)
saveRDS(links, path); rm(path)

# set FDR < 0.01
fdr.thresh <- 0.01
links$adj <- p.adjust(links$pvalue, method = "BH")
links.pos <- links[links$score > 0 & links$adj < fdr.thresh]
links.neg <- links[links$score < 0 & links$adj < fdr.thresh]

# create data frames (gene, peak, coef, pval, adj)
columns <- c("gene", "peak", "coef", "pval", "adj")

pos.df <- as.data.frame(mcols(links.pos))
pos.df <- pos.df[, c(2, 3, 1, 5, 6)]
colnames(pos.df) <- columns

neg.df <- as.data.frame(mcols(links.neg))
neg.df <- neg.df[, c(2, 3, 1, 5, 6)]
colnames(neg.df) <- columns

# save to csv files
file <- paste0("linkpeaks.pos.", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(pos.df, path); rm(path)

file <- paste0("linkpeaks.neg.", fdr.thresh, ".csv")
path <- file.path(dir.out, file)
write.csv(neg.df, path); rm(path)
