# load packages
library(Seurat)
library(Signac)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(ggvenn)
library(nbpMatching)
library(gridBase)
library(Sushi)
library(BiocParallel)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2020)
library(ggplot2)
library(AnnotationHub)
library(rtracklayer)
library(Biostrings)
library(pheatmap)
require(RColorBrewer)

# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"

# download data into the source data directory
# http://hgdownload.cse.ucsc.edu/goldenpath/mm9/liftOver/mm9ToMm10.over.chain.gz
# https://hgdownload-test.gi.ucsc.edu/goldenPath/rn4/vsMm10/rn4.mm10.all.chain.gz
# https://hgdownload.soe.ucsc.edu/goldenPath/rn5/vsMm10/rn5.mm10.all.chain.gz
# olig2
# GSE42454_RAW/GSM1040156_Olig2_ChIP-seq_OPC.bed.gz
# GSE42454_RAW/GSM1040157_Olig2_ChIP-seq_iOL.bed.gz
# GSE42454_RAW/GSM1040158_Olig2_ChIP-seq_mOL.bed.gz
# neurog2
# GSE63620_RAW/GSM1553880_Neurog2.bed.gz
# eomes
# GSE63620_RAW/GSM1553879_Tbr2.bed.gz
# tbr1
# GSM1833461_mm9.tbr1.idr01.MACS2.bed.gz

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "genes.rds"
path <- file.path(dir.out, file)
genes <- readRDS(path); rm(path)

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "trio.string.list.0.01.rds"
path <- file.path(dir.out, file)
trio.string.list <- readRDS(path); rm(path)

# import chain files
file <- "mm9ToMm10.over.chain"
path <- file.path(dir.in, file)
mm9tomm10 <- import.chain(path)

file <- "rn4.mm10.all.chain"
path <- file.path(dir.in, file)
rn4tomm10 <- import.chain(path); rm(path)

file <- "rn5.mm10.all.chain"
path <- file.path(dir.in, file)
rn5tomm10 <- import.chain(path); rm(path)

# get file names for ChIP-seq data
files <- list.files(dir.in, pattern = ".bed.gz$", full.names = TRUE)

# create a list of GRanges object for olig2
files.olig2 <- files[grep("Olig2", files)]
peaks.olig2.rn4.list <- lapply(files.olig2, import)
names(peaks.olig2.rn4.list) <- c("opc", "iol", "mol") # the order is kept
peaks.olig2.list <- lapply(peaks.olig2.rn4.list, function(x) unlist(liftOver(x, rn4tomm10)))

# create a GRanges object for neurog2
file <- "GSM1553880_Neurog2.bed.gz"
path <- file.path(dir.in, file)
peaks.neurog2.mm9 <- import(path); rm(path)
peaks.neurog2 <- unlist(liftOver(peaks.neurog2.mm9, mm9tomm10))
strand(peaks.neurog2) <- "*"

# create a GRanges object for eomes
file <- "GSM1553879_Tbr2.bed.gz"
path <- file.path(dir.in, file)
peaks.eomes.mm9 <- import(path); rm(path)
peaks.eomes <- unlist(liftOver(peaks.eomes.mm9, mm9tomm10))
strand(peaks.eomes) <- "*"

# create a GRanges object for tbr1
file <- "GSM1833461_mm9.tbr1.idr01.MACS2.bed.gz"
path <- file.path(dir.in, file)
peaks.tbr1.mm9 <- import(path); rm(path)
peaks.tbr1 <- unlist(liftOver(peaks.tbr1.mm9, mm9tomm10))

chip.gr.list <- list(
	olig2.opc = peaks.olig2.list$opc,
	olig2.iol = peaks.olig2.list$iol,
	olig2.mol = peaks.olig2.list$mol,
	neurog2 = peaks.neurog2,
	eomes = peaks.eomes,
	tbr1 = peaks.tbr1
)

file <- "chip.gr.list.rds"
path <- file.path(dir.out, file)
saveRDS(chip.gr.list, path); rm(path)

# convert the strings back to data frames
hit.list <- lapply(
	X = trio.string.list,
	FUN = convertTrioStringToDF)

TF.names <- c(rep("Olig2", 3), "Neurog2", "Eomes", "Tbr1")

# compute q, m, n, and k
res.list.list <- list()
for (i in 1:length(chip.gr.list)) {
	print(i)
  chip.gr <- chip.gr.list[[i]]
  data.name <- names(chip.gr.list)[i]
  TF.name <- TF.names[i]
  res.list <- list()
  for (j in 1:length(hit.list)) {
  	# cat(j, "\t")
  	hits <- hit.list[[j]]
	  res <- getOverlapsAcrossGenes(
	  	genes = genes,
	  	TF.name = TF.name,
	  	xymats.list = xymats.list,
	  	hits = hits,
	  	subject = chip.gr
	  )
	  res.list[[j]] <- res
  }
  names(res.list) <- names(hit.list)
  res.list.list[[i]] <- res.list
  # cat("\n")
}
names(res.list.list) <- names(chip.gr.list)
file <- "res.chip.seq.list.list.0.01.rds"
path <- file.path(dir.out, file)
saveRDS(res.list.list, path); rm(path)

# write q, m, n, k, fold enrichment, p-value, and -log(p-value) to files
for (i in 1:length(res.list.list)) {
	x <- data.frame(t(sapply(res.list.list[[i]], colSums)))
  x$fe <- (x$q/x$k)/(x$m/(x$m+x$n))
  x$pval <- apply(x, 1, getHypergeometricPvalue)
  x$log <- -log10(x$pval)
  x <- x[c(seq(1, nrow(x), 2), seq(2, nrow(x), 2)), ]
  outfile <- paste0("res.chip.seq.0.01.", names(res.list.list)[i], ".csv")
  path <- file.path(dir.out, outfile)
  write.csv(x, path, quote = FALSE); rm(path)
}

enrichment.score.df <- data.frame(t(sapply(res.list.list, getEnrichmentScore)))
enrichment.score.pos <- enrichment.score.df[, seq(1, ncol(enrichment.score.df), 2)]
file <- "enrichment.score.pos.0.01.csv"
path <- file.path(dir.out, file)
write.csv(enrichment.score.pos, path); rm(path)

# generate heatmaps
color <- colorRampPalette(brewer.pal(n = 9, name = "Greens"))(100)

# file <- "enrichment.score.pos.0.01.csv"
# path <- file.path(dir.out, file)
# enrichment.score.pos <- read.csv(path, row.names = 1); rm(path)

x1 <- t(enrichment.score.pos)
x1 <- x1[c(1, 8), ]
limit <- ceiling(max(x1, na.rm = TRUE)/10)*10
breaks <- seq(0, limit, limit/100)

labels.row <- c(
	"Interaction",
	"TRIPOD union"
)

labels.col <- c(
	"Olig2 (OPC)",
	"Olig2 (iOL)",
	"Olig2 (mOL)",
	"Neurog2",
	"Eomes",
	"Tbr1"
)

file <- "tf_chipseq_enrichment_pos_0.01.pdf"
path <- file.path(dir.fig, file)
x <- pheatmap(
	x1,
  color = color,
	breaks = breaks,
  cluster_rows = F,
	cluster_cols = F,
  labels_row = labels.row,
	labels_col = labels.col,
  # show_rownames = F,
	# show_colnames=F,
  border_color = NA,
	fontsize = 6,
  fontsize_row = 7.5,
  fontsize_col = 6,
	display_numbers = TRUE,
	legend = FALSE
)
path <- file.path(dir.fig, file)
pdf(path, width = 2.5, height = 1)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()

# dev.print(pdf, path, width = 2.5, height = 1); rm(path)
# graphics.off()

palette <- color
main <- expression(-log[10](~italic(P)~value))

# generate  color legend
file <- "legend_chipseq_horizontal.pdf"
path <- file.path(dir.fig, file)
pdf(path, width = 1, height = 0.6); rm(path)
par(mar = c(1, 0.5, 1, 0.5))
par(mgp = c(2, -0.2, 0))
# par(mar = c(1.5, 1, 1, 1.5), mgp = c(2, 0.3, 0))
plot(
	c(-1, 1),
	c(0, 1),
  main = main,
  font.main = 1,
  cex.main = 0.7,
  type = "n",
  xlab = "",
  ylab = "",
  axes = FALSE,
  xaxs = "r",
  yaxs = "r",
  frame.plot = FALSE
)
n <- length(palette)
bin <- 2/n
for (i in 1:n){
	rect(
		xleft = -1 + (i-1)*bin,
    ybottom = 0,
    xright = -1 + i*bin,
    ytop = 1,
    col = palette[i],
    border = palette[i])
}
axis(1, at = c(-1, 0, 1), labels = c(0, 10, 20), las = 1, cex.axis = 0.7, lwd = 0)
dev.off()
