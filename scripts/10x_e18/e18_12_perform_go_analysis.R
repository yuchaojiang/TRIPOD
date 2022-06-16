# set directories
dir.in <- "data"
dir.out <- "output"
dir.fig <- "figures"
dir.r <- "functions"
dir.david <- "david"

# load packages
library(RColorBrewer)
library(pheatmap)

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
hit.list <- readRDS(path); rm(path)

# select hits from the interaction model and TRIPOD
# with significantly positive coefficients
hit.list <- hit.list[seq(1, length(hit.list), 2)]
names(hit.list) <- gsub("[.]pos$", "", names(hit.list))
hit.list <- hit.list[-(1:3)]

# get neurogenesis and gliogenesis TFs
tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")
tfs <- c(tfs.neuro, tfs.glio)

# get the number of target genes for all TFs and all models
# get lists of target genes for all TFs and all models
gene.list.list <- lapply(
	X = hit.list,
	FUN = function(x) lapply(
		X = tfs,
		FUN = function(y) unique(x$gene[x$TF == y])))

# get the number of target genes
num.genes <- sapply(
  X = gene.list.list,
	FUN = function(x) sapply(x, length)
)
rownames(num.genes) <- tfs
file <- "num.target.genes.csv"
path <- file.path(dir.out, file)
write.csv(num.genes, path); rm(path)

# take the union of hits from all models
union.list <- list()
for (i in 1:length(tfs)) {
	tmp <- c()
	for (j in 1:length(gene.list.list)) {
		tmp <- c(tmp, gene.list.list[[j]][[i]])
	}
	union.list[[i]] <- unique(tmp)
}
names(union.list) <- tfs
num.genes.union <- sapply(union.list, length)

# create an input txt file for DAVID
max.num <- max(sapply(union.list, length))
tmp.list <- lapply(union.list, function(x) c(x, rep(NA, max.num - length(x))))
tmp.matrix <- do.call("cbind", tmp.list)
colnames(tmp.matrix) <- tfs
file <- paste0("david_input_genes.txt")
path <- file.path(dir.out, file)
write.table(tmp.matrix, path, sep = "\t", quote = FALSE, row.names = FALSE); rm(path)

# read in david results
fdr.thresh <- 0.05
sig.union.list <- list()
david.union.list <- list()
files <- list.files(dir.david, pattern = ".txt", full.names = TRUE)
names <- list.files(dir.david, pattern = ".txt")
names <- gsub(".txt", "", names)
david.union.list <- lapply(files, read.delim, header = TRUE)
names(david.union.list) <- names
sig.union.list <- lapply(david.union.list,
	function(x) x[x$FDR < fdr.thresh, c("Term", "Count", "FDR")])

file <- "david.sig.union.list.pos.rds"
path <- file.path(dir.out, file)
saveRDS(sig.union.list, path); rm(path)
file <- "david.union.list.pos.rds"
path <- file.path(dir.out, file)
saveRDS(david.union.list, path); rm(path)

num.david.union <- sapply(david.union.list, nrow)
num.sig.union <- sapply(sig.union.list, nrow)
num.union <- merge(num.david.union, num.sig.union, by = 0, all = TRUE)
num.union <- merge(num.genes.union, num.union, by.x = 0, by.y = 1, all = TRUE)
rownames(num.union) <- num.union[, 1]
num.union <- num.union[, -1]
num.union <- num.union[tfs, ]
colnames(num.union) <- c("gene", "david", "sig")

file <- "david.num.union.csv"
path <- file.path(dir.out, file)
write.csv(num.union, path); rm(path)

terms <- c()
for (i in 1:length(sig.union.list)) {
  tmp <- sig.union.list[[i]]$Term
  terms <- c(terms, tmp)
}
terms.table <- table(terms)
sorted.table <- sort(terms.table, decreasing = TRUE)
sorted.table[sorted.table >= 5]
# manually select seemingly most "informative" terms
selected.terms <- c(
  "GO:0006357~regulation of transcription from RNA polymerase II promoter",
  "GO:0045944~positive regulation of transcription from RNA polymerase II promoter",
  "GO:0000122~negative regulation of transcription from RNA polymerase II promoter",
  "GO:0048663~neuron fate commitment",
  "GO:0050767~regulation of neurogenesis",
  "GO:0030182~neuron differentiation",
  "GO:0021915~neural tube development",
  "GO:0021861~forebrain radial glial cell differentiation",
  "GO:0030900~forebrain development",
  "GO:0048708~astrocyte differentiation",
  "GO:0048709~oligodendrocyte differentiation",
  "GO:0001755~neural crest cell migration",
  "GO:0008285~negative regulation of cell proliferation",
  "GO:0045665~negative regulation of neuron differentiation",
  "GO:0035019~somatic stem cell population maintenance",
  "GO:0002052~positive regulation of neuroblast proliferation",
  "GO:0007049~cell cycle",
  "GO:0051301~cell division",
  "GO:0008284~positive regulation of cell proliferation"
)

merged <- merge(selected.terms,
	david.union.list[[1]][, c("Term", "Count", "FDR")], by = 1, all.x = TRUE)
for (i in 2:length(david.union.list)) {
  merged <- merge(merged, david.union.list[[i]][, c("Term", "Count", "FDR")],
  	by = 1, all.x = TRUE)
}
rownames(merged) <- merged[, 1]
merged <- merged[, -1]
merged <- merged[selected.terms, ]
count <- merged[, seq(1, ncol(merged), 2)]
fdr <- merged[, seq(2, ncol(merged), 2)]
colnames(count) <- colnames(fdr) <- names(david.union.list)
count <- count[, tfs]
fdr <- fdr[, tfs]
log.fdr <- apply(fdr, 2, function(x) -log(x, base = 10))

file <- "david.count.union.csv"
path <- file.path(dir.out, file)
write.csv(count, path); rm(path)

file <- "david.log.fdr.union.csv"
path <- file.path(dir.out, file)
write.csv(log.fdr, path); rm(path)

# generate a heatmap
color <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100)

limit <- ceiling(max(log.fdr, na.rm = TRUE))
breaks <- seq(0, limit, limit/100)

labels.row <- sapply(rownames(log.fdr), function(x) unlist(strsplit(x, "~"))[2])
names(labels.row) <- NULL

v <- c(log.fdr)
number.color <- rep("gray10", length(v))
number.color[v > 8] <- "gray90"

file <- "david_log_fdr_union_no_legend.pdf"
path <- file.path(dir.fig, file)
x <- pheatmap(
	log.fdr,
  color = color,
	breaks = breaks,
  cluster_rows = F,
	cluster_cols = F,
  labels_row = labels.row,
	# labels_col = labels.col,
  # show_rownames = F,
	# show_colnames = F,
  border_color = NA,
	fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
	display_numbers = count,
	number_color = number.color,
	legend = FALSE
)

path <- file.path(dir.fig, file)
pdf(path, width = 4.75, height = 4)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
# dev.print(pdf, path, width = 4.75, height = 4); rm(path)
# graphics.off()

# generate a simple heatmap
simple.terms <- c(
  "GO:0006357~regulation of transcription from RNA polymerase II promoter",
  "GO:0048709~oligodendrocyte differentiation",
  "GO:0045665~negative regulation of neuron differentiation",
  "GO:0035019~somatic stem cell population maintenance",
  "GO:0007049~cell cycle"
)

limit <- ceiling(max(log.fdr[simple.terms, ], na.rm = TRUE))
# limit <- ceiling(max(log.fdr, na.rm = TRUE)/10)*10
limit <- 12
breaks <- seq(0, limit, limit/100)

labels.row <- sapply(rownames(log.fdr[simple.terms, ]), function(x) unlist(strsplit(x, "~"))[2])
names(labels.row) <- NULL

v <- c(log.fdr[simple.terms, ])
number.color <- rep("gray10", length(v))
number.color[v > 8] <- "gray90"

file <- "david_log_fdr_union_simple_no_legend.pdf"
path <- file.path(dir.fig.go, file)
x <- pheatmap(
	log.fdr[simple.terms, ],
  color = color,
	breaks = breaks,
  cluster_rows = F,
	cluster_cols = F,
  labels_row = labels.row,
	# labels_col = labels.col,
  # show_rownames = F,
	# show_colnames = F,
  border_color = NA,
	fontsize = 6,
  fontsize_row = 6,
  fontsize_col = 6,
	display_numbers = count[simple.terms, ],
	number_color = number.color,
	legend = FALSE
)
path <- file.path(dir.fig, file)
pdf(path, width = 4.75, height = 1.6)
grid::grid.newpage()
grid::grid.draw(x$gtable)
dev.off()
# dev.print(pdf, path, width = 4.75, height = 1.6); rm(path)
# graphics.off()

# generate color legend
palette <- color
main <- expression(-log[10](~italic(P)~value))
file <- "legend_go_horizontal.pdf"
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
axis(1, at = c(-1, 0, 1), labels = c(0, 6, 12), las = 1, cex.axis = 0.7, lwd = 0)
dev.off()
