# set directories
dir.in <- "source_data"
dir.out <- "derived_data"
dir.fig <- "figures"
dir.r <- "functions"
dir.david <- "david"

# load packages
# library(Seurat)
# library(Signac)
library(GenomicRanges)
library(TFBSTools)
library(JASPAR2020)
library(BSgenome.Mmusculus.UCSC.mm10)

# source functions
scripts <- list.files(dir.r, full.names = T, pattern = ".R$")
tmp <- lapply(scripts, source)
rm(tmp)

# read in data
file <- "transcripts.gr.rds"
path <- file.path(dir.out, file)
transcripts.gr <- readRDS(path); rm(path)

file <- "motifxTF.rds"
path <- file.path(dir.out, file)
motifxTF <- readRDS(path); rm(path)

file <- "xymats.list.rds"
path <- file.path(dir.out, file)
xymats.list <- readRDS(path); rm(path)

file <- "hit.list.0.01.rds"
path <- file.path(dir.out, file)
hit.list <- readRDS(path); rm(path)

file <- "chip.gr.list.rds"
path <- file.path(dir.out, file)
chip.gr.list <- readRDS(path); rm(path)

### examine whether neurogenesis TFs may regulate gliogenesis TFs (crosstalk) ###

# select hits with significantly positive coefficients
hit.list <- hit.list[seq(1, length(hit.list), 2)]

# select hits from TRIPOD
names(hit.list) <- gsub("[.]pos$", "", names(hit.list))
hit.list <- hit.list[-(1:4)]

# get neurogenesis and gliogenesis TFs
# tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")

# extract Neurog2, Eomes, and Tbr1 data
chip.gr.list <- chip.gr.list[match(c("neurog2", "eomes", "tbr1"), names(chip.gr.list))]

TF.names <- c("Neurog2", "Eomes", "Tbr1")

ol.list <- list()
for (i in 1:length(tfs.glio)) {
    gene.name <- tfs.glio[i]
	xymats <- xymats.list[[gene.name]]
	vicinity.size <- c(2e5, 2e5)
	peak.gr.g <- xymats$peak.gr.g
	peakxmotif.g <- xymats$peakxmotif.g
	ol.matrix <- matrix(data = NA, nrow = length(peak.gr.g), ncol = length(chip.gr.list))
    for (j in 1:length(chip.gr.list)) {
        TF.name <- TF.names[j]
    chip.gr <- chip.gr.list[[j]]
    motif.list <- getMotifSites(
	    gene.name = gene.name,
      TF.name = TF.name,
	    vicinity.size = vicinity.size,
	    transcripts = transcripts.gr,
	    motifxTF = motifxTF,
	    BSgenome = BSgenome.Mmusculus.UCSC.mm10
    )
    motif.gr <- motif.list$sites
    vicinity.gr <- motif.list$vicinity
    intersect.gr <- GenomicRanges::intersect(chip.gr, motif.gr, ignore.strand = TRUE)
    ol <- as.numeric(overlapsAny(peak.gr.g, intersect.gr, ignore.strand = TRUE))
                ol.matrix[, j] <- ol
    }
	colnames(ol.matrix) <- TF.names
	rownames(ol.matrix) <- rownames(peakxmotif.g)
	ol.list[[i]] <- ol.matrix
}
names(ol.list) <- tfs.glio

file <- "crosstalk.ol.list.rds"
path <- file.path(dir.out, file)
saveRDS(ol.list, path); rm(path)

# get indicators as to whether a given ATAC peak is linked to a target gene via a given TF
sig.list <- list()
for (i in 1:length(tfs.glio)) {
	gene.name <- tfs.glio[i]
	xymats <- xymats.list[[gene.name]]
  peak.gr <- xymats$peak.gr
  peakxmotif.g <- xymats$peakxmotif.g
  mat <- matrix(0, nrow = nrow(peakxmotif.g), ncol = length(TF.names))
  rownames(mat) <- rownames(peakxmotif.g)
  colnames(mat) <- (TF.names)
	for (j in 1:length(TF.names)) {
		TF.name <- TF.names[j]
		x.list <- sapply(hit.list, function(x) x$peak[x$gene == gene.name & x$TF == TF.name])
    x <- do.call("c", x.list)
    indicator <- rep(0, nrow(peakxmotif.g))
    indicator[!is.na(match(rownames(peakxmotif.g), x))] <- 1
    mat[, j] <- indicator
	}
  sig.list[[i]] <- mat
}
names(sig.list) <- tfs.glio

file <- "crosstalk.sig.list.rds"
path <- file.path(dir.out, file)
saveRDS(sig.list, path); rm(path)

# combine the list of matrices
combined.list <- mapply(cbind, ol.list, sig.list)
selected.list <- lapply(combined.list, function(x) x[rowSums(x) > 0, , drop = FALSE])
selected.list <- lapply(selected.list, function(x) x[, c(1, 4, 2, 5, 3, 6), drop = FALSE])
selected.list <- selected.list[sapply(selected.list, function(x) nrow(x) > 0)]
validated.list <- lapply(
	selected.list,
	function(x) x[x[, 1] * x[, 2] == 1 | x[, 3] * x[, 4] == 1 | x[, 5] * x[, 6] == 1, ,
		drop = FALSE]
)
validated.list <- validated.list[sapply(validated.list, function(x) sum(x) > 0)]

file <- "crosstalk.validated.list.rds"
path <- file.path(dir.out, file)
saveRDS(validated.list, path); rm(path)

validated <- data.frame()
for (i in 1:length(validated.list)) {
	gene.name <- names(validated.list)[i]
  x <- validated.list[[i]]
  y <- as.data.frame(matrix(c(t(x)), ncol = 2, byrow = TRUE))
  y <- cbind(rep(colnames(x)[seq(1, ncol(x), 2)], nrow(x)), y)
  y <- cbind(rep(rownames(x), each = ncol(x)/2), y)
  y <- cbind(rep(gene.name, nrow(y)), y)
  colnames(y) <- c("gene", "peak", "TF", "chipseq", "model")
  validated <- rbind(validated, y)
}
validated <- validated[rowSums(validated[, 4:5]) == 2, 1:3]

file <- "crosstlak.validated.csv"
path <- file.path(dir.out, file)
write.csv(validated, path); rm(path)

### validate known pathways ###

# get matrix representing the known regulatory relationships
known.matrix <- rbind(
	c("Pax6", "Neurog2"),
	c("Neurog2", "Eomes"),
	c("Eomes", "Neurod1"),
  c("Neurod1", "Tbr1"),
  c("Olig2", "Sox10"),
  c("Olig2", "Nkx2-2"),
	c("Sox9", "Nfia")
)
known.sig.list <- list()
for (i in 1:nrow(known.matrix)) {
	TF.name <- known.matrix[i, 1]
	gene.name <- known.matrix[i, 2]
	xymats <- xymats.list[[gene.name]]
	peakxmotif.g <- xymats$peakxmotif.g
	x.list <- sapply(hit.list, function(x) x$peak[x$gene == gene.name & x$TF == TF.name])
  x <- do.call("c", x.list)
  indicator <- rep(0, nrow(peakxmotif.g))
  indicator[!is.na(match(rownames(peakxmotif.g), x))] <- 1
  names(indicator) <- rownames(peakxmotif.g)
  known.sig.list[[i]] <- indicator
}
names(known.sig.list) <- paste(known.matrix[, 1], known.matrix[, 2], sep = "_")

file <- "known.sig.list.rds"
path <- file.path(dir.out, file)
saveRDS(known.sig.list, path); rm(path)

sig.df <- data.frame()
for (i in 1:length(known.sig.list)) {
	pair <- unlist(strsplit(names(known.sig.list)[i], "_"))
	gene.name <- pair[2]
	TF.name <- pair[1]
	peaks <- names(known.sig.list[[i]])
	indicators <- known.sig.list[[i]]
	df <- data.frame(
		gene = rep(gene.name, length(peaks)),
		TF = rep(TF.name, length(peaks)),
		peak = peaks,
		indicator = indicators
	)
	sig.df <- rbind(sig.df, df)
}
rownames(sig.df) <- NULL

all.df <- sig.df
file <- "known.all.csv"
path <- file.path(dir.out, file)
write.csv(all.df, path); rm(path)

sig.df <- sig.df[sig.df$indicator == 1, 1:3]
file <- "known.sig.csv"
path <- file.path(dir.out, file)
write.csv(sig.df, path); rm(path)
