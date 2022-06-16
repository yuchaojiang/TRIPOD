# load packages
library(GenomicRanges)
library(nbpMatching)
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
file <- "transcripts.gr.rds"
path <- file.path(dir.out, file)
transcripts.gr <- readRDS(path); rm(path)

file <- "peaks.gr.rds"
path <- file.path(dir.out, file)
peak.gr <- readRDS(path); rm(path)

file <- "motifxTF.rds"
path <- file.path(dir.out, file)
motifxTF <- readRDS(path); rm(path)

file <- "peakxmotif.rds"
path <- file.path(dir.out, file)
peakxmotif <- readRDS(path); rm(path)

file <- "metacell.rna.rds"
path <- file.path(dir.out, file)
wnn.rna <- readRDS(path); rm(path)

file <- "metacell.peak.rds"
path <- file.path(dir.out, file)
wnn.peak <- readRDS(path); rm(path)

file <- "hvg.rds"
path <- file.path(dir.out, file)
hvg <- readRDS(path); rm(path)

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
wnn.celltype <- readRDS(path); rm(path)

file <- "metacell.celltype.col.rds"
path <- file.path(dir.out, file)
wnn.celltype.col <- readRDS(path); rm(path)

# get a list of target genes of interest
tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")
genes <- unique(c(hvg[1:1000], tfs.neuro, tfs.glio))
file <- "genes.rds"
path <- file.path(dir.out, file)
saveRDS(genes, path); rm(path)

# set window size
ext.upstream <- ext.downstream <- 2e5

# get XY matrices
system.time(
xymats.list <- bplapply(
  genes,
  getXYMatrices,
  ext.upstream = ext.upstream,
  transcripts.gr = transcripts.gr,
  peak.gr = peak.gr,
  wnn.rna = wnn.rna,
  wnn.peak = wnn.peak,
  peakxmotif = peakxmotif,
  motifxTF = motifxTF,
  wnn.celltype = wnn.celltype,
  wnn.celltype.col = wnn.celltype.col
)
)
names(xymats.list) <- genes
file <- "xymats.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats.list, path); rm(path)

# fit model 1a (marginal)
system.time(
xymats1.list <- bplapply(
  xymats.list,
  fitModel,
  modelName = "1a"
)
)
names(xymats1.list) <- genes
file <- "xymats1.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats1.list, path); rm(path)

# fit model 4a (interaction)
system.time(
xymats4.list <- bplapply(
  xymats.list,
  fitModel,
  modelName = "4a"
)
)
names(xymats4.list) <- genes
file <- "xymats4.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats4.list, path); rm(path)

# fit model 5a (TRIPOD)
system.time(
xymats5X.list <- bplapply(
  xymats.list,
  fitModel,
  modelName = "5a",
  match.by = "Xt",
  do.capValues = TRUE,
  VERBOSE = FALSE
)
)
names(xymats5X.list) <- genes
file <- "xymats5X.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats5X.list, path); rm(path)

system.time(
xymats5Y.list <- bplapply(
  xymats.list,
  fitModel,
  modelName = "5a",
  match.by = "Yj",
  do.capValues = TRUE,
  VERBOSE = FALSE
)
)
names(xymats5Y.list) <- genes
file <- "xymats5Y.list.rds"
path <- file.path(dir.out, file)
saveRDS(xymats5Y.list, path); rm(path)
