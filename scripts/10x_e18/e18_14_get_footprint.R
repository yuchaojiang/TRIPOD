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
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

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

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
wnn.celltype <- readRDS(path); rm(path)

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
wnn.celltype.col <- readRDS(path); rm(path)

file <- "color.map.rds"
path <- file.path(dir.out, file)
col.map <- readRDS(path); rm(path)

file <- "e18.metacell.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

# update the motif object
DefaultAssay(e18) <- "ATAC"

# create a motif object
pfm.set <- getMatrixSet(
	x = JASPAR2020,
	opts = list(species = 9606, all_versions = FALSE)) # human
# replace the olig2 motif with what is supported by the published literature
pfm.olig2 <- PFMatrix(
	ID = "PFM.Olig2",
	name = "OLIG2", # name = "olig2",
  # matrixClass = "Zipper-Type",
	strand = "+",
	bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  tags = list(
  	species = "10090",
    tax_group = "vertebrates"),
    profileMatrix = matrix(
    	c(0L,  1L, 0L,  0L,  0L,  0L, # A
        1L, 0L,  0L, 1L,  0L,  0L, # C
        0L,  0L,  1L,  0L, 0L,  1L, # G
        0L,  0L,  0L,  0L,  1L, 0L), # T
        byrow = TRUE, nrow = 4,
        dimnames=list(c("A", "C", "G", "T")))
)
pfm.set <- pfm.set[-which(sapply(pfm.set, function(x) x@name == "OLIG2"))]
pfm.set$PFM.Olig2 <- pfm.olig2

file <- "pfm.set.rds"
path <- file.path(dir.out, file)
saveRDS(pfm.set, path); rm(path)

motif.positions <- matchMotifs(
  pwms = pfm.set,
  subject = granges(e18),
  out = "positions",
  genome = "mm10"
)
motif.matrix <- CreateMotifMatrix(features = granges(e18),
	pwm = pfm.set, genome = "mm10", use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix,
	pwm = pfm.set, positions = motif.positions)

file <- "motif.positions.rds"
path <- file.path(dir.out, file)
saveRDS(motif.positions, path); rm(path)

file <- "motif.matrix.rds"
path <- file.path(dir.out, file)
saveRDS(motif.matrix, path); rm(path)

file <- "motif.object.rds"
path <- file.path(dir.out, file)
saveRDS(motif.object, path); rm(path)

# focus on neurogenesis and gliogenesis TFs
tfs.neuro <- c("Pax6", "Neurog2", "Eomes", "Neurod1", "Tbr1")
tfs.glio <- c("Olig2", "Sox10", "Nkx2-2", "Sox9", "Nfia", "Ascl1")
tfs <- c(tfs.neuro, tfs.glio)

e18 <- SetAssayData(e18, assay = "ATAC",
	slot = "motifs", new.data = motif.object)
TF.names <- tfs
for (i in 1:length(TF.names)) {
	TF.name <- toupper(TF.names[i])
  e18 <- Footprint(
	  object = e18,
	  genome = BSgenome.Mmusculus.UCSC.mm10,
    assay = "ATAC",
	  motif.name = TF.name
  )
}

file <- "e18.footprint.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)
