# load packages
library(Seurat)
library(Signac)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(GenomeInfoDb)
library(patchwork)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
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
file <- "e18.qc.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

# get motif data
pfm.set <- getMatrixSet(
	x = JASPAR2020,
	opts = list(species = 9606, all_versions = FALSE))

pfm.olig2 <- PFMatrix(
	ID = "PFM.Olig2",
	name = "Olig2",
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

motif.matrix <- CreateMotifMatrix(
  features = granges(e18),
  pwm = pfm.set,
  genome = "mm10",
  use.counts = FALSE
)

motif.object <- CreateMotifObject(data = motif.matrix, pwm = pfm.set)
e18 <- SetAssayData(e18, assay = "ATAC", slot = "motifs",
	new.data = motif.object)
e18 <- RunChromVAR(
  object = e18,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

file <- "e18.chromvar.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)
