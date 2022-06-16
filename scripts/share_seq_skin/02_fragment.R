
dir.out <- "./"
dir.in <- "./"

library(Seurat)
library(Signac)

## add a fragment object
# read in the skin seurat object
file <- "skin.rda"
path <- paste(dir.in, file, sep="")
load(path)

# read in the barcode mapping
file <- 'celltype_skin.txt'
path <- paste(dir.in, file, sep="")
ct <- read.table(path, header = T, sep='\t')
dim(ct)
# [1] 34774     3

# rename cells
head(x = colnames(x = skin))
# [1] "R1.01.R2.01.R3.06.P1.55" "R1.01.R2.03.R3.68.P1.55"
# [3] "R1.01.R2.05.R3.15.P1.53" "R1.01.R2.05.R3.40.P1.55"
# [5] "R1.01.R2.05.R3.49.P1.55" "R1.01.R2.06.R3.14.P1.55"
length(colnames(skin))
# [1] 29308
table(!is.na(match(colnames(skin), ct$rna.bc)))
# 
# TRUE
# 29308
atac.bc <- ct$atac.bc[match(colnames(skin), ct$rna.bc)]
# atac.bc <- ct$atac.bc[!is.na(match(ct$rna.bc, colnames(skin)))]
skin <- RenameCells(object = skin, new.names = atac.bc)
head(x = colnames(x = skin))
# [1] "R1.01.R2.01.R3.06.P1.07" "R1.01.R2.03.R3.68.P1.07"
# [3] "R1.01.R2.05.R3.15.P1.05" "R1.01.R2.05.R3.40.P1.07"
# [5] "R1.01.R2.05.R3.49.P1.07" "R1.01.R2.06.R3.14.P1.07"

# create a fragment object
file <- "skin.atac.fragments.tsv.gz"
path <- paste(dir.in, file, sep="")
frags <- CreateFragmentObject(path)
DefaultAssay(skin) <- "ATAC"
# data <- GetAssayData(skin, slot = "data")
# chrom_assay <- skin[["ATAC"]] 
# Fragments(chrom_assay) 
Fragments(skin) <- frags
class(Fragments(skin)[[1]])
# [1] "Fragment"
# attr(,"package")
# [1] "Signac"

file <- "skin.frag.rds"
path <- paste(dir.in, file, sep="")
if (!file.exists(path)) saveRDS(skin.frag, file=path)

# # sanity check (see Figure 2H in Ma et al. (2020))
# # Gapdh
# fig <- paste("covplot_Gapdh.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Gapdh",
#   annotation = T,
#   peaks = T,
#   extend.downstream = 10000
# )
# dev.off()
# 
# # Krt5
# fig <- paste("covplot_Krt5.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Krt5",
#   annotation = T,
#   peaks = T,
#   extend.downstream = 10000
# )
# dev.off()
# 
# # Krt15
# fig <- paste("covplot_Krt15.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Krt15",
#   annotation = T,
#   peaks = T,
#   extend.upstream = 10000
# )
# dev.off()
# 
# # Lgr5
# fig <- paste("covplot_Lgr5.pdf", sep="")
# path <- paste(dir.fig, fig, sep="")
# pdf(path, height=8, width=12)
# CoveragePlot(
#   object = skin,
#   region = "Lgr5",
#   annotation = T,
#   peaks = T,
#   extend.upstream = 10000
# )
# dev.off()
