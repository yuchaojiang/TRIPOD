# load packages
library(Seurat)
library(Signac)

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
file <- "e18.intersect.rds"
path <- file.path(dir.out, file)
e18 <- readRDS(path); rm(path)

# get rna meta-cells by clustering cells on basis of their scRNA-seq profiles
DefaultAssay(e18) <- "SCT"
e18 <- FindNeighbors(e18, reduction = "pca", dims = 1:30)

# optimize clustering resolution
resolutions <- seq(10, 35, 5)
num.clusters <- rep(NA, length(resolutions))
for (i in 1:length(resolutions)){
  res <- resolutions[i]
  e18 <- FindClusters(e18, resolution = res, verbose = F)
  snn.name <- paste("SCT_snn_res.", res, sep="")
  num.clusters[i] <- length(levels(e18@meta.data[[snn.name]]))
}
snn.names <- paste("SCT_snn_res.", resolutions, sep = "")
min.cells <- sapply(snn.names, function(x) min(table(e18@meta.data[[x]])))
num.cluster.df <- data.frame(resolution = resolutions,
                             num_clusters = num.clusters,
                             min_cells = min.cells)

# write the result to a csv file
file <- "num_clusters.csv"
path <- file.path(dir.out, file)
write.csv(num.cluster.df, path); rm(path)

res <- 15
e18 <- FindClusters(e18, graph.name = "SCT_snn",
    	resolution = res, verbose = FALSE)

# reorder the factor levels
tmp <- as.character(e18$SCT_snn_res.15)
levels.15 <- as.character(sort(as.numeric(levels(e18$SCT_snn_res.15))))
tmp <- factor(tmp, levels = levels.15)
e18$SCT_snn_res.15 <- tmp
e18@meta.data$seurat_clusters <- e18@meta.data$SCT_snn_res.15
e18@meta.data$SCT_snn_res.15 <- NULL

# get metacell matrices
metacell.rna <- matrix(nrow = length(levels(e18$seurat_clusters)),
	ncol = nrow(e18@assays$RNA))
rownames(metacell.rna) <- paste("metacell_", levels(e18$seurat_clusters), sep = "")
colnames(metacell.rna) <- rownames(e18@assays$RNA)
for(i in 1:nrow(metacell.rna)){
  metacell.rna[i,] <- apply(e18@assays$RNA@counts[, e18$seurat_clusters == (i-1)], 1, sum)
  # adjust for library size for each metacell
  metacell.rna[i,] <- metacell.rna[i,]/sum(metacell.rna[i,])*10^6
}

metacell.peak <- matrix(nrow = length(levels(e18$seurat_clusters)),
	ncol = nrow(e18@assays$ATAC))
rownames(metacell.peak) <- paste("metacell_", levels(e18$seurat_clusters), sep = "")
colnames(metacell.peak) <- rownames(e18@assays$ATAC)
for(i in 1:nrow(metacell.peak)){
  metacell.peak[i,] <- apply(e18@assays$ATAC@counts[, e18$seurat_clusters == (i-1)], 1, sum)
  # adjust for library size for each metacell
  metacell.peak[i,] <- metacell.peak[i,]/sum(metacell.peak[i,])*10^6
}

# remove metacell clusters with fewer than 20 cells
remove <- as.vector(table(e18$seurat_clusters)) < 20
names.remove <- names(table(e18$seurat_clusters))[remove]
metacell.rna <- metacell.rna[!remove, ]
metacell.peak <- metacell.peak[!remove, ]

file <- "metacell.rna.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.rna, path); rm(path)

file <- "metacell.peak.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.peak, path); rm(path)

# remove cells in the removed metacell clusters from the Seurat object
e18 <- e18[, !(e18$seurat_clusters %in% names.remove)]
e18$seurat_clusters <- droplevels(e18$seurat_clusters)

file <- "e18.metacell.rds"
path <- file.path(dir.out, file)
saveRDS(e18, path); rm(path)

# assigning cell types and colors for visualization
temp <- table(e18$celltype, e18$seurat_clusters)
metacell.celltype <- rep(NA, nrow(metacell.rna))
for(i in 1:length(metacell.celltype)){
  temp.i_1 <- temp[, colnames(temp) == as.character(i-1)]
  metacell.celltype[i] <- names(temp.i_1)[which.max(temp.i_1)]
}

# get the corresponding color for each cell type from Seurat
p <- Seurat::DimPlot(e18, reduction = "umap.rna", label = TRUE, group.by = "celltype")
pbuild <- ggplot2::ggplot_build(p) # use ggplot_build to deconstruct the ggplot object
pdata <- pbuild$data[[1]] # this is to get the color palette by Seurat
pdata <- cbind(e18$celltype, pdata)

metacell.celltype.col <- rep(NA, length(metacell.celltype))
for(i in 1:length(metacell.celltype)){
  metacell.celltype.col[i] <- pdata$colour[min(which(pdata$`e18$celltype` == metacell.celltype[i]))]
}

file <- "metacell.celltype.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.celltype, path); rm(path)

file <- "metacell.celltype.col.rds"
path <- file.path(dir.out, file)
saveRDS(metacell.celltype.col, path); rm(path)

col.map <- data.frame(celltype = metacell.celltype, color = metacell.celltype.col)
col.map <- unique(col.map)
rownames(col.map) <- col.map$celltype
col.map <- col.map[levels(e18$seurat_clusters), ]
rownames(col.map) <- NULL

file <- "color.map.rds"
path <- file.path(dir.out, file)
saveRDS(col.map, path); rm(path)

# obtain top 3000 highly variable genes based on SCT
DefaultAssay(e18) <- "SCT"
hvg <- VariableFeatures(e18)

file <- "hvg.rds"
path <- file.path(dir.out, file)
saveRDS(hvg, path); rm(path)
