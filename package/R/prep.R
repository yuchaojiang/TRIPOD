#' Get objects required for model fitting
#'
#' This function takes a Seurat object and returns a list object containing
#' elements that are required for model fitting.
#'
#' @param object a Seurat object.
#' @param chr a character string specifying chromosomes. This would be
#' `paste0("chr", 1:22)` for the autosomes in human.
#'
#' @return This function returns a list object with the following elements.
#' \describe{
#'   \item{transcripts.gr}{a GRanges object containing chromosomal coordinates
#'   of annotated protein containing genes.}
#'   \item{peaks.gr}{a GRanges object containing chromosomal coordinates
#'   of the ATAC peak regions in the multiomic data.}
#'   \item{motifxTF}{a matrix containing mapping between TFs and their binding motifs.}
#'   \item{peakxmotif}{a binary sparse matrix containing indicators
#'   as to whether a given motif is present in a given ATAC peak region.}
#' }
#' @import GenomicRanges Seurat Signac
#' @export
#'
getObjectsForModelFit <- function(object,
                                  chr){
    DefaultAssay(object) <- "ATAC"
    transcripts.gr <- Signac:::CollapseToLongestTranscript(ranges = Annotation(object))
    transcripts.gr <- transcripts.gr[transcripts.gr$gene_biotype == "protein_coding"]
    transcripts.gr <- transcripts.gr[as.character(seqnames(transcripts.gr)) %in% chr]
    transcripts.gr <- sort(transcripts.gr)
    peaks.gr <- object@assays$ATAC@ranges
    motifxTF <- unlist(object@assays$ATAC@motifs@motif.names)
    motifxTF <- cbind(names(motifxTF), motifxTF)
    colnames(motifxTF) <- c("motif", "TF")
    peakxmotif <- object@assays$ATAC@motifs@data
    
    # keep motifs when genes encoding corresponding TFs are present in the RNA-seq data
    sel <- motifxTF[, 2]%in%toupper(rownames(object@assays$RNA))
    peakxmotif <- peakxmotif[, sel]
    motifxTF <- motifxTF[sel, ]
  
    # Change the TF name to be consistent with genes (esp. for the mouse genome)
    motifxTF[, 2] <- rownames(object@assays$RNA)[match(motifxTF[, 2], toupper(rownames(object@assays$RNA)))]
	
    # perform further data filtering to ensure that the same set of genes is included in all objects
    genes <- intersect(transcripts.gr$gene_name, rownames(object@assays$SCT))
    DefaultAssay(object) <- "RNA"
    transcripts.gr <- transcripts.gr[match(genes, transcripts.gr$gene_name)]
    peakxmotif <- peakxmotif[, motifxTF[, 2] %in% genes]
    motifxTF <- motifxTF[motifxTF[, 2] %in% genes, ]

    result <- list(transcripts.gr = transcripts.gr,
                   peaks.gr = peaks.gr,
                   motifxTF = motifxTF,
                   peakxmotif = peakxmotif
                   )
    return(result)
}

#' Filter a Seurat object
#'
#' This function filters a Seurat object and an output
#' from {\code{\link{getObjectsForModelFit}}}
#' to ensure that
#' TF binding motifs are kept only when genes encoding corresponding TFs are
#' present in the RNA-seq data and that the same set of genes is included
#' in all objects that are necessary for model fitting.
#'
#' @param object a Seurat object.
#' @param tripod.object a list object returned by {\code{\link{getObjectsForModelFit}}}
#'
#' @return This function returns a Seurat object.
#' @import Seurat Signac
#' @export
filterSeuratObject <- function(object, tripod.object){
    genes <- tripod.object$transcripts.gr$gene_name
    motifxTF <- tripod.object$motifxTF
    object@assays$RNA <- subset(object@assays$RNA,
                                features = match(genes, rownames(object@assays$RNA)))
    object@assays$SCT <- subset(object@assays$SCT,
                                features = match(genes, rownames(object@assays$SCT)))
    object@assays$chromvar <- subset(object@assays$chromvar,
                                     features = match(motifxTF[, 1], rownames(object@assays$chromvar)))
    return(object)
}

#' Process a Seurat object
#'
#' This is a wrapper function to perform normalization, dimension reduction,
#' nearest-neighbor graph construction, and UMAP embedding for the RNA and
#' ATAC modalities. It also calculates a weighted nearest neighbor (WNN) graph.
#'
#' @param object a Seurat object.
#' @param dim.rna a integer vector passed to {\code{\link{RunUMAP}}},
#' {\code{\link{FindNeighbors}}}, and {\code{\link{FindMultiModalNeighbors}}}
#' for the RNA modality.
#' @param dim.atac a integer vector passed to {\code{\link{RunUMAP}}},
#' {\code{\link{FindNeighbors}}}, and {\code{\link{FindMultiModalNeighbors}}}
#' for the ATAC modality.
#' @param verbose a logical variable.
#'
#' @return This function returns a Seurat object.
#' @import dplyr Seurat Signac
#' @export
processSeuratObject <- function(object, dim.rna = 1:50, dim.atac = 2:50,
	verbose = TRUE){
    # RNA analysis
    DefaultAssay(object) <- "RNA"
    object <- SCTransform(object, verbose = FALSE) %>%
	RunPCA(verbose = verbose) %>%
	RunUMAP(dims = dim.rna, reduction.name = "umap.rna", reduction.key = "rnaUMAP_",
		verbose = verbose) %>%
        FindNeighbors(reduction = "pca", dims = dim.rna)

    # ATAC analysis
    # exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(object) <- "ATAC"
    object <- RunTFIDF(object, verbose = verbose) %>%
	FindTopFeatures(min.cutoff = "q0", verbose = verbose) %>%
	RunSVD(verbose = verbose) %>%
	RunUMAP(reduction = "lsi", dims = dim.atac,
		reduction.name = "umap.atac", reduction.key = "atacUMAP_",
		verbose = verbose) %>%
        FindNeighbors(reduction = "lsi", dims = dim.atac)

  # calculate a WNN graph
  object <- FindMultiModalNeighbors(object, reduction.list = list("pca", "lsi"),
  	dims.list = list(dim.rna, dim.atac), verbose = verbose) %>%
  	RunUMAP(nn.name = "weighted.nn", reduction.name = "wnn.umap",
  		reduction.key = "wnnUMAP_", verbose = verbose)

  return(object)
}

#' Get a metacell matrix
#'
#' This function takes a Seurat object as an input and returns a matrix
#' containing normalized RNA expression or chromatin accessibility per
#' metacell.
#'
#' @param object a Seurat object.
#' @param cluster.name a character string specifying the cell clusters on which
#' metacells are based. This must be one of the column names in the meta data
#' of the Seurat object.
#' @param assay a character string specifying the assay. This must be one of
#' "RNA", "SCT", and "ATAC".
#' @param slot a character string specifying the slot. This must be one of
#' "counts" or "data".
#'
#' @return a matrix
#' @import Seurat Matrix GenomicRanges
#' @export
getMetacellMatrix <- function(object, cluster.name, assay, slot) {
    if (all(assay != c("RNA", "ATAC"))) {
        stop("The assay argument should be either \"RNA\" or \"ATAC\".")
    }
    if (!(cluster.name %in% colnames(object@meta.data))) {
        stop("The cluster.name argument does not match the meta data.")
    }
    clusters <- object@meta.data[[cluster.name]]
    clust.levels <- levels(clusters)
    assay.matrix <- GetAssay(object = object, assay = assay) %>% slot(slot)
    nrow <- length(clust.levels)
    ncol <- nrow(assay.matrix)
    metacell <- matrix(nrow = nrow, ncol = ncol)
    rownames(metacell) <- paste0("metacell_", clust.levels)
    colnames(metacell) <- rownames(assay.matrix)
    for (i in 1:nrow(metacell)) {
        metacell[i, ] <- apply(assay.matrix[, clusters ==
            as.character(i - 1)], 1, sum)
        metacell[i, ] <- metacell[i, ]/sum(metacell[i, ]) * 10^6
    }
    return(metacell)
}

#' Get metacell matrices
#'
#' This is a wrapper function of {\code{\link{getMetacellMatrix}}}.
#'
#' @inheritParams getMetacellMatrix
#' @param slot.rna a character string specifying the slot for the RNA modality.
#' This must be one of "counts" and "data".
#' @param slot.atac a character string specifying the slot for the ATAC modality.
#' This must be one of "counts" and "data."
#' @param min.num an integer to represent the minimum number of cells. Metacells
#' containing single cells fewer than this threshold will be removed.
#'
#' @return This function returns a list object containing the following elements.
#' \item{rna}{a matrix.}
#' \item{atac}{a matrix.}
#'
#' @import Seurat
#' @export
getMetacellMatrices <- function(object,
                                cluster.name = "seurat_clusters",
                                slot.rna = "counts",
                                slot.atac = "counts",
                                min.num = 0) {
    metacell.rna <- getMetacellMatrix(
        object = object,
        cluster.name = cluster.name,
        assay = "RNA",
        slot = slot.rna
    )
    metacell.peak <- getMetacellMatrix(
        object = object,
        cluster.name = cluster.name,
        assay = "ATAC",
        slot = slot.atac
    )
    remove <- as.vector(table(object@meta.data[[cluster.name]])) <
        min.num
    names.remove <- names(table(object@meta.data[[cluster.name]]))[remove]
    metacell.rna <- metacell.rna[!remove, ]
    metacell.peak <- metacell.peak[!remove, ]
    result <- list(rna = metacell.rna, peak = metacell.peak)
    return(result)
}

#' Remove small metacells
#'
#' This functions removes cells corresponding to metacells with cells fewer than
#' a threshold from a Seurat object. {\code{\link{getMetacellMatrices}}} must be
#' run before running this function.
#'
#' @inheritParams getMetacellMatrix
#' @param min.num an integer to represent the minimum number of cells. Metacells
#' containing single cells fewer than this threshold will be removed.
#' This must be set to the same value as the input argument to
#' {\code{\link{getMetacellMatrices}}}.
#'
#' @return A Seurat object.
#'
#' @import Seurat
#' @export
removeSmallMetacells <- function(
	object,
	min.num
) {
	remove <- as.vector(table(object$seurat_clusters)) < min.num
  names.remove <- names(table(object$seurat_clusters))[remove]
  # remove cells in the removed metacell clusters from the Seurat object
  object <- object[, !(object$seurat_clusters %in% names.remove)]
  object$seurat_clusters <- droplevels(object$seurat_clusters)
  return(object)
}

#' Optimize cluster resolutions
#'
#' This function takes a Seurat object as input and returns a data frame
#' containing cluster numbers at specified resolutions.
#'
#' @param object a Seurat object.
#' @param assay.name a character string specifying the assay name. This must
#' be one of "SCT", "ATAC", or "WNN".
#' @param graph.name a character string specifying the graph name. This must be
#' one of the graphs stored in the Seurat object.
#' @param algorithm an integer representing a clustering algorithm. See
#' {\code{\link{FindClusters}}} in the Seurat package.
#' @param resolutions a numerical vector specifying resolutions to be examined.
#' @param min.num an integer setting a threshold for the minimum number of
#' cells per cluster.
#' @param ... further arguments to be passed to {\code{\link{FindClusters}}}.
#'
#' @return a data frame with three columns:
#' \describe{
#'   \item{resolusion}{the resolution examined}
#'   \item{num_clusters}{the number of clusters}
#'   \item{num_below}{the number of clusters with fewer than
#'   the threshold number of single cells}
#' }
#' @import Seurat
#' @export
optimizeResolution <- function(
	object, assay.name, graph.name = NULL,
	algorithm = NULL, resolutions = seq(10, 35, 5), min.num, ...
) {
	if (all(assay.name != c("SCT", "ATAC", "WNN"))) {
    stop('The assay.name argument must be one of "SCT", "ATAC", or "WNN".')
	}
  # object <- FindNeighbors(object, dims = 1:30)
  num.clusters <- num.below <- rep(NA, length(resolutions))
  for (i in 1:length(resolutions)){
    res <- resolutions[i]
    if (assay.name %in% c("SCT", "ATAC")) {
    	# DefaultAssay(object) <- assay.name
    	if (is.null(graph.name)) graph.name <- paste0(assay.name, "_snn")
      if (is.null(algorithm)) algorithm <- 1
      # snn.name <- paste0(assay.name, "_snn_res.", res)
    } else if (assay.name == "WNN") {
    	if (is.null(graph.name)) graph.name <- "wsnn"
      if (is.null(algorithm)) algorithm <- 3
      # snn.name <- paste0("wsnn_res.", res)
    }
    object <- FindClusters(object, graph.name = graph.name,
    	resolution = res, algorithm = algorithm, verbose = FALSE, ...)
    meta.col.name <- paste0(graph.name, "_res.", res)
    num.clusters[i] <- length(levels(object@meta.data[[meta.col.name]]))
    num.below[i] <- sum(table(object@meta.data[[meta.col.name]]) < min.num)
  }
  results <- data.frame(resolution = resolutions,
                        num_clusters = num.clusters,
                        num_below = num.below)
  return(results)
}

#' Get metacell assignment
#'
#' This function takes an Seurat object as an input and
#' returns the object with a meta data column, `seurat_clusters`.
#' containing metacell assignment.
#'
#' @param object a Seurat object.
#' @param resolution a numeric variable specifying the resolution for clustering.
#' @param graph.name a character string specifying the graph name. This must be
#' one of the graphs stored in the Seurat object.
#' @param ... further arguments passed to {\code{\link{FindClusters}}}.
#'
#' @return This function returns a Seurat object.
#' @export
#' @import Seurat
getClusters <- function(object, resolution, graph.name, ...){
    object <- FindClusters(object, graph.name = graph.name,
    	resolution = resolution, ...)

    # reorder the factor levels
    column <- paste0(graph.name, "_res.", resolution)
        if (!(column %in% colnames(object@meta.data))) column <- "seurat_clusters"
    tmp <- as.character(object@meta.data[[column]])
    levels.clust <- as.character(sort(as.numeric(levels(object@meta.data[[column]]))))
    tmp <- factor(tmp, levels = levels.clust)
    object$seurat_clusters <- tmp
    if (column != "seurat_clusters") object[[column]] <- NULL

    return(object)
}

#' Get cell types for metacells
#'
#' @param object a Seurat object.
#' @param celltype.col.name a character string specifying the name of the
#' meta data column in the Seurat object containing cell types.
#' @param cluster.col.name a character string specifying the name of the
#' meta data column in the Seurat object containing cluster numbers.
#'
#' @return a character vector
#' @export
getCellTypeForMetacell <- function(
	object,
	celltype.col.name = "celltype",
	cluster.col.name = "seurat_clusters"
) {
  tmp <- table(unlist(object[[celltype.col.name]]),
  	unlist(object[[cluster.col.name]]))
  metacell.celltype <- rep(NA, ncol(tmp))
  for (i in 1:length(metacell.celltype)) {
    tmp.i_1 <- tmp[, colnames(tmp) == as.character(i-1)]
    metacell.celltype[i] <- names(tmp.i_1)[which.max(tmp.i_1)]
  }
  return(metacell.celltype)
}

#' Get a mapping between cell types and colors for single cells
#'
#' @param object a Seurat object.
#' @param reduction a character string specifying a dimension reduction method.
#' @param celltype.col.name a character string specifying the name of the meta
#' data column in the Seurat object containing cell types.
#'
#' @return a data frame
#' @export
getColorsForSingleCell <- function(
	object, reduction, celltype.col.name
) {
	# get the corresponding color for each cell type from Seurat
	p <- Seurat::DimPlot(object, reduction = reduction, label = TRUE,
		group.by = celltype.col.name)
  # use ggplot_build to deconstruct the ggplot object
  pbuild <- ggplot2::ggplot_build(p)
  # get the color palette by Seurat
  pdata <- pbuild$data[[1]]
  pdata <- cbind(object[[celltype.col.name]], pdata)
  sc.color.map <- pdata[, 1:2]
  colnames(sc.color.map) <- c("celltype", "color")
  return(sc.color.map)
}

#' Get a mapping between cell types and colors for meta cells
#'
#' @param metacell.celltype a character string
#' @param sc.color.map a data frame
#'
#' @return a data frame
#' @export
getColorsForMetacell <- function(metacell.celltype, sc.color.map) {
  metacell.celltype.col <- rep(NA, length(metacell.celltype))
  for (i in 1:length(metacell.celltype)) {
    metacell.celltype.col[i] <-
    	sc.color.map$color[min(which(sc.color.map$celltype == metacell.celltype[i]))]
  }
  metacell.color.map <- data.frame(
  	celltype = metacell.celltype,
	  color = metacell.celltype.col
  )
}

#' Get a mapping between cell types and colors
#'
#' @param ordered.celltypes a character string.
#' @param metacell.color.map a data frame.
#'
#' @return a data frame
#' @export
getColorsForCellType <- function(ordered.celltypes, metacell.color.map) {
  color.map <- unique(metacell.color.map)
  color.map <- color.map[match(ordered.celltypes, color.map$celltype), ]
  rownames(color.map) <- NULL
  return(color.map)
}

#' Get color coding based on cell types
#'
#' This function takes a Seurat object as input and
#' returns data frames containing mappings of colors to
#' single cells, metacells, and cell types.
#'
#' @inheritParams getCellTypeForMetacell
#' @inheritParams getColorsForSingleCell
#' @inheritParams getColorsForMetacell
#' @inheritParams getColorsForCellType
#'
#' @return a list of data frames
#' @export
getColors <- function(
	object, reduction = "umap.rna",
	celltype.col.name = "celltype", cluster.col.name = "seurat_clusters"
) {
	metacell.celltype <- getCellTypeForMetacell(
		object = object,
		celltype.col.name = celltype.col.name,
		cluster.col.name = cluster.col.name
	)
	sc.color.map <- getColorsForSingleCell(
		object = object, reduction = reduction,
		celltype.col.name = celltype.col.name
	)
	metacell.color.map <- getColorsForMetacell(
	  metacell.celltype = metacell.celltype, sc.color.map = sc.color.map)
  # metacell.celltype.col <- metacell.color.map$color
	celltype.color.map <- getColorsForCellType(
		ordered.celltypes = levels(droplevels(object@meta.data[[celltype.col.name]])),
		metacell.color.map = metacell.color.map)
	result <- list(sc = sc.color.map,
		metacell = metacell.color.map,
		celltype = celltype.color.map)
	return(result)
}
