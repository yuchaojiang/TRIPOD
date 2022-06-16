#' Get objects required for model fitting
#'
#' This function takes a Seurat object and returns a list object containing
#' elements that are required for model fitting.
#' It allows for a motif corresponding to multiple TFs as wlle as
#' a TF corresponding to multiple motifs.
#'
#' @param object a Seurat object.
#' @param chr a character string specifying chromosomes. This would be
#' `paste0("chr", 1:22)` for the autosomes in human.
#' @param convert a logical variable to indicate whether the gene names need
#' to be converted. This needs to be set to `TRUE`
#' when a human motif dataset is used for analyzing mouse data.The current
#' version only supports conversion from human to mouse.
#'
#'
#' @return This function returns a list object with the following elements.
#' \item{transcripts.gr}{a GRanges object containing chromosomal coordinates
#' of annotated protein containing genes.}
#' \item{peaks.gr}{a GRanges object containing chromosomal coordinates
#' of the ATAC peak regions in the multiomic data.}
#' \item{motif.TF.map}{a data frame containing mapping between TFs and their binding motifs.}
#' \item{peakxmotif}{a binary sparse matrix containing indicators
#' as to whether a given motif is present in a given ATAC peak region.}
#'
#' @import GenomicRanges Seurat Signac
#' @export
#'
getObjectsForModelFit2 <- function(object,
                                   chr,
	                           convert = FALSE) {
    DefaultAssay(object) <- "ATAC"
    transcripts.gr <- Signac:::CollapseToLongestTranscript(ranges = Annotation(object))
    transcripts.gr <- transcripts.gr[transcripts.gr$gene_biotype == "protein_coding"]
    transcripts.gr <- transcripts.gr[as.character(seqnames(transcripts.gr)) %in% chr]
    transcripts.gr <- sort(transcripts.gr)

    # get ATAC peak regions
    peaks.gr <- object@assays$ATAC@ranges

    # get peak by motif, overlap of ATAC peaks and motifs (binary)
    peakxmotif <- object@assays$ATAC@motifs@data

    # get motif.TF.map
    motif.names <- unlist(object@assays$ATAC@motifs@motif.names)
    motif.TF.map <- data.frame(motif = names(motif.names), TF = motif.names)
    rownames(motif.TF.map) <- NULL

    # copy the original TF notations
    motif.TF.map$var <- motif.TF.map$TF

    # delete variant numbers (e.g., "BACH2(var.2)").
    motif.TF.map$TF <- gsub("[(]var[.].[)]", "", motif.TF.map$TF)

    # note that some rows in the data frames contain double colon-separated lists of TF gene symbols
    # expand the data frame so that each row contains one gene symbol
    # to do this, convert the data frame to a list,
    # in which each element corresponds to one row of the original data frame first
    # then expand each of the one-row data frames so that each row contains one target gene
    # then combine the data frames in the list

    expandRow <- function(x) {
        tfs <- unlist(strsplit(x[[2]], "::"))
        n <- length(tfs)
        d <- data.frame(motif = rep(x[[1]], n),
                        TF = tfs,
                        var = rep(x[[2]], n))
        d
    }

    motif.TF.list <- split(motif.TF.map, seq(nrow(motif.TF.map)))
    expanded.list <- lapply(motif.TF.list, expandRow)
    motif.TF.map <- as.data.frame(do.call("rbind", expanded.list))

    # the mouse has gene symbol as Brca1 v.s. the gene/TF symbol as BRCA1 in human
    if (convert) {
        motif.TF.map$TF <- stringr::str_to_title(tolower(motif.TF.map$TF))
    }

    # perform further data filtering to ensure that the same set of genes is included in all objects
    genes <- intersect(transcripts.gr$gene_name, rownames(object@assays$SCT))
    transcripts.gr <- transcripts.gr[match(genes, transcripts.gr$gene_name)]
    motif.TF.map <- motif.TF.map[motif.TF.map$TF %in% genes, ]
    motifs <- unique(motif.TF.map$motif)
    motifs <- intersect(motifs, colnames(peakxmotif))
    peakxmotif <- peakxmotif[, motifs]

    motif.TF.map <- motif.TF.map[motif.TF.map$motif %in% motifs, ]

    result <- list(transcripts.gr = transcripts.gr,
                   peaks.gr = peaks.gr,
                   motif.TF.map = motif.TF.map,
                   peakxmotif = peakxmotif)
    result
}

#' Filter a Seurat object
#'
#' This function filters a Seurat object and an output
#' from {\code{\link{getObjectsForModelFit2}}}
#' to ensure that
#' TF binding motifs are kept only when genes encoding corresponding TFs are
#' present in the RNA-seq data and that the same set of genes is included
#' in all objects that are necessary for model fitting.
#'
#' @param object a Seurat object.
#' @param tripod.object a list object returned by {\code{\link{getObjectsForModelFit2}}}
#'
#' @return This function returns a Seurat object.
#' @import Seurat Signac
#' @export
filterSeuratObject2 <- function(object, tripod.object){
    genes <- tripod.object$transcripts.gr$gene_name
    peakxmotif <- tripod.object$peakxmotif
    motifs <- colnames(peakxmotif)
    object@assays$RNA <- subset(object@assays$RNA,
                                features = match(genes, rownames(object@assays$RNA)))
    object@assays$SCT <- subset(object@assays$SCT,
                                features = match(genes, rownames(object@assays$SCT)))
    object@assays$chromvar <- subset(object@assays$chromvar,
                                     features = match(motifs, rownames(object@assays$chromvar)))
    object
}

