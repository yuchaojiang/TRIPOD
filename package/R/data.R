#' Transcripts
#'
#' A GRanges object containing 14476 protein_coding genes detected
#' in single-cell multiome data of mouse embryonic brain at day 18.
#'
#' @format a GRanges object with 14476 elements:
#' \describe{
#'   \item{gene_name}{gene symbols}
#'   \item{gene_biotype}{gene types}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets}
"transcripts.gr"

#' ATAC peaks
#'
#' A GRanges object containing 139083 ATAC peaks identified
#' in single-cell multiome data of mouse embryonic brain at day 18.
#'
#' @format a GRanges object with 139083 elements
#' @source \url{https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets}
"peaks.gr"

#' TFs and their binding motifs
#'
#' TF binding motif data for 351 TFs based on the JASPAR2020 database.
#'
#' @format a matrix object with 351 rows and 2 columns:
#' \describe{
#'   \item{motif}{JASPAR motif IDs}
#'   \item{TF}{gene types}
#' }
"motifxTF"

#' The presence and absence of TF binding motifs in ATAC peaks
#'
#' A binary sparse matrix indicating whether TF binding motifs are present
#' in ATAC peaks identified in single-cell multiome data of mouse embryonic brain
#' at day 18.
#'
#' @format a dgCMatrix matrix with 139083 rows and 351 columns:
#' \describe{
#'   \item{row}{ATAC peak chromosomal coordinates}
#'   \item{column}{JASPAR motif IDs}
#' }
"peakxmotif"

#' RNA expression in metacells
#'
#' Normalized RNA expression levels of 14476 genes in 84 metacells identified
#' in single-cell multiome data of mouse embryonic brain at day 18.
#'
#' @format a matrix with 84 rows and 14476 columns:
#' \describe{
#'   \item{row}{metacell indices}
#'   \item{column}{gene symbols}
#' }
"metacell.rna"

#' Chromatin accessibility in metacells
#'
#' Normalized Chromatin accessibility of 139083 ATAC peaks in 84 metacells.
#' identified in single-cell multiome data of mouse embryonic brain at day 18.
#'
#' @format a matrix object with 84 rows and 139083 columns:
#' \describe{
#'   \item{row}{metacell indices}
#'   \item{column}{ATAC peak chromosomal coordinates}
#' }
"metacell.peak"

#' Cell types of metacells
#'
#' Assigned cell types of metacells identified in single-cell multiome data
#' of mouse embryonic brain at day 18.
#'
#' @format a character vector with 84 elements
"metacell.celltype"

#' Colors assigned to metacells
#'
#' Cell type colors of metacells identified in single-cell multiome data
#' of mouse embryonic brain at day 18.
#'
#' @format a character vector with 84 elements
"metacell.celltype.col"

#' Cell-type labels of single cells obtained by SAVERCAT
#'
#' Cell-type labels of single-cells multiome data of mouse embryonic brain
#' at day 18 obtained using SAVERCAT.
#'
#' @format a character vector with 4881 elements
"label.savercat"

#' Cell-type labels of single cells obtained by Seurat
#'
#' Cell-type labels of single-cells multiome data of mouse embryonic brain
#' at day 18 obtained using Seurat.
#'
#' @format a character vector with 4881 elements
"label.seurat"

#' Transcripts in PBMC
#'
#' A GRanges object containing 14508 protein_coding genes detected
#' in single-cell multiome data of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a GRanges object with 14508 elements:
#' \describe{
#'   \item{gene_name}{gene symbols}
#'   \item{gene_biotype}{gene types}
#' }
#' @source \url{https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets}
"pbmc.transcripts.gr"

#' ATAC peaks in PBMC
#'
#' A GRanges object containing 106056 ATAC peaks identified
#' in single-cell multiome data of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a GRanges object with 106056 elements
#' @source \url{https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets}
"pbmc.peaks.gr"

#' TFs and their binding motifs for analyzing the PBMC data
#'
#' TF binding motif data for 342 TFs based on the JASPAR2020 database.
#'
#' @format a matrix object with 342 rows and 2 columns:
#' \describe{
#'   \item{motif}{JASPAR motif IDs}
#'   \item{TF}{gene types}
#' }
"pbmc.motifxTF"

#' The presence and absence of TF binding motifs in ATAC peaks in PBMC
#'
#' A binary sparse matrix indicating whether TF binding motifs are present
#' in ATAC peaks identified in single-cell multiome data of
#' human peripheral blood mononuclear cells (PBMC).
#'
#' @format a dgCMatrix matrix with 106056 rows and 342 columns:
#' \describe{
#'   \item{row}{ATAC peak chromosomal coordinates}
#'   \item{column}{JASPAR motif IDs}
#' }
"pbmc.peakxmotif"

#' RNA expression in metacells in PBMC
#'
#' Normalized RNA expression levels of 14508 genes in 80 metacells identified
#' in single-cell multiome data of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a matrix with 80 rows and 14508 columns:
#' \describe{
#'   \item{row}{metacell indices}
#'   \item{column}{gene symbols}
#' }
"pbmc.metacell.rna"

#' Chromatin accessibility in metacells in PBMC
#'
#' Normalized Chromatin accessibility of 106056 ATAC peaks in 80 metacells.
#' identified in single-cell multiome data of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a matrix object with 80 rows and 106056 columns:
#' \describe{
#'   \item{row}{metacell indices}
#'   \item{column}{ATAC peak chromosomal coordinates}
#' }
"pbmc.metacell.peak"

#' Cell types of metacells in PBMC
#'
#' Assigned cell types of metacells identified in single-cell multiome data
#' of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a character vector with 80 elements
"pbmc.metacell.celltype"

#' Colors assigned to metacells in PBMC
#'
#' Cell type colors of metacells identified in single-cell multiome data
#' of human peripheral blood mononuclear cells (PBMC).
#'
#' @format a character vector with 80 elements
"pbmc.metacell.celltype.col"
