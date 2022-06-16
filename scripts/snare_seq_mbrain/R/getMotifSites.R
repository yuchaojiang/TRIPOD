# getMotifSites.R

#' Get motif genomic locations
#'
#' @param gene.name a character string
#' @param TF.name a character string
#' @param vicinity.size a integer vector with two elements. The first and second elements represent sizes of upstream and downstream regions around TSS, respectively.
#' @param transcripts a GRanges object containing annotated protein-coding genes
#' @param motifxTF a matrix with two columns. The first and second columns contain motif IDs and TF names, respectively.
#' @param BSgenome a BSgenome object
#' 
#' @import GenomicRanges
#' @import TFBSTools
#' @import JASPAR2020
#' @import Biostrings
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @return a list with two elements
#' @name getMotifSites
#' @rdname getMotifSites 
#' @export
#'
#' @examples
getMotifSites <- function(gene.name,
                          TF.name,
                          vicinity.size,
                          transcripts,
                          motifxTF,
                          BSgenome){
  # get the motif name
  motif.name<- motifxTF[motifxTF[, 2] == TF.name, 1]
  # get pfm
  pfm <- getMatrixByID(JASPAR2020, ID = motif.name)
  # convert pfm to pwm
  pwm <- toPWM(pfm)
  # get chr, start, end, and strand
  transcripts <- transcripts[transcripts$gene_name == gene.name]
  chr <- seqnames(transcripts)
  # get TSS
  TSS.position <- ifelse(strand(transcripts) == "+", start(transcripts), end(transcripts))
  TSS <- GRanges(seqnames = seqnames(transcripts),
                 ranges = IRanges(start = TSS.position, width = 1),
                 strand = strand(transcripts),
                 gene_name = transcripts$gene_name)
  # get the vicinity region
  vicinity.gr <- getRegion(vicinity.size, TSS = TSS)
  origin <- start(vicinity.gr)
  # get the sequence of the vicinity region
  string <- getSeq(BSgenome, vicinity.gr)
  # get motif sites
  sites <- searchSeq(pwm, string, min.score = "80%", strand = "*")
  # convert to a data frame
  sites.df <- writeGFF3(sites)
  sites.gr <- GRanges(seqnames = rep(chr, nrow(sites.df)),
                      ranges = IRanges(start = sites.df$start + origin, # GFF3 is 1-based
                                       end = sites.df$end + origin), 
                      strand = sites.df$strand, 
                      TF_name = rep(TF.name, nrow(sites.df)))
  results <- list(vicinity = vicinity.gr, sites = sites.gr)
  return(results)
}
