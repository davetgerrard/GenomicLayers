#' Subset a BSgenome object
#' 
#' Create a BSgenome object with some sequences (chromosomes) removed. 
#' This is useful when wanting to work without mitochondria or unplaced chromosomes.
#' This function was posted by Hervé Pagès on the
#'  [bioconductor support pages around 2016](https://support.bioconductor.org/p/83588/) 
#'  but I cannot find a version within his BSgenome package.
#' 
#' @param genome  a BSgenome object
#' @param seqnames  names of chromosomes to keep
#'
#' @return \code{"genome"}
#' 
#'  
#'
#' @import BSgenome
#' 
#' @examples
#' 
#' library(BSgenome.Scerevisiae.UCSC.sacCer3)
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' seqinfo(genome)
#' seqnames(genome)
#' sequences_to_keep <- setdiff(seqnames(genome), "chrM")   # nuclear chroms only
#' genomeNuc <- keepBSgenomeSequences(genome, sequences_to_keep)
#' seqinfo(genomeNuc)
#' 
#' @export

# hack from https://support.bioconductor.org/p/83588/  to keep only some chromosomes within a BSgenome object.
keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

