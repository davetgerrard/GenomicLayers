# could set some signature methods for this to use createLayerSet for both BSgenomes and arbitrary string collections.
# this itself could be a class...
# need to tidy up naming of layerSet and layerList.
#' Create a genome-wide layerList object
#'
#' Convenience wrapper to quickly create a layerList from a BSgenome object.
#'
#' @param genome a BSgenome objects
#' @param n.layers    1
#' @param layer.names  paste("LAYER",1:n.layers, sep=".")
#' @param verbose    Give more output [FALSE]
#'
#' @return \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding.BSgenome}}
#'
#' @examples
#' require(Biostrings)
#' require(BSgenome.Scerevisiae.UCSC.sacCer3)  # a relatively small genome
#'
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3 
#' scLayerSet <- createLayerSet.BSgenome(genome, n.layers=1)
#'
#' @export
createLayerSet.BSgenome <- function(genome, n.layers=1, layer.names=paste("LAYER",1:n.layers, sep="."), verbose=FALSE) {
  stopifnot(class(genome) == "BSgenome")
  
  #Need to use Genome wide ranges objects to do intersections.
  
  if(verbose) print(paste("creating layerSet with", n.layers, "layers on", bsgenomeName(genome)))
  layerSet <- list(LAYER.0 = genome)
  for(i in 1:n.layers) {
    layerSet[[layer.names[i]]] <- GRanges(seqinfo=seqinfo(genome))
  }
  
  layerList <- list(layerSet=layerSet, history=NULL)
  
  return(layerList)
  
}