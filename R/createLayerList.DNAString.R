#' Create a LayerList object from a DNA sequence
#'
#' Create a LayerList object from a DNAString object
#'
#' @param seq a (\code{"DNAString"} object
#' @param layerNames a character vector of names for the layers. Should be unique 
#' @param nLayers  How many layers to make. Can be specified _instead_ of _layerNames_
#' @param verbose Output extra information
#'
#' @return \code{"LayerList"}
#'
#' @examples
#' library(Biostrings)
#' mySeq <- DNAString("TGACATCGTCTATCGATCG")
#' createLayerList.DNAstring(seq=mySeq, nLayers=1)
#'
#' @import GenomicRanges
#' @import Biostrings
#' 
#'
#' @export
createLayerList.DNAstring <- function(seq, layerNames=NULL, nLayers=length(layerNames), verbose=TRUE)  {
 
  stopifnot(class(seq) == "DNAString")
  stopifnot(!any(layerNames=="LAYER.0"))
  if(is.null(layerNames)  & nLayers > 0)  layerNames <- paste0("LAYER.", 1:nLayers)
  #stopifnot(length(layerNames) == nLayers)
  
  layerSet <- list(LAYER.0 = seq)
  for(thisLayer in layerNames)  {
    layerSet[[thisLayer]] <- IRanges()  
  }

  
  layerList <- list(layerSet=layerSet, history=NULL)  # add some metadata
  return(layerList)
}

# testSeq <- DNAString("TACGTAGCCTGTGATCGTACTAGCGAGTGCAG")
# badSeq <- "ATGTSGDTASGDTGT"
# createLayerList.DNAstring(seq=testSeq)     # allow LayerLists with no layers above the sequence.
# createLayerList.DNAstring(seq=badSeq)
# createLayerList.DNAstring(seq=testSeq, layerNames=c("CpG_island","PRC","H3K27me3"))
# createLayerList.DNAstring(seq=testSeq, nLayers=5)
# createLayerList.DNAstring(seq=testSeq, nLayers=0)

