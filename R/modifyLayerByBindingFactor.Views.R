#' Modify a layerSet object
#'
#' Modify a layerSet object according to a bindingFactor and a pre-calculated set of hits.
#'
#' @param layerSet method to do something to (\code{"hsv"} or \code{"cluster"})
#' @param hits a \code{"Views"} object representing the places the bindingFactor can bind.
#' @param bindingFactor  The \code{"bindingFactor"}
#' @param verbose you get the idea
#'
#' @return \code{"layerSet"}
#'
#' @examples
#' x <- 1   # great!
#'
#' @import GenomicRanges
#' @import Biostrings
#' 
#'
#' @export
modifyLayerByBindingFactor.Views <- function(layerSet, hits, bindingFactor, verbose=FALSE) {
  require(Biostrings)
  newLayerSet <- layerSet
  for(thisLayer in names(bindingFactor$mods))  {


    #seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
    thisState <- bindingFactor$mods[[thisLayer]]$state

    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    hits <- resize(hits, width=stateWidth, fix="center")    # adjust the width
    # if there is an offset.method, use it. Otherwise just use the offset. 
    # Written so that function is applied for each one of hits independently.
    if(!is.null(bindingFactor$mods[[thisLayer]]$offset.method )) {
      thisOffset <- replicate(length(hits),bindingFactor$mods[[thisLayer]]$offset.method(x=bindingFactor$mods[[thisLayer]]$offset))
    } else {
      thisOffset <- bindingFactor$mods[[thisLayer]]$offset
    }
    hits <- shift(hits, shift=thisOffset)                 # move up- or downstream

    # restrict hits to range


    newLayerSet[[thisLayer]] <- switch(thisState,
                                       "1" = reduce(union(newLayerSet[[thisLayer]], hits)),
                                       "0" = setdiff(newLayerSet[[thisLayer]], hits),
                                       stop(paste("unknown state", thisState)))


  }
  return(newLayerSet)
}
