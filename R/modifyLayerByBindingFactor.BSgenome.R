# Layers are now Views or Iranges objects. marks (1) are contiguous ranges, absence of marks (0) are gaps between.
# hits is now a Views object
# TODO, edit to alter layerset IN-PLACE i.e. modifyLayerByBindingFactor.Views(layerSet, position.vec,bindingFactor)
# function
# Modify all hits genome-wide on a layerSet built upon a BSgenome object.
# Parameters:-
# layerSet,
# hits,
# bindingFactor,
# verbose=FALSE
#' modifyLayerByBindingFactor.BSgenome
#'
#'  Modify all hits genome-wide on a layerSet built upon a BSgenome object.
#'
#' @param layerSet A layerset object
#' @param hits  The set of hits (GRanges) determined by \code{\link{matchBindingFactor}} or \code{\link{matchBindingFactor.BSgenome}}
#' @param bindingFactor   An individual bindingFactor object
#' @param verbose  Give more output
#'
#' @return \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding.BSgenome}}  \code{\link{matchBindingFactor}}  \code{\link{matchBindingFactor.BSgenome}}
#'
#' @examples
#' x <- 1   # great!
#'
#' @export
modifyLayerByBindingFactor.BSgenome <- function(layerSet, hits, bindingFactor, verbose=FALSE) {
  require(Biostrings)
  require(GenomicRanges)
  newLayerSet <- layerSet
  for(thisLayer in names(bindingFactor$mods))  {


    #seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
    thisState <- bindingFactor$mods[[thisLayer]]$state

    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    hits <- resize(hits, width=stateWidth, fix="center")    # adjust the width
    thisOffset <- bindingFactor$mods[[thisLayer]]$offset
    hits <- shift(hits, shift=thisOffset)                 # move up- or downstream

    # restrict hits to range


    newLayerSet$layerSet[[thisLayer]] <- switch(as.character(thisState),
                                       "1" = reduce(union(newLayerSet$layerSet[[thisLayer]], hits), ignore.strand=TRUE),
                                       "0" = setdiff(newLayerSet$layerSet[[thisLayer]], hits, ignore.strand=TRUE),
                                       stop(paste("unknown state", thisState)))


  }
  return(newLayerSet)
}
