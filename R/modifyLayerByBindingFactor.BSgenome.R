# Layers are now Views or Iranges objects. marks (1) are contiguous ranges, absence of marks (0) are gaps between.
# hits is now a Views object
# TODO, edit to alter layerset IN-PLACE i.e. modifyLayerByBindingFactor.Views(layerSet, position.vec,bindingFactor)
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
#' @seealso \code{\link{createLayerSet.BSgenome}} \code{\link{createBindingFactor.DNA_motif}} \code{\link{runLayerBinding.BSgenome}}  \code{\link{matchBindingFactor}}  \code{\link{matchBindingFactor.BSgenome}}
#'
#' @examples
#' require(Biostrings)
#' require(BSgenome.Scerevisiae.UCSC.sacCer3)
#' 
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
#' 
#' 
#' scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
#' 
#' testFactor3 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
#'                                              mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))
#' 
#' listOfHits <- matchBindingFactor.BSgenome(layerSet=scLayerSet, bindingFactor=testFactor3)
#' 
#' moddedLayerSet <- modifyLayerByBindingFactor.BSgenome(layerSet= scLayerSet, hits=listOfHits, bindingFactor=testFactor3)
#'
#' as.data.frame(lapply(moddedLayerSet$layerSet, length))  # should be hits on LAYER.0 and LAYER.4. Not LAYER.2 as these are set to 0.
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
