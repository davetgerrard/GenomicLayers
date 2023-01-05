#' Generate a cache of binding factor within a LayerList
#'
#' Generate and attach a cache of binding factor hit positions for named binding factors on certain levels.
#' This is an internal function called from runLayerBinding.BSgenome 
#' In turn it uses \code{\link{matchBindingFactor.BSgenome}} to generate hits.
#'
#' @param layerList a \code{"Layerlist"} object containing a layerSet based on a BSgenome and other meta-data
#' @param factorSet a \code{"list"} of \code{"bindingFactor"} objects
#' @param cache.layers which named layers to cache hits on (default NULL)
#' @param verbose give more output
#'
#' @return \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding.BSgenome}} \code{\link{matchBindingFactor.BSgenome}}
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
#' generateHitsCache(layerList=scLayerSet, factorSet=list(testFactor3))
#' 
#' @export
generateHitsCache <- function(layerList, factorSet, cache.layers=NULL, verbose=FALSE)  {
  
  newlayerList <- layerList
  if(is.null(newlayerList$cache))  newlayerList$cache <- list()
  for(thisFactor in names(factorSet)) {
    if(is.null(newlayerList$cache[[thisFactor]]))  newlayerList$cache[[thisFactor]] <- list()
    for(thisLayer in cache.layers) {
      newlayerList$cache[[thisFactor]][[thisLayer]]   <- matchBindingFactor.BSgenome(layerSet = newlayerList, bindingFactor = factorSet[[thisFactor]], match.layers = thisLayer, verbose=verbose)
    }
  }
  return(newlayerList)
  
}