#' Generate a cache of binding factor within a LayerList
#'
#' Generate and attach a cache of binding factor hit positions for named binding factors on certain levels.
#'
#' @param layerList a \code{"Layerlist"} object containing a layerSet and other meta-data
#' @param factorSet a \code{"list"} of \code{"bindingFactor"} objects
#' @param cache.layers which named layers to cache hits on (default NULL)
#' @param verbose give more output
#'
#' @return \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding}}
#'
#' @examples
#' x <- 1   # great!
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