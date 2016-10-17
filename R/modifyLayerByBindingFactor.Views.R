# Layers are now Views or Iranges objects. marks (1) are contiguous ranges, absence of marks (0) are gaps between.
# hits is now a Views object
# TODO, edit to alter layerset IN-PLACE i.e. modifyLayerByBindingFactor.Views(layerSet, position.vec,bindingFactor)
# function
# description
# Parameters:-
# layerSet,
# hits,
# bindingFactor,
# verbose=FALSE
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
