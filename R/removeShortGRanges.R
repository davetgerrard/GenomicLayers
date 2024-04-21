

removeShortGRanges <- function(x, minSize = 1)  {
 require(GenomicRanges)
  
  stopifnot(class(x) == "GRanges")
  
  #return()
  keepIndex <- width(x)  >= minSize
  return(x[keepIndex])
}

#thisLayer <- "H3K4me1"

#hist(width(removeShortGRanges(x=nucLayerSet_2$layerSet[[thisLayer]], minSize = 147)), breaks=50)
#range(width(removeShortGRanges(x=nucLayerSet_2$layerSet[[thisLayer]], minSize = 147)))
