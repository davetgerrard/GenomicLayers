

removeGRangesBySize <- function(x, maxSize=-1, minSize = 1, verbose=FALSE)  {
 require(GenomicRanges)
  
  stopifnot(class(x) == "GRanges")
  
  aboveMin <- width(x)  >= minSize
  if(maxSize == -1)  {
    withinMax <- rep(TRUE, length(x))
  } else {
    withinMax <- width(x)  <= maxSize
  }
  #return()
  keepIndex <- withinMax & aboveMin
  if(verbose)  {
    print(paste("Removing features:" , sum(!keepIndex) , "/", length(x)))
  }
  
  return(x[keepIndex])
}

#thisLayer <- "H3K4me1"

#hist(width(removeShortGRanges(x=nucLayerSet_2$layerSet[[thisLayer]], minSize = 147)), breaks=50)
#range(width(removeShortGRanges(x=nucLayerSet_2$layerSet[[thisLayer]], minSize = 147)))
