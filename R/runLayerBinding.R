#' Run the layer binding simulation
#'
#' Description text here
#'
#' @param layerList a \code{"Layerlist"} object containing a layerSet and other meta-data
#' @param factorSet a \code{"list"} of \code{"bindingFactor"} objects
#' @param iterations how many changes to make  NEED TO RENAME THIS
#' @param bindingFactorFreqs the relative proportions of each \code{"bindingFactor"} in \code{"factorSet"}
#' @param watch.function have this function execute during each iteration e.g. print something
#' @param collect.stats collect a table of stats each iteration
#' @param target.layer NOT IMPLEMENTED
#' @param verbose give more output
#'
#' @return A list containing a \code{"LayerSet"} and meta-data, a.k.a. a  \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding.BSgenome}} \code{\link{matchBindingFactor}} \code{\link{modifyLayerByBindingFactor.Views}}
#' 
#' @examples
#' require(Biostrings)     # hopefully this should be available?!
#' 
#' basicLayerList <- createLayerList.DNAstring(seq=DNAString("ACGTTGAAGT"), layerNames=c( "LAYER.1"))
#' 
#' simpleBindingFactor <- createBindingFactor.layer_region(name="simpleBF")
#' 
#' runLayerBinding(basicLayerList, factorSet= list(simpleBF = simpleBindingFactor))  # not working.
#' 
#' # The above find a match to a single DNA pattern in the sequence and is not different from using matchPattern()
#' # Here, we model two layers on the sequence and a second binding factor that only responds to changes 
#' # that have already been made to LAYER.1
#' 
#' twoLayerList <- createLayerList.DNAstring(seq=DNAString("ACGTTGCCATAAACGTTGCCATAAGTGT"), 
#'                layerNames=c( "LAYER.1", "LAYER.2"))
#' 
#' bfList <- list( DNA_A = createBindingFactor.DNA_motif(name="DNA_A",patternString = "CAT" ) ,
#'                 LAYER_1_2 = createBindingFactor.layer_region(name="LAYER_1_2", 
#'                                                              profile.layers =c("LAYER.1"), profile.marks = c(1), 
#'                                                              mod.layers = c("LAYER.2"), mod.marks = c(1)) )
#' 
#' 
#' 
#' runLayerBinding(twoLayerList, factorSet= bfList)  
#' runLayerBinding(twoLayerList, factorSet= rev(bfList))     # should not produce hits on LAYER.2
#' 
#'
#' @export
runLayerBinding <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE, ...)  {
  if(verbose) print(paste(Sys.time(), "runLayerBinding pos 1", sep=" "))
  #bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  bindingOrder <- names(factorSet)  # JUST USE EACH FACTOR ONCE, IN ORDER GIVEN
  newLayerList <- layerList
  seqRange <- c(start(layerList$layerSet[['LAYER.0']])[1], end(layerList$layerSet[['LAYER.0']])[1])
  if(collect.stats) {
    stats.table <- data.frame()
  } else {
    stats.table <- NULL
  }
  max.hits <- ceiling(iterations/length(factorSet))  # TODO could tailor this to be different for each factor.
  for(thisBF in bindingOrder)  {
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast thisBF =", thisBF, factorSet[[thisBF]]$profile$LAYER.0$pattern, sep=" "))
    theseHits <- matchBindingFactor(newLayerList$layerSet, factorSet[[thisBF]], verbose=verbose, ...)
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits =", length(theseHits), sep=" "))
	#print(length(theseHits))
    if(length(theseHits) < 1) {
      if(verbose) print(paste(thisBF, "0 hits, skipping"))
      next ;
    }

    # how many of the potential hits to mark?
    # iterations/n.factors (rounded up).
    # need to check if hits object is a single hit spanning whole chrom (e.g. for DNA_region matching 'N' and nothing else)
    # TODO: if so, need to make pseudo-hits that can be sampled.
    if(identical(theseHits, IRanges(start=seqRange[1], end=seqRange[2]))) {
      if(verbose) print("Hits match whole chromosome")
      hit.width <- factorSet[[thisBF]]$profile$LAYER.0$length
      starts <- sample(1:(seqRange[2]-hit.width), min(max.hits,seqRange[2]-hit.width), replace=FALSE)
      hits.sample <- IRanges(starts, starts+hit.width)
    } else {
      hits.sample <- theseHits[sample(1:length(theseHits) ,min(length(theseHits),max.hits))]   # now multiple
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits.used =", length(hits.sample), sep=" "))
    }
    if(verbose) {print(paste(thisBF, length(theseHits), length(hits.sample)))
    } else { cat( ".")}
	#thisHitPosition <- start(hits.sample) + floor(width(hits.sample)/2)

    newLayerList$layerSet <-  modifyLayerByBindingFactor.Views(newLayerList$layerSet, hits=hits.sample, bindingFactor=factorSet[[thisBF]], verbose=verbose)
    watch.function(x= newLayerList$layerSet, ...)
    if(collect.stats) {
      coverages <- unlist(lapply(newLayerList$layerSet[-1], FUN=function(x) sum(width(x))))
      block.counts <- unlist(lapply(newLayerList$layerSet[-1], length))
      #hit.counts <-
      thisRow <- data.frame(bf=thisBF,hits=length(hits.sample) , target.coverage=sum(width(hits.sample)))
      thisRow[, paste("Coverage.", names(coverages), sep="")] <- coverages
      thisRow[, paste("nBlocks.", names(block.counts), sep="")] <- block.counts
      stats.table <- rbind(stats.table, thisRow)
    }
    #print(as.numeric(letterFrequency(newLayerSet$LAYER.1, letters= "1")))
  }
  newLayerList$history <- stats.table
  #print(letterFrequency(newLayerSet$LAYER.1, letters= "1"))   # how many of layer.1 were set to 1
  if(verbose) {print(paste(Sys.time(), "runLayerBinding.fast pos 2", sep=" ")) } else {print("")}
  return(newLayerList)
}