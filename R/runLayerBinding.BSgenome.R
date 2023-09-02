#' Run the layer binding simulation on a BSgenome object
#'
#' Simulate layer binding by the binding factor on the layerList
#'
#' @param layerList a \code{"Layerlist"} object containing a layerSet and other meta-data
#' @param factorSet a \code{"list"} of \code{"bindingFactor"} objects
#' @param iterations deprecated - do not use
#' @param bf.abundances 10000  the quantities of each \code{"bindingFactor"} in \code{"factorSet"}
#' @param watch.function have this function execute during each iteration e.g. print something
#' @param collect.stats FALSE collect a table of stats each iteration
#' @param keep.stats TRUE  whether to retain data contained in $history.
#' @param target.layer NOT IMPLEMENTED
#' @param cache.layers store hits for these layers if they are immutable (e.g. LAYER.0 sequence) default= "LAYER.0". Set to NULL to prevent caching
#' @param verbose give more output
#'
#' @return \code{"LayerList"}
#'
#' @seealso \code{\link{runLayerBinding}} \code{\link{matchBindingFactor.BSgenome}} \code{\link{modifyLayerByBindingFactor.BSgenome}} \code{\link{generateHitsCache}} \code{\link{createLayerSet.BSgenome}} \code{\link{createBindingFactor.DNA_motif}}
#'
#' @examples
#' # test to show whole genome layer binding 
#' # three factors are created, the third can only hit after the first has modified.
#' 
#' require(Biostrings)
#' require(BSgenome.Scerevisiae.UCSC.sacCer3)
#' 
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
#' 
#' 
#' scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
#' 
#' testFactor2 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA")
#' 
#' testFactor3 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
#'                                              mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))
#' 
#' # check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
#' testFactor4 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
#'                                              mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))
#' 
#'  
#' # now can match things genome wide. Need to run layerBinding and modification.
#' 
#' # need to have a factorSet, a list of bindingFactors
#' 
#' testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)
#' 
#'  
#' # with the above configuration, there are 41 possible sites across the genome, setting bf.abundances=30, restricts the number that are marked, so the number of potential sites reduces.
#' modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, bf.abundances=30)
#'
#' @export
runLayerBinding.BSgenome <- function(layerList, factorSet, iterations=1, bf.abundances=rep(10000, length(factorSet)), 
                                     watch.function=function(x){}, collect.stats=FALSE, keep.stats=TRUE, 
                                     target.layer=2, cache.layers="LAYER.0", verbose=FALSE, ...)  {
  if(verbose) print(paste(Sys.time(), "runLayerBinding pos 1", sep=" "))
  #bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  
  bindingOrder <- names(factorSet)  # JUST USE EACH FACTOR ONCE, IN ORDER GIVEN
  if(length(bf.abundances)==1) {   # if a single value is given, assign it to all binding factors. If only one bf, does nothing.
    bf.abundances <- rep(bf.abundances, length(factorSet))
  }
  names(bf.abundances) <- bindingOrder
  newLayerList <- layerList
  #seqRange <- c(start(layerList$layerSet[['LAYER.0']])[1], end(layerList$layerSet[['LAYER.0']])[1])
  genome <- layerList$layerSet[[1]]
  genome.sl <- seqlengths(genome)
  genome.starts <- rep(1, length(genome))
  names(genome.starts) <- seqnames(genome)
  genome.ends <- genome.sl
  
  if(collect.stats) {
    stats.table <- data.frame()
  } else {
    stats.table <- NULL
  }
  #max.hits <- ceiling(iterations/length(factorSet))  # TODO could tailor this to be different for each factor.
  
  if(!is.null(cache.layers)  & is.null(newLayerList$cache)) {  # generate a new cache of hits for each factor on the specified layers.
    if(verbose) print("Generating hits cache")
    newLayerList <- generateHitsCache(newLayerList, factorSet=factorSet, cache.layers=cache.layers, verbose=verbose)
  }
  for(thisBF in bindingOrder)  {
    max.hits <- bf.abundances[thisBF]
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast thisBF =", thisBF, factorSet[[thisBF]]$profile$LAYER.0$pattern, sep=" "))
    theseHits <- matchBindingFactor.BSgenome(newLayerList, factorSet[[thisBF]], cache.layers=cache.layers, verbose=verbose, ...)
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
    # 2016-09-05  this section commented out while testing DNA_motif hits against whole genome. May need to rework for other motif types.
    #if(identical(theseHits, IRanges(start=seqRange[1], end=seqRange[2]))) {
    #  if(verbose) print("Hits match whole chromosome")
    #  hit.width <- factorSet[[thisBF]]$profile$LAYER.0$length
    #  starts <- sample(1:(seqRange[2]-hit.width), min(max.hits,seqRange[2]-hit.width), replace=FALSE)
    #  hits.sample <- IRanges(starts, starts+hit.width)
    #} else {
      hits.sample <- theseHits[sample(1:length(theseHits) ,min(length(theseHits),max.hits))]   # now multiple
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits.used =", length(hits.sample), sep=" "))
    #}
    if(verbose) {print(paste(thisBF, length(theseHits), length(hits.sample)))
    } else { cat( ".")}
	#thisHitPosition <- start(hits.sample) + floor(width(hits.sample)/2)

    newLayerList <-  modifyLayerByBindingFactor.BSgenome(newLayerList, hits=hits.sample, bindingFactor=factorSet[[thisBF]], verbose=verbose)
    watch.function(x= newLayerList$layerSet, ...)
    stats.table <-  newLayerList$history  # copy the existing stats table if there is one.  
    if(collect.stats) {
      coverages <- unlist(lapply(newLayerList$layerSet[-1], FUN=function(x) sum(width(x))))
      block.counts <- unlist(lapply(newLayerList$layerSet[-1], length))
      #hit.counts <-
      thisRow <- data.frame(time=Sys.time(),bf=thisBF,hits=length(hits.sample) , target.coverage=sum(width(hits.sample)))
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