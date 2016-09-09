# runLayerBinding.BSgenome
# as runLayerBinding but on whole genome
# Parameters:-
# layerList,
# factorSet,
# iterations=1,
# bindingFactorFreqs=rep(1, length(factorSet)),
# watch.function=function(x){},
# collect.stats=FALSE,
# target.layer=2,
# verbose=FALSE
runLayerBinding.BSgenome <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE, ...)  {
  if(verbose) print(paste(Sys.time(), "runLayerBinding pos 1", sep=" "))
  #bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  
  bindingOrder <- names(factorSet)  # JUST USE EACH FACTOR ONCE, IN ORDER GIVEN
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
  max.hits <- ceiling(iterations/length(factorSet))  # TODO could tailor this to be different for each factor.
  for(thisBF in bindingOrder)  {
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast thisBF =", thisBF, factorSet[[thisBF]]$profile$LAYER.0$pattern, sep=" "))
    theseHits <- matchBindingFactor.BSgenome(newLayerList, factorSet[[thisBF]], verbose=verbose)
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