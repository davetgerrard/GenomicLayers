# HUMAN TEST RUN -----------------------
library(Biostrings)
library(GenomicLayers)
library(BSgenome.Hsapiens.UCSC.hg38)


#genome <- BSgenome.Hsapiens.UCSC.hg38 

hsLayerSet <- createLayerSet.BSgenome(genome=BSgenome.Hsapiens.UCSC.hg38, n.layers = 5, verbose=TRUE)

bf.motif <- createBindingFactor.DNA_motif("tf3", patternString="TGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                          mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = hsLayerSet, bindingFactor = bf.motif)

# non-sequence bindingfactor with offset 
# want a binding factor that creates mods either up or downstream of binding site. New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}
bf.spread4 <- createBindingFactor.layer_region("offset", patternLength = 5, profile.layers="LAYER.4", profile.marks=1, 
                                               mod.layers="LAYER.4", mod.marks = 1,  stateWidth = 7, offset = 150, offset.method=upDownFunc)

results5 <- modifyLayerByBindingFactor.BSgenome(layerSet=hsLayerSet, hits=results3, bindingFactor=bf.spread4)
# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list( tf3=bf.motif, offset=bf.spread4)
n.runs <- 10

resultsList <- list()
timeList <- list()
for(use.cache in c(FALSE, TRUE)) {
  cacheLayers <- NULL
  if(use.cache)  cacheLayers <- "LAYER.0"
  finalLayer <- hsLayerSet
  i <- 1
  statsTrace <- data.frame()
  coverTrace <- data.frame()
  thisTime <- system.time(
    while(i <= n.runs) {
      print(i)
      finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=testFS, verbose=TRUE, cache.layers = cacheLayers, iterations=30, collect.stats = TRUE)
      thisRow <- cbind(i, finalLayer$history)
      statsTrace <- rbind(statsTrace , thisRow)
      covL <- lapply(coverage(finalLayer$layerSet$LAYER.4), sum)   # get sum of marked regions on specific chromosome
      thisRow <- cbind(i, as.data.frame(covL))
      coverTrace <- rbind(coverTrace, thisRow)
      if( i %% 100 == 0)  {   # export tables every 100th iteration
        write.table(statsTrace, file=paste0("simpleSpread.", n.runs, ".cache.", use.cache, ".statsTrace.tab"), sep="\t", quote=F, row.names=F)
        write.table(coverTrace, file=paste0("simpleSpread.", n.runs, ".cache.", use.cache, ".coverTrace.tab"), sep="\t", quote=F, row.names=F)
      }
      i <- i +1
    }
  )
  timeList <- c(timeList, list(thisTime))
  resultsList <- c(resultsList, list(statsTrace))
}
timeList
