

# sandbox to trial caching of genome wide hits and re-use them.

#DESIGN:
#      store hits object within factorSet, within bindingFactor, within layerSet, or separately?
#     Want it to be automatic and behind the scenes, each sim may begin by generating cache (though pre-loading would be good)
#     Multiple layers or just on sequence?
#           cache.layers=NULL    set to "LAYER.0" , could additionally have use.cache=FALSE/TRUE but this is somewhat redundant
#    Meta-data of cache: genome or object run against? Does layerSet currently retain genome?
#     This works for meta-data of BSgenomes:-    
#         scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
#         bsgenomeName(scLayerSet$layerSet[["LAYER.0"]])
#   Some binding factors will be cached/cachable, others will not.
# Might be tidier to link cache to the layerSet. A set of named hits based on names of bindingFactors
#  Better check for replicated binding factors and warn. 




#IMPACTS :  runLayerBinding.BSgenome, matchBindingFactor.BSgenome
#             ?  createBindingFactor(s), matchBindingFactor(?)

library(Biostrings)
library(GenomicLayers)
library(BSgenome.Scerevisiae.UCSC.sacCer3)




# test matchBindingFactor.BSgenome()
pattern.long <- "TGAAACR"


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)



bf.motif <- createBindingFactor.DNA_motif("tf3", patternString="TGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                          mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf.motif)

resultsNoCache <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf.motif, cache.layers = "LAYER.0", verbose=TRUE)

resultsSingleLayer <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf.motif, match.layers = "LAYER.0", verbose=TRUE)
# 


# see R/generateHitsCache.R
# generateHitsCache <- function(layerList, factorSet, cache.layers=NULL, verbose=FALSE)  {
#   
#   newlayerList <- layerList
#   if(is.null(newlayerList$cache))  newlayerList$cache <- list()
#   for(thisFactor in names(factorSet)) {
#     if(is.null(newlayerList$cache[[thisFactor]]))  newlayerList$cache[[thisFactor]] <- list()
#     for(thisLayer in cache.layers) {
#       newlayerList$cache[[thisFactor]][[thisLayer]]   <- matchBindingFactor.BSgenome(layerSet = newlayerList, bindingFactor = factorSet[[thisFactor]], match.layers = thisLayer, verbose=verbose)
#     }
#   }
#   return(newlayerList)
#   
# }

cachedLayerList <- generateHitsCache(scLayerSet, factorSet = list(tf3=bf.motif), cache.layers = "LAYER.0", verbose=TRUE)
 
# now try matchBindingFactor but use the cache for the selected layer.
system.time(
  resultsNoCache <- matchBindingFactor.BSgenome(layerSet = cachedLayerList, bindingFactor = bf.motif, cache.layers = NULL, verbose=TRUE)
)
system.time(
  resultsWithCache <- matchBindingFactor.BSgenome(layerSet = cachedLayerList, bindingFactor = bf.motif, cache.layers = "LAYER.0", verbose=TRUE)
)

# try with multiple bindingFactors
cachedLayerList <- generateHitsCache(scLayerSet, factorSet = list(tf2=bf.motif, tf4=bf.motif), cache.layers = "LAYER.0", verbose=TRUE)
resultsWithWrongCache <- matchBindingFactor.BSgenome(layerSet = cachedLayerList, bindingFactor = bf.motif, cache.layers = "LAYER.0", verbose=TRUE)


# now try runLayerBinding.BSgenome, using a cache (starting with none)

testFS <- list( bf.motif=bf.motif, tf3=bf.motif)   # these confuse names!

# to run through once with two factors, caching is slower (because the caching is currently a separate step.
system.time(
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, cache.layers= "LAYER.0", verbose=TRUE)
)
system.time(
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, cache.layers= NULL, verbose=TRUE)
)

# try iterations using a cache


n.runs <- 10


finalLayer1 <- scLayerSet
i <- 1
system.time(
while(i <= n.runs) {
  finalLayer1 <- runLayerBinding.BSgenome(layerList=finalLayer1, factorSet=testFS, verbose=TRUE, cache.layers= "LAYER.0", iterations=30, collect.stats = TRUE)
  i <- i +1
}
)

finalLayer2 <- scLayerSet
i <- 1
system.time(
  while(i <= n.runs) {
    finalLayer2 <- runLayerBinding.BSgenome(layerList=finalLayer2, factorSet=testFS, verbose=TRUE, cache.layers= NULL, iterations=30, collect.stats = TRUE)
    i <- i +1
  }
)


stopifnot(FALSE)



# create a layerSet to use as a target (does it have a name?)

# create or load a bindingFactor.

# runLayerBinding









# 
# pattern.core <- "GAAAC"
 pattern.long <- "TGAAACR"  # R = puRine (A,G)
# IUPAC_CODE_MAP

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
# genome
# seqnames(genome)
# summary(genome)
# organism(genome)
# provider(genome) 
# bsgenomeName(genome)  # - could use this for simple comparison 

tf.hits <- vmatchPattern(pattern.long, genome, fixed=F) 


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
bsgenomeName(scLayerSet$layerSet[["LAYER.0"]])

testFactor <- createRandomBindingFactor(name="test.1", layerSet=scLayerSet$layerSet, type="DNA_motif")


results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor)



# need to prescribe a simple binding factor, made a simple function to create one.

testFactor2 <- createBindingFactor.DNA_motif("tf2", patternString="TGGGCTA")  # ACTGGGCTA does not hit all chromosomes

results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

bf.motif <- createBindingFactor.DNA_motif("tf3", patternString="TGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                          mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = bf.motif)




# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.DNA_motif("tf4", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(1,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results4 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor4)


# non-sequence bindingfactor with offset 
# want a binding factor that creates mods either up or downstream of binding site. New parameter 'offset.method' to accept function as argument
upDownFunc <- function(x)  {
  return(sample(c(x, -x), 1))
  
}
bf.spread4 <- createBindingFactor.layer_region("offset", patternLength = 5, profile.layers="LAYER.4", profile.marks=1, 
                                               mod.layers="LAYER.4", mod.marks = 1,  stateWidth = 7, offset = 150, offset.method=upDownFunc)

results5 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=bf.spread4)
# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list( bf.motif=bf.motif, bf.spread4=bf.spread4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)
# 2016-09-05 

# with the above configuration, there are 41 possible sites across the genome, setting iterations=30, restricts the number that are marked, so the number of potential sites reduces.
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)

n.runs <- 5

finalLayer <- scLayerSet
i <- 1
statsTrace <- data.frame()
coverTrace <- data.frame()
while(i <= n.runs) {
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)
  thisRow <- cbind(i, finalLayer$history)
  statsTrace <- rbind(statsTrace , thisRow)
  covL <- lapply(coverage(finalLayer$layerSet$LAYER.4), sum)   # get sum of marked regions on specific chromosome
  thisRow <- cbind(i, as.data.frame(covL))
  coverTrace <- rbind(coverTrace, thisRow)
  if( i %% 100 == 0)  {   # export tables every 100th iteration
    write.table(statsTrace, file=paste0("simpleSpread.", n.runs, ".statsTrace.tab"), sep="\t", quote=F, row.names=F)
    write.table(coverTrace, file=paste0("simpleSpread.", n.runs, ".coverTrace.tab"), sep="\t", quote=F, row.names=F)
  }
  i <- i +1
}

