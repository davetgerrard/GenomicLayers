

library(Biostrings)
library(GenomicLayers)


pattern.core <- "GAAAC"
pattern.long <- "TGAAACR"  # R = puRine (A,G)
IUPAC_CODE_MAP

library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
genome
seqnames(genome)
summary(genome)
organism(genome)
provider(genome) 


tf.hits <- vmatchPattern(pattern.long, genome, fixed=F) 


scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)


testFactor <- createRandomBindingFactor(name="test.1", layerSet=scLayerSet$layerSet, type="DNA_motif")


results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor)



# need to prescribe a simple binding factor, made a simple function to create one.

testFactor2 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")

results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

testFactor3 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor3)


# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(1,0), 
                                                            mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor4)


# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)
# 2016-09-05 

# with the above configuration, there are 41 possible sites across the genome, setting iterations=30, restricts the number that are marked, so the number of potential sites reduces.
modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)

finalLayer <- scLayerSet
i <- 1
statsTrace <- data.frame()
coverTrace <- data.frame()
while(i <= 5) {
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=testFS, verbose=TRUE, iterations=30, collect.stats = TRUE)
  thisRow <- cbind(i, finalLayer$history)
  statsTrace <- rbind(statsTrace , thisRow)
  covL <- lapply(coverage(finalLayer$layerSet$LAYER.4), sum)   # get sum of marked regions on specific chromosome
  coverTrace <- rbind(coverTrace, cbind(i, as.data.frame(covL)))
  i <- i +1
}

# TODO convert the above into a test.


# will this run on human genome?  No, need much more memory to do human genome pattern matching.

library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19

matchPattern("TTTCCCTAATC", genome, fixed=F)

tf.hits <- vmatchPattern("TTTCCCTAATC", genome, fixed=F) 
library(BSgenome.Scerevisiae.UCSC.sacCer3)
genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
genome
seqnames(genome)
summary(genome)
organism(genome)
provider(genome) 


tf.hits <- vmatchPattern(pattern.long, genome, fixed=F) 


scLayerSet <- createLayerSet.BSgenome(genome=BSgenome.Scerevisiae.UCSC.sacCer3, n.layers = 5, verbose=TRUE)

#genome <- BSgenome.Hsapiens.UCSC.hg19
genome


matchBindingFactor.BSgenome(scLayerSet, bf.CpGisland)


  "(CG.{0,20}){9}CG"

  gregexpr("(CG.{0,20}){9}CG",  scLayerSet$layerSet[["LAYER.0"]])
  gregexpr("(CG.{0,20}){9}CG",  scLayerSet$layerSet[["LAYER.0"]][[1]])
  bsapply(scLayerSet$layerSet[["LAYER.0"]], FUN=gregexpr, pattern="(CG.{0,20}){9}CG")
  bsapply(BSgenome.Scerevisiae.UCSC.sacCer3, FUN=gregexpr, pattern="(CG.{0,20}){9}CG")
  bsapply(genome, FUN=gregexpr, pattern="(CG.{0,20}){9}CG")
  
  
  bsParams <- new("BSParams", X=genome, FUN=gregexpr)
  grepResultBS <- bsapply(bsParams, pattern="(CG.{0,20}){9}CG")
  
  bsGrepGR <- GRanges()
  for(chromName in names(grepResultBS)) {
    grepResult <- grepResultBS[[chromName]]
  #grepResult <- gregexpr(bindingFactor$profile[[thisLayer]]$pattern, layerSet[[thisLayer]][win.starts[i]: win.ends[i]])
  if(grepResult[[1]][1] == -1 ) {  # no grep hits
    win.hits <- IRanges() 
  } else {
    win.hits <- GRanges(chromName, IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE)))
  }
    bsGrepGR <- c(bsGrepGR, win.hits)
  }
  
  grepResult <- gregexpr(bindingFactor$profile[[thisLayer]]$pattern, layerSet[[thisLayer]][win.starts[i]: win.ends[i]])
  
  
  
  
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=XFS, verbose=TRUE, iterations=30, collect.stats = TRUE)
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=XFS[1], verbose=TRUE, iterations=30, collect.stats = TRUE)
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=XFS[2], verbose=TRUE, iterations=30, collect.stats = TRUE)
  finalLayer <- runLayerBinding.BSgenome(layerList=finalLayer, factorSet=XFS[3], verbose=TRUE, iterations=30, collect.stats = TRUE)
  
  matchBindingFactor.BSgenome(finalLayer, bf.CpGisland)
  matchBindingFactor.BSgenome(finalLayer, testFactor2)
  