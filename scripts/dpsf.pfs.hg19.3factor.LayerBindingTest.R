# Script to run whole genome layer binding on hg19 - may need high memory computer
# three factors are created, the third can only hit after the first has modified.


require(Biostrings)
#setwd('C:/Users/Dave/HalfStarted/predictFromSequence/')
#source('C:/Users/Dave/Dropbox/Temp/predictFromSequence.functions.R')
#source('scripts/predictFromSequence.functions.R')
source('scripts/pfs.functions.R')

library(BSgenome.Hsapiens.UCSC.hg19)



genome <- BSgenome.Hsapiens.UCSC.hg19   # for convenience



scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)



testFactor2 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")

#results2 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

testFactor3 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                             mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

#results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor3)
#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)

# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                             mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))

#results4 <- matchBindingFactor.BSgenome(layerSet = mfLayer, bindingFactor = testFactor4)
#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=mfLayer, hits=results4, bindingFactor=testFactor4)

# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


#modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)


# with the above configuration, there are 41 possible sites across the genome, setting iterations=30, restricts the number that are marked, so the number of potential sites reduces.
system.time(
  modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, iterations=5000)
)
modTest

gc()








