# test to show whole genome layer binding 
# three factors are created, the third can only hit after the first has modified.

require(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience



scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)



testFactor2 <- createBindingFactor.DNA_consensus("testFactor2", patternString="ACTGGGCTA")

results2 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

testFactor3 <- createBindingFactor.DNA_consensus("testFactor3", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                                 mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor3)
#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)

# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor4 <- createBindingFactor.layer_region("testFactor4", patternLength=5,
                                                profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))

results4 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor4)
#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=mfLayer, hits=results4, bindingFactor=testFactor4)

# test that a bindingFactor with profile layer not in the layerSet will not match.
testFactor5 <- createBindingFactor.layer_region("testFactor5", patternLength=5,
                                                profile.layers = c("LAYER.3", "LAYER.400"), profile.marks = c(0,1), 
                                                mod.layers = c("LAYER.1", "LAYER.2"), mod.marks=c(0,1))


# test that a bindingFactor with mod layer not in the layerSet will not match.
testFactor6 <- createBindingFactor.layer_region("testFactor6", patternLength=5,
                                                profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c("LAYER.1", "LAYER.200"), mod.marks=c(0,1))

#results6 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor6) # should be error


# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


#modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)


# with the above configuration, there are 41 possible sites across the genome, setting bf.abundances=30, restricts the number that are marked, so the number of potential sites reduces.
#modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, bf.abundances=30)



test_that("Binding factor are compatible with layerset", {
  expect_s4_class(results2, "GRanges")
  expect_true(length(results2)  == 41 )
  expect_error( results5 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor5))
  expect_error( results6 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor6)) 
  
})