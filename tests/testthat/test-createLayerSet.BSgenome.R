
# test createLayerSet()
# three factors are created, the third can only hit after the first has modified.

require(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience



scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)

#length(scLayerSet)
#length(scLayerSet[['layerSet']])

test_that("Genome loaded as LAYER.0", {
  expect_true(class(scLayerSet$layerSet[[1]])  == "BSgenome")
  
})

