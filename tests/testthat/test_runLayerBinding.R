
require(Biostrings)
data(yeastSEQCHR1)   # included in Biostrings


motifPattern <- "TGGCCTC"


hits <- matchPattern(motifPattern, yeastSEQCHR1)

modList <- list()
modList[['LAYER.1']] <- list(state=1, stateWidth=nchar(motifPattern), offset=0, align="centre") 

bf <-   list(name="testBf.motif", type="DNA_motif",
             profile=list(LAYER.0=list(pattern=DNAString(motifPattern) , mismatch.rate=0, length=nchar(motifPattern))),
             mods= modList)


factorSetTest <- list(bf.test = bf)     # binding factors must be in a list

n.layers <- 5
layerSet.5 <- list(LAYER.0 = DNAString(yeastSEQCHR1))
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}
layerList.5 <- list(layerSet=layerSet.5, history=NULL)


matches <- matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor=factorSetTest[['bf.test']], verbose=FALSE)




modLayerSet <- runLayerBinding(layerList.5, factorSet=factorSetTest)


test_that("Modification matches pattern length", {
  expect_equal(nchar(motifPattern),  width(modLayerSet$layerSet[["LAYER.1"]]))

})


