
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


#matches <- matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor=factorSetTest[['bf.test']], verbose=FALSE)

factorTypes <- c("DNA_motif", "DNA_region","layer_region","layer_island")

set.seed(12345)

factorSetRandom <- list()
for(thisType in factorTypes)  {
  thisName <- paste0("testBf.", thisType)
  factorSetRandom[[thisName]] <- createRandomBindingFactor(thisName, layerSet = layerList.5$layerSet, type = thisType)
}
#modLayerSet <- runLayerBinding(layerList.5, factorSet=factorSetTest)

(sum.table <- print.bfSet(factorSetRandom))

names(factorSetRandom)
all(as.character(sum.table$types)  == factorTypes)


test_that("Factors have expected types", {
  expect_true(all(as.character(sum.table$types)  == factorTypes))

})

test_that("Number of factors correct", {
  expect_equal(length(factorSetRandom), length(factorTypes))
  
})
