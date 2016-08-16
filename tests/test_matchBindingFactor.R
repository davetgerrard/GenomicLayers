
require(Biostrings)
data(yeastSEQCHR1)   # included in Biostrings


motifPattern <- "CCACACACC"


matchPattern(motifPattern, yeastSEQCHR1)

modList <- list()
modList[['LAYER.1']] <- list(state=1, stateWidth=nchar(motifPattern), offset=0, align="centre") 

bf <-   list(name="testBf.motif", type="DNA_motif",
             profile=list(LAYER.0=list(pattern=DNAString(motifPattern) , mismatch.rate=0, length=nchar(motifPattern))),
             mods= modList)


n.layers <- 5
layerSet.5 <- list(LAYER.0 = DNAString(yeastSEQCHR1))
for(i in 1:n.layers) {
  layerSet.5[[paste('LAYER.', i , sep="")]] <- IRanges()    # use IRanges to store state of layers. TODO limit to chrom length
}
layerList.5 <- list(layerSet=layerSet.5, history=NULL)


matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor=bf, verbose=TRUE)


test_that("trigonometric functions match identities", {
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})



thisLayer <- "LAYER.0"
round(bf$profile[[thisLayer]]$mismatch.rate * bf$profile[[thisLayer]]$length)
