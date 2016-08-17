
require(Biostrings)
data(yeastSEQCHR1)   # included in Biostrings


motifPattern <- "TGGCCTC"


hits <- matchPattern(motifPattern, yeastSEQCHR1)

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


matches <- matchBindingFactor(layerSet=layerList.5$layerSet, bindingFactor=bf, verbose=FALSE)


test_that("Factor matches equal pattern matches equal eleven", {
  expect_equal(length(hits), 11)
  expect_equal(length(matches), 11)
  expect_equal(length(hits), length(matches))
})


