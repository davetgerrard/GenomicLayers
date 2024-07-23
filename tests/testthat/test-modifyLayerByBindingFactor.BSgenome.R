# Important to use an actual BSgenome object in tests
# However, tests should be more deterministic and rely on known motifs
#   Some layers could be pre-populated with specific ranges to test conditional matching.

require(Biostrings)
library(BSgenome.Scerevisiae.UCSC.sacCer3)

genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience

scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)


testFactor2 <- createBindingFactor.DNA_consensus("testFactor2", patternString="ACTGGGCTA")  # 41 matches across SacCer3

results2 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor2)

# a factor to change some layer to match.
testFactor3 <- createBindingFactor.DNA_consensus("testFactor3", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
                                                 mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(1,1))

results3 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor3)
#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)

scLayerSet3 <- modifyLayerByBindingFactor.BSgenome(layerSet = scLayerSet, hits= results3, bindingFactor=testFactor3)

# changes to LAYER.2 and LAYER.4 should match results3 but need to be ordered to be "identical"


test_that("modified layer matches simple hits with same bf", {
  expect_true(identical(results3[order(results3)], scLayerSet3$layerSet$LAYER.2[order(scLayerSet3$layerSet$LAYER.2)])
)
})

# modify hits to increase them in size
testFactor4 <- createBindingFactor.layer_region("testFactor4", patternLength=9,
                                                profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c( "LAYER.3"), mod.marks=c(1),
                                                stateWidth = 500)
results4 <- matchBindingFactor.BSgenome(layerSet = scLayerSet3, bindingFactor = testFactor4)
scLayerSet4 <- modifyLayerByBindingFactor.BSgenome(layerSet = scLayerSet3, hits= results4, bindingFactor=testFactor4)

testFactor5 <- createBindingFactor.layer_region("testFactor5", patternLength=9,
                                                profile.layers = c("LAYER.1", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c( "LAYER.1"), mod.marks=c(1),
                                                stateWidth = 5)
results5 <- matchBindingFactor.BSgenome(layerSet = scLayerSet3, bindingFactor = testFactor5)
scLayerSet5 <- modifyLayerByBindingFactor.BSgenome(layerSet = scLayerSet3, hits= results5, bindingFactor=testFactor5)


test_that("binding factor can change size of regions during modification", {
  expect_true(all(width(scLayerSet4$layerSet$LAYER.3) == 500))
  expect_true(all(width(scLayerSet5$layerSet$LAYER.1) == 5))
})

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

# check that a profile looking for 1 will not find any.  N.B this WILL bind AFTER testFactor2
testFactor7 <- createBindingFactor.layer_region("testFactor7", patternLength=5,
                                                profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c( "LAYER.5"), mod.marks=c(1),
                                                stateWidth = 500)

results7 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor7)


# now can match things genome wide. Need to run layerBinding and modification.

# need to have a factorSet, a list of bindingFactors

testFS <- list(testFactor2=testFactor2, testFactor3=testFactor3, testFactor4=testFactor4)


# also need for modifyLayerByBindingFactor.Views to work on BSgenome and hits

#mfLayer <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results3, bindingFactor=testFactor3)


#modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE)


# with the above configuration, there are 41 possible sites across the genome, setting bf.abundances=30, restricts the number that are marked, so the number of potential sites reduces.
#modTest <- runLayerBinding.BSgenome(layerList=scLayerSet, factorSet=testFS, verbose=TRUE, bf.abundances=30)

# there should be results2 entries
# chrII 579821-579829      *
# chrII 682352-682360
# add a single region to LAYER.5 such that an overlap with the results2 produces exactly one valid hit.
scLayerSet$layerSet$LAYER.5 <- GRanges(seqnames = "chrII", ranges=IRanges(start=579823, end=682366), seqinfo = seqinfo(genome))

testFactor8 <- createBindingFactor.DNA_consensus("testFactor8", patternString="ACTGGGCTA", profile.layers = c( "LAYER.5"), profile.marks = c(1), 
                                                 mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))

results8 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor8)
width(results8)

# set LAYER.3 to be the results from results2 only larger
# use this to test if a layer-region type is correctly selecting by size.
scLayerSet$layerSet$LAYER.3  <- resize(results2, width=50)

width(intersect(scLayerSet$layerSet$LAYER.3, scLayerSet$layerSet$LAYER.5))  # 48 and 15


# puts some random blocks down on LAYER.5
testFactor10 <- createBindingFactor.layer_region("testFactor10", patternLength=40,
                                                 profile.layers = c("LAYER.3", "LAYER.5"), profile.marks = c(1,1), 
                                                 mod.layers = c("LAYER.2"), mod.marks=c(1))
results10 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor10)
width(results10 )
# ideally, none of these should be less than patternLength.

scLayerSet$layerSet$LAYER.4 <- GRanges(seqnames = "chrIV", ranges = IRanges(start=c(1,300, 500), end=c(100, 400, 2000)), seqinfo = seqinfo(genome))
# leaves two gaps ofs 100bp and 200bp
scLayerSet$layerSet$LAYER.1 <- GRanges(seqnames = "chrIV", ranges = IRanges(start=c(1), end=c(2000)), seqinfo = seqinfo(genome))
# block out the first 2000bp, this can be used to test patterns just within the gaps on LAYER.4
# a bindingFactor that should insert into the first gap but not the second.
testFactor11 <- createBindingFactor.layer_region("testFactor11", patternLength=180,
                                                 profile.layers = c("LAYER.1", "LAYER.4"), profile.marks = c(1,0), 
                                                 mod.layers = c("LAYER.4"), mod.marks=c(1), stateWidth=147)
results11 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor11)
# correctly only returning potential matches in the larger gap. 

# what to do about hits that would overlap if two or more applied?
# If valid hits and dispersed across genome and n.hits >> abundance, not much of a problem
# however, as n.hits reduces (e.g. nucleosome model) and abundance stays high, 
#  these will overlap within one iteration and clash when use runLayerBinding() and modifyLayerByBindingFactor()
mfLayer11 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet, hits=results11[c(2,5,10)], bindingFactor=testFactor11)
width(mfLayer11$layerSet$LAYER.4)    # ideally the second regions should now be width 147 but is greater than this due to multiple overlapping hits.
# how would changing the sampling method better represent biology?
# want to avoid nucleoomes of > 147bp.   
# but also may want to model sink effect such that the availability of sites drives BFs onto particular sites. 
# easiest option seems to be to keep the cleaner functions that remove anything != 147 on certain layers each round. 
# Somewhat unsatisfactory and forces the BF ordering to be carefully considered (unless eventyally added as option to runLayerBinding). 

test_that("Binding factor are compatible with layerset", {
  expect_s4_class(results2, "GRanges")
  expect_true(length(results2)  == 41 )
  expect_error( results5 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor5))
  expect_error( results6 <- matchBindingFactor.BSgenome(layerSet = scLayerSet, bindingFactor = testFactor6)) 
})
  
test_that("intersects are correct", {
  expect_true(length(results8)  == 1 )
  expect_true(min(width(results10)) == 40)
  })


# SECTION testing the offset function

upDownFuncRnorm <- function(n, offset.mean, offset.sd)  {
  y <- round(rnorm(n, mean=offset.mean, sd=offset.sd))
  z <- sample(c(1, -1), length(y), replace=T)  # random vector of 1,-1  to negate half the values
  
  return(round(y*z))
}

# mods layer 3 with NO offset
testFactor40 <- createBindingFactor.layer_region("testFactor40", patternLength=9,
                                                 profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                 mod.layers = c( "LAYER.3"), mod.marks=c(1)
                                                 )                                                
results40 <- matchBindingFactor.BSgenome(layerSet = scLayerSet3, bindingFactor = testFactor40)
scLayerSet340 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet3, hits=results40, bindingFactor=testFactor40)

# 41 is same as 40 but offsets results by -8 (to give 1bp overlap)
testFactor41 <- createBindingFactor.layer_region("testFactor41", patternLength=9,
                                                 profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                 mod.layers = c( "LAYER.3"), mod.marks=c(1),
                                                 offset = -8)  
results41 <- matchBindingFactor.BSgenome(layerSet = scLayerSet3, bindingFactor = testFactor41)
scLayerSet341 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet3, hits=results41, bindingFactor=testFactor41)

ov340.341 <- intersect(scLayerSet340$layerSet$LAYER.3, scLayerSet341$layerSet$LAYER.3)
table(width(ov340.341))  # should be a count table showing all values are '1'
sumDiff340_341 <- sum(start(scLayerSet340$layerSet$LAYER.3) - end(scLayerSet341$layerSet$LAYER.3))

test_that("fixed offsets are working", {
  expect_true(sumDiff340_341  == 0 )
})

# next is a demonstration of using a function to dynamically generate offset according to a distribution
# as it involves randomisation and quite variable results, I'm not sure a unit test is appropriate.
testFactor43 <- createBindingFactor.layer_region("testFactor43", patternLength=9,
                                                profile.layers = c("LAYER.3", "LAYER.4"), profile.marks = c(0,1), 
                                                mod.layers = c( "LAYER.3"), mod.marks=c(1),
                                                offset.method=upDownFuncRnorm, 
                                                offset.params=list(offset.mean=c(5,20), offset.sd=c(1:4))                                                )
results43 <- matchBindingFactor.BSgenome(layerSet = scLayerSet3, bindingFactor = testFactor43)

scLayerSet343 <- modifyLayerByBindingFactor.BSgenome(layerSet=scLayerSet3, hits=results43, bindingFactor=testFactor43, verbose=T)

diff340_343 <- start(scLayerSet340$layerSet$LAYER.3) - start(scLayerSet343$layerSet$LAYER.3)
sum(diff340_343)
    