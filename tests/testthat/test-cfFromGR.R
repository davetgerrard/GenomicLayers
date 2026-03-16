

# test cfFromGR.R function
library(GenomicRanges)
# create a mock Seqinfo object with chromosome names and lengths based on SacCer3
test_Si  <- Seqinfo(seqnames=c("chrI", "chrII"),
                        seqlengths=c(230218, 813184),
                        isCircular=c(FALSE, FALSE),
                        genome="test_genome")


# create a set of perfectly regular GenomicRanges regions. 
regEnd_I <- seq(from=1000, to = seqlengths(test_Si)[1], by=1000 )
regEnd_II <- seq(from=1000, to = seqlengths(test_Si)[2], by=1000 )
gr_reg_I <- GRanges(seqnames = "chrI", ranges = IRanges(end=regEnd_I, width=500), seqinfo=test_Si)
gr_reg_II <- GRanges(seqnames = "chrII", ranges = IRanges(end=regEnd_II, width=500), seqinfo=test_Si)
# need to specify size of chromosomes using seqinfo
gr_reg <- c(gr_reg_I, gr_reg_II)
#width(gr_reg)
gr_regGaps <- gaps(gr_reg, ignore.strand=T )

# runs some checks

#sum(countOverlaps(gr_reg, gr_regGaps) )

gr_fullChroms <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=seqlengths(test_Si)), seqinfo=test_Si)
# a bad prediction might predict ALL the chromosomes is in a given state (or none of it). 

# now create a GR with exactly half the chromosome in covered in one single range. 
gr_halfChrom <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=floor(seqlengths(test_Si)/2)), seqinfo=test_Si)

gr_noSeqinfo <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=seqlengths(test_Si)))

# for method = "features", subject must have seqinfo
#cfFromGR(query = gr_noSeqinfo, subject = gr_noSeqinfo)   # this throws an error

# for method = "bases", either subject must have seqlengths > 0 or genomeSize is > 0

test_that("Missing seqinfo in features method throws an error.", {  # this test could be faster if limited to one chromosome.
    expect_error(cfFromGR(query = gr_noSeqinfo, subject = gr_noSeqinfo))
})


# #TODO:  deduce a proper test from the below code. Note the 0.2 result. 
# testSI <- Seqinfo(seqnames="chr1",
#                   seqlengths=c(3000),
#                   isCircular=c(FALSE),
#                   genome="test_genome")
# testSubj_GR <- GRanges(seqnames= "chr1", ranges=IRanges(start=1000, width=1000),  seqinfo=testSI)
# testQ_GR <- GRanges(seqnames= "chr1", ranges=IRanges(start=c(300, 600, 1300, 1600, 2300, 2600), width=100),  seqinfo=testSI)
# 
# # there is one feature and two gaps. The queries cover a combined 200bp within each 1000 gap or feature. 
# cfFromGR(query=testQ_GR, subject=testSubj_GR)  # default minoverlap
# cfFromGR(query=testQ_GR, subject=testSubj_GR, minPropOverlap = 0.1) # calls all gaps and regions as positive. [CORRECT]
# cfFromGR(query=testQ_GR, subject=testSubj_GR, minPropOverlap = 0.19) # calls the one region but also both gaps. [CORRECT]
# cfFromGR.features(query=testQ_GR, subject=testSubj_GR, minPropOverlap = 0.2)  # calls one gap...  [FALSE]
# cfFromGR(query=testQ_GR, subject=testSubj_GR, minPropOverlap = 0.21 )  # Fails to call anything (correct here with minPropoOverlap). [CORRECT]
# cfFromGR(query=testQ_GR, subject=testSubj_GR, minPropOverlap = 0.3) # Fails to call anything (correct here with minPropoOverlap). [CORRECT]
# 



# 
# 
# cfFromGR(query = gr_halfChrom, subject = gr_reg) # approx equal in each class
# cfFromGR(query = gr_halfChrom, subject = gr_reg, method="bases", ) 
# 
# 
# cfFromGR(gr_reg, gr_reg, verbose=T) # perfect overlap
# cfFromGR(gr_reg, gr_reg, verbose=T)
# 
# 
# cfFromGR(gr_halfChrom, gr_halfChrom, verbose=T)  ### perfect overlap between two very large regions. 
# cfFromGR(gr_halfChrom, gr_halfChrom, verbose=T)
# 
# # zero negative predictions 
# cfFromGR(gr_fullChroms, gr_reg)
# cfFromGR(gr_fullChroms, gr_reg, method="bases", verbose=TRUE)
# 
# # perfect missing
# cfFromGR(gr_regGaps, gr_reg) 
# cfFromGR(gr_regGaps, gr_reg, method="bases", verbose=TRUE) 
# 
# # subject is a very large domain. 
# cfFromGR(gr_regGaps, gr_halfChrom )  # misses one because overlap < minPropOverlap
# # will call pretty much any overlap as a positive. 
# cfFromGR(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001) 