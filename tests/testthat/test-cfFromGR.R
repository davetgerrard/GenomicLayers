

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