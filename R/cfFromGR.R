#' Calculate the confusion matrix from two GRanges objects. 
#'
#' Returns a confusion matrix 
#' listing TP, FP, FN, TN. 
#'
#' @param query a foreground set of regions (perhaps predictions) \code{"GRanges"}
#' @param subject  a background or reference set of regions \code{"GRanges"}
#' @param minPropOverlap a value between 0.0 and 1.0 specifying the proportion of each feature (or gap) in subject covered by query to be counted.
#' @param
#' @param
#' @param
#' @param verbose set to TRUE for more output
#'
#' @return \code{"matrix"}
#' 
#' @seealso 
#'
#' @import Biostrings
#' 
#' @examples
#' library(GenomicRanges)
#' # create a mock Seqinfo object with chromosome names and lengths based on SacCer3
#' test_Si <- x <- Seqinfo(seqnames=c("chrI", "chrII"),
#'                         seqlengths=c(230218, 813184),
#'                         isCircular=c(FALSE, FALSE),
#'                         genome="test_genome")
#' 
#' 
#' 
#' # create a set of perfectly regular regions. 
#' regEnd_I <- seq(from=1000, to = seqlengths(test_Si)[1], by=1000 )
#' regEnd_II <- seq(from=1000, to = seqlengths(test_Si)[2], by=1000 )
#' gr_reg_I <- GRanges(seqnames = "chrI", ranges = IRanges(end=regEnd_I, width=500), seqinfo=test_Si)
#' gr_reg_II <- GRanges(seqnames = "chrII", ranges = IRanges(end=regEnd_II, width=500), seqinfo=test_Si)
#' # need to specify size of chromosomes using seqinfo, borrow this from BSgenome.Scerevisiae.UCSC.sacCer3
#' (gr_reg <- c(gr_reg_I, gr_reg_II))
#' width(gr_reg)
#' gr_regGaps <- gaps(gr_reg, ignore.strand=T )
#' 
#' # runs some checks
#' countOverlaps(gr_reg, gr_regGaps)  # should all be zero.  Perfect opposites. 
#' sum(countOverlaps(gr_reg, gr_regGaps) )
#' 
#' 
#' gr_fullChroms <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=seqlengths(test_Si)), seqinfo=test_Si)
#' # a bad prediction might predict ALL the chromosomes is in a given state (or none of it). 
#' 
#' # now create a GR with exactly half the chromosome in covered in one single range. 
#' gr_halfChrom <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=floor(seqlengths(test_Si)/2)), seqinfo=test_Si)
#' 
#' 
#' cfFromGR(gr_halfChrom, gr_reg)
#' 
#' cfFromGR(gr_reg, gr_reg)  # perfect overlap
#' 
#' cfFromGR(gr_halfChrom, gr_halfChrom)  ### perfect overlap between two very large regions. 
#' 
#' cfFromGR(gr_fullChroms, gr_reg)   # zero negative predictions 
#' cfFromGR(gr_regGaps, gr_reg)    # perfect missing
#' 
#' cfFromGR(gr_halfChrom, gr_regGaps) 
#' cfFromGR(gr_regGaps, gr_halfChrom ) 
#' cfFromGR(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001)   # will call pretty much any overlap as a positive. 
#' 
#' 
#' 
#' @export
cfFromGR <- function(query, subject, method=c("features", "bases"), minPropOverlap=0.5, maxgap = -1L, minoverlap = 0L, genomeSize=0L, verbose=FALSE) {
    
    method <- method[1]  # if none given, default to "features"
    # remove strand/ In future we might use this but we need a set of non-overlapping strand-less ranges.
    strand(query) <- '*'
    strand(subject) <- '*'
    # reduce each set to a non-overlapping set
    query <- reduce(query)
    subject <- reduce(subject)
    
     
    confMat <- switch(method[1], 
                    features =  cfFromGR.features(query=query, subject=subject, minPropOverlap=minPropOverlap, verbose=verbose),
                    bases = cfFromGR.bases(query=query, subject=subject, maxgap = maxgap, minoverlap = minoverlap, genomeSize=genomeSize, verbose=verbose)
    )
    if(verbose) {
        print(outMat)
        print( matrix(c("TP", "FP", "FN", "TN"), ncol=2, byrow = T))  # to check layout.
        
    }
    return(confMat)
}


library(GenomicRanges)
# create a mock Seqinfo object with chromosome names and lengths based on SacCer3
test_Si <- x <- Seqinfo(seqnames=c("chrI", "chrII"),
                        seqlengths=c(230218, 813184),
                        isCircular=c(FALSE, FALSE),
                        genome="test_genome")



# create a set of perfectly regular regions. 
regEnd_I <- seq(from=1000, to = seqlengths(test_Si)[1], by=1000 )
regEnd_II <- seq(from=1000, to = seqlengths(test_Si)[2], by=1000 )
gr_reg_I <- GRanges(seqnames = "chrI", ranges = IRanges(end=regEnd_I, width=500), seqinfo=test_Si)
gr_reg_II <- GRanges(seqnames = "chrII", ranges = IRanges(end=regEnd_II, width=500), seqinfo=test_Si)
# need to specify size of chromosomes using seqinfo, borrow this from BSgenome.Scerevisiae.UCSC.sacCer3
(gr_reg <- c(gr_reg_I, gr_reg_II))
width(gr_reg)
gr_regGaps <- gaps(gr_reg, ignore.strand=T )

# runs some checks
countOverlaps(gr_reg, gr_regGaps)  # should all be zero.  Perfect opposites. 
sum(countOverlaps(gr_reg, gr_regGaps) )


gr_fullChroms <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=seqlengths(test_Si)), seqinfo=test_Si)
# a bad prediction might predict ALL the chromosomes is in a given state (or none of it). 

# now create a GR with exactly half the chromosome in covered in one single range. 
gr_halfChrom <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=floor(seqlengths(test_Si)/2)), seqinfo=test_Si)


cfFromGR.features(gr_halfChrom, gr_reg)
cfFromGR(query = gr_halfChrom, subject = gr_reg)

cfFromGR.features(gr_reg, gr_reg)  # perfect overlap
cfFromGR(gr_reg, gr_reg, verbose=T) 

cfFromGR.features(gr_halfChrom, gr_halfChrom)  ### perfect overlap between two very large regions. 
cfFromGR(gr_halfChrom, gr_halfChrom, verbose=T)

cfFromGR.features(gr_fullChroms, gr_reg)   # zero negative predictions 
cfFromGR(gr_fullChroms, gr_reg)
cfFromGR.bases(query=gr_fullChroms, subject=gr_reg)
cfFromGR(gr_fullChroms, gr_reg, method="bases", verbose=TRUE)

cfFromGR.features(gr_regGaps, gr_reg)    # perfect missing
cfFromGR(gr_regGaps, gr_reg) 
cfFromGR(gr_regGaps, gr_reg, method="bases", verbose=TRUE) 

cfFromGR.features(gr_halfChrom, gr_regGaps) 
cfFromGR.features(gr_regGaps, gr_halfChrom ) 
cfFromGR(gr_regGaps, gr_halfChrom ) 
cfFromGR.features(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001)   # will call pretty much any overlap as a positive. 
cfFromGR(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001) 




