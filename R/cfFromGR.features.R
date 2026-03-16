#' Calculate the confusion matrix from counts of features shared between two GRanges objects. 
#'
#' Counts the number of features (TP) or gaps (FP) in a subject (reference) GRanges object
#' that are covered by the features defined in the query.  Returns a confusion matrix 
#' listing TP, FP, FN, TN. 
#'
#' @param query a foreground set of regions (perhaps predictions) \code{"GRanges"}
#' @param subject  a background or reference set of regions \code{"GRanges"}
#' @param minPropOverlap a value between 0.0 and 1.0 specifying the proportion of each feature (or gap) in subject covered by query to be counted.
#' @param verbose set to TRUE for more output
#'
#' @return \code{"matrix"}
#' 
#' @seealso \code{\link{cfFromGR}} 
#'
#' @import Biostrings
#' @import IRanges
#' 
#' @export
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
#' cfFromGR.features(gr_halfChrom, gr_reg)
#' 
#' cfFromGR.features(gr_reg, gr_reg)  # perfect overlap
#' 
#' cfFromGR.features(gr_halfChrom, gr_halfChrom)  ### perfect overlap between two very large regions. 
#' 
#' cfFromGR.features(gr_fullChroms, gr_reg)   # zero negative predictions 
#' cfFromGR.features(gr_regGaps, gr_reg)    # perfect missing
#' 
#' cfFromGR.features(gr_halfChrom, gr_regGaps) 
#' cfFromGR.features(gr_regGaps, gr_halfChrom ) 
#' cfFromGR.features(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001)   # will call pretty much any overlap as a positive. 
#' 
#' 
#' 
cfFromGR.features <- function(query, subject, minPropOverlap=0.5, verbose=FALSE)  {
    
    if(verbose)  print("Calculating confusion matrix using number of features")
    
    require(GenomicRanges)
    stopifnot("Missing seqlengths for subject. Define using seqinfo()" = any(!is.na(seqlengths(subject))))  # check that the subject has seqlengths for all chroms. Required to calculate gaps. 
    
    # by default, reports counts of "subject" features that are TP, FN, FP, TN. 
    
    gap_sub <- gaps(subject, ignore.strand=T )
    
    compInt <- intersect(subject, query)   # calculate the intesection ranges. Need this to know the width. 
    # we want to know how many of the query ranges are covered 
    
    ov <- findOverlaps(compInt, subject)    # get indices for the original regions that DO have overlaps. (should be sames as overlapsAny())
    
    # need to sum the overlapping sections for each of the original subjects.
    # This makes the function quite a bit slower. 
    bySubject <- by(compInt, INDICES=subjectHits(ov), FUN=function(x) sum(width(x)))
    hitInd <- as.integer(names(bySubject))  # the indices in the original subject GR
    hitSpans <-  as.integer(bySubject)    # the summed intersection lengths. 
    TP_count <- length(which(hitSpans > width(subject)[hitInd] * minPropOverlap))
    
    # create index of original features that have an overlap that is sufficiently long. 
    #true_hits <- width(compInt[queryHits(ov)])  >=  width(subject[subjectHits(ov)]) *minPropOverlap
    #TP_count <- sum(true_hits)
    # the above two lines removed as this does not account for multiple query hits on a long subject summing to greater than minPropOverlap of the subject
    FN_count <- length(subject) - TP_count   # number of subject features missed
    
    # now repeat above using gaps between features. 
    ### Does this need to use " 1 - minPropOverlap"?
    (compInt <- intersect(gap_sub, query))   # calculate the intesection ranges. Need this to know the width. 
    # we want to know how many of the query ranges are covered 
    
    ov <- findOverlaps(compInt, gap_sub)    # get indices for the original gaps that DO have overlaps. (should be sames as overlapsAny())
    
    # need to sum the overlapping sections for each of the original gaps
    bySubject <- by(compInt, INDICES=subjectHits(ov), FUN=function(x) sum(width(x)))
    hitInd <- as.integer(names(bySubject))  # the indices in the gaps GR
    hitSpans <-  as.integer(bySubject)    # the summed intersection lengths. 
    FP_count <- length(which(hitSpans > width(gap_sub)[hitInd] * minPropOverlap))
    
    
    # create index of original features that have an overlap that is sufficiently long.
    #false_hits <- width(compInt[queryHits(ov)])  >=  width(gap_sub[subjectHits(ov)]) *minPropOverlap
    #FP_count <- sum(false_hits)
    TN_count <- length(gap_sub) - FP_count
    
    # should return confusion matrix in format expected by Caret package. 
    # create a matrix in the correct format
    # confMatrix is
    #         reference
    # predicted TRUE FALSE
    # TRUE      TP   FP
    # FALSE     FN   TN
    confMat <- matrix(c(TP_count, FP_count, FN_count, TN_count), ncol=2, byrow=T)
    dimnames(confMat) <- list(query=c("TRUE", "FALSE"), subject=c("TRUE", "FALSE"))
    
    # check we have counted each feature and gap in the subject exactly once. 
    stopifnot( TP_count + FN_count + FP_count + TN_count == length(subject)  + length(gap_sub))
    stopifnot(!any(confMat < 0))
    
    return(confMat)
    
}

# 
# 
# # EXAMPLES
# 
# library(GenomicRanges)
# # create a mock Seqinfo object with chromosome names and lengths based on SacCer3
# test_Si <- x <- Seqinfo(seqnames=c("chrI", "chrII"),
#                         seqlengths=c(230218, 813184),
#                         isCircular=c(FALSE, FALSE),
#                         genome="test_genome")
# 
# 
# 
# # create a set of perfectly regular regions. 
# regEnd_I <- seq(from=1000, to = seqlengths(test_Si)[1], by=1000 )
# regEnd_II <- seq(from=1000, to = seqlengths(test_Si)[2], by=1000 )
# gr_reg_I <- GRanges(seqnames = "chrI", ranges = IRanges(end=regEnd_I, width=500), seqinfo=test_Si)
# gr_reg_II <- GRanges(seqnames = "chrII", ranges = IRanges(end=regEnd_II, width=500), seqinfo=test_Si)
# # need to specify size of chromosomes using seqinfo, borrow this from BSgenome.Scerevisiae.UCSC.sacCer3
# (gr_reg <- c(gr_reg_I, gr_reg_II))
# width(gr_reg)
# gr_regGaps <- gaps(gr_reg, ignore.strand=T )
# 
# # runs some checks
# countOverlaps(gr_reg, gr_regGaps)  # should all be zero.  Perfect opposites. 
# sum(countOverlaps(gr_reg, gr_regGaps) )
# 
# 
# gr_fullChroms <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=seqlengths(test_Si)), seqinfo=test_Si)
# # a bad prediction might predict ALL the chromosomes is in a given state (or none of it). 
# 
# # now create a GR with exactly half the chromosome in covered in one single range. 
# gr_halfChrom <- GRanges(seqnames=seqnames(test_Si), ranges=IRanges(start=rep(1, length(seqnames(test_Si))), end=floor(seqlengths(test_Si)/2)), seqinfo=test_Si)
# 
# 
# cfFromGR.features(gr_halfChrom, gr_reg)
# 
# cfFromGR.features(gr_reg, gr_reg)  # perfect overlap
# 
# cfFromGR.features(gr_halfChrom, gr_halfChrom)  ### perfect overlap between two very large regions. 
# 
# cfFromGR.features(gr_fullChroms, gr_reg)   # zero negative predictions 
# cfFromGR.features(gr_regGaps, gr_reg)    # perfect missing
# 
# cfFromGR.features(gr_halfChrom, gr_regGaps) 
# cfFromGR.features(gr_regGaps, gr_halfChrom ) 
# cfFromGR.features(gr_regGaps, gr_halfChrom , minPropOverlap = 0.001)   # will call pretty much any overlap as a positive. 

