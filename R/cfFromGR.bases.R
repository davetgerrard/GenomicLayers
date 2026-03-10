

## this version uses total size of overlaps (in bp). 
cfFromGR.bases <- function(query, subject,   maxgap = -1L, minoverlap = 0L, genomeSize=0L, verbose=FALSE)  {
 
  if(verbose)  print("Calculating confusion matrix using number of base-pairs")
  
  # check for consistent genomeSize
    # GRanges object may have no seqinfo/seqlengths with sum == 0
    q_len <- sum(seqlengths(query))
    s_len <- sum(seqlengths(subject))
    
 stopifnot("subject and query must have same genomeSize or both be zero" = q_len == s_len)
 
 if(genomeSize == 0)  genomeSize <- q_len
 stopifnot("genomeSize or seqlengths of provided GR ranges must be greater than 0" = genomeSize > 0)
      
  # created with Jeremy George (UG project student 2022-23)
  #require(GenomicRanges)
  # create a table suitable for a chi-square test
  #genomeSize the sum of the lengths of all sequences in the genome.  Only include sequences
  #         that made it through filtering
  #maxgap  (not used, for future use)
  #minoverlap (not used, for future use)
  #query first set of GRanges. If comparing query to true data, list the prediction as query
  #subject Second set of GRanges. If comparing query to true data, list the true data as subject
  
  
  # remove strand/ In future we might use this but we need a set of non-overlapping strand-less ranges.
  strand(query) <- '*'
  strand(subject) <- '*'
  # reduce each set to a non-overlapping set
  query <- reduce(query)
  subject <- reduce(subject)
  
  # count True-positives
  TP <- sum(width(intersect(query, subject)))
  
  #count false-positives
  FP <- (sum(width(query)))-TP
  
  # count false-negative
  FN <- (sum(width(subject)))-TP
  # count true-negatives
  TN <- genomeSize-(TP+FP+FN)
  
  # create a matrix in the correct format
  # confMatrix is
  #         subject
  # query TRUE FALSE
  # TRUE      TP   FP
  # FALSE     FN   TN
  outMat <- matrix( c(TP, FP, FN, TN),ncol=2, byrow=T) 
  dimnames(outMat) <- list(query=c("TRUE", "FALSE"), subject=c("TRUE", "FALSE"))
  

  
  # check for negative values . This can happen if genomeSize is wrong.
  stopifnot(!any(outMat < 0))   # should produce FALSE if all above 0 so use ! to make true. 
  
  return(outMat)
}


# # examples taken from GenomicRanges package  ?findOverlaps
# 
# gr <- GRanges(
#   seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#   ranges=IRanges(1:10, width=10:1, names=head(letters,10)),
#   strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#   score=1:10,
#   GC=seq(1, 0, length=10)
# )
# gr1 <- GRanges(seqnames="chr2", ranges=IRanges(4:3, 6),
#                strand="+", score=5:4, GC=0.45)
# gr2 <- GRanges(seqnames=c("chr1", "chr1"),
#                ranges=IRanges(c(7,13), width=3),
#                strand=c("+", "-"), score=3:4, GC=c(0.3, 0.5))
# gr3 <- GRanges(seqnames=c("chr1", "chr2"),
#                ranges=IRanges(c(1, 4), c(3, 9)),
#                strand=c("-", "-"), score=c(6L, 2L), GC=c(0.4, 0.1))
# 
# cfFromGR.bases(query=gr, subject=gr1, genomeSize = 40)
# cfFromGR.bases(query=gr1, subject=gr, genomeSize = 40)
# cfFromGR.bases(query=gr1, subject=gr2, genomeSize = 30)   # no overlap
# cfFromGR.bases(query=gr2, subject=gr3, genomeSize = 20)   # no overlap
# cfFromGR.bases(query=gr1, subject=gr3, genomeSize = 20)  
# cfFromGR.bases(query=gr3, subject=gr1, genomeSize = 20)  # check that results are reversed when query/subject reversed
# 
# # If using BSgenome, take genomeSize from sum of seqlengths(genome)
# 
# tprFromConfMat <- function(confMatrix) {   #  true positive rate, TPR a.k.a.sensitivity
#   # confMatrix is
#   #  TP FN
#   #  FP TN
#   # true positive rate = TP/P =  TP/(TP + FN) 
#   return(confMatrix[1,1]/(confMatrix[1,1] + confMatrix[1,2]))
# }
# 
# thisCf <- cfFromGR.bases(gr1, gr3, genomeSize = 20)  
# tprFromConfMat(thisCf)
