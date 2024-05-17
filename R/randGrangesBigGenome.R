#' randGrangesBigGenome
#'
#'  Generate random non-overlapping GRanges features within parameters of size and distance
#'
#' @param genome A BSgenome object
#' @param n  NULL
#' @param sizeFunc NULL  must be a function that generates a vector of positive integers
#   and one parameter must be "n" corresponding to the number of integers returned e.g.  gapFunc = function(value, n) rep(x=value, times=n)
#' @param gapFunc NULL must be a function that generates a vector of positive integers
#   and one parameter must be "n" corresponding to the number of integers returned
#' @param argsSizeFunc NULL a list giving named parameter values passed to sizeFunc
#' @param argsGapFunc NULL a list giving named parameter values passed to gapFunc
#' @param minGapSize 1
#' @param minRange 1 
#' @param verbose FALSE  Give more output
#'
#' @return \code{"GRanges"} object with seqinfo matching genome
#' 
#'
#' @examples
#'
#' library(BSgenome.Scerevisiae.UCSC.sacCer3)
#' 
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3
#' 
#' 
#' 
#' 
#' test <- randGrangesBigGenome(genome=genome, 
#'                              sizeFunc = function(value, n) rep(x=value, times=n), 
#'                              gapFunc = function(value, n) rep(x=value, times=n), 
#'                              argsSizeFunc = list(value=147, n=5),
#'                              argsGapFunc = list(value=30, n=7), verbose=T)
#' # fixed nuclesomes of 147bp with random gaps following a poisson.
#' testGR <- randGrangesBigGenome(genome=genome, 
#'                                sizeFunc = function(value, n) rep(x=value, times=n), 
#'                                gapFunc = function(n, value) rpois(n = n , lambda = value), 
#'                                argsSizeFunc = list(value=147, n=5),
#'                                argsGapFunc = list(value=40, n=7))
#' 
#' 
#' 
#' 
#' @export
# 
# could not see an answer here https://www.biostars.org/p/225520/
#sizeFunc and gapFunc each must be a function that generates a vector of positive integers
#   and one parameter must be "n" corresponding to the number of integers returned
# gapFunc = function(value, n) rep(x=value, times=n)
# this version an attempt to generate features one chrom at a time and then rejoin them.
randGrangesBigGenome <- function(genome, n=NULL, sizeFunc=NULL, gapFunc=NULL, 
                        argsSizeFunc, argsGapFunc, minGapSize=1, minRange=1, verbose=FALSE)  {
  
  # what would happen if gaps of size zero were allowed?
  # would we wnat items to be merged?  (that would probably break the even-ness 
  # of the chromosomal allocation if a merged region spanned a chrom boundary
  
  # must only determine two of the three inputs n, featureSize, gapSize because 
  # knowing two, determines the third
  
  # strategy:  create a linear set of features and gaps longer than the chromosome

  chromNames <- names(seqlengths(genome))
  genomeLength <- sum(seqlengths(genome)) 
  
  # decisions to make on whether it is preferable to start a chrom with a gap or a feature
  # should warn/error if any gaps or features are larger than half the size of the smallest chrom
  
  #test passed funcs
  
  local_argsSizeFunc <- argsSizeFunc
  local_argsSizeFunc['n']  <- 100
  #print(do.call(sizeFunc, args=local_argsSizeFunc))  # should output 100 values
  local_argsGapFunc <- argsGapFunc
  local_argsGapFunc['n']  <- 100
  #print(do.call(gapFunc, args=local_argsGapFunc))  # should output 100 values
  
  sampleWidths <- do.call(sizeFunc, args=local_argsSizeFunc)  # feature widths
  sampleGaps <- do.call(gapFunc, args=local_argsGapFunc)  # gaps between features (+1 to get to end of chrom)
  
  #sampleWidths <- rpois(300, 50)  # feature widths
  #sampleGaps <- rpois(length(fWidths) +1, 400)  # gaps between features (+1 to get to end of chrom)
  
  m_sw <- mean(sampleWidths)
  m_sg <- mean(sampleGaps)
  
  allGR <- GRanges(seqinfo=seqinfo(genome))
  
  
  for(thisChrom in chromNames) {
    if (verbose) print(thisChrom)
    chromLength <- seqlengths(genome)[thisChrom]
    
    #genomeLength/m_sw
    #genomeLength/m_sg
    minPairs <- chromLength/(m_sw + m_sg)
    if (verbose) print(paste("Chromosome size: ", chromLength, "Mean feature width: ", round(m_sw,2), "Mean gap width", round(m_sg,2),  "Calculated number of feature/gaps:", round(minPairs,2)))
    
    stopifnot( "requires more features than half genome length" = minPairs < chromLength/2)
    stopifnot("fewer than one feature per chromosome" = minPairs > 1)
    
    excessPairs <- ceiling(minPairs + 1)  # add extra 
    local_argsGapFunc['n']  <- excessPairs
    local_argsSizeFunc['n']  <- excessPairs
    gaps <- do.call(gapFunc, args=local_argsGapFunc)
    widths <-   do.call(sizeFunc, args=local_argsSizeFunc)
   
    # set gaps or widths less than 1 to one.  Should not be necessary. 
    gaps[gaps < 1]  <- 1
    widths[widths < 1]  <- 1   
    #gaps <- rpois(excessPairs , 400)
    #widths <-   rpois(excessPairs, 50)
    #interleave the pairs into a single vector so that cumulative length of all can be cumsum
    interleaved <- integer(length=excessPairs*2)
    interleaved[seq(by=2, from=1, length.out=excessPairs)]  <- gaps   # gaps first
    interleaved[seq(by=2, from=2, length.out=excessPairs)]  <- widths
    interleavedCumsum <- cumsum(interleaved)
    
    # extract start and end points (relative to start of whole genome)
    ends <- interleavedCumsum[seq(by=2, from=2, length.out=excessPairs)]
    starts <- (ends - widths) + 1
    # check results using alternate method
    #head(interleavedCumsum[seq(by=2, from=1, length.out=excessPairs)]+1)  # alt calculation.
    #head(gaps)
    #head(starts)
    #hist((ends - starts))   # should match original widths
    #hist(widths)
    
    # now need to assign to chroms, taking care with overlap at ends. 
    # also need to decide whether to start with a feature or a gap.
    # could randomise that depending on proportions of genome covered by features to gaps e.g. if equal sized, 50:50
    
    #chromCums <- cumsum(seqlengths(genome))
    
    #thisChromCumValue <- 0  # to subtract from cumulative chrom positions, zero for first chrom.
    
    # for loop starts here
    #for(i in 1:length(chromCums)) {
    #thisChrom <- names(chromCums)[i]
    thisChromEndValue <- chromLength
    
    chromFirstIndex <- 1  
    
    
    chromFinalIndex <- max(which(ends < chromLength))  # which is the highest end value within this chrom
    
    chromGR <- GRanges(seqnames=thisChrom, seqinfo=seqinfo(genome), 
                       IRanges(start=starts[1:chromFinalIndex],
                               end=ends[1:chromFinalIndex]))
    #hist(width(chromGR))
    # now set the index along to use the next set of values.
    #chromFirstIndex <-chromFinalIndex +1
    
    allGR <- c(allGR, chromGR)
  }
 
  #hist(width(allGR))
  
  return(allGR)
  
}



