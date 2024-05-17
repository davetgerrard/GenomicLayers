

#  DEFUNCT see randGRangesBigGenome

# Generate random non-overlapping GRanges features within parameters of size and distance
# could not see an answer here https://www.biostars.org/p/225520/


# this version works for SacCer but not hg38 - does not like size of vector required to 
#  cover whole chromosomes. Also does not like integer overflow with cumsum
randGranges <- function(genome, n=NULL, sizeFunc=NULL, gapFunc=NULL, 
                        argsSizeFunc, argsGapFunc, minGapSize=1, minRange=1)  {
  
  # what would happen if gaps of size zero were allowed?
  # would we wnat items to be merged?  (that would probably break the even-ness 
  # of the chromosomal allocation if a merged region spanned a chrom boundary
  
  # must only determine two of the three inputs n, featureSize, gapSize because 
  # knowing two, determines the third
  
  # strategy:  create a linear set of features and gaps longer than the genome.
  # use cumsum to re-map these back onto the chromosomes
  genomeLength <- sum(seqlengths(genome)) 
  
  # decisions to make on whether it is preferable to start a chrom with a gap or a feature
  # should warn/error if any gaps or features are larger than half the size of the smallest chrom
  
  #test passed funcs
  print(do.call(sizeFunc, args=argsSizeFunc))
  
  sampleWidths <- rpois(300, 50)  # feature widths
  sampleGaps <- rpois(length(fWidths) +1, 400)  # gaps between features (+1 to get to end of chrom)
  
  m_sw <- mean(sampleWidths)
  m_sg <- mean(sampleGaps)
  
  #genomeLength/m_sw
  #genomeLength/m_sg
  minPairs <- genomeLength/(m_sw + m_sg)
  print(paste("Genome size: ", genomeLength, "Mean feature width: ", round(m_sw,2), "Mean gap width", round(m_sg,2),  "Calculated number of feature/gaps:", round(minPairs,2)))
  
  stopifnot( "requires more features than half genome length" = minPairs < genomeLength/2)
  stopifnot("fewer than one feature per chromosome" = minPairs > length(seqlengths(genome)))
  
  excessPairs <- ceiling(minPairs + length(seqlengths(genome)))  # add extra such that their is at least one pair extra per chrom
  gaps <- rpois(excessPairs , 400)
  widths <-   rpois(excessPairs, 50)
  #interleave the pairs into a single vector so that cumulative length of all can be cumsum
  interleaved <- integer(length=excessPairs*2)
  interleaved[seq(by=2, from=1, length.out=excessPairs)]  <- gaps   # gaps first
  interleaved[seq(by=2, from=2, length.out=excessPairs)]  <- widths
  interleavedCumsum <- cumsum(interleaved)
  
  # extract start and end points (relative to start of whole genome)
  ends <- interleavedCumsum[seq(by=2, from=2, length.out=length(interleavedCumsum))]
  starts <- (ends - widths) + 1
  # check results using alternate method
  #head(interleavedCumsum[seq(by=2, from=1, length.out=length(interleavedCumsum))] + 1)  # alt calculation.
  #head(gaps)
  #head(starts)
  #hist((ends - starts))   # should match original widths
  #hist(widths)
  
  # now need to assign to chroms, taking care with overlap at ends. 
  # also need to decide whether to start with a feature or a gap.
  # could randomise that depending on proportions of genome covered by features to gaps e.g. if equal sized, 50:50
  
  chromCums <- cumsum(seqlengths(genome))
  allGR <- GRanges(seqinfo=seqinfo(genome))
  chromFirstIndex <- 1  # set the first one to use the first feature. 
  thisChromCumValue <- 0  # to subtract from cumulative chrom positions, zero for first chrom.
  
  # for loop starts here
  for(i in 1:length(chromCums)) {
    thisChrom <- names(chromCums)[i]
    thisChromEndValue <- chromCums[thisChrom]
    
    chromFirstIndex <- min(which(starts > thisChromCumValue))  # which is the highest end value within this chrom
    
    
    chromFinalIndex <- max(which(ends < thisChromEndValue))  # which is the highest end value within this chrom
    
    chromGR <- GRanges(seqnames=thisChrom, seqinfo=seqinfo(genome), 
                       IRanges(start=starts[chromFirstIndex:chromFinalIndex] - thisChromCumValue,
                               end=ends[chromFirstIndex:chromFinalIndex] - thisChromCumValue))
    #hist(width(chromGR))
    # now set the index along to use the next set of values.
    #chromFirstIndex <-chromFinalIndex +1
    thisChromCumValue <- thisChromEndValue
    allGR <- c(allGR, chromGR)
  }
 
  #hist(width(allGR))
  
  return(allGR)
  
}
# 
# library(BSgenome.Scerevisiae.UCSC.sacCer3)
# library(BSgenome.Hsapiens.UCSC.hg38)
# genome <- BSgenome.Hsapiens.UCSC.hg38
# library(GenomicLayers)
# sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
# #sequences_to_keep <- "chr17"
# genomeNuc <- keepBSgenomeSequences(genome, sequences_to_keep)
# genome <- BSgenome.Scerevisiae.UCSC.sacCer3
# 
# randGranges(genome=genome)
# randGranges(genome=genomeNuc)
# 
# myFunc <- function(sizeFunc, gapFunc, argsSizeFunc, argsGapFunc) {
#   
#   print(do.call(sizeFunc, args=argsSizeFunc))
#   print(do.call(gapFunc, args=argsGapFunc))
# }
# 
# 
# 
# myFunc( sizeFunc = function(x, times) rep(x, times), 
#         gapFunc = function(x, times) rep(x, times), 
#         argsSizeFunc=list(x=10, times=5),
#         argsGapFunc = list(x=20, times=7))
# 

