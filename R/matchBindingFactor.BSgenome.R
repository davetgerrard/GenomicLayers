# matchBindingFactor.BSgenome  
# attempt to match a binding factor over a whole genome
# Parameters:-
# layerSet,
# bindingFactor,
# clusterGap=10,
# max.window=10000000,
# verbose=FALSE
matchBindingFactor.BSgenome <- function(layerSet, bindingFactor, clusterGap=10, max.window=10000000, verbose=FALSE)  {
  
  require(Biostrings)
  stopifnot( class(layerSet$layerSet[[1]]) == "BSgenome")
  genome <- layerSet$layerSet[[1]]
  genome.sl <- seqlengths(genome)
  genome.starts <- rep(1, length(genome))
  names(genome.starts) <- seqnames(genome)
  genome.ends <- genome.sl
  #seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
  #max.window <- min(max.window, seqRange[2])
  hitList <- list()
  #validHits <-
  for(thisLayer in names(bindingFactor$profile)) {
    thisPattern <- bindingFactor$profile[[thisLayer]]$pattern
    max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
    if(thisLayer == "LAYER.0") {
      patternLength <- bindingFactor$profile[[thisLayer]]$length
      if(bindingFactor$type == "layer_region"  || bindingFactor$type == "layer_island" || (bindingFactor$type == "DNA_region" && length(grep("^N", bindingFactor$profile$LAYER.0$pattern)) >0 ) ) {     # lAYER.0 does not matter
        hitList[[thisLayer]] <- GRanges(seqnames(genome), IRanges(start=1, end=seqlengths(genome)),seqinfo=seqinfo(genome))
      } else {
        #if(verbose) print(paste("Sequence of length ", seqRange[2], ", using ",length(win.starts) ,"windows of length", max.window))
        all.hits <- vmatchPattern(thisPattern, genome, fixed=F) 
        hitList[[thisLayer]] <- reduce(all.hits, ignore.strand=TRUE)   # perhaps make strand-specific later.
        #hitList[[thisLayer]] <-  as(matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]], fixed=FALSE, max.mismatch= max.mismatches), "IRanges") # allows matching with IUPAC codes
      }
      #validHits <- hitList[[thisLayer]]
    } else {
      # with binary patterns, take ranges that have value 0 (gaps) 1 (IRanges)
      patternLength <- bindingFactor$profile[[thisLayer]]$length
      pattern <- bindingFactor$profile[[thisLayer]]$pattern
      if(bindingFactor$type == "layer_island") {
        pos <- layerSet[[thisLayer]][width(layerSet[[thisLayer]]) >= patternLength]
        neg <- gaps(layerSet[[thisLayer]])[width(gaps(layerSet[[thisLayer]]))>=patternLength]
        if(pattern==1) {
          island.index <- which(countOverlaps(pos,neg, maxgap = 1) > 1)
          hitList[[thisLayer]] <- pos[island.index]
          #print(paste("Layer island hits", length(island.index), "pos", pattern))
        } else {
          island.index <- which(countOverlaps(neg,pos, maxgap = 1) > 1)
          hitList[[thisLayer]] <- neg[island.index]
          #print(paste("Layer island hits", length(island.index), "neg", pattern))
        }
        #countOverlaps(pos,neg, maxgap = 1)   # values of 2 are positive regions overlapping two negative regions by 1 bp (i.e. on either side).

        # would need to be inverted for negative islands.
      } else {
        these.hits <- layerSet$layerSet[[thisLayer]][width(layerSet$layerSet[[thisLayer]]) >= patternLength]
        these.gaps <- gaps(layerSet$layerSet[[thisLayer]])[width(gaps(layerSet$layerSet[[thisLayer]]))>=patternLength]   # TODO include length of feature..
        if(pattern == 1) {
          hitList[[thisLayer]] <-these.hits
        } else {
          if (pattern == 0) {
            hitList[[thisLayer]] <-these.gaps
          } else{
            stop(paste("Unrecognised pattern", pattern))
          }
        }

        # hitList[[thisLayer]] <- ifelse(pattern == 1, these.hits, these.gaps) # this threw weird error: Error in NSBS(i, x, exact = exact, upperBoundIsStrict = !allow.append) :  subscript contains NAs or out-of-bounds indices
      }
          #validHits <- union(validHits, )
    }
    # trim the hitList to be within bounds for the sequence.
    #print(paste(thisLayer, class(hitList[[thisLayer]])))
    hitList[[thisLayer]] <- as(hitList[[thisLayer]], "GRanges")
    if(length(hitList[[thisLayer]]) > 0)   hitList[[thisLayer]] <- restrict(hitList[[thisLayer]] , start=genome.starts, end=genome.ends)
    # remove those shorter than patternLength (those overlapping the edges.
    hitList[[thisLayer]] <- hitList[[thisLayer]][width(hitList[[thisLayer]]) >= patternLength]

  }

  if(any(lapply(hitList, length) == 0))  {   # one or more of the patterns were not matched.
    validHits <- IRanges()
  } else{
    validHits <- hitList[[1]]
    for(i in 1:length(hitList)) {
      #overlaps <- findOverlaps(validHits, hitList[[i]])
      #validHits <- validHits[unique(queryHits(overlaps))]   # temp value to return
      validHits <- intersect(validHits, hitList[[i]])
    }

    #overlaps <- findOverlaps(hitList[[1]], hitList[[2]])
    #validHits <- hitList[[1]][unique(queryHits(overlaps))]   # temp value to return

  }
  # intersect hits to get proper valid hits.
  # Not sure how best to do this.
  # Don't expect matches to align perfectly
  # so just allow midpoints to within clusterGap?
  # --------+++++-----------------------
  # ----------++++++--------------------
  # ------------------------------------    # some parts will have no match


  # ?what to return
  return(validHits)

}