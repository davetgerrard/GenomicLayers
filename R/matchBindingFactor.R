#' Find matches for a binding factor on a layer set
#'
#' Generate a list of matches for a binding factor against a layerSet object. 
#'
#' @param layerSet method to do something to (\code{"hsv"} or \code{"cluster"})
#' @param bindingFactor description of that param
#' @param clusterGap  =10 you get the idea
#' @param max.window    =10000000 you get the idea
#' @param verbose you get the idea
#'
#' @return \code{"hits"}
#'
#' @examples
#' x <- 1   # great!
#'
#' @import GenomicRanges
#' 
#'
#' @export
matchBindingFactor <- function(layerSet, bindingFactor, clusterGap=10, max.window=10000000, verbose=FALSE)  {
  require(Biostrings)
  seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
  max.window <- min(max.window, seqRange[2])
  hitList <- list()
  #validHits <-
  for(thisLayer in names(bindingFactor$profile)) {
    thisPattern <- bindingFactor$profile[[thisLayer]]$pattern
    max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
    if(thisLayer == "LAYER.0") {
      patternLength <- bindingFactor$profile[[thisLayer]]$length
      max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
      if(bindingFactor$type == "layer_region"  || bindingFactor$type == "layer_island" || (length(bindingFactor$profile$LAYER.0$pattern)==0) || (bindingFactor$type == "DNA_region" && length(grep("^N", bindingFactor$profile$LAYER.0$pattern)) >0 ) ) {     # lAYER.0 does not matter
        hitList[[thisLayer]] <- IRanges(start=seqRange[1], end=seqRange[2])
      } else {
        
        # currently DNA matches using fixed length string with IUPAC codes (adapt later for pwm matching)
        #TODO add in WINDOWING for long searches here.
        # for very long sequences >100k, need to break the sequence into sections, get results and concatenate them.
        win.starts <- seq(1, seqRange[2]-patternLength, by=max.window-patternLength)
        win.ends <- c(seq(max.window, seqRange[2], by=max.window), seqRange[2])
        if(length(win.starts) ==1)  win.ends <- win.ends[1]
        stopifnot(length(win.starts) == length(win.ends))
        if(verbose) print(paste("Sequence of length ", seqRange[2], ", using ",length(win.starts) ,"windows of length", max.window))
        all.hits <- IRanges()
        for(i in 1:length(win.starts)) {
          if(bindingFactor$type == "DNA_regexp" ) {
            grepResult <- gregexpr(bindingFactor$profile[[thisLayer]]$pattern, layerSet[[thisLayer]][win.starts[i]: win.ends[i]])
            if(grepResult[[1]][1] == -1 ) {  # no grep hits
              win.hits <- IRanges() 
            } else {
              win.hits <- IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE))
            }
          } else  {
            win.hits <-  as(matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]][win.starts[i]: win.ends[i]], fixed='subject', max.mismatch= max.mismatches), "IRanges")   # fixed='subject' ignores NNNs in subject (e.g. telomeres). See ?`lowlevel-matching` for more information.
          }
          win.hits <- shift(win.hits , win.starts[i] - 1)
          #print(win.hits)
          all.hits <- c(all.hits, win.hits)
        }
        hitList[[thisLayer]] <- reduce(all.hits)
        #hitList[[thisLayer]] <-  as(matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]], fixed=FALSE, max.mismatch= max.mismatches), "IRanges") # allows matching with IUPAC codes
        
      }
      #validHits <- hitList[[thisLayer]]
    } else {      # Not LAYER.0
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
        these.hits <- layerSet[[thisLayer]][width(layerSet[[thisLayer]]) >= patternLength]
        these.gaps <- gaps(layerSet[[thisLayer]])[width(gaps(layerSet[[thisLayer]]))>=patternLength]   # TODO include length of feature..
        if(length(these.hits) < 1) these.gaps <- IRanges(start=seqRange[1], end=seqRange[2])   # if no hits, whole chrom is a gap.
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
    }   # end of other LAYER names
    # trim the hitList to be within bounds for the sequence.
    #print(paste(thisLayer, class(hitList[[thisLayer]])))
    hitList[[thisLayer]] <- as(hitList[[thisLayer]], "IRanges")
    hitList[[thisLayer]] <- restrict(hitList[[thisLayer]] , start=seqRange[1], end=seqRange[2])
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