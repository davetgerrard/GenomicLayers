#' Find matches for a binding factor on a layer set containing a BSgenome sequence
#'
#' Generate a list of matches for a binding factor against a layerSet object containing a BSgenome sequence. 
#' For this (.BSgenome) form, the \code{"Hits"} object returned will be a \code{"GRanges"} object
#'
#' @param layerSet the \code{"layerSet"} target
#' @param bindingFactor the \code{"bindingFactor"} to match
#' @param match.layers  restrict matches to only these named layers (default: all layers in names(bindingFactor$profile))
#' @param clusterGap  =10 NOT IMPLEMENTED
#' @param max.window    =10000000 on less powerful computers, break up the search into windows of this size.
#' @param cache.layers NULL which named layers to cache hits on (default NULL). Only use for fixed matches (e.g. to genome sequence "LAYER.0")
#' @param verbose output more information to the screen
#'
#' @return \code{"GRanges"}
#' 
#' @seealso \code{\link{runLayerBinding}} \code{\link{runLayerBinding.BSgenome}} \code{\link{modifyLayerByBindingFactor.BSgenome}}
#'
#' @examples
#' require(Biostrings)
#' require(BSgenome.Scerevisiae.UCSC.sacCer3)
#' 
#' genome <- BSgenome.Scerevisiae.UCSC.sacCer3   # for convenience
#' 
#' 
#' scLayerSet <- createLayerSet.BSgenome(genome=genome, n.layers = 5, verbose=TRUE)
#' 
#' testFactor3 <- createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA", profile.layers = c("LAYER.1", "LAYER.3"), profile.marks = c(0,0), 
#'                                              mod.layers = c("LAYER.2", "LAYER.4"), mod.marks=c(0,1))
#' 
#' listOfHits <- matchBindingFactor.BSgenome(layerSet=scLayerSet, bindingFactor=testFactor3)
#'
#'
#' @import GenomicRanges
#' 
#' @export
matchBindingFactor.BSgenome <- function(layerSet, bindingFactor, match.layers=names(bindingFactor$profile),
                                        clusterGap=10, max.window=10000000,
                                        cache.layers=NULL, verbose=FALSE)  {
  
  require(Biostrings)
  # check input  
  stopifnot(exprs = {
        "LayerSet does not have a BSgenome in LAYER.0" = class(layerSet$layerSet[[1]]) == "BSgenome"
        "binding factor profile.layers matches layers not present in layerSet" = all(names(bindingFactor$profile) %in% names(layerSet$layerSet))
        "binding factor mod.layers list matches layers not present in layerSet" = all(names(bindingFactor$mods) %in% names(layerSet$layerSet))
    })

  # capture data from the genome
  genome <- layerSet$layerSet[[1]]
  genome.sl <- seqlengths(genome)
  genome.starts <- rep(1, length(genome))
  names(genome.starts) <- seqnames(genome)
  genome.ends <- genome.sl
  if(is.null(match.layers))  return( GRanges(seqinfo=seqinfo(genome)))   # in the case that the BF has no profile.
  #seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
  #max.window <- min(max.window, seqRange[2])
  hitList <- list() # list to store hits for each layer
  #validHits <-
  for(thisLayer in match.layers) {
    #thisPattern <- bindingFactor$profile[[thisLayer]]$pattern   # content will vary depending on bf type and layer. Could be regular expression?
    # calculate as integer the length of mismatch to tolerate.
    #max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
    if(!(thisLayer %in% cache.layers) | is.null(layerSet$cache[[bindingFactor$name]][[thisLayer]])) {  
      # this layer has no cache of hits for this factor or is not to be cached
      if(verbose) print(paste0("Finding matches for ",  bindingFactor$name, " on ", thisLayer)) 
      if(thisLayer == "LAYER.0") {
        patternLength <- bindingFactor$profile[[thisLayer]]$length
        if(bindingFactor$type == "layer_region"  || bindingFactor$type == "layer_island" || (bindingFactor$type == "DNA_region" && length(grep("^N", bindingFactor$profile$LAYER.0$pattern)) >0 ) ) {     # lAYER.0 does not matter
          hitList[[thisLayer]] <- GRanges(seqnames(genome), IRanges(start=1, end=seqlengths(genome)),seqinfo=seqinfo(genome))
        } else {    # bf type uses sequence
          if (bindingFactor$type == "DNA_regexp" )  {   # a regular expression on DNA e.g. "(CG.{0,20}){9}CG"
           all.hits <- matchRegExp.BSgenome(genome=genome, forRegExp = bindingFactor$profile[[thisLayer]]$forRegExp, revRegExp = bindingFactor$profile[[thisLayer]]$revRegExp)
            #  bsParams <- new("BSParams", X=genome, FUN=gregexpr)  # set up params for using bsapply
            # grepResultBS <- bsapply(bsParams, pattern=thisPattern)  # run gregexpr over each chromosome separately
            # all.hits <- GRanges(seqinfo=seqinfo(genome))   # empty GRanges to collate the results from different chroms
            # for(chromName in names(grepResultBS)) {  # cycle through chroms and convert matches to GRanges
            #   grepResult <- grepResultBS[[chromName]]
            #   if(grepResult[[1]][1] == -1 ) {  # no grep hits, do nothing
            #     #win.hits <- IRanges() 
            #   } else {
            #     win.hits <- GRanges(chromName, IRanges(start= as.integer(grepResult[[1]]), width = attr(grepResult[[1]], which="match.length", exact=TRUE)),seqinfo=seqinfo(genome))
            #   }
            #   all.hits <- c(all.hits, win.hits)
            # }
          } else {         # not a DNA_regexp
            #if(verbose) print(paste("Sequence of length ", seqRange[2], ", using ",length(win.starts) ,"windows of length", max.window))
            if(bindingFactor$type == "DNA_consensus")  {
              if(verbose) print("DNA_consensus class, using vmatchPattern()")
              # type DNA_consensus, use Biostrings vmatchPattern
              thisPattern <- bindingFactor$profile[[thisLayer]]$pattern   # content will vary depending on bf type and layer. Could be regular expression?
              # calculate as integer the length of mismatch to tolerate.
              max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
            all.hits <- vmatchPattern(thisPattern, genome, 
                                      fixed=bindingFactor$profile[[thisLayer]]$fixed,
                                      max.mismatch = bindingFactor$profile[[thisLayer]]$max.mismatch,
                                      min.mismatch = bindingFactor$profile[[thisLayer]]$min.mismatch, 
                                      with.indels= bindingFactor$profile[[thisLayer]]$with.indels, 
                                      algorithm =bindingFactor$profile[[thisLayer]]$algorithm) 
            } else { 
              if(bindingFactor$type == "DNA_motif") {
                # type DNA_motif, for now use Biostrings vmatchPattern
                # TODO change this to accept and use pwm
                if(verbose) print("DNA_motif class, using matchPWM()")
                all.hits <- matchPWM(pwm=bindingFactor$profile[[thisLayer]]$pattern,
                                     subject=genome,
                                    min.score=bindingFactor$profile[[thisLayer]]$min.score,
                                     with.score=bindingFactor$profile[[thisLayer]]$with.score)
                
                #
              } else {
              # unknown type of binding factor?
              stop(paste("Unknown binding factor type",bindingFactor$type ))
              }
            }
          }
          ###### ToDO make strand specific..  #######
          strand(all.hits) <- "*" # v3   - trying not to use reduce() and therefore keep adjacent hits distinct
          hitList[[thisLayer]] <- all.hits   # v3   - trying not to use reduce() and therefore keep adjacent hits distinct
          # v2 #hitList[[thisLayer]] <- reduce(all.hits, ignore.strand=TRUE)   # perhaps make strand-specific later.
          # v1 #hitList[[thisLayer]] <-  as(matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]], fixed=FALSE, max.mismatch= max.mismatches), "IRanges") # allows matching with IUPAC codes
        }   
        #validHits <- hitList[[thisLayer]]
      } else {
        # NOT LAYER.0
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
          # not a layer_island
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
      
    }  else  {  
      # attempt to use cached layers
      if(verbose) print(paste0("Using cached hits for " , bindingFactor$name, " on layer ", thisLayer))
      if(is.null(layerSet$cache[[bindingFactor$name]][[thisLayer]])) {
        stop(paste("No cache available for",  bindingFactor$name, "on layer", thisLayer))
      } else {
        hitList[[thisLayer]]  <- layerSet$cache[[bindingFactor$name]][[thisLayer]]
      }
      
      
    }
  }  # end of thisLayer loop 
  
  # all layers have been matched or their cached results added to hitList. 
  # now need to intersect the results. 
  # first test if any required layers are empty as this will mean there are no hits overall. 
  # Remember that profiles looking for '0' will be inverse matches.
  if(any(lapply(hitList, length) == 0))  {   # one or more of the patterns were not matched.
    validHits <- GRanges()    #####TODO change IRanges  to GRAnges()  to fix the intersect below.
  } else{
    # 2-step strategy to combine hits. 
    # Need to keep sequence based matches distinct (even if adjacent or overlapping) as these
    # are all equally valid binding sites and represent availability of such sites.
    # currently, GRanges::intersect()  causes adjacent or overlapping ranges to be reduced. e.g. 1-3 & 4-6 becomes 1-6.
    # Step 1. Combine all layers EXCEPT LAYER.0 (sequence) using intersect.
    # Step 2. Subset the LAYER.0 hits to those overlapped within the output from step 1
    ## N.B. check this works when LAYER.0 does not matter. Or run separate layer binding in this case? 
    #       this would be true for bf type "layer_region"  but what if a regular expression matches whole chrom?
    validHits  <- GRanges(seqnames(genome), IRanges(start=1, end=seqlengths(genome)),seqinfo=seqinfo(genome))  # start with whole genome valid
    for(thisLayer in names(hitList))  {
      if(thisLayer == "LAYER.0") {    # sequence layer don't yet combine these.
      } else {
        validHits <- intersect(validHits, hitList[[thisLayer]]) 
      }
    }
    # now add back in the sequence layer matches, keeping each distinct if that is required.
    if(bindingFactor$type == "DNA_regexp" | bindingFactor$type == "DNA_consensus" | bindingFactor$type == "DNA_motif" ) {   # LAYER.0 matches are important and should be kept distinct.
      # subset seq based LAYER.0 hits. Must be entirely "within" the valid hits generated by other layers.
      validHits <- hitList[["LAYER.0"]][overlapsAny(hitList[["LAYER.0"]], validHits, type="within")]  
    } else {  # LAYER.0 matches are not important and can be intersected and reduced
      # generate slidingWindowHits of length patternLength within every valid hit.
      validHits <- unlist(slidingWindows(validHits, width=patternLength, step=1))
      #slidingWindows can truncate the last window if it does not fit
      validHits <- validHits[width(validHits) >= patternLength]
      #validHits <- intersect(validHits, hitList[["LAYER.0"]]) 
    }
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