# December 2015 
# Decided that base objects are not use-able for whole chromosomes.
# Need to re-write using hits objects

# function 
# description
# Parameters:-
# factorSet
#createBindingFactor <-  function()  {
#  bindingFactor <- list(name=name, type=type, 
#                        profile=profileList, 
 #                       mods=modList)
#  
#  return(bindingFactor)
#}

# createRandomBindingFactor 
# description
# Parameters:-
# name, 
# layerSet, 
# type=c("DNA_motif", "DNA_region","layer_region","layer_island"),
# test.layer0.binding=FALSE, 
# test.mismatch.rate=.1 , 
# max.pattern.tries=1000, 
# min.DM.length=2, 
# min.DR.length=10
createRandomBindingFactor <- function(name, layerSet, type=c("DNA_motif", "DNA_region","layer_region","layer_island"),
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  if(!test.layer0.binding) max.pattern.tries <- 1  # only try one pattern.
  
  if(type == "DNA_motif") {
    
    for(i in 1:max.pattern.tries) {
      patternLength <- max(min.DM.length,rnbinom(1, 50, mu= 12))  # patterns must be at least of length 1
      pattern <- paste(sample(names(IUPAC_CODE_MAP), patternLength, replace=T), collapse="")
      max.mismatches <- round(test.mismatch.rate * patternLength )
      if(test.layer0.binding)  {
        matches.length <- length( matchPattern(pattern, layerSet[['LAYER.0']], max.mismatch=max.mismatches, fixed="subject" ))
        
      }else {
        matches.length <- 1
      }
      if(matches.length > 0)  {
        # print(i)
        break; 
      }
      
    }
    if(test.layer0.binding & verbose) print(paste(pattern , "matches training sequence" , matches.length, "times"))    
    
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0, length=patternLength))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    #openPattern <- paste(rep("1", patternLength  ), collapse="")
    #closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      profileList[[thisLayer]] <- list(pattern=sample(c(0,1),1), mismatch.rate=0.1, length=patternLength)
    }
    
    modList <- list()
    n.modPatterns <- sample(1:(length(layerSet)-1), 1, prob=1/(1:(length(layerSet)-1)), replace=T)   # higher numbers successively less likely, zero not allowed
    for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      modState <- as.character(sample(0:1,1))
      modList[[thisLayer]] <- list(state=modState, stateWidth=patternLength, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
    }
  }   
  
  if(type == "DNA_region") {
    
    for(i in 1:max.pattern.tries) {
      patternLength <- max(10, rnbinom(min.DR.length, 10, mu= 50))  # min 10 long
      pattern <- paste(rep(sample(names(IUPAC_CODE_MAP), 1), patternLength), collapse="")
      max.mismatches <- round(test.mismatch.rate * patternLength )
      if(test.layer0.binding)  {
        matches.length <- length( matchPattern(pattern, layerSet[['LAYER.0']], max.mismatch=max.mismatches, fixed=FALSE ))
      }else {
        matches.length <- 1
      }
      if(matches.length > 0)  {
        # print(i)
        break; 
      }
      
    }
    if(test.layer0.binding & verbose) print(paste(pattern , "matches training sequence" , matches.length, "times"))    
    
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1, length=patternLength))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    #openPattern <- paste(rep("1", patternLength  ), collapse="")
    #closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      profileList[[thisLayer]] <- list(pattern=sample(c(0,1),1), mismatch.rate=0.1, length=patternLength)
    }
    
    modList <- list()
    n.modPatterns <- sample(1:(length(layerSet)-1), 1, prob=1/(1:(length(layerSet)-1)), replace=T)   # higher numbers successively less likely, zero not allowed
    for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      modState <- as.character(sample(0:1,1))
      modList[[thisLayer]] <- list(state=modState, stateWidth=patternLength, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
    }
  } 
  
  if(type == "layer_region") {
    patternLength <- max(1,rnbinom(1, 10, mu= 200))
    #pattern <- paste(rep("N", patternLength), collapse="")  # caused too much memory usage.
    pattern <- "N"
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1, length=patternLength))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    #openPattern <- paste(rep("1", patternLength  ), collapse="")
    #closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      profileList[[thisLayer]] <- list(pattern=sample(c(0,1),1), mismatch.rate=0.1, length=patternLength)
    }
    
    modList <- list()
    n.modPatterns <- sample(1:(length(layerSet)-1), 1, prob=1/(1:(length(layerSet)-1)), replace=T)   # higher numbers successively less likely, zero not allowed
    for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      modState <- as.character(sample(0:1,1))
      modList[[thisLayer]] <- list(state=modState, stateWidth=patternLength, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
    }
  } 
  
  if(type == "layer_island") {
    islandLength <- max(1,rnbinom(1, 10, mu= 100))
    shoulderLength <- max(1,rnbinom(1, 10, mu= 100))
    patternLength <- islandLength + (2*islandLength)
    #pattern <- paste(rep("N", patternLength), collapse="")  # caused too much memory usage.
    pattern <- "N"
    
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1, length=patternLength))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    #openPattern <- paste(rep(c("0","1","0"), times=c(shoulderLength, islandLength, shoulderLength)  ), collapse="")
    #closedPattern <- paste(rep(c("1","0","1"), times=c(shoulderLength, islandLength, shoulderLength)  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      profileList[[thisLayer]] <- list(pattern=sample(c(0,1),1), mismatch.rate=0.1, length=islandLength, shoulder=shoulderLength)
    }
    
    modList <- list()
    n.modPatterns <- sample(1:(length(layerSet)-1), 1, prob=1/(1:(length(layerSet)-1)), replace=T)   # higher numbers successively less likely, zero not allowed
    for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      modState <- as.character(sample(0:1,1))
      modList[[thisLayer]] <- list(state=modState, stateWidth=islandLength, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
    }
  } 
  
  bindingFactor <- list(name=name, type=type, 
                        profile=profileList, 
                        mods=modList)
  
  return(bindingFactor)
}


# print.bfSet 
# description
# Parameters:-
# factorSet 
print.bfSet <- function(factorSet) {
  print("A list of binding factors:-")
  bf.names <- names(factorSet)
  types <- lapply(factorSet, FUN=function(x) x$type)
  dna.patterns <- unlist(lapply(factorSet, FUN=function(x) as.character(x$profile$LAYER.0$pattern)))
  pattern.lengths <- unlist(lapply(factorSet, FUN=function(x) as.character(x$profile$LAYER.0$length)))
  layer.patterns <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$profile))), paste, collapse=","))
  layer.mods <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$mods))), paste, collapse=","))
  as.data.frame(cbind(bf.names, types,pattern.lengths, dna.patterns,  layer.patterns, layer.mods) )
}

# matchBindingFactor 
# description
# Parameters:-
# layerSet, 
# bindingFactor, 
# clusterGap=10, 
# max.window=10000000, 
# verbose=FALSE
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
      if(bindingFactor$type == "layer_region"  || bindingFactor$type == "layer_island" || (bindingFactor$type == "DNA_region" && length(grep("^N", bindingFactor$profile$LAYER.0$pattern)) >0 ) ) {     # lAYER.0 does not matter
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
          win.hits <-  as(matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]][win.starts[i]: win.ends[i]], fixed='subject', max.mismatch= max.mismatches), "IRanges")   # fixed='subject' ignores NNNs in subject (e.g. telomeres). See ?`lowlevel-matching` for more information.
          win.hits <- shift(win.hits , win.starts[i] - 1)
          #print(win.hits)
          all.hits <- c(all.hits, win.hits)
        }
        hitList[[thisLayer]] <- reduce(all.hits)
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
        these.hits <- layerSet[[thisLayer]][width(layerSet[[thisLayer]]) >= patternLength]
        these.gaps <- gaps(layerSet[[thisLayer]])[width(gaps(layerSet[[thisLayer]]))>=patternLength]   # TODO include length of feature..
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
      overlaps <- findOverlaps(validHits, hitList[[i]])
      validHits <- validHits[unique(overlaps@queryHits)]   # temp value to return
    }
    
    #overlaps <- findOverlaps(hitList[[1]], hitList[[2]])
    #validHits <- hitList[[1]][unique(overlaps@queryHits)]   # temp value to return
    
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



# Layers are now Views or Iranges objects. marks (1) are contiguous ranges, absence of marks (0) are gaps between.
# hits is now a Views object
# TODO, edit to alter layerset IN-PLACE i.e. modifyLayerByBindingFactor.Views(layerSet, position.vec,bindingFactor)
# function 
# description
# Parameters:-
# layerSet, 
# hits, 
# bindingFactor, 
# verbose=FALSE
modifyLayerByBindingFactor.Views <- function(layerSet, hits, bindingFactor, verbose=FALSE) {
  require(Biostrings)
  newLayerSet <- layerSet
  for(thisLayer in names(bindingFactor$mods))  {
    
    
    seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
    thisState <- bindingFactor$mods[[thisLayer]]$state
    
    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    hits <- resize(hits, width=stateWidth, fix="center")    # adjust the width
    thisOffset <- bindingFactor$mods[[thisLayer]]$offset
    hits <- shift(hits, shift=thisOffset)                 # move up- or downstream
    
    # restrict hits to range
    
    
    newLayerSet[[thisLayer]] <- switch(thisState, 
                                       "1" = reduce(union(newLayerSet[[thisLayer]], hits)), 
                                       "0" = setdiff(newLayerSet[[thisLayer]], hits),
                                       stop(paste("unknown state", thisState)))


  }
  return(newLayerSet)
}


#runLayerBinding <- function() {

# runLayerBinding 
# description
# Parameters:-
# layerList, 
# factorSet, 
# iterations=1, 
# bindingFactorFreqs=rep(1, length(factorSet)), 
# watch.function=function(x){}, 
# collect.stats=FALSE, 
# target.layer=2, 
# verbose=FALSE
runLayerBinding <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE, ...)  {
  if(verbose) print(paste(Sys.time(), "runLayerBinding pos 1", sep=" "))
  #bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  bindingOrder <- names(factorSet)  # JUST USE EACH FACTOR ONCE, IN ORDER GIVEN
  newLayerList <- layerList
  seqRange <- c(start(layerList$layerSet[['LAYER.0']])[1], end(layerList$layerSet[['LAYER.0']])[1])
  if(collect.stats) {
    stats.table <- data.frame()
  } else {
    stats.table <- NULL
  }  
  max.hits <- ceiling(iterations/length(factorSet))  # TODO could tailor this to be different for each factor.
  for(thisBF in bindingOrder)  {
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast thisBF =", thisBF, factorSet[[thisBF]]$profile$LAYER.0$pattern, sep=" "))
    theseHits <- matchBindingFactor(newLayerList$layerSet, factorSet[[thisBF]], verbose=verbose)  
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits =", length(theseHits), sep=" "))
	#print(length(theseHits))
    if(length(theseHits) < 1) {
      if(verbose) print(paste(thisBF, "0 hits, skipping"))
      next ;
    }
    
    # how many of the potential hits to mark? 
    # iterations/n.factors (rounded up).
    # need to check if hits object is a single hit spanning whole chrom (e.g. for DNA_region matching 'N' and nothing else)
    # TODO: if so, need to make pseudo-hits that can be sampled.
    if(identical(theseHits, IRanges(start=seqRange[1], end=seqRange[2]))) {
      if(verbose) print("Hits match whole chromosome")
      hit.width <- factorSet[[thisBF]]$profile$LAYER.0$length
      starts <- sample(1:(seqRange[2]-hit.width), min(max.hits,seqRange[2]-hit.width), replace=FALSE)
      hits.sample <- IRanges(starts, starts+hit.width)
    } else {
      hits.sample <- theseHits[sample(1:length(theseHits) ,min(length(theseHits),max.hits))]   # now multiple
    #if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits.used =", length(hits.sample), sep=" "))
    }
    if(verbose) {cat(paste(thisBF, length(theseHits), length(hits.sample)))
    } else { cat( ".")}
	#thisHitPosition <- start(hits.sample) + floor(width(hits.sample)/2)
    
    newLayerList$layerSet <-  modifyLayerByBindingFactor.Views(newLayerList$layerSet, hits=hits.sample, bindingFactor=factorSet[[thisBF]], verbose=verbose)
    watch.function(x= newLayerList$layerSet, ...)
    if(collect.stats) {
      coverages <- unlist(lapply(newLayerList$layerSet[-1], FUN=function(x) sum(width(x))))
      block.counts <- unlist(lapply(newLayerList$layerSet[-1], length))
      #hit.counts <- 
      thisRow <- data.frame(bf=thisBF,hits=length(hits.sample) , target.coverage=sum(width(hits.sample)))
      thisRow[, paste("Coverage.", names(coverages), sep="")] <- coverages
      thisRow[, paste("nBlocks.", names(block.counts), sep="")] <- block.counts
      stats.table <- rbind(stats.table, thisRow)
    }
    #print(as.numeric(letterFrequency(newLayerSet$LAYER.1, letters= "1")))
  }
  newLayerList$history <- stats.table
  #print(letterFrequency(newLayerSet$LAYER.1, letters= "1"))   # how many of layer.1 were set to 1
  if(verbose) {print(paste(Sys.time(), "runLayerBinding.fast pos 2", sep=" ")) } else {print("")}
  return(newLayerList)
}

# mutateFactorSet 
# description
# Parameters:-
# factorSet, 
# layerSet , 
# mut_type="subRandomFactor", 
# n.muts=1, 
# verbose=FALSE, 
# test.layer0.binding=FALSE,  
# type_list = c("subRandomFactor" , "duplicate", "switch"),   "insert", "delete", "translocate", "move"
# prob=rep(1/length(type_list), length(type_list))
# fix.set.size= TRUE,   Keep the factor set the same size. 
# name.prefix=""    give new factors a new name beginning with this name 
#
# TODO re-write this so that a mixed set of mutation types can be generated.
mutateFactorSet <- function(factorSet, layerSet , mut_type="subRandomFactor", n.muts=1, verbose=FALSE, test.layer0.binding=FALSE,  
                            type_list = c("subRandomFactor" , "duplicate", "switch"), prob=rep(1/length(type_list), length(type_list)), fix.set.size= TRUE, name.prefix=""){
  
  newFactorSet <- factorSet
  if(mut_type == "random")  {
   mut_type <- sample( type_list, 1, prob=prob)
  }

  
  if(mut_type== "duplicate")  {
    for(j in 1:n.muts)  {
      copy <- sample(1:length(newFactorSet), 1)  # which factor to duplicate
      index <- sample(1:length(newFactorSet), 1)  # what position to insert
      position  <- sample(c("before", "after"), 1)
      newName <- paste("dup", copy, sep=".")
      newFactor <- newFactorSet[[copy]]
      newFactor$name <-  newName
      # first join the newFactor 
      append.index <- ifelse(position=="after", index, index-1)   # subtract one if "before"
      newFactorSet <- append(newFactorSet, list(newFactor), after=append.index)
      names(newFactorSet)[append.index+1] <- newName

      if(verbose)  {
        print(paste("Factor", copy,  "duplicated to ", position , "position", index ))
      }
    }
    
  }
  
  if(mut_type== "insert")  {
    for(j in 1:n.muts)  {
    index <- sample(1:length(newFactorSet), 1)  # what position to insert
    position  <- sample(c("before", "after"), 1)
    #for(i in index) {
    newName <- paste("ins", index, sep=".")
    newType <- sample(c("DNA_motif", "DNA_region","layer_region","layer_island"),1, replace=T)
    newFactor <- createRandomBindingFactor(newName,layerSet, type=newType, test.layer0.binding=test.layer0.binding, test.mismatch.rate=.1, verbose=verbose ) 
    #if(verbose) print(newFactor)
    # first join the newFactor 
    append.index <- ifelse(position=="after", index, index-1)   # subtract one if "before"
    #fS.l <- length(factorSet)
    #append.vec <- append(1:fS.l, fS.l+1, after=append.index)
    newFactorSet <- append(newFactorSet, list(newFactor), after=append.index)
    names(newFactorSet)[append.index+1] <- newName
    #newFactorSet <- newFactorSet[c(append.vec)]  # reorder
      
    if(verbose)  {
        print(paste("Factor inserted", position , "position", index ))
    }
    #}
    }
  }
  
  if(mut_type== "delete")  {
    index <- sample(1:length(factorSet), n.muts)
    #rem.names <-   # when rework factor set, will need to remove the corresponding abundances (or only keep those we keep.)
    newFactorSet <- newFactorSet[-c(index)]
    if(verbose)  {
      print(paste("Factor(s)", paste(index, collapse=",") ,"deleted"))
    }   
    
  }  
  
  if(mut_type== "subRandomFactor")  {
    index <- sample(1:length(factorSet), n.muts)
    for(i in index) {
      thisName <- names(factorSet)[i]
      newName <- ifelse(name.prefix == "", thisName, paste(name.prefix, i, sep="."))
      newFactorSet[[thisName]] <-  createRandomBindingFactor(newName,layerSet, type=factorSet[[i]]$type, test.layer0.binding=test.layer0.binding, test.mismatch.rate=.1, verbose=verbose ) 
      if(verbose)  {
        print(paste("Factor", i ,"substituted"))
      }
    }
    
  }
  
  if(mut_type =="switch")  {  # switch a factor with another from the set
    index <- sample(1:length(factorSet), n.muts)
    for(i in index) {
      partner <- sample(setdiff(1:length(factorSet), i), 1)
      new.index <- 1:length(factorSet)
      new.index[c(i, partner)] <- c(partner, i)
      newFactorSet <- newFactorSet[new.index]
      if(verbose)  {
        print(paste("Factors", i ,"and", partner, "switched"))
      }
    }
    
  }
  
  
  return(newFactorSet)
}



# optimiseFactorSet 
# description
# Parameters:-
# layerList, 
# factorSet, 
# testing.function, 
# target.layer, 
# target.vec, 
# n.iter=10, 
# target.score=1, 
# mut.rate=0.1, 
# modsPerCycle=100,
# test.layer0.binding=FALSE, 
# method="fast", 
# logFile="",
# logCycle=10, record scores every 'logCycle' iterations
# maxNoChange=n.iter, 
# verbose=FALSE, 
# use.parallel=FALSE, 
# n.cores=1
# max.factors=50    # max number of factors allowed. (when reached, removes "insert" and "duplicate" mut_types
optimiseFactorSet <- function(layerList, factorSet, testing.function, target.layer, target.vec, n.iter=10, target.score=1, mut.rate=0.1, modsPerCycle=100,
                              test.layer0.binding=FALSE, method="fast", logFile="",logCycle=10, maxNoChange=n.iter, verbose=FALSE, use.parallel=FALSE, n.cores=1, 
                              max.factors=50, mut_type_list = c("subRandomFactor" ,  "switch", "duplicate" , "insert", "delete"))  {
  
  if(logFile != "")  write.table(cbind("iter", "best.score"), row.names=F, col.names=F, sep="\t", quote=F,file=logFile)
  
  target.count <- length(target.vec)
  target.coverage <- sum(width(target.vec))
  chrom.size <- nchar(layerList$layerSet$LAYER.0)
  currentFactorSet <- factorSet
  scores.vector <- numeric()
  print("Calculating initial scores")
  if(method =="fast") {
    print("using fast algorithm")
    initialModLayer <- runLayerBinding(layerList=layerList, factorSet = factorSet, iterations=modsPerCycle, verbose=verbose)
  } else {
    print("using slow algorithm")
    initialModLayer <- runLayerBinding(layerList=layerList, factorSet = factorSet, iterations=modsPerCycle, verbose=verbose)
  }
  
  initialScore <- testing.function(initialModLayer, targetLayer=target.layer, target.vec=target.vec)
  if(is.na(initialScore )) {
    stop("initialScore was NA")
    
  }
  
  currentBestScore <- initialScore
  if(logFile != "") write.table(cbind(0, currentBestScore), row.names=F, col.names=F, quote=F,sep="\t",file=logFile, append=TRUE)
  iSinceLastImprovement <- 0  # use to break from following loop if very many rounds without improvement
 
  if(use.parallel) {
    # attempt at parallel runs
    require(parallel)
    mc <- getOption("cl.cores", n.cores)
    if(verbose) print(paste("Cores to use:" ,mc))
    
    cl <- makeCluster(mc)
    # TEMP need to export all functions to the nodes but don't know what they are.
    func.vec <- sapply(ls(.GlobalEnv), FUN =function(x) class(get(x)))
    func.names <- names(func.vec[func.vec == "function"])
    if(verbose) print(func.names)
    # TODO package code to avoid above
    clusterExport(cl, c( func.names, "layerList", 
                         "modsPerCycle", "verbose", "target.layer","target.vec"), envir = environment())
    n.tries <- 10
    
  }
  
  
  for(i in 1:n.iter) {
    # write to log file if appropriate
    if((i %% logCycle ) == 0  & logFile != "") {
      write.table(cbind(i, currentBestScore), row.names=F, col.names=F, quote=F,sep="\t",file=logFile, append=TRUE)
      # could also save the current best factor set at this point.
      save(currentFactorSet, file=paste(dirname(logFile), "/", "currentFactorSet.",i,".Rdata",sep=""))
    }
    if(iSinceLastImprovement > maxNoChange)  {
      print(paste("No improvement in ", iSinceLastImprovement, " iterations, exiting!", sep=""))
      break ; 
    }
    
    if(length(currentFactorSet) >= max.factors)  {   # stop the factor set growing larger than this
      use_mut_types <- setdiff(mut_type_list, c("insert", "duplicate"))
      if(verbose)  print(paste("Factor set limited to max.factors: ", max.factors))
    } else {
      use_mut_types <- mut_type_list
      
    }
    
    better <- FALSE
    if(use.parallel) {
      # create n.cores new factorSet by mutating the currentSet.
      factorSet.list <- replicate(mc, mutateFactorSet(currentFactorSet, layerList$layerSet,type_list = use_mut_types, mut_type= "random", n.muts=floor(length(currentFactorSet) * mut.rate), verbose=verbose, test.layer0.binding=T, name.prefix=paste("m",i,sep="")), simplify=FALSE)
      
      # list(...) added when cluster not picking up layerList
      #n.tries <- 10  # VERY WEIRDLY , the below command sometimes fails on first call but works on second (or later).  TODO: fix this!
      tries=0
      while(tries <= n.tries)  {   
        worked <- try(mod.list <- parLapply(cl, factorSet.list, fun=function(x) runLayerBinding(layerList=layerList, factorSet = x, iterations=modsPerCycle, verbose=verbose)), TRUE)
        tries <- tries +1
        if(class(worked) != "try-error") {
          if(verbose) print(paste("parallel mod list tries: " , tries))
          break;
        }
      }
      if(verbose) print(paste(tries, "of", n.tries, "attempted"))
      #clusterExport(cl, "mod.list", envir = environment())
      par.scores <- parSapply(cl, mod.list, FUN=function(x) testing.function(layerList=x, targetLayer=target.layer, target.vec=target.vec))
      print(paste("x", mc, sep=""))
      
      
      best.index <- which.max(par.scores)
      newFactorSet <- factorSet.list[[best.index]]
      newModLayer <- mod.list[[best.index]]
      newScore <- par.scores[best.index]
    } else {
      newFactorSet <- mutateFactorSet(currentFactorSet, layerList$layerSet,n.muts=floor(length(currentFactorSet) * mut.rate),type_list = use_mut_types, mut_type= "random",  verbose=verbose, test.layer0.binding=test.layer0.binding, name.prefix=paste("m",i,sep=""))
      if(method =="fast") {
        #print("using fast algorithm")
        newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
      } else {
        #print("using slow algorithm")
        newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
      }
      newScore <- testing.function(newModLayer, targetLayer=target.layer, target.vec=target.vec)
      
    }
    
    n.regions <- length(newModLayer$layerSet[[target.layer]])
    coverage <- sum(width(newModLayer$layerSet[[target.layer]]))
    regionsWithHit <- sum(overlapsAny(newModLayer$layerSet[[target.layer]], target.vec))
    targetsWithHit <- sum(overlapsAny(target.vec, newModLayer$layerSet[[target.layer]]))                      
    
    
    scores.vector[i] <- newScore
    print(paste("Round", i, ". Facors:", length(newFactorSet),". Marks on target layer:", n.regions, ", Coverage:", coverage, ", Regions with a hit:", regionsWithHit,
                ", Targets Hit:", targetsWithHit, ", Chrom size:", chrom.size, ", Target count:", target.count, ", Target coverage:", target.coverage))
    print(paste("Round", i, ". OldScore", currentBestScore, "NewScore", newScore, "Better?"))
    if(is.na(newScore)) {
      print("Newscore = NA, skipping to next")
      next ;
    }
    if(newScore > currentBestScore) {
      better <- TRUE
      currentBestScore <- newScore
      currentFactorSet <- newFactorSet
      iSinceLastImprovement <- 0 
      print("Yes!")
    } else {
      iSinceLastImprovement <- iSinceLastImprovement + 1 
      print("No!")
    }
    #print(paste("Round", i, "oldScore", currentBestScore, "newScore", newScore, "Better?",better))
    
  }
  if(use.parallel)   stopCluster(cl)
  currentFactorSet$optimScores <- scores.vector
  return(currentFactorSet)
}



# create a matrix of nucleotide frequencies in the LAYER.0 patterns of the factorSet 
af.factorSet <- function(factorSet) {
  for(i in 1:length(factorSet)) {
    
    if(i == 1)  {
      af.table <-  alphabetFrequency(factorSet[[1]]$profile[['LAYER.0']]$pattern)
    }else {
      af.table <- rbind(af.table,alphabetFrequency( factorSet[[i]]$profile[['LAYER.0']]$pattern))
      
      
    }
    
    
  }
  return(af.table)
}
# colSums(af.factorSet(result[1:length(result)-1]))  # quite useful to see composition of bfs

plot.layerBindingHistory <- function(modLayerSet ) {
  
  coverage.cols <-  sort(grep("Coverage", names(modLayerSet$history), value=T))
  blockCount.cols <-  sort(grep("nBlocks", names(modLayerSet$history), value=T))
  par(mfrow=c(1,2))
  matplot(1:nrow(modLayerSet$history),modLayerSet$history[,coverage.cols], ylab="Coverage", xlab="Binding factor")
  matplot(1:nrow(modLayerSet$history),modLayerSet$history[,blockCount.cols], ylab="Number of blocks", xlab="Binding factor")
  
}

#plot.layerBindingHistory(modLayerSet)

# show the profile and modification specificities for a set of binding factors. 
plot.factorSet <- function(factorSet)  {
  
  layerMap <- data.frame()
  
  profileLayers <- character()
  modLayers <- character()
  
  for(i in 1:length(factorSet)) {
    profileLayers <- sort(union(profileLayers, names(factorSet[[i]]$profile)))
    modLayers <-  sort(union(modLayers, names(factorSet[[i]]$mods)))
  }
  # could just set modLayers <- profileLayers ?
  
  profile.layerMap <- matrix(NA, ncol=length(profileLayers), nrow=length(factorSet), dimnames=list(factor=1:length(factorSet), layer=profileLayers ))
  mod.layerMap <- matrix(NA, ncol=length(modLayers), nrow=length(factorSet), dimnames=list(factor=1:length(factorSet), layer=modLayers ))
  
  for(i in 1:length(factorSet)) {
    
    
    
    for(thisLayer in names(factorSet[[i]]$profile)) {
      if(thisLayer == "LAYER.0") {
        profile.layerMap[i,thisLayer] <- 1
      } else {
        profile.layerMap[i,thisLayer] <- factorSet[[i]]$profile[[thisLayer]]$pattern
      }
      
    }
    
    
    for(thisLayer in names(factorSet[[i]]$mods)) {
      #mod.layerMap[i,thisLayer] <- TRUE  
      mod.layerMap[i,thisLayer] <- as.integer(factorSet[[i]]$mods[[thisLayer]]$state)
      
    }
    
  }
  
  par(mfrow=c(1,2))
  image(t(profile.layerMap), main="profile", axes=F, xlab="Layer") ; box()
  mtext(row.names(profile.layerMap), side=2, at=seq(0,1, length.out=nrow(profile.layerMap)), las=1)
  mtext(0:(ncol(profile.layerMap)-1), side=1, at=seq(0,1, length.out=ncol(profile.layerMap)), las=1)
  image(t(mod.layerMap), main="mods", axes=F, xlab="Layer"); box()
  mtext(row.names(mod.layerMap), side=2, at=seq(0,1, length.out=nrow(mod.layerMap)), las=1)
  mtext(1:ncol(mod.layerMap), side=1, at=seq(0,1, length.out=ncol(mod.layerMap)), las=1)
}
# plot.factorSet(factorSetRandom)


# plot multiple score measures. Green for good, red for bad.
plot.score.hits <- function(query, target, methods=c("acc", "tpr","tnr", "ppv", "npv", "fpr", "fnr"), 
                            colours=structure(c("green", "green", "green", "green", "green", "red", "red"), 
                                              names=c("acc", "tpr","tnr", "ppv", "npv", "fpr", "fnr")))  {
  scores <- numeric()
  for(i  in 1:length(methods)) {
    thisMethod <- methods[i]
    scores[thisMethod] <- score.hits(query=query, target = target, method = get(methods[i]))
  }
  #return(scores)
  barplot(scores, ylim=c(0,1), col=colours[methods])
}

# plot.score.hits(query=modLayerSet$layerSet$LAYER.5, target = tss.IR)