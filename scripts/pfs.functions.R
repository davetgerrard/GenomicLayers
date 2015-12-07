# December 2015 
# Decided that base objects are not use-able for whole chromosomes.
# Need to re-write using hits objects

createBindingFactor <-  function()  {
  bindingFactor <- list(name=name, type=type, 
                        profile=profileList, 
                        mods=modList)
  
  return(bindingFactor)
}

# name, layerSet, 
# type=c("DNA_motif", "DNA_region","layer_region","layer_island"),
# test.layer0.binding=FALSE, 
# test.mismatch.rate=.1 , 
# max.pattern.tries=1000, 
# min.DM.length=2, 
# min.DR.length=10
createRandomBindingFactor <- function(name, layerSet, type=c("DNA_motif", "DNA_region","layer_region","layer_island"),
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, min.DM.length=2, min.DR.length=10) {
  
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
    if(test.layer0.binding) print(paste(pattern , "matches training sequence" , matches.length, "times"))    
    
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
    if(test.layer0.binding) print(paste(pattern , "matches training sequence" , matches.length, "times"))    
    
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


#
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
runLayerBinding <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE)  {
  if(verbose) print(paste(Sys.time(), "runLayerBinding.fast pos 1", sep=" "))
  #bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  bindingOrder <- names(factorSet)  # JUST USE EACH FACTOR ONCE, IN ORDER GIVEN
  newLayerList <- layerList
  
  if(collect.stats) {
    stats.table <- data.frame()
  } else {
    stats.table <- NULL
  }  
  max.hits <- ceiling(iterations/length(factorSet))  # TODO could tailor this to be different for each factor.
  for(thisBF in bindingOrder)  {
    if(verbose) print(paste(Sys.time(), "runLayerBinding.fast thisBF =", thisBF, factorSet[[thisBF]]$profile$LAYER.0$pattern, sep=" "))
    theseHits <- matchBindingFactor(newLayerList$layerSet, factorSet[[thisBF]], verbose=verbose)  
    if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits =", length(theseHits), sep=" "))
    #print(length(theseHits))
    if(length(theseHits) < 1) { next ;}
    
    # how many of the potential hits to mark? 
    # iterations/n.factors (rounded up).
    
    hits.sample <- theseHits[sample(1:length(theseHits) ,min(length(theseHits),max.hits))]   # now multiple
    if(verbose) print(paste(Sys.time(), "runLayerBinding.fast n.hits.used =", length(hits.sample), sep=" "))
    #thisHitPosition <- start(hits.sample) + floor(width(hits.sample)/2)
    
    newLayerList$layerSet <-  modifyLayerByBindingFactor.Views(newLayerList$layerSet, hits=hits.sample, bindingFactor=factorSet[[thisBF]], verbose=verbose)
    watch.function(newLayerList$layerSet)
    if(collect.stats) {
      thisRow <- data.frame(bf=thisBF,hits=length(hits.sample) , target.coverage=sum(width(hits.sample)))
      stats.table <- rbind(stats.table, thisRow)
    }
    #print(as.numeric(letterFrequency(newLayerSet$LAYER.1, letters= "1")))
  }
  newLayerList$history <- stats.table
  #print(letterFrequency(newLayerSet$LAYER.1, letters= "1"))   # how many of layer.1 were set to 1
  if(verbose) print(paste(Sys.time(), "runLayerBinding.fast pos 2", sep=" "))
  return(newLayerList)
}

# name.prefix # give new factors a new name beginning with this name 
mutateFactorSet <- function(factorSet, layerSet , mut_type="subRandomFactor", n.muts=1, verbose=FALSE, test.layer0.binding=FALSE,  name.prefix=""){
  newFactorSet <- factorSet
  if(mut_type== "subRandomFactor")  {
    index <- sample(1:length(factorSet), n.muts)
    for(i in index) {
      thisName <- names(factorSet)[i]
      newName <- ifelse(name.prefix == "", thisName, paste(name.prefix, i, sep="."))
      newFactorSet[[thisName]] <-  createRandomBindingFactor(newName,layerSet, type=factorSet[[i]]$type, test.layer0.binding=test.layer0.binding, test.mismatch.rate=.1 ) 
      if(verbose)  {
        print(paste("Factor", i ,"substituted"))
      }
    }
    
  }
  
  return(newFactorSet)
}

#}

# logFile
# logCycle record scores every 'logCycle' iterations
optimiseFactorSet <- function(layerList, factorSet, testing.function, target.layer, target.vec, n.iter=10, target.score=1, mut.rate=0.1, modsPerCycle=100,
                              test.layer0.binding=FALSE, method="fast", logFile="",logCycle=10, maxNoChange=n.iter, verbose=FALSE)  {
  
  if(logFile != "")  write.table(cbind("iter", "best.score"), row.names=F, col.names=F, sep="\t", quote=F,file=logFile)
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
    better <- FALSE
    newFactorSet <- mutateFactorSet(currentFactorSet, layerList$layerSet,n.muts=floor(length(currentFactorSet) * mut.rate), verbose=T, test.layer0.binding=test.layer0.binding, name.prefix=paste("m",i,sep=""))
    if(method =="fast") {
      #print("using fast algorithm")
      newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
    } else {
      #print("using slow algorithm")
      newModLayer <- runLayerBinding(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
    }
    newScore <- test_function(newModLayer, targetLayer=target.layer, target.vec=target.vec)
    scores.vector[i] <- newScore
    print(paste("Round", i, "oldScore", currentBestScore, "newScore", newScore, "Better?"))
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
  currentFactorSet$optimScores <- scores.vector
  return(currentFactorSet)
}




