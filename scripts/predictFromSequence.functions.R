matchBindingFactor <- function(layerSet, bindingFactor, clusterGap=10)  {
  require(Biostrings)
  hitList <- list()
  #validHits <- 
  for(thisLayer in names(bindingFactor$profile)) {
    thisPattern <- bindingFactor$profile[[thisLayer]]$pattern
    max.mismatches <- round(bindingFactor$profile[[thisLayer]]$mismatch.rate * nchar(thisPattern))
    if(thisLayer == "LAYER.0") {  
      
      # currently DNA matches using fixed length string with IUPAC codes (adapt later for pwm matching)
      hitList[[thisLayer]] <-  matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]], fixed=FALSE, max.mismatch= max.mismatches) # allowas matching with IUPAC codes
      #validHits <- hitList[[thisLayer]]
    } else {
      hitList[[thisLayer]] <-  matchPattern(bindingFactor$profile[[thisLayer]]$pattern,layerSet[[thisLayer]], fixed=TRUE )   # TODO: how to do fixed=F for BStrings? 
      #validHits <- union(validHits, )
    }
  }
  
  if(any(lapply(hitList, length) == 0))  {   # one or more of the patterns were not matched.
    validHits <- list()
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

# utility function to ensure value is within range
truncateValueWithinRange <- function(range, value)  {
  newValue <- max(value, min(range))
  newValue <- min(newValue, max(range))
  return(newValue)
}

# How to modify the layers
modifyLayer <- function(layerSet, position, layer, state="1", offset=0, align="centre")  {
  layerSet[[layer]][position]  <- state
  return(layerSet)
}

modifyLayerByBindingFactor <- function(layerSet, position, bindingFactor, verbose=FALSE) {
  
  newLayerSet <- layerSet
  for(thisLayer in names(bindingFactor$mods))  {
    thisState <- bindingFactor$mods[[thisLayer]]$state
    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    thisOffset <- bindingFactor$mods[[thisLayer]]$offset
    thisStart <- position - floor(stateWidth/2) + thisOffset   # using midpoint of state
    thisEnd <- thisStart + stateWidth
    seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
    thisStart <- truncateValueWithinRange(seqRange, thisStart)
    thisEnd <- truncateValueWithinRange(seqRange, thisEnd)
    # TODO: design decisions on whether to encode modifications as strings or dimensions: e.g. '111111111' or paste(rep('1', 8), collapse="")
    #     currently using the former
    # DECISION FORCED by having to deal with negative ranges after offset. 
    # need to have state as single value
    #thisPosition  <- position + thisOffset
    # TODO: add code for using centre/left/right for adjustment
    newLayerSet <- modifyLayer(layerSet=newLayerSet, position=thisStart:thisEnd, layer=thisLayer, state =thisState)
  }
  return(newLayerSet)
}

# as modifyLayerByBindingFactor but modify in multiple places with the same patterning.
modifyLayerByBindingFactor.multiHits <- function(layerSet, position.vec, bindingFactor, verbose=FALSE) {
  
  newLayerSet <- layerSet
  for(thisLayer in names(bindingFactor$mods))  {
    seqRange <- c(start(layerSet[['LAYER.0']])[1], end(layerSet[['LAYER.0']])[1])
    thisState <- bindingFactor$mods[[thisLayer]]$state
    stateWidth <- bindingFactor$mods[[thisLayer]]$stateWidth
    thisOffset <- bindingFactor$mods[[thisLayer]]$offset
    thisStart.vec <- position.vec - floor(stateWidth/2) + thisOffset   # using midpoint of state
    mod.vec <- integer()
    for(thisStart in thisStart.vec) {
      thisEnd <- thisStart + stateWidth
      thisStart <- truncateValueWithinRange(seqRange, thisStart)
      thisEnd <- truncateValueWithinRange(seqRange, thisEnd)
      mod.vec <- c(mod.vec, thisStart:thisEnd)
    }
    # TODO: design decisions on whether to encode modifications as strings or dimensions: e.g. '111111111' or paste(rep('1', 8), collapse="")
    #     currently using the former
    # DECISION FORCED by having to deal with negative ranges after offset. 
    # need to have state as single value
    #thisPosition  <- position + thisOffset
    # TODO: add code for using centre/left/right for adjustment
    newLayerSet <- modifyLayer(layerSet=newLayerSet, position=mod.vec, layer=thisLayer, state =thisState)
  }
  return(newLayerSet)
}

# might want this to store some metrics along the way. 
# perhaps the layerSet class should carry some of these?
runLayerBinding <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE)  {
  
  if(verbose) print(paste(Sys.time, "runLayerBinding pos 1", sep="\t"))
  bindingOrder <- sample(names(factorSet), size=iterations,prob=bindingFactorFreqs, replace=T)
  newLayerList <- layerList
  
   if(collect.stats) {
	stats.table <- data.frame()
	} else {
	stats.table <- NULL
	}  

  for(thisBF in bindingOrder)  {
    if(verbose) print(paste(Sys.time, "runLayerBinding thisBF =", thisBF, sep="\t"))
    theseHits <- matchBindingFactor(newLayerList$layerSet, factorSet[[thisBF]])  
    #print(length(theseHits))
    if(length(theseHits) < 1) { next ;}
    # might it be better/faster to just determine all the potential hits once?
    # NO, because the order determines whether some features can bind later on.
    #choose one hit
    thisHit <- theseHits[sample(1:length(theseHits) , 1)]
    thisHitPosition <- start(thisHit) + floor(width(thisHit)/2)
    
    newLayerList$layerSet <-  modifyLayerByBindingFactor(newLayerList$layerSet, position=thisHitPosition, bindingFactor=factorSet[[thisBF]])
	watch.function(newLayerList$layerSet)
     if(collect.stats) {
	thisRow <- data.frame(bf=thisBF,position=thisHitPosition , target.coverage=as.numeric(letterFrequency(newLayerList$layerSet[[target.layer]], letters="1")))
	stats.table <- rbind(stats.table, thisRow)
	}
    #print(as.numeric(letterFrequency(newLayerSet$LAYER.1, letters= "1")))
  }
	newLayerList$history <- stats.table
  #print(letterFrequency(newLayerSet$LAYER.1, letters= "1"))   # how many of layer.1 were set to 1
  if(verbose) print(paste(Sys.time, "runLayerBinding pos 2", sep="\t"))
  return(newLayerList)
}


# might want this to store some metrics along the way. 
# perhaps the layerSet class should carry some of these?
runLayerBinding.fast <- function(layerList, factorSet, iterations=1, bindingFactorFreqs=rep(1, length(factorSet)), watch.function=function(x){}, collect.stats=FALSE, target.layer=2, verbose=FALSE)  {
  if(verbose) print(paste(Sys.time, "runLayerBinding.fast pos 1", sep="\t"))
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
    if(verbose) print(paste(Sys.time, "runLayerBinding.fast thisBF =", thisBF, sep="\t"))
    theseHits <- matchBindingFactor(newLayerList$layerSet, factorSet[[thisBF]])  
    #print(length(theseHits))
    if(length(theseHits) < 1) { next ;}

    # how many of the potential hits to mark? 
    # iterations/n.factors (rounded up).
    
    thisHit <- theseHits[sample(1:length(theseHits) ,min(length(theseHits),max.hits))]   # now multiple
    thisHitPosition <- start(thisHit) + floor(width(thisHit)/2)
    
    newLayerList$layerSet <-  modifyLayerByBindingFactor.multiHits(newLayerList$layerSet, position=thisHitPosition, bindingFactor=factorSet[[thisBF]])
    watch.function(newLayerList$layerSet)
    if(collect.stats) {
      thisRow <- data.frame(bf=thisBF,position=thisHitPosition , target.coverage=as.numeric(letterFrequency(newLayerList$layerSet[[target.layer]], letters="1")))
      stats.table <- rbind(stats.table, thisRow)
    }
    #print(as.numeric(letterFrequency(newLayerSet$LAYER.1, letters= "1")))
  }
  newLayerList$history <- stats.table
  #print(letterFrequency(newLayerSet$LAYER.1, letters= "1"))   # how many of layer.1 were set to 1
  if(verbose) print(paste(Sys.time, "runLayerBinding.fast pos 2", sep="\t"))
  return(newLayerList)
}

print.bfSet <- function(factorSet) {
  print("A list of binding factors:-")
  bf.names <- names(factorSet)
  types <- lapply(factorSet, FUN=function(x) x$type)
  dna.patterns <- unlist(lapply(factorSet, FUN=function(x) as.character(x$profile$LAYER.0$pattern)))
  layer.patterns <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$profile))), paste, collapse=","))
  layer.mods <- unlist(lapply(lapply(factorSet, FUN=function(x) as.character(names(x$mods))), paste, collapse=","))
  as.data.frame(cbind(bf.names, types,dna.patterns, layer.patterns, layer.mods) )
}


#createRandomDNABindingPattern <- function(length, homogeneity, alphabet)
#createRandomBinaryBindingPattern <- function(length, homogeneity, alphabet)
    

# how to classify different types of patterns?
# e.g.  CCCCCCCCCC versus GCTGAGCAAATGCG  

createRandomBindingFactor <- function(name, layerSet, type=c("DNA_motif", "DNA_region","layer_region","layer_island"),
	 test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, min.DM.length=2, min.DR.length=10) {
  
   if(!test.layer0.binding) max.pattern.tries <- 1  # only try one pattern.
  
  if(type == "DNA_motif") {
    
    for(i in 1:max.pattern.tries) {
      patternLength <- max(min.DM.length,rnbinom(1, 50, mu= 12))  # patterns must be at least of length 1
      pattern <- paste(sample(names(IUPAC_CODE_MAP), patternLength, replace=T), collapse="")
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

    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    openPattern <- paste(rep("1", patternLength  ), collapse="")
    closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      pattern <- unlist(sample(c(openPattern,closedPattern),1))
      profileList[[thisLayer]] <- list(pattern=BString(pattern), mismatch.rate=0.1)
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

    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    openPattern <- paste(rep("1", patternLength  ), collapse="")
    closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      pattern <- unlist(sample(c(openPattern,closedPattern),1))
      profileList[[thisLayer]] <- list(pattern=BString(pattern), mismatch.rate=0.1)
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
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    openPattern <- paste(rep("1", patternLength  ), collapse="")
    closedPattern <- paste(rep("0", patternLength  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      pattern <- unlist(sample(c(openPattern,closedPattern),1))
      profileList[[thisLayer]] <- list(pattern=BString(pattern), mismatch.rate=0.1)
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
    
    profileList <- list(LAYER.0=list(pattern=DNAString(pattern) , mismatch.rate=0.1))
    
    n.layerPatterns <- sample(0:(length(layerSet)-1), 1, prob=1/(1:length(layerSet)), replace=T)   # higher numbers successively less likely
    openPattern <- paste(rep(c("0","1","0"), times=c(shoulderLength, islandLength, shoulderLength)  ), collapse="")
    closedPattern <- paste(rep(c("1","0","1"), times=c(shoulderLength, islandLength, shoulderLength)  ), collapse="")
    for(thisLayer in sample(names(layerSet)[-1], n.layerPatterns, replace=F)) {
      pattern <- unlist(sample(c(openPattern,closedPattern),1))
      profileList[[thisLayer]] <- list(pattern=BString(pattern), mismatch.rate=0.1)
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

#createRandomBindingFactor("testRandom", layerSet, type="DNA_motif", test.layer0.binding=FALSE, test.mismatch.rate=.1 ) 


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
    initialModLayer <- runLayerBinding.fast(layerList=layerList, factorSet = factorSet, iterations=modsPerCycle, verbose=verbose)
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
      newModLayer <- runLayerBinding.fast(layerList=layerList, factorSet = newFactorSet, iterations=modsPerCycle, verbose=verbose)
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


# measure correlation between tow binary layers

cor.layers <- function(x,y,method="pearson")  {
  require(Biostrings)
  layer.vector.x  <- as.numeric(strsplit(as.character(x),"")[[1]]) 
  layer.vector.y  <- as.numeric(strsplit(as.character(y),"")[[1]]) 
  stopifnot(length(layer.vector.x) == length(layer.vector.y))
  cor(layer.vector.x, layer.vector.y)
}


# create table from a factor set.
table.factorSet <- function(factorSet, classes="",patternTrimLength=20)  {
  #type.table <- data.frame(factor.name=names(factorSet))
  type.table <- data.frame()
  for(i in 1:length(factorSet)) {
    thisFactor <- names(factorSet)[i]
    thisType <- factorSet[[thisFactor]]$type
    thisName <- factorSet[[thisFactor]]$name
    names(factorSet[[thisFactor]]$profile)
    for(thisLayer in names(factorSet[[thisFactor]]$profile)) {
      mismatchRate <- factorSet[[thisFactor]]$profile[[thisLayer]]$mismatch.rate
      raw.pattern <- factorSet[[thisFactor]]$profile[[thisLayer]]$pattern
      if(nchar(raw.pattern) > patternTrimLength) {
        pattern <- paste(as.character(raw.pattern[1:patternTrimLength]), "...", sep="")
      } else {
        pattern <- as.character(raw.pattern)
      }
      thisRow <- data.frame(id=thisFactor, name=thisName,type=thisType, group="target",layer=thisLayer, pattern=pattern, mm.rate=mismatchRate, state=NA, width=nchar(raw.pattern), offset=NA, align=NA )
      type.table <- rbind(type.table, thisRow)
    }
    
    for(thisLayer in names(factorSet[[thisFactor]]$mods)) {
      thisState <- factorSet[[thisFactor]]$mods[[thisLayer]]$state
      thisStateWidth <- factorSet[[thisFactor]]$mods[[thisLayer]]$stateWidth
      thisOffset <- factorSet[[thisFactor]]$mods[[thisLayer]]$offset
      thisAlign <- factorSet[[thisFactor]]$mods[[thisLayer]]$align
      thisRow <- data.frame(id=thisFactor, name=thisName,type=thisType, group="mod",layer=thisLayer, pattern=NA, mm.rate=NA, state=thisState, width=thisStateWidth, offset=thisOffset, align=thisAlign )
      type.table <- rbind(type.table, thisRow)
    }
  }
  return(type.table)
}

# following table.factorSet(), create a matrix with one column per target/mod and layer
# can be used to classify factor, group similar and to produce eulerGrids.
type.patterns <- function(df, id.col=1) {
  layers <- unique(df[,'layer'])
  groups <- unique(df[,'group'])
  colNames <- character()
  for(thisGroup in groups) {
    for(thisLayer in layers) {
      colNames <- c(colNames, paste(thisGroup, thisLayer, sep="_"))
    }
  }
  pattern.table <- matrix(0, nrow=length(unique(df[,id.col])), ncol=length(colNames), dimnames=list(unique(df[,id.col]), colNames))
  for(thisGroup in groups) {
    for(thisLayer in layers) {
      factors <- as.character(subset(df, group==thisGroup & layer==thisLayer)[,id.col])
      pattern.table[factors,paste(thisGroup, thisLayer, sep="_")] <- 1
    }
  }
  return(pattern.table)
}
