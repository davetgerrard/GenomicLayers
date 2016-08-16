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
