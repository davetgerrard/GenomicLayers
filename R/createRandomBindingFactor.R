#' Randomly generate a binding factor
#'
#' Create a new binding factor that \emph{may}  also require 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#' Need to specify the type of factor. 
#'
#' @param name give the binding factor a name
#' @param type  [c("DNA_motif", "DNA_region","layer_region","layer_island")]  to differentiate from other types
#' @param patternString   ["N"] put motif here (using IUPAC codes for degenerate bases)
#' @param patternLength   length of pattern to be matched [nchar(patternString)]
#' @param stateWidth the width of pattern to recognise on other layers
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate  proportion of mismatches to tolerate when testing [.1]
#' @param max.pattern.tries  1000
#' @param min.DM.length 2
#' @param min.DR.lengt 10
#' @param verbose set to TRUE for more output
#'
#' @return \code{"bindingFactor"}
#' 
#' @seealso \code{\link{runLayerBinding}} \code{\link{createBindingFactor.DNA_regexp}} 
#'
#' @import Biostrings
#' 
#' @examples
#' rBf <- createRandomBindingFactor(name="DNA_A",type="DNA_motif",patternString = "CAT" )
#'
#' @export

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
