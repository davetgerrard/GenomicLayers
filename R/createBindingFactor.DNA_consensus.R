
#' Create a binding factor object to match a given DNA consensus
#'
#' Create a new binding factor based on a DNA consensus that \emph{may}  also require 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#' Makes use of biostrings function \code{\link{vmatchPattern}}
#' Can use IUPAC codes and allow mismatches 
#'
#' @param name give the binding factor a name
#' @param type  "DNA_consensus"  to differentiate from other types
#' @param patternString   "N" put consensus here (using IUPAC codes for degenerate bases)
#' @param patternLength   length of pattern to be matched [nchar(patternString)]
#' @param stateWidth the width of pattern to recognise on other layers
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate  proportion of mismatches to tolerate when testing [.1]
#' @param max.pattern.tries  NA
#' @param min.DM.length NA
#' @param min.DR.lengt NA
#' @param max.mismatch  0   see \code{\link{vmatchPattern}}
#' @param min.mismatch  0   see \code{\link{vmatchPattern}}
#' @param with.indels  FALSE  see \code{\link{vmatchPattern}}
#' @param fixed  TRUE  see \code{\link{vmatchPattern}}
#' @param algorithm  "auto"  see \code{\link{vmatchPattern}}
#' @param verbose set to TRUE for more output
#'
#' @return \code{"bindingFactor"}
#' 
#' @seealso \code{\link{runLayerBinding}} \code{\link{createBindingFactor.DNA_regexp}} 
#'
#' @import Biostrings
#' 
#' @examples
#' DNA_A <- createBindingFactor.DNA_consensus(name="DNA_A",patternString = "CAT" )
#'
#' gc1-2 <- createBindingFactor.DNA_consensus(name="gc1-2",
#'                           patternString="SSSSSSSSSSSSSSSS",
#'                           max.mismatch= 2, 
#'                           min.mismatch= 1,
#'                           fixed="subject",
#'                           profile.layers = NULL,
#'                           profile.marks=NULL,
#'                           mod.layers="LAYER.1", mod.marks = 1)
#' 
#' @export
createBindingFactor.DNA_consensus <- function(name,  type="DNA_consensus", patternString="N",
                                          patternLength = nchar(patternString), stateWidth=patternLength,
                                          profile.layers=NULL,profile.marks=NULL,
                                          mod.layers=NULL,mod.marks=NULL,
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, 
                                      min.DM.length=2, min.DR.length=10, verbose=FALSE,max.mismatch=0, min.mismatch=0,
                                      with.indels=FALSE, fixed=TRUE,
                                      algorithm="auto") {

  # check input  
  stopifnot(exprs = {
    "profile.layers has non-unique names" = length(profile.layers) == length(unique(profile.layers))
    "mod.layers has non-unique names" = length(mod.layers) == length(unique(mod.layers))
    })  

  # create a list to store the profile that would constitute a match. 
  # parameters to pass to Biostrings::vmatchPattern()
  profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=0, 
                                   length=patternLength, max.mismatch=max.mismatch, min.mismatch=min.mismatch,
                                   with.indels=with.indels, fixed=fixed, algorithm=algorithm))
  
  
  if(length(profile.layers) >0) {  # there are layers to match beyond the sequence layer
    stopifnot("profile.marks does not match length of profile.layers" = length(profile.layers) == length(profile.marks))  
    for(i in 1:length(profile.layers)) {
      thisLayer <- profile.layers[i]
      profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
    }
  }
  # now create a second list of intended modifications.
  modList <- list()
  if(length(mod.layers) >0) { 
    stopifnot("mod.marks does not match length of mod.layers" = length(mod.layers) == length(mod.marks)) 
    for(i in 1:length(mod.layers)) {
      #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
      thisLayer <- mod.layers[i]
      modState <- mod.marks[i]
      modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=0, align="centre")   #  make stateWidth independent of patternLength
    }
  }
  
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



#createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA")
#createBindingFactor.DNA_consensus("test", patternString="ACTGGGCTA" , profile.layers=c("bob", "bob"))  # should give error
