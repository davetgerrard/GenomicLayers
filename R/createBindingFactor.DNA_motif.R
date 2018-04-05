
#' Create a binding factor object to match a given DNA motif
#'
#' Create a new binding factor based on a DNA motif that \emph{may}  also require 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#'
#' @param name give the binding factor a name
#' @param type  ["DNA_motif"]  to differentiate from other types
#' @param patternString   ["N"] put motif here (using IUPAC codes for degenerate bases)
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
#' @param verbose set to TRUE for more output
#'
#' @return \code{"bindingFactor"}
#' 
#' @seealso \code{\link{runLayerBinding}} \code{\link{createBindingFactor.DNA_regexp}}
#'
#' @import Biostrings
#' 
#' @examples
#' DNA_A <- createBindingFactor.DNA_motif(name="DNA_A",patternString = "CAT" )
#'
#' @export
createBindingFactor.DNA_motif <- function(name,  type="DNA_motif", patternString="N",
                                          patternLength = nchar(patternString), stateWidth=patternLength,
                                          profile.layers="LAYER.1",profile.marks=0,
                                          mod.layers="LAYER.1",mod.marks=1,
                                      test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000, 
                                      min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  
  profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=0, length=patternLength))
  
  
  if(length(profile.layers) >0) {
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
  }
  }
  modList <- list()
  for(i in 1:length(mod.layers)) {
  #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
    thisLayer <- mod.layers[i]
    modState <- mod.marks[i]
    modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=0, align="centre")   #  make stateWidth independent of patternLength
  }
  
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



#createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")
