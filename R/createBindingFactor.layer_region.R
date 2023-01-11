#' Create a binding factor object to match a given pattern of layers, typically ignoring DNA
#'
#' Create a new binding factor based on a simple pattern of 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#'
#' @param name give the binding factor a name
#' @param type "layer_region"  to differentiate from other types
#' @param patternString  = NOT USED in this case
#' @param patternLength   [= 1] length of pattern
#' @param stateWidth the width of pattern to recognise on other layers, default is same as patternLength
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param offset 0  integer value to indicate relative distance from pattern to apply modifications
#' @param offset.method NULL   a function to apply to apply offset.
#' @param align "centre"
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate proportion of mismatches to tolerate when testing
#' @param max.pattern.tries  NA
#' @param min.DM.length NA
#' @param min.DR.length NA
#' @param verbose set to TRUE for more output
#'
#' @return \code{"hits"}
#'
#' @examples
#' bf.LR <- createBindingFactor.layer_region("layerBf", type="layer_region", 
#'          patternLength = 1, patternString = "N",
#'           stateWidth = patternLength, profile.layers = "LAYER.1",
#'           profile.marks = 0, mod.layers = "LAYER.1", mod.marks = 1)
#'
#' @export
createBindingFactor.layer_region <- function(name,  type="layer_region", patternLength=1, patternString="N",
                                             mismatch.rate=0, stateWidth=patternLength,
                                             profile.layers="LAYER.1",  profile.marks=0,
                                             mod.layers="LAYER.1",mod.marks=1,
                                             offset=0, offset.method=NULL, align="centre",
                                            test.layer0.binding=FALSE, test.mismatch.rate=.1 , 
                                            max.pattern.tries=1000, min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  #patternLength <- nchar(patternString)
  #profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=0, length=patternLength))
  profileList <- list(LAYER.0=list(pattern=DNAString(patternString) , mismatch.rate=mismatch.rate, length=patternLength))
  
  
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=mismatch.rate, length=patternLength)
  }
  
  modList <- list()
  for(i in 1:length(mod.layers)) {
  #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
    thisLayer <- mod.layers[i]
    modState <- mod.marks[i]
    modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=offset, offset.method=offset.method,align=align)   # TODO make stateWidth independent of patternLength
  }
  
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



#createBindingFactor.DNA_motif("test", patternString="ACTGGGCTA")
