#' Create a binding factor object to match a given pattern of layers
#'
#' Create a new binding factor based on a simple pattern of 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#'
#' @param name give the binding factor a name
#' @param type "DNA_regexp"  to differentiate from other types
#' @param patternString  = NOT USED in this case
#' @param patternLength   [= 0] length of pattern
#' @param stateWidth the width of pattern to recognise on other layers
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate % mismatches to tolerate when testing
#' @param max.pattern.tries  NA
#' @param min.DM.length NA
#' @param min.DR.lengt NA
#' @param verbose set to TRUE for more output
#'
#' @return \code{"hits"}
#'
#' @seealso \code{\link{runLayerBinding}} \code{\link{createBindingFactor.DNA_motif}}
#'
#' @examples
#' LAYER_1_2 <-  createBindingFactor.layer_region(name="LAYER_1_2", 
#'                                   profile.layers =c("LAYER.1"), profile.marks = c(1), 
#'                                   mod.layers = c("LAYER.2"), mod.marks = c(1)) 
#'
#' @export
createBindingFactor.DNA_regexp <- function(name,  type="DNA_regexp", patternString="N",patternLength=0, 
					profile.layers=NULL,profile.marks=NULL,
                                        mod.layers=NULL,mod.marks=NULL, stateWidth=patternLength,
                                      	test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000,
                                      	min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  # patternLength will be variable for regular expressions. Need separate parameter for modLength and may become a vector or list with different lengths for each layer.
  #patternLength <- nchar(patternString)
  #TODO sort out how to define patternLength or re-write other functions to accomodate variable patternLength
  profileList <- list(LAYER.0=list(pattern=patternString , mismatch.rate=0, length=patternLength))
  
  
  if(length(profile.layers) >0) {
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
  }
  }
  modList <- list()
  if(length(mod.layers) >0) {
  for(i in 1:length(mod.layers)) {
  #for(thisLayer in sample(names(layerSet)[-1], n.modPatterns, replace=F)) {
    thisLayer <- mod.layers[i]
    modState <- mod.marks[i]
    modList[[thisLayer]] <- list(state=modState, stateWidth=stateWidth, offset=0, align="centre")   # TODO make stateWidth independent of patternLength
  }
  }
  bindingFactor <- list(name=name, type=type,
                        profile=profileList,
                        mods=modList)
  
  return(bindingFactor)
  
}



createBindingFactor.DNA_regexp("test", patternString="ACTGGGCTA")
