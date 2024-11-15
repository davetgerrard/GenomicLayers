#' Create a binding factor object to match a given pattern of layers
#'
#' Create a new binding factor based on a simple pattern of 
#' marks on others layers and \emph{may} (when used) set marks on other layers.
#'
#' @param name give the binding factor a name
#' @param type "DNA_regexp"  to differentiate from other types
#' @param forRegExp  = regular expression to match DNA sequence
#' @param revRegExp  = reverse complement of forRegExp
#' @param patternLength   [= 0] an approximate length of match but may not be the actual matched length, which may vary for regular expression matches
#' @param stateWidth the width of pattern to recognise on other layers
#' @param profile.layers a vector of named layers to set as a match
#' @param profile.marks a vector of 0/1 to match the layers in profile.layers
#' @param mod.layers a vector of named layers to alter on a match
#' @param mod.marks a vector of 0/1 to set on the mod.layers
#' @param offset 0  integer value to indicate relative distance from pattern to apply modifications. Very simple. 
#' @param offset.method NULL   a \code{function} to apply to apply offset. MUST have parameter "n" that is used internally to represent the number of hits. 
#' @param offset.params NULL  a \code{list} of named parameters to pass to offset.method function
#' @param test.layer0.binding when creating, test if the DNA sequence has a match.
#' @param test.mismatch.rate proportion of mismatches to tolerate when testing
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
#' simpleBF <- createBindingFactor.DNA_regexp("test", patternString="ACTGGGCTA")
#' 
#' # this regular expression finds 4 CpGs with 0-4 bases between them
#' CGI<- createBindingFactor.DNA_regexp("CGI", patternString="(CG.{0,4}){3}CG", 
#'                           patternLength=20, mod.layers = "CpG_island",
#'                           mod.marks=1, stateWidth=20)
#'                           
#' RAP1 <- createBindingFactor.DNA_regexp("RAP1", forRegExp="GGTGT(.{0,3})GGTGT",
#'                                     revRegExp="ACACC(.{0,3})ACACC",patternLength=20, 
#'                                     mod.layers = "RAP1_bound", mod.marks=1, stateWidth=20)                         
#'
#' @export
createBindingFactor.DNA_regexp <- function(name,  type="DNA_regexp", forRegExp="N",
                                           revRegExp="N", patternLength=0, 
					                              profile.layers=NULL,profile.marks=NULL,
                                        mod.layers=NULL,mod.marks=NULL, 
					                              offset=0, offset.method=NULL, offset.params=NULL, stateWidth=patternLength,
                                      	test.layer0.binding=FALSE, test.mismatch.rate=.1 , max.pattern.tries=1000,
                                      	min.DM.length=2, min.DR.length=10, verbose=FALSE) {
  
  
  # check input  
  stopifnot(exprs = {
    "profile.layers has non-unique names" = length(profile.layers) == length(unique(profile.layers))
    "mod.layers has non-unique names" = length(mod.layers) == length(unique(mod.layers))
  })  
  
  # patternLength will be variable for regular expressions. Need separate parameter for modLength and may become a vector or list with different lengths for each layer.
  #patternLength <- nchar(patternString)
  #TODO sort out how to define patternLength or re-write other functions to accomodate variable patternLength
  profileList <- list(LAYER.0=list(forRegExp=forRegExp , revRegExp=revRegExp, mismatch.rate=0, length=patternLength))
  
  
  if(length(profile.layers) >0) {
    stopifnot("profile.marks does not match length of profile.layers" = length(profile.layers) == length(profile.marks))  
  for(i in 1:length(profile.layers)) {
    thisLayer <- profile.layers[i]
    profileList[[thisLayer]] <- list(pattern=profile.marks[i], mismatch.rate=0.1, length=patternLength)
  }
  }
  modList <- list()
  if(length(mod.layers) >0) {
    stopifnot("mod.marks does not match length of mod.layers" = length(mod.layers) == length(mod.marks))
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



#createBindingFactor.DNA_regexp("test", forRegExp="GGTGT(.{0,3})GGTGT",revRegExp="ACACC(.{0,3})ACACC", patternLength = 10)
