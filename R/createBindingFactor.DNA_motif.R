
#' Find matches for a binding factor on a layer set
#'
#' Generate a list of matches for a binding factor against a layerSet object. 
#'
#' @param name method to do something to (\code{"hsv"} or \code{"cluster"})
#' @param type description of that param
#' @param patternString  =10 you get the idea
#' @param patternLength    =10000000 you get the idea
#' @param stateWidth you get the idea
#' @param profile.layers you get the idea
#' @param profile.marks you get the idea
#' @param mod.layers you get the idea
#' @param mod.marks you get the idea
#' @param test.layer0.binding you get the idea
#' @param test.mismatch.rate you get the idea
#' @param max.pattern.tries you get the idea
#' @param min.DM.length you get the idea
#' @param min.DR.lengt you get the idea
#' @param verbose you get the idea
#'
#' @return \code{"hits"}
#'
#' @examples
#' x <- 1   # great!
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
